# test/integration/srcg_option_equivalent.jl
#
# C mVMC 1.4.0 SR-CG option contract (Julia-mVMC#95): first SR-CG step of the
# hand-authored NSRCG fixtures with NSRCGAbortOnFail / NSRCGFallback / NSRCG 2
# (DiagScale) / RescaleSmat variants, compared with C 1.4.0 references generated
# by tools/regenerate_c_reference.sh (provenance: reference/C_1.4.0_PROVENANCE.txt).
#
# Exact: C exit (abort vs continue), the SR-CG warning/error lines, and the
# integer zvo_SRinfo.dat fields (Npara Msize optCut diagCut imax info), i.e. the
# discrete CG status, fallback and cut decisions. Numerical: explicit tolerances.

using Test
using MVMCOptimizers

const SRCG_REF = joinpath(@__DIR__, "reference")
# First zvo_out row: computed before the SR step; same budget as runtests.jl.
const SRCG_ZVO_TOL = 1e-10
const SRCG_ZVO_LOOSE = 1e-9
# S diagonal extrema: sample averages of O*O (no solver involved), relative.
const SRCG_SDIAG_RTOL = 1e-9
# CG increments (absRmax, updated parameters) after a truncated or stopped CG are
# sensitive to BLAS/reduction order (see runtests.jl NSRCG_PARAM_TOL); keep the
# same coarse e2e budget. Observed differences are printed with @info.
const SRCG_PARAM_TOL = 1e-2

parse_pairs(path) = (v = parse.(Float64, split(read(path, String)));
                     [ComplexF64(v[2i-1], v[2i]) for i in 1:div(length(v), 2)])

function srinfo_row(path)
    lines = filter(l -> !startswith(l, "#"), readlines(path))
    @assert length(lines) == 1 "expected one SRinfo row in $path"
    t = split(lines[1])
    ints = parse.(Int, t[[1, 2, 3, 4, 8, 9]])
    floats = parse.(Float64, t[5:7])
    return ints, floats
end

function run_variant(model, variant)
    refdir = joinpath(SRCG_REF, model)
    vdir = joinpath(refdir, "variants", variant)
    work = mktempdir()
    cp(joinpath(refdir, "inputs"), joinpath(work, "inputs"))
    cp(joinpath(vdir, "modpara.def"), joinpath(work, "inputs", "modpara.def"); force = true)
    out = joinpath(work, "output")
    errlog = joinpath(work, "stderr.txt")
    thrown = nothing
    open(errlog, "w") do io
        redirect_stderr(io) do
            try
                MVMCOptimizers.run_para_opt_from_namelist(
                    joinpath(work, "inputs", "namelist.def");
                    nsteps = 1, nsmp = 1, mode = :real, output_dir = out)
            catch e
                thrown = e
            end
        end
    end
    sr_lines = filter(l -> occursin(r"^(warning|error|remark): .*SR-", l), readlines(errlog))
    return (; vdir, out, thrown, sr_lines)
end

@testset "C 1.4.0 SR-CG option variants (first step)" begin
    for (model, variants) in (
        ("heisenberg_chain_real_nsrcg",
         ("abort_default", "fallback", "diag_scale", "nsrcg2_fallback", "rescale_smat")),
        ("hubbard_chain_real_nsrcg", ("plain", "diag_scale", "rescale_smat")),
    ), variant in variants
        @testset "$model/$variant" begin
            r = run_variant(model, variant)
            c_exit = parse(Int, strip(read(joinpath(r.vdir, "exit_status.txt"), String)))
            c_lines = filter(l -> occursin("SR-", l), readlines(joinpath(r.vdir, "stderr.txt")))
            @test r.sr_lines == c_lines
            if c_exit != 0
                # C MPI_Abort: Julia throws with the same final message.
                @test r.thrown isa ErrorException
                @test r.thrown !== nothing && r.thrown.msg == last(c_lines)
                @test !isfile(joinpath(r.out, "zvo_SRinfo.dat")) ||
                      count(!startswith("#"), readlines(joinpath(r.out, "zvo_SRinfo.dat"))) == 0
                continue
            end
            @test r.thrown === nothing
            r.thrown === nothing || continue
            c_int, c_float = srinfo_row(joinpath(r.vdir, "zvo_SRinfo.dat"))
            j_int, j_float = srinfo_row(joinpath(r.out, "zvo_SRinfo.dat"))
            @test j_int == c_int
            @test isapprox(j_float[1], c_float[1]; rtol = SRCG_SDIAG_RTOL)
            @test isapprox(j_float[2], c_float[2]; rtol = SRCG_SDIAG_RTOL, atol = 1e-14)
            @test abs(j_float[3] - c_float[3]) <= SRCG_PARAM_TOL
            @test readlines(joinpath(r.out, "zvo_SRinfo.dat"))[1] ==
                  readlines(joinpath(r.vdir, "zvo_SRinfo.dat"))[1]
            c_row = parse.(Float64, split(readline(joinpath(r.vdir, "zvo_out_first1.dat"))))
            j_row = parse.(Float64, split(readline(joinpath(r.out, "zvo_out.dat"))))
            @test length(j_row) == length(c_row)
            for j in eachindex(c_row)
                @test abs(j_row[j] - c_row[j]) <= (j in (3, 4) ? SRCG_ZVO_LOOSE : SRCG_ZVO_TOL)
            end
            c_par = parse_pairs(joinpath(r.vdir, "zqp_opt_1step.dat"))
            j_par = parse_pairs(joinpath(r.out, "zqp_opt.dat"))
            @test length(j_par) == length(c_par)
            if length(j_par) == length(c_par)
                maxdiff = maximum(abs.(j_par[3:end] .- c_par[3:end]))
                @info "SR-CG variant parameter difference" model variant maxdiff absRmax_diff = j_float[3] - c_float[3]
                @test maxdiff <= SRCG_PARAM_TOL
            end
        end
    end
end
