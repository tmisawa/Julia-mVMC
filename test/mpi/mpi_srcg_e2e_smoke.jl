# mpiexec 配下で NSRCG=1 の vmc_para_opt! end-to-end path を C rank2
# reference と比較する worker。
# 使い方: julia --project=<workspace> test/mpi/mpi_srcg_e2e_smoke.jl <output_dir>
using Test
using MPI
using MVMCOptimizers

const refdir = joinpath(@__DIR__, "..", "integration", "reference",
                        "heisenberg_chain_real_nsrcg")
const fixture = joinpath(refdir, "inputs", "namelist.def")
const outdir = ARGS[1]
const expected_ncols = 6
const tol_default = 1e-10
const tol_loose = 1e-9
const loose_cols = (3, 4)
# Same tolerance class as the serial NSRCG=1 parameter gate: truncated SR-CG is
# sensitive to BLAS/reduction order even when the energy row remains tight.
const nsrcg_param_tol = 1e-2

function parse_zvo_row(path)
    row = parse.(Float64, split(strip(readline(path))))
    @test length(row) == expected_ncols
    return row
end

function parse_complex_pairs(path)
    vals = parse.(Float64, split(read(path, String)))
    @test iseven(length(vals))
    return [ComplexF64(vals[2*i-1], vals[2*i]) for i = 1:div(length(vals), 2)]
end

function parse_orbital_parameter_indices(path)
    idxs = Int[]
    for line in eachline(path)
        parts = split(strip(line))
        length(parts) == 3 || continue
        parsed = tryparse.(Int, parts)
        any(isnothing, parsed) && continue
        push!(idxs, parsed[3] + 1)  # C idx -> Julia 1-based parameter slot
    end
    return idxs
end

try
    result = MVMCOptimizers.run_para_opt_from_namelist(
        fixture; nsteps = 1, nsmp = 1, mode = :real, output_dir = outdir)

    if isempty(result.zvo_first_n)
        @test isnan(result.final_energy_per_site)
        println("srcg-e2e worker: non-root rank ok")
    else
        @test result.status == 0
        @test length(result.zvo_first_n) == 1

        ours_vals = parse.(Float64, split(strip(result.zvo_first_n[1])))
        ref_vals = parse_zvo_row(joinpath(refdir, "zvo_out_mpi2_first1.dat"))
        @test length(ours_vals) == expected_ncols
        for j in eachindex(ours_vals)
            tol = j in loose_cols ? tol_loose : tol_default
            @test abs(ours_vals[j] - ref_vals[j]) <= tol
        end

        @test strip.(readlines(joinpath(outdir, "zvo_SRinfo.dat"))) ==
              strip.(readlines(joinpath(refdir, "zvo_SRinfo_mpi2_1step.dat")))

        c_params = parse_complex_pairs(joinpath(refdir, "zqp_opt_mpi2_1step.dat"))[3:end]
        julia_pairs = parse_complex_pairs(joinpath(outdir, "zqp_opt.dat"))
        orbital_idxs = parse_orbital_parameter_indices(joinpath(refdir, "inputs", "orbitalidx.def"))

        @test length(c_params) == 14
        @test length(julia_pairs) == 2 + length(orbital_idxs)
        @test maximum(abs.(julia_pairs[1:2] .- c_params[1:2])) <= tol_default
        orbital_maxdiff = maximum(
            abs(julia_pairs[2+k] - c_params[2+orbital_idxs[k]]) for k in eachindex(orbital_idxs)
        )
        @test orbital_maxdiff <= nsrcg_param_tol

        println("srcg-e2e worker: root rank ok")
    end
catch err
    @error "mpi_srcg_e2e_smoke worker failed" exception = (err, catch_backtrace())
    if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD, 1)
    end
    rethrow()
end
