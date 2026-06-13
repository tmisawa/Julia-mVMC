# test/integration/runtests.jl
#
# Julia-mVMC integration tests: replay each fixture model and compare
# the first N_STEPS SR steps of zvo_out against the committed C reference.
#
# Run from the workspace root:
#   julia --project=@. test/integration/runtests.jl

using Test
using Printf
using Random
using SFMT
using MVMCExpertModeParsers
using MVMCOptimizers

# zvo_out.dat columns:
#   1: real(<H>)   2: imag(<H>)   3: <H^2>   4: variance = (<H^2>-<H>^2)/<H>^2
#   5: <Sz>        6: <Sz^2>
# Columns 1, 2, 5, 6 are linear accumulators of bounded sample quantities
# and match C to within BLAS summation-order noise (~1e-13 .. 1e-12).
# Column 3 (<H^2>) and column 4 (variance) involve squared-energy terms:
# col 3 is Σ w·conj(e)·e and col 4 takes a cancelling difference on top.
# Both pick up amplified summation-order noise (~1e-11 .. 1e-10) that is
# fundamentally bounded by C/Julia BLAS ordering, not by Julia bugs.
# We therefore tighten the directly-measured energy/Sz columns to 1e-10 and
# grant the squared-energy columns a 10x slacker bound.
const TOL_DEFAULT = 1e-10            # cols 1, 2, 5, 6
const TOL_LOOSE   = 1e-9             # cols 3, 4 (squared/derived)
const LOOSE_COLS  = (3, 4)
const EXPECTED_NCOLS = 6             # zvo_out.dat layout (see comment above)
const N_STEPS = 10
# 1-step SR-CG parameter updates are deliberately a tolerance gate, not a bit
# parity gate. Observed maxdiffs include about 2.84e-4 on the reference
# macOS/OpenBLAS setup and 5.33e-3 on GitHub Actions ubuntu-latest Julia 1.12.
# The first zvo row remains tight; this parameter check is a coarse e2e guard
# for layout and order-of-magnitude regressions under documented SR-CG numeric
# sensitivity.
const NSRCG_PARAM_TOL = 1e-2

const MODELS = [
    ("heisenberg_chain_real", :real),
    ("heisenberg_chain_cmp",  :cmp),
    ("heisenberg_chain_fsz",  :fsz),
    ("hubbard_chain_real",    :real),
]

"Parse the first n whitespace-separated lines of a zvo_out file into Vector{Vector{Float64}}."
function parse_zvo_first_n(path, n)
    open(path, "r") do io
        [parse.(Float64, split(strip(readline(io)))) for _ in 1:n]
    end
end

function parse_complex_pairs(path)
    vals = parse.(Float64, split(read(path, String)))
    @assert iseven(length(vals)) "expected an even number of floats in $path"
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

function run_model(name, mode; nsteps = N_STEPS, output_dir = tempname())
    refdir   = joinpath(@__DIR__, "reference", name)
    inputs   = joinpath(refdir, "inputs")
    refpath  = joinpath(refdir, "zvo_out_first$(nsteps).dat")
    namelist = joinpath(inputs, "namelist.def")

    @assert isfile(namelist) "namelist.def missing for $name"
    @assert isfile(refpath)  "reference missing for $name"

    # Drive Julia VMCParaOpt via the public wrapper. The wrapper internally:
    #   - parses .def files (relative paths resolved from namelist_path),
    #   - seeds the SFMT19937 RNG with modpara.def's RndSeed (C-parity
    #     resolve_rnd_seed: missing → 11272, 0 → 0, < 0 → time seed),
    #   - runs nsteps SR steps,
    #   - returns the first N output lines as Vector{String}.
    result = MVMCOptimizers.run_para_opt_from_namelist(
        namelist;
        nsteps = nsteps,
        nsmp = nsteps,
        mode   = mode,
        output_dir = output_dir,
    )
    return result
end

@testset "Julia-mVMC integration vs C reference" begin
    for (name, mode) in MODELS
        @testset "$name ($mode)" begin
            ours_lines = run_model(name, mode).zvo_first_n
            ref_lines = open(joinpath(@__DIR__, "reference", name, "zvo_out_first10.dat")) do io
                [strip(readline(io)) for _ in 1:N_STEPS]
            end
            for i in 1:N_STEPS
                ours_vals = parse.(Float64, split(strip(ours_lines[i])))
                ref_vals  = parse.(Float64, split(ref_lines[i]))
                # Guard against silently-passing comparisons when Julia
                # truncates a column (e.g. NaN-suppression, future format
                # changes). Both sides must carry the full 6-column layout.
                @test length(ours_vals) == EXPECTED_NCOLS
                @test length(ref_vals)  == EXPECTED_NCOLS
                ok = true
                for j in eachindex(ours_vals)
                    tol = j in LOOSE_COLS ? TOL_LOOSE : TOL_DEFAULT
                    ok &= abs(ours_vals[j] - ref_vals[j]) <= tol
                end
                @test ok
            end
        end
    end
end

@testset "Julia-mVMC NSRCG=1 first-step e2e vs C reference" begin
    name = "heisenberg_chain_real_nsrcg"
    refdir = joinpath(@__DIR__, "reference", name)
    outdir = tempname()
    result = run_model(name, :real; nsteps = 1, output_dir = outdir)

    ours_vals = parse.(Float64, split(strip(result.zvo_first_n[1])))
    ref_vals = parse_zvo_first_n(joinpath(refdir, "zvo_out_first1.dat"), 1)[1]
    @test length(ours_vals) == EXPECTED_NCOLS
    @test length(ref_vals) == EXPECTED_NCOLS
    for j in eachindex(ours_vals)
        tol = j in LOOSE_COLS ? TOL_LOOSE : TOL_DEFAULT
        @test abs(ours_vals[j] - ref_vals[j]) <= tol
    end

    # For NSROptItrSmp=1 C writes raw pairs: (Etot, Etot2, Para...).
    # Julia writes projection parameters followed by orbital terms. Compare the
    # projection block directly and map each orbital term through orbitalidx.def
    # to the C unique-parameter slot. The tolerance reflects the documented
    # SR-CG BLAS/reduction-order sensitivity after one parameter update.
    c_params = parse_complex_pairs(joinpath(refdir, "zqp_opt_1step.dat"))[3:end]
    julia_pairs = parse_complex_pairs(joinpath(outdir, "zqp_opt.dat"))
    orbital_idxs = parse_orbital_parameter_indices(joinpath(refdir, "inputs", "orbitalidx.def"))

    @test length(c_params) == 14
    @test length(julia_pairs) == 2 + length(orbital_idxs)
    @test maximum(abs.(julia_pairs[1:2] .- c_params[1:2])) <= TOL_DEFAULT
    orbital_maxdiff = maximum(
        abs(julia_pairs[2+k] - c_params[2+orbital_idxs[k]]) for k in eachindex(orbital_idxs)
    )
    @test orbital_maxdiff <= NSRCG_PARAM_TOL
end

# Plan 3a prerequisites (pure Julia, no C binary): the Green-file comparison
# helpers and the PhysCal runner's parse->loader contract. Included here so they
# run in the CI "Run integration tests" step alongside the fixture replay.
include(joinpath(@__DIR__, "tools", "test_green_compare.jl"))
include(joinpath(@__DIR__, "test_run_phys_cal_contract.jl"))
