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

function run_model(name, mode)
    refdir   = joinpath(@__DIR__, "reference", name)
    inputs   = joinpath(refdir, "inputs")
    refpath  = joinpath(refdir, "zvo_out_first10.dat")
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
        nsteps = N_STEPS,
        nsmp = N_STEPS,
        mode   = mode,
    )
    return result.zvo_first_n
end

@testset "Julia-mVMC integration vs C reference" begin
    for (name, mode) in MODELS
        @testset "$name ($mode)" begin
            ours_lines = run_model(name, mode)
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

# Plan 3a prerequisites (pure Julia, no C binary): the Green-file comparison
# helpers and the PhysCal runner's parse->loader contract. Included here so they
# run in the CI "Run integration tests" step alongside the fixture replay.
include(joinpath(@__DIR__, "tools", "test_green_compare.jl"))
include(joinpath(@__DIR__, "test_run_phys_cal_contract.jl"))
