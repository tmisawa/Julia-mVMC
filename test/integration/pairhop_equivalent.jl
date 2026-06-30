# test/integration/pairhop_equivalent.jl
#
# PairHop parity gate against C-mVMC. The fixtures are one-step Hubbard
# para-opt runs derived from the existing real and FSZ Hubbard references, with
# one PairHop input row. C expands that row to both (i,j) and (j,i); the parser
# and Hamiltonian path must reproduce the same first zvo_out row.

using Test
using MVMCOptimizers

const PAIRHOP_TOL_DEFAULT = 1e-10
const PAIRHOP_TOL_LOOSE = 1e-9
const PAIRHOP_LOOSE_COLS = (3, 4)
const PAIRHOP_EXPECTED_NCOLS = 6

const PAIRHOP_MODELS = [
    (fixture = "hubbard_chain_pairhop_real", mode = :real),
    (fixture = "hubbard_chain_pairhop_fsz", mode = :fsz),
]

function parse_pairhop_zvo_first1(path::AbstractString)
    return parse.(Float64, split(strip(readline(path))))
end

function run_pairhop_model(model)
    refdir = joinpath(@__DIR__, "reference", model.fixture)
    namelist = joinpath(refdir, "inputs", "namelist.def")
    refpath = joinpath(refdir, "zvo_out_first1.dat")

    @assert isfile(namelist) "namelist.def missing for $(model.fixture)"
    @assert isfile(refpath) "zvo_out_first1.dat missing for $(model.fixture)"

    result = run_para_opt_from_namelist(
        namelist;
        nsteps = 1,
        nsmp = 1,
        mode = model.mode,
        output_dir = tempname(),
    )
    return (actual = parse.(Float64, split(strip(result.zvo_first_n[1]))),
            expected = parse_pairhop_zvo_first1(refpath))
end

@testset "Julia-mVMC PairHop e2e (C reference)" begin
    for model in PAIRHOP_MODELS
        @testset "$(model.fixture) ($(model.mode))" begin
            run = run_pairhop_model(model)
            @test length(run.actual) == PAIRHOP_EXPECTED_NCOLS
            @test length(run.expected) == PAIRHOP_EXPECTED_NCOLS
            for idx in eachindex(run.actual)
                tol = idx in PAIRHOP_LOOSE_COLS ? PAIRHOP_TOL_LOOSE : PAIRHOP_TOL_DEFAULT
                @test abs(run.actual[idx] - run.expected[idx]) <= tol
            end
        end
    end
end
