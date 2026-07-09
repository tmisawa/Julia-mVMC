# test/integration/lanczos_equivalent.jl
#
# Julia-mVMC Full Lanczos PhysCal end-to-end integration gate.
#
# For each committed system, this runner:
#   1. loads fixed variational parameters from the committed C `zqp_opt.dat`,
#   2. runs `run_phys_cal_from_namelist` on a committed PhysCal input set with
#      NLanczosMode=2 into a fresh tempdir, and
#   3. compares the produced Lanczos files against committed C-mVMC references:
#        - zvo_ls_out_001.dat   R1 energy, norm, and alpha
#        - zvo_ls_qqqq_001.dat  16 flattened QQQQ moments
#        - zvo_ls_cisajs_001.dat
#        - zvo_ls_cisajscktalt_001.dat
#        - zvo_ls_cisajscktaltex_001.dat for TwoBodyGEx fixtures
#
# Run all systems:
#   julia --project=@. test/integration/lanczos_equivalent.jl
# Run a subset (by fixture dir name or C model name):
#   JULIA_MVMC_LANCZOS_MODELS=hubbard_chain_lanczos \
#     julia --project=@. test/integration/lanczos_equivalent.jl

using Test
using MVMCExpertModeParsers
using MVMCOptimizers

const LANCZOS_MODEL_FILTER_ENV = "JULIA_MVMC_LANCZOS_MODELS"
const LANCZOS_LS_OUT_TOL = 1e-8
const LANCZOS_QQQQ_TOL = 1e-10
const LANCZOS_GREEN_TOL = 1e-8

const LANCZOS_MODELS = [
    (
        fixture = "hubbard_chain_lanczos",
        c_model = "HubbardChainLanczos",
        mode = :real,
        n_para = 31,
        has_two_body_ex = false,
    ),
    (
        fixture = "spin_chain_lanczos",
        c_model = "SpinChainLanczos",
        mode = :real,
        n_para = 34,
        has_two_body_ex = false,
    ),
    (
        fixture = "hubbard_chain_real",
        c_model = "HubbardChain",
        mode = :real,
        n_para = 19,
        has_two_body_ex = true,
    ),
]

function selected_lanczos_models()
    raw = strip(get(ENV, LANCZOS_MODEL_FILTER_ENV, ""))
    isempty(raw) && return LANCZOS_MODELS
    requested = Set(strip.(split(raw, ",")))
    return filter(LANCZOS_MODELS) do model
        model.fixture in requested || model.c_model in requested
    end
end

parse_float_vector(path::AbstractString) = parse.(Float64, split(read(path, String)))

function assert_close_vector(actual_path, expected_path; nvalues = nothing, atol)
    actual = parse_float_vector(actual_path)
    expected = parse_float_vector(expected_path)

    if nvalues !== nothing
        @test length(actual) == nvalues
        @test length(expected) == nvalues
    end
    @test length(actual) == length(expected)
    @test maximum(abs.(actual .- expected)) <= atol
end

function run_lanczos_model(model)
    refdir = joinpath(@__DIR__, "reference", model.fixture, "physcal_ref")
    inputs = joinpath(refdir, "inputs")
    expected = joinpath(refdir, "expected")
    namelist = joinpath(inputs, "namelist.def")
    opt_para = joinpath(refdir, "zqp_opt.dat")

    @assert isfile(namelist) "Lanczos namelist.def missing for $(model.fixture)"
    @assert isfile(opt_para) "Lanczos zqp_opt.dat missing for $(model.fixture)"
    @assert isdir(expected) "Lanczos expected/ missing for $(model.fixture)"

    mktempdir() do dir
        result = run_phys_cal_from_namelist(
            namelist; opt_para = opt_para, mode = model.mode, output_dir = dir,
        )

        @test result.status == 0
        @test result.n_para_consumed == model.n_para

        actual_ls_out = joinpath(dir, "zvo_ls_out_001.dat")
        expected_ls_out = joinpath(expected, "zvo_ls_out_001.dat")
        actual_qqqq = joinpath(dir, "zvo_ls_qqqq_001.dat")
        expected_qqqq = joinpath(expected, "zvo_ls_qqqq_001.dat")

        @test isfile(actual_ls_out)
        @test isfile(actual_qqqq)
        assert_close_vector(
            actual_ls_out, expected_ls_out;
            nvalues = 3, atol = LANCZOS_LS_OUT_TOL,
        )
        assert_close_vector(
            actual_qqqq, expected_qqqq;
            nvalues = 16, atol = LANCZOS_QQQQ_TOL,
        )

        for file in ("zvo_ls_cisajs_001.dat", "zvo_ls_cisajscktalt_001.dat")
            actual_path = joinpath(dir, file)
            expected_path = joinpath(expected, file)
            @test isfile(actual_path)
            @test isfile(expected_path)
            assert_close_vector(actual_path, expected_path; atol = LANCZOS_GREEN_TOL)
        end

        if model.has_two_body_ex
            file = "zvo_ls_cisajscktaltex_001.dat"
            actual_path = joinpath(dir, file)
            expected_path = joinpath(expected, file)
            @test isfile(actual_path)
            @test isfile(expected_path)
            assert_close_vector(actual_path, expected_path; atol = LANCZOS_GREEN_TOL)
        end
    end
end

@testset "Julia-mVMC Full Lanczos PhysCal e2e (C reference)" begin
    models = selected_lanczos_models()
    @test !isempty(models)

    for model in models
        @testset "$(model.fixture) ($(model.c_model), $(model.mode))" begin
            run_lanczos_model(model)
        end
    end
end
