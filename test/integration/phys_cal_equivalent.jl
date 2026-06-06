# test/integration/phys_cal_equivalent.jl
#
# Julia-mVMC PhysCal end-to-end integration gate (Plan 3b).
#
# For each committed system, this runner:
#   1. loads the fixed variational parameters from the committed C `zqp_opt.dat`
#      (read_opt_para_file!, asserting the consumed count),
#   2. runs `run_phys_cal_from_namelist` on the committed PhysCal input set
#      (NVMCCalMode=1 + greentwoex.def) into a fresh tempdir, and
#   3. compares the three produced Green files against the committed C-mVMC
#      references with the per-quantity tolerances of tools/green_compare.jl:
#        - one-body  zvo_cisajs         rtol 1e-10 (linear accumulator)
#        - direct    zvo_cisajscktalt   rtol 1e-9  (4-operator product)
#        - factored  zvo_cisajscktaltex rtol 1e-9  (Sum w*g*conj(g))
#
# Unlike runtests.jl (first-10 SR steps) and ctest_equivalent.jl (final opt
# summary), this gate exercises the NVMCCalMode=1 PhysCal path and the factored
# two-body Green (TwoBodyGEx). References are read-only; the run writes only into
# the tempdir.
#
# Run all systems:
#   julia --project=@. test/integration/phys_cal_equivalent.jl
# Run a subset (by fixture dir name or C model name):
#   JULIA_MVMC_PHYS_CAL_MODELS=hubbard_chain_real \
#     julia --project=@. test/integration/phys_cal_equivalent.jl

using Test
using MVMCExpertModeParsers
using MVMCOptimizers

include(joinpath(@__DIR__, "tools", "green_compare.jl"))

const PHYS_CAL_MODEL_FILTER_ENV = "JULIA_MVMC_PHYS_CAL_MODELS"

# Systems carrying a committed reference/<fixture>/physcal_ref/ fixture.
# n_para is the expected read_opt_para_file! consumption
# (= NProj + NOrbitalIdx for these NRBM=0 fixtures); asserting it catches a
# silently-truncated zqp_opt.dat. It equals (token_count - 6) / 3.
const PHYS_CAL_MODELS = [
    (fixture = "heisenberg_chain_real", c_model = "HeisenbergChain",     mode = :real, n_para = 14),
    (fixture = "heisenberg_chain_cmp",  c_model = "HeisenbergChain_cmp", mode = :cmp,  n_para = 14),
    (fixture = "hubbard_chain_real",    c_model = "HubbardChain",        mode = :real, n_para = 19),
    (fixture = "hubbard_chain_dh_real", c_model = "HubbardChain_DH",     mode = :real, n_para = 35),
    (fixture = "kondo_chain_real",      c_model = "KondoChain",          mode = :real, n_para = 76),
]

function selected_phys_cal_models()
    raw = strip(get(ENV, PHYS_CAL_MODEL_FILTER_ENV, ""))
    isempty(raw) && return PHYS_CAL_MODELS
    requested = Set(strip.(split(raw, ",")))
    return filter(PHYS_CAL_MODELS) do model
        model.fixture in requested || model.c_model in requested
    end
end

# The single `<prefix>_NNN.dat` produced for a one-sample run (NDataQtySmp=1).
# Returns "" when absent so the green_compare helpers report a clean "missing".
# More than one match is a hard error: a stray extra `zvo_*_NNN.dat` must not let
# the gate silently compare only one sample and pass. The trailing underscore
# disambiguates the three nested prefixes
# (zvo_cisajs / zvo_cisajscktalt / zvo_cisajscktaltex).
function find_green_file(dir::AbstractString, prefix::AbstractString)
    isdir(dir) || return ""
    matches = filter(readdir(dir)) do name
        startswith(name, prefix * "_") && endswith(name, ".dat")
    end
    isempty(matches) && return ""
    length(matches) == 1 ||
        error("find_green_file: expected exactly one $(prefix)_*.dat in $dir, " *
              "found $(length(matches)): $(sort(matches))")
    return joinpath(dir, only(matches))
end

function run_phys_cal_model(model)
    refdir = joinpath(@__DIR__, "reference", model.fixture, "physcal_ref")
    inputs = joinpath(refdir, "inputs")
    expected = joinpath(refdir, "expected")
    namelist = joinpath(inputs, "namelist.def")
    opt_para = joinpath(refdir, "zqp_opt.dat")

    @assert isfile(namelist) "physcal namelist.def missing for $(model.fixture)"
    @assert isfile(opt_para) "physcal zqp_opt.dat missing for $(model.fixture)"
    @assert isdir(expected) "physcal expected/ missing for $(model.fixture)"

    mktempdir() do dir
        result = run_phys_cal_from_namelist(
            namelist; opt_para = opt_para, mode = model.mode, output_dir = dir,
        )

        @test result.status == 0
        @test result.n_para_consumed == model.n_para

        comparisons = [
            ("one-body zvo_cisajs", compare_green_one(
                find_green_file(dir, "zvo_cisajs"),
                find_green_file(expected, "zvo_cisajs"))),
            ("direct zvo_cisajscktalt", compare_green_two_dc(
                find_green_file(dir, "zvo_cisajscktalt"),
                find_green_file(expected, "zvo_cisajscktalt"))),
            ("factored zvo_cisajscktaltex", compare_green_factored(
                find_green_file(dir, "zvo_cisajscktaltex"),
                find_green_file(expected, "zvo_cisajscktaltex"))),
        ]
        for (label, r) in comparisons
            if r.ok
                @info "PhysCal Green match" fixture = model.fixture quantity = label exact = r.exact fallback = r.fallback_used max_abs = r.max_abs_err n = r.n_values
            else
                @info "PhysCal Green mismatch" fixture = model.fixture quantity = label detail = r.detail
            end
            @test r.ok
        end
    end
end

@testset "Julia-mVMC PhysCal e2e (C reference)" begin
    models = selected_phys_cal_models()
    @test !isempty(models)

    for model in models
        @testset "$(model.fixture) ($(model.c_model), $(model.mode))" begin
            run_phys_cal_model(model)
        end
    end
end
