# test/integration/ctest_equivalent.jl
#
# Julia-mVMC ctest-equivalent integration tests.
#
# This runner mirrors C-mVMC's standard `test/python/runtest.py` criterion:
# compare two final optimisation summary values against `ref_mean.dat` and
# `ref_std.dat`, failing only when both
#
#   abs(diff) >= 3 * ref_std && abs(diff) >= 1e-8
#
# hold. Unlike `test/integration/runtests.jl`, this is not a first-10-step
# bit-level comparison.
#
# Run all supported standard models:
#   julia --project=@. test/integration/ctest_equivalent.jl
#
# Run a subset by fixture name or C model name:
#   JULIA_MVMC_CTEST_MODELS=HeisenbergChain,hubbard_chain_cmp \
#     julia --project=@. test/integration/ctest_equivalent.jl

using Test
using MVMCExpertModeParsers
using MVMCOptimizers

include(joinpath(@__DIR__, "ctest_models.jl"))

const CTEST_ABS_FLOOR = 1e-8
const CTEST_MODEL_FILTER_ENV = "JULIA_MVMC_CTEST_MODELS"

function parse_numeric_file(path::AbstractString)
    values = Float64[]
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        append!(values, parse.(Float64, split(stripped)))
    end
    return values
end

function selected_ctest_models()
    raw_filter = strip(get(ENV, CTEST_MODEL_FILTER_ENV, ""))
    isempty(raw_filter) && return CTEST_STANDARD_MODELS

    requested = Set(strip.(split(raw_filter, ",")))
    return filter(CTEST_STANDARD_MODELS) do model
        model.fixture in requested || model.c_model in requested
    end
end

function effective_counts(model, namelist::AbstractString)
    parsed = MVMCExpertModeParsers.parse_expert_mode_files(String(namelist))
    nsteps = model.nsteps_override === nothing ?
        parsed.modpara.nsr_opt_itr_step :
        Int(model.nsteps_override)
    nsmp = model.nsmp_override === nothing ?
        parsed.modpara.nsr_opt_itr_smp :
        Int(model.nsmp_override)
    return (nsteps = nsteps, nsmp = nsmp)
end

function ctest_passes(calculated::Float64, expected::Float64, sigma::Float64)
    diff = abs(calculated - expected)
    return !(diff >= 3 * sigma && diff >= CTEST_ABS_FLOOR)
end

function run_ctest_model(model)
    refdir = joinpath(@__DIR__, "reference", model.fixture)
    inputs = joinpath(refdir, "inputs")
    ctest_ref = joinpath(refdir, "ctest_ref")
    namelist = joinpath(inputs, "namelist.def")
    ref_mean_path = joinpath(ctest_ref, "ref_mean.dat")
    ref_std_path = joinpath(ctest_ref, "ref_std.dat")

    @assert isfile(namelist) "namelist.def missing for $(model.fixture)"
    @assert isfile(ref_mean_path) "ref_mean.dat missing for $(model.fixture)"
    @assert isfile(ref_std_path) "ref_std.dat missing for $(model.fixture)"

    counts = effective_counts(model, namelist)
    result = MVMCOptimizers.run_para_opt_from_namelist(
        namelist;
        nsteps = counts.nsteps,
        nsmp = counts.nsmp,
        mode = model.mode,
    )

    return (
        result = result,
        ref_mean = parse_numeric_file(ref_mean_path)[1:2],
        ref_std = parse_numeric_file(ref_std_path)[1:2],
        nsteps = counts.nsteps,
        nsmp = counts.nsmp,
    )
end

@testset "Julia-mVMC C ctest-equivalent integration" begin
    models = selected_ctest_models()
    @test !isempty(models)

    for model in models
        @testset "$(model.fixture) ($(model.c_model), $(model.mode))" begin
            run = run_ctest_model(model)

            @test run.result.status == 0
            @test run.result.effective_nsteps == run.nsteps
            @test run.result.effective_nsmp == run.nsmp

            for idx in 1:2
                calculated = run.result.ctest_values[idx]
                expected = run.ref_mean[idx]
                sigma = run.ref_std[idx]
                ok = ctest_passes(calculated, expected, sigma)
                if !ok
                    @info "ctest comparison failed" fixture=model.fixture column=idx calculated expected sigma diff=abs(calculated - expected)
                end
                @test ok
            end
        end
    end
end
