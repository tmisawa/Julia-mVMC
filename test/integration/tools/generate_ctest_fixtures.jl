#!/usr/bin/env julia

# Generate Julia-mVMC ctest-equivalent fixtures from a built C-mVMC tree.
#
# Expected C-side preparation:
#   cd private-mVMC/mVMC/build/test/python
#   python3 runtest.py HeisenbergChain
#   python3 runtest.py HubbardChain
#   ...
#
# Then from Julia-mVMC:
#   julia --project=@. test/integration/tools/generate_ctest_fixtures.jl \
#     --c-test-dir ../private-mVMC/mVMC/build/test/python

include(joinpath(@__DIR__, "..", "ctest_models.jl"))

const DEFAULT_C_TEST_DIR = abspath(joinpath(
    @__DIR__,
    "..",
    "..",
    "..",
    "..",
    "private-mVMC",
    "mVMC",
    "build",
    "test",
    "python",
))

function usage()
    println("""
    usage: julia --project=@. test/integration/tools/generate_ctest_fixtures.jl [options]

    options:
      --c-test-dir PATH   C-mVMC build/test/python directory, or mVMC root.
                          Default: $DEFAULT_C_TEST_DIR
      --models LIST       Comma-separated fixture names or C model names.
                          Default: all CTEST_STANDARD_MODELS.
      --help              Show this message.
    """)
end

function parse_args(args)
    c_test_dir = DEFAULT_C_TEST_DIR
    model_filter = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help"
            usage()
            exit(0)
        elseif arg == "--c-test-dir"
            i += 1
            i <= length(args) || error("--c-test-dir requires a value")
            c_test_dir = args[i]
        elseif arg == "--models"
            i += 1
            i <= length(args) || error("--models requires a value")
            model_filter = args[i]
        else
            error("unknown argument: $arg")
        end
        i += 1
    end
    return (c_test_dir = c_test_dir, model_filter = model_filter)
end

function normalize_c_test_dir(path::AbstractString)
    expanded = abspath(expanduser(path))
    candidates = [
        expanded,
        joinpath(expanded, "build", "test", "python"),
        joinpath(expanded, "test", "python"),
    ]
    for candidate in candidates
        if isfile(joinpath(candidate, "runtest.py")) && isdir(joinpath(candidate, "data"))
            return candidate
        end
    end
    error("could not find C-mVMC test/python directory from: $path")
end

function selected_models(model_filter)
    model_filter === nothing && return CTEST_STANDARD_MODELS
    requested = Set(strip.(split(model_filter, ",")))
    models = filter(CTEST_STANDARD_MODELS) do model
        model.fixture in requested || model.c_model in requested
    end
    isempty(models) && error("no CTEST_STANDARD_MODELS matched --models=$model_filter")
    return models
end

function require_file(path::AbstractString)
    isfile(path) || error("required file missing: $path")
    return path
end

function require_dir(path::AbstractString)
    isdir(path) || error("required directory missing: $path")
    return path
end

function copy_file(src::AbstractString, dst::AbstractString)
    mkpath(dirname(dst))
    cp(src, dst; force = true)
    return dst
end

function copy_work_defs(work_dir::AbstractString, inputs_dir::AbstractString)
    def_files = sort(filter(path -> endswith(path, ".def"), readdir(work_dir; join = true)))
    isempty(def_files) && error("no generated .def files found in $work_dir")
    for src in def_files
        copy_file(src, joinpath(inputs_dir, basename(src)))
    end
    return length(def_files)
end

function copy_optional_initial_files(data_dir::AbstractString, inputs_dir::AbstractString)
    copied = String[]
    for name in ("initial.def", "zqp_opt.dat")
        src = joinpath(data_dir, name)
        if isfile(src)
            copy_file(src, joinpath(inputs_dir, name))
            push!(copied, name)
        end
    end
    return copied
end

function generate_fixture(model, c_test_dir::AbstractString)
    work_dir = require_dir(joinpath(c_test_dir, "work", model.c_model))
    data_dir = require_dir(joinpath(c_test_dir, "data", model.c_model))
    source_ref_dir = require_dir(joinpath(data_dir, "ref"))

    fixture_dir = joinpath(@__DIR__, "..", "reference", model.fixture)
    inputs_dir = joinpath(fixture_dir, "inputs")
    ctest_ref_dir = joinpath(fixture_dir, "ctest_ref")

    n_defs = copy_work_defs(work_dir, inputs_dir)
    copied_initial = copy_optional_initial_files(data_dir, inputs_dir)

    copy_file(require_file(joinpath(source_ref_dir, "ref_mean.dat")), joinpath(ctest_ref_dir, "ref_mean.dat"))
    copy_file(require_file(joinpath(source_ref_dir, "ref_std.dat")), joinpath(ctest_ref_dir, "ref_std.dat"))

    return (
        fixture = model.fixture,
        c_model = model.c_model,
        n_defs = n_defs,
        copied_initial = copied_initial,
    )
end

function main(args = ARGS)
    parsed = parse_args(args)
    c_test_dir = normalize_c_test_dir(parsed.c_test_dir)
    models = selected_models(parsed.model_filter)

    println("C test dir: $c_test_dir")
    println("Generating $(length(models)) ctest fixture(s)")

    for model in models
        result = generate_fixture(model, c_test_dir)
        optional = isempty(result.copied_initial) ? "none" : join(result.copied_initial, ", ")
        println("  $(result.fixture) <= $(result.c_model): $(result.n_defs) .def files, optional inputs: $optional")
    end
end

main()
