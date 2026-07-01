# mpiexec worker for standard-projection NQPFull > 1 with NSplitSize.
# Usage:
#   julia --project=<workspace> test/mpi/mpi_nsplit_standard_projection_smoke.jl <fixture> <mode> <nsplit> <nstore> <output_dir>
using MVMCOptimizers
using MPI
using Logging

if get(ENV, "JULIA_MVMC_SMOKE_LOG_STDOUT", "0") == "1"
    global_logger(ConsoleLogger(stdout, Logging.Info))
end

length(ARGS) == 5 ||
    error("usage: mpi_nsplit_standard_projection_smoke.jl <fixture> <real|cmp|fsz> <nsplit> <nstore> <output_dir>")

const fixture_name = ARGS[1]
const mode = Symbol(ARGS[2])
const nsplit = parse(Int, ARGS[3])
const nstore = parse(Int, ARGS[4])
const outdir = ARGS[5]
const smoke_nsp_gauss_leg = get(ENV, "JULIA_MVMC_SMOKE_NSPGAUSSLEG", "")

mode in (:real, :cmp, :fsz) || error("mode must be real, cmp, or fsz; got $mode")
nsplit >= 1 || error("nsplit must be >= 1; got $nsplit")
nstore in (0, 1) || error("nstore must be 0 or 1; got $nstore")

const src_inputs = joinpath(@__DIR__, "..", "integration", "reference", fixture_name, "inputs")

function replace_modpara_values!(path::AbstractString, replacements)
    lines = readlines(path)
    open(path, "w") do io
        for line in lines
            stripped = strip(line)
            if isempty(stripped) || startswith(stripped, "-")
                println(io, line)
                continue
            end
            key = first(split(stripped))
            if haskey(replacements, key)
                replacement = replacements[key]
                replacement == "" ? println(io, line) : println(io, rpad(key, 15), replacement)
            else
                println(io, line)
            end
        end
    end
end

function ensure_modpara_value!(path::AbstractString, key::AbstractString, value::AbstractString)
    isempty(value) && return
    text = read(path, String)
    occursin(Regex("(?m)^" * key * "\\s+"), text) && return
    open(path, "a") do io
        println(io, rpad(key, 15), value)
    end
end

function prepare_inputs()
    workdir = mktempdir()
    inputs = joinpath(workdir, "inputs")
    cp(src_inputs, inputs)

    modpara_path = joinpath(inputs, "modpara.def")
    replace_modpara_values!(
        modpara_path,
        Dict(
            "NSPGaussLeg" => smoke_nsp_gauss_leg,
            "NSROptItrStep" => "1",
            "NSROptItrSmp" => "1",
            "NVMCWarmUp" => "1",
            "NVMCSample" => "4",
            "NSplitSize" => string(nsplit),
            "NStore" => string(nstore),
            "NSRCG" => "0",
        ),
    )
    ensure_modpara_value!(modpara_path, "NSPGaussLeg", smoke_nsp_gauss_leg)
    return joinpath(inputs, "namelist.def")
end

try
    namelist = prepare_inputs()
    result = MVMCOptimizers.run_para_opt_from_namelist(
        namelist;
        nsteps = 1,
        nsmp = 1,
        mode = mode,
        output_dir = outdir,
    )

    label = "nsplit-standard-projection worker: $(fixture_name) nsplit=$(nsplit) nstore=$(nstore)"
    if isempty(result.zvo_first_n)
        @assert isnan(result.final_energy_per_site)
        println("$label non-root rank ok")
    else
        @assert result.status == 0
        @assert length(result.zvo_first_n) == 1
        println("$label root rank ok")
    end
catch err
    @error "mpi_nsplit_standard_projection_smoke worker failed" exception = (err, catch_backtrace())
    if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD, 1)
    end
    rethrow()
end
