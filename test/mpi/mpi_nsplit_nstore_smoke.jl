# mpiexec worker for the direct-SR NSplitSize/NStoreO smoke matrix.
# Usage:
#   julia --project=<workspace> test/mpi/mpi_nsplit_nstore_smoke.jl <fixture> <mode> <nsplit> <nstore> <output_dir>
#
# The worker intentionally canonicalizes the copied fixture to NQPFull=1
# (NSPGaussLeg=1, NMPTrans=1, identity qptrans). This is a self-consistency
# gate for comm1 sample splitting and NStore handling, not a C-reference physics
# fixture for the original momentum/projector configuration.
using MVMCOptimizers
using MPI
using Logging

if get(ENV, "JULIA_MVMC_SMOKE_LOG_STDOUT", "0") == "1"
    global_logger(ConsoleLogger(stdout, Logging.Info))
end

length(ARGS) == 5 ||
    error("usage: mpi_nsplit_nstore_smoke.jl <fixture> <real|cmp|fsz> <nsplit> <nstore> <output_dir>")

const fixture_name = ARGS[1]
const mode = Symbol(ARGS[2])
const nsplit = parse(Int, ARGS[3])
const nstore = parse(Int, ARGS[4])
const outdir = ARGS[5]

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
                println(io, rpad(key, 15), replacements[key])
            else
                println(io, line)
            end
        end
    end
end

function read_nsite(modpara_path::AbstractString)
    for line in eachline(modpara_path)
        parts = split(strip(line))
        if length(parts) >= 2 && parts[1] == "Nsite"
            return parse(Int, parts[2])
        end
    end
    error("Nsite not found in $modpara_path")
end

function write_identity_qptrans!(path::AbstractString, nsite::Int)
    open(path, "w") do io
        println(io, "=============================================")
        println(io, "NQPTrans          1")
        println(io, "=============================================")
        println(io, "======== TrIdx_TrWeight_and_TrIdx_i_xi ======")
        println(io, "=============================================")
        println(io, "0    1.00000")
        for site = 0:(nsite-1)
            println(io, lpad(0, 5), lpad(site, 7), lpad(site, 7), lpad(1, 6))
        end
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
            "NSPGaussLeg" => "1",
            "NMPTrans" => "1",
            "NSROptItrStep" => "1",
            "NSROptItrSmp" => "1",
            "NVMCWarmUp" => "1",
            "NVMCSample" => "4",
            "NSplitSize" => string(nsplit),
            "NStore" => string(nstore),
            "NSRCG" => "0",
        ),
    )
    write_identity_qptrans!(joinpath(inputs, "qptransidx.def"), read_nsite(modpara_path))
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

    label = "nsplit-nstore worker: $(fixture_name) nsplit=$(nsplit) nstore=$(nstore)"
    if isempty(result.zvo_first_n)
        @assert isnan(result.final_energy_per_site)
        println("$label non-root rank ok")
    else
        @assert result.status == 0
        @assert length(result.zvo_first_n) == 1
        println("$label root rank ok")
    end
catch err
    @error "mpi_nsplit_nstore_smoke worker failed" exception = (err, catch_backtrace())
    if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD, 1)
    end
    rethrow()
end
