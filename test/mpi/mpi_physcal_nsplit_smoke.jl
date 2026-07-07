# mpiexec worker for VMCPhysCal NSplitSize self-consistency smokes.
# Usage:
#   julia --project=<workspace> test/mpi/mpi_physcal_nsplit_smoke.jl <fixture> <mode> <nsplit> <output_dir>
using MVMCOptimizers
using MPI
using Logging

if get(ENV, "JULIA_MVMC_SMOKE_LOG_STDOUT", "0") == "1"
    global_logger(ConsoleLogger(stdout, Logging.Info))
end

length(ARGS) == 4 ||
    error("usage: mpi_physcal_nsplit_smoke.jl <fixture> <real|cmp> <nsplit> <output_dir>")

const fixture_name = ARGS[1]
const mode = Symbol(ARGS[2])
const nsplit = parse(Int, ARGS[3])
const outdir = ARGS[4]

mode in (:real, :cmp) || error("mode must be real or cmp; got $mode")
nsplit >= 1 || error("nsplit must be >= 1; got $nsplit")

const refdir = joinpath(
    @__DIR__,
    "..",
    "integration",
    "reference",
    fixture_name,
    "physcal_ref",
)
const src_inputs = joinpath(refdir, "inputs")
const opt_para = joinpath(refdir, "zqp_opt.dat")
const expected_files = (
    "zvo_out.dat",
    "zvo_var.dat",
    "zvo_cisajs_001.dat",
    "zvo_cisajscktalt_001.dat",
    "zvo_cisajscktaltex_001.dat",
)

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

function prepare_inputs()
    workdir = mktempdir()
    inputs = joinpath(workdir, "inputs")
    cp(src_inputs, inputs)
    replace_modpara_values!(
        joinpath(inputs, "modpara.def"),
        Dict(
            "NDataQtySmp" => "1",
            "NVMCWarmUp" => "1",
            "NVMCSample" => "4",
            "NSplitSize" => string(nsplit),
            "NLanczosMode" => "0",
        ),
    )
    return joinpath(inputs, "namelist.def")
end

try
    namelist = prepare_inputs()
    result = MVMCOptimizers.run_phys_cal_from_namelist(
        namelist;
        opt_para = opt_para,
        mode = mode,
        output_dir = outdir,
    )

    rank = MPI.Initialized() && !MPI.Finalized() ? MPI.Comm_rank(MPI.COMM_WORLD) : 0
    label = "physcal-nsplit worker: $(fixture_name) nsplit=$(nsplit)"
    @assert result.status == 0

    if rank == 0
        for name in expected_files
            @assert isfile(joinpath(outdir, name)) "missing PhysCal output: $name"
        end
        println("$label root rank ok")
    else
        println("$label non-root rank ok")
    end
catch err
    @error "mpi_physcal_nsplit_smoke worker failed" exception = (err, catch_backtrace())
    if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD, 1)
    end
    rethrow()
end
