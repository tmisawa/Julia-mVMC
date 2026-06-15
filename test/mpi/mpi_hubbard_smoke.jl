# mpiexec 配下で Hubbard chain の itinerant para-opt path を 1 回実行する worker。
# 使い方: julia --project=<workspace> test/mpi/mpi_hubbard_smoke.jl <output_dir>
using MVMCOptimizers
using MPI
using Logging

if get(ENV, "JULIA_MVMC_SMOKE_LOG_STDOUT", "0") == "1"
    global_logger(ConsoleLogger(stdout, Logging.Info))
end

const fixture = joinpath(@__DIR__, "..", "integration", "reference",
                         "hubbard_chain_real", "inputs", "namelist.def")
const outdir = ARGS[1]

try
    result = MVMCOptimizers.run_para_opt_from_namelist(
        fixture; nsteps = 4, nsmp = 4, mode = :real, output_dir = outdir)

    if isempty(result.zvo_first_n)
        @assert isnan(result.final_energy_per_site)
        println("hubbard worker: non-root rank ok")
    else
        @assert result.status == 0
        @assert length(result.zvo_first_n) == 4
        println("hubbard worker: root rank ok (E/site = $(result.final_energy_per_site))")
    end
catch err
    @error "mpi_hubbard_smoke worker failed" exception = (err, catch_backtrace())
    if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD, 1)
    end
    rethrow()
end
