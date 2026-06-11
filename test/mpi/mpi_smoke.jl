# mpiexec 配下（または serial）で run_para_opt_from_namelist を 1 回実行する worker。
# 使い方: julia --project=<workspace> test/mpi/mpi_smoke.jl <output_dir>
using MVMCOptimizers
using MPI
using Logging

if get(ENV, "JULIA_MVMC_SMOKE_LOG_STDOUT", "0") == "1"
    global_logger(ConsoleLogger(stdout, Logging.Info))
end

const fixture = joinpath(@__DIR__, "..", "integration", "reference",
                         "heisenberg_chain_real", "inputs", "namelist.def")
const outdir = ARGS[1]

# top-level abort guard（design §6.1、plan review F6）: 片 rank だけが exception で
# 死ぬと他 rank が collective / barrier で hang するため、worker 全体を try/catch で
# 包み、MPI 初期化済みなら MPI.Abort で全 rank を fail-fast させる。
try
    result = MVMCOptimizers.run_para_opt_from_namelist(
        fixture; nsteps = 4, nsmp = 4, mode = :real, output_dir = outdir)

    # rank の役割確認（serial 実行では常に output rank）。
    if isempty(result.zvo_first_n)
        # 非 rank0: minimal result（readback なし）であること。
        @assert isnan(result.final_energy_per_site)
        println("worker: non-root rank ok")
    else
        @assert result.status == 0
        @assert length(result.zvo_first_n) == 4
        println("worker: root rank ok (E/site = $(result.final_energy_per_site))")
    end
catch err
    @error "mpi_smoke worker failed" exception = (err, catch_backtrace())
    if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD, 1)   # 他 rank を hang させない
    end
    rethrow()
end
