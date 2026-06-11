# mpiexec 配下で WeightAverageWE(comm0) の数値 contract を直接検査する worker。
# 使い方: julia --project=<workspace> test/mpi/mpi_weight_average_smoke.jl
using MVMCOptimizers
using MPI
using Logging

if get(ENV, "JULIA_MVMC_SMOKE_LOG_STDOUT", "0") == "1"
    global_logger(ConsoleLogger(stdout, Logging.Info))
end

try
    ctx = MVMCOptimizers.build_parallel_context(1)
    state = MVMCOptimizers.VMCOptimizationState(1, 1, 0, 1, 1, 1, true, false)

    # Rank-local sums. After allreduce, Wc=sum(1:size), and all average fields
    # normalize to the known constants below.
    local_weight = ComplexF64(ctx.rank0 + 1, 0.0)
    state.energy.wc = local_weight
    state.energy.etot = 10.0 * local_weight
    state.energy.etot2 = 100.0 * local_weight
    state.energy.sztot = 2.0 * local_weight
    state.energy.sztot2 = 3.0 * local_weight

    MVMCOptimizers.weight_average_we!(ctx, state)

    expected_wc = ComplexF64(ctx.size0 * (ctx.size0 + 1) / 2, 0.0)
    @assert state.energy.wc == expected_wc
    @assert state.energy.etot ≈ 10.0 + 0.0im
    @assert state.energy.etot2 ≈ 100.0 + 0.0im
    @assert state.energy.sztot ≈ 2.0 + 0.0im
    @assert state.energy.sztot2 ≈ 3.0 + 0.0im

    if get(ENV, "JULIA_MVMC_SMOKE_TINY_WC", "0") == "1"
        state.energy.wc = 0.0 + 0.0im
        state.energy.etot = 0.0 + 0.0im
        MVMCOptimizers.weight_average_we!(ctx, state)
    end

    if ctx.rank0 == 0
        println("weight-average worker: root rank ok")
    else
        println("weight-average worker: non-root rank ok")
    end
catch err
    @error "mpi_weight_average_smoke worker failed" exception = (err, catch_backtrace())
    if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD, 1)
    end
    rethrow()
end
