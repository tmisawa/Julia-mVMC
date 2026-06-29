using Test
using MPI
using MVMCOptimizers

try
    ctx = MVMCOptimizers.build_parallel_context(1)
    ctx.size0 == 2 || error("mpi_srcg_operate_smoke.jl expects mpiexec -n 2")

    ws = MVMCOptimizers.CGWorkspace(2, 2, false)
    if ctx.rank0 == 0
        ws.stcOs_real .= [1.0 2.0; 3.0 4.0]
        x = [0.5, -1.0]
    else
        ws.stcOs_real .= [5.0 6.0; 7.0 8.0]
        x = [999.0, 999.0]
    end
    ws.stcO .= [0.25, -0.5]
    ws.sdiag .= [2.0, 3.0]

    z = zeros(2)
    MVMCOptimizers.operate_by_s!(
        z,
        x,
        ws,
        2,
        2,
        0.25,
        0.1,
        false,
        ctx,
    )

    @test x ≈ [0.5, -1.0] atol = 0.0
    @test z ≈ [-15.30625, -22.7375] atol = 1e-12

    ws_complex = MVMCOptimizers.CGWorkspace(2, 2, true)
    if ctx.rank0 == 0
        ws_complex.stcOs_real .= [1.0 2.0; 3.0 4.0]
        ws_complex.stcOs_imag .= [0.5 -1.0; 1.5 0.25]
        x_complex = [0.5, -1.0]
    else
        ws_complex.stcOs_real .= [5.0 6.0; 7.0 8.0]
        ws_complex.stcOs_imag .= [-0.5 1.25; 0.75 -1.5]
        x_complex = [999.0, 999.0]
    end
    ws_complex.stcO .= [0.25, -0.5]
    ws_complex.sdiag .= [2.0, 3.0]

    z_complex = zeros(2)
    MVMCOptimizers.operate_by_s!(
        z_complex,
        x_complex,
        ws_complex,
        2,
        2,
        0.25,
        0.1,
        true,
        ctx,
    )

    @test x_complex ≈ [0.5, -1.0] atol = 0.0
    @test z_complex ≈ [-14.4859375, -24.2375] atol = 1e-12

    if MVMCOptimizers.is_output_rank(ctx)
        println("srcg-operate worker: root rank ok")
    else
        println("srcg-operate worker: non-root rank ok")
    end
catch err
    @error "mpi_srcg_operate_smoke worker failed" exception = (err, catch_backtrace())
    if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD, 1)
    end
    rethrow()
end
