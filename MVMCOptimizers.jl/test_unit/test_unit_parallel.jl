# Unit tests for src/parallel.jl (serial paths only; no mpiexec needed).
using Test
using MVMCOptimizers
using MVMCOptimizers: split_loop, split_range

@testset "split_loop matches C SplitLoop (mVMC/src/mVMC/splitloop.c:33-63)" begin
    # (loop_length, size) => [(ist, ien) for rank 0..size-1]   (0-based, half-open)
    cases = Dict(
        (10, 2) => [(0, 5), (5, 10)],                       # divisible
        (10, 3) => [(0, 3), (3, 6), (6, 10)],               # non-divisible (imod=1)
        (3, 2)  => [(0, 1), (1, 3)],                        # F8: size0=3, NSplitSize=2 相当
        (5, 4)  => [(0, 1), (1, 2), (2, 3), (3, 5)],        # non-divisible (imod=1)
        (4, 4)  => [(0, 1), (1, 2), (2, 3), (3, 4)],        # size == loop
        (2, 4)  => [(0, 1), (1, 2), (2, 2), (2, 2)],        # A8: 空 range (size > loop)
        (1, 4)  => [(0, 1), (1, 1), (1, 1), (1, 1)],        # A8: loop=1, size=4
        (0, 2)  => [(0, 0), (0, 0)],                        # A8: loop=0 → 全 rank 空 (C と同じ)
    )
    for ((n, sz), expected) in cases
        for r in 0:(sz - 1)
            @test split_loop(n, r, sz) == expected[r + 1]
        end
    end
end

@testset "split_range is the 1-based Julia view of split_loop" begin
    @test split_range(10, 0, 3) == 1:3
    @test split_range(10, 2, 3) == 7:10
    @test isempty(split_range(2, 3, 4))     # 空 range → empty Julia range
end

@testset "group1 assignment (C vmcmain.c:239)" begin
    # size0=3, NSplitSize=2: 最後の group が小さい (F8, warning のみで継続)
    @test [div(r, 2) for r in 0:2] == [0, 0, 1]
end

using MVMCOptimizers: ParallelContext, serial_context, is_output_rank,
                      mpi_env_detected, resolve_mpi_mode

@testset "serial_context" begin
    ctx = serial_context()
    @test ctx.is_mpi == false
    @test ctx.comm0 === nothing && ctx.comm1 === nothing && ctx.comm2 === nothing
    @test (ctx.rank0, ctx.size0) == (0, 1)
    @test (ctx.rank1, ctx.size1) == (0, 1)
    @test (ctx.rank2, ctx.size2) == (0, 1)
    @test ctx.group1 == 0
    @test is_output_rank(ctx)
end

@testset "JULIA_MVMC_MPI policy (spec §4.1, F12+A7)" begin
    clean = ("JULIA_MVMC_MPI" => nothing, "OMPI_COMM_WORLD_SIZE" => nothing,
             "PMI_SIZE" => nothing, "PMI_RANK" => nothing)
    withenv(clean...) do
        @test !mpi_env_detected()
        @test resolve_mpi_mode() === :serial                       # auto + 未検出
    end
    withenv(clean..., "OMPI_COMM_WORLD_SIZE" => "2") do
        @test mpi_env_detected()
        @test resolve_mpi_mode() === :mpi                          # auto + 検出
    end
    withenv(clean..., "PMI_SIZE" => "2") do
        @test resolve_mpi_mode() === :mpi                          # MPICH hydra
    end
    withenv(clean..., "JULIA_MVMC_MPI" => "1") do
        @test resolve_mpi_mode() === :mpi                          # =1 は常に MPI 必須
    end
    withenv(clean..., "JULIA_MVMC_MPI" => "0") do
        @test resolve_mpi_mode() === :serial                       # =0 + 未検出 → serial
    end
    withenv(clean..., "JULIA_MVMC_MPI" => "0", "PMI_RANK" => "0") do
        @test resolve_mpi_mode() === :mpi_guarded_serial           # =0 + 検出 → guarded
    end
    withenv(clean..., "JULIA_MVMC_MPI" => "yes") do
        @test_throws ErrorException resolve_mpi_mode()             # 不正値は明示 error
    end
end

using MVMCOptimizers: bcast!, bcast_scalar, allreduce_sum!, reduce_sum_to_root!,
                      barrier, reduce_counter!, _chunk_ranges

@testset "serial wrappers are no-ops" begin
    ctx = serial_context()
    v = ComplexF64[1.0 + 2.0im, 3.0]
    @test bcast!(ctx, v) === v && v == ComplexF64[1.0 + 2.0im, 3.0]
    @test bcast_scalar(ctx, 42) == 42
    @test allreduce_sum!(ctx, v) === v && v == ComplexF64[1.0 + 2.0im, 3.0]
    @test reduce_sum_to_root!(ctx, v) === v
    @test barrier(ctx) === nothing
    c = collect(1:11)
    @test reduce_counter!(ctx, c) === c && c == collect(1:11)   # F10: 不変
end

@testset "_chunk_ranges (C SafeMpiAllReduce, safempi.c:29 D_MpiSendMax)" begin
    @test _chunk_ranges(0, 4) == UnitRange{Int}[]
    @test _chunk_ranges(3, 4) == [1:3]
    @test _chunk_ranges(8, 4) == [1:4, 5:8]
    @test _chunk_ranges(9, 4) == [1:4, 5:8, 9:9]
end

using MVMCOptimizers: resolve_rnd_seed

@testset "resolve_rnd_seed C parity (spec §5-1, A5)" begin
    ctx = serial_context()
    @test resolve_rnd_seed(ctx, 11272, nothing) == 11272   # missing → parser default
    @test resolve_rnd_seed(ctx, 0, nothing) == 0           # RndSeed==0 → seed 0 (C parity)
    @test resolve_rnd_seed(ctx, 123, nothing) == 123       # 正値
    t = resolve_rnd_seed(ctx, -1, nothing)                 # 負値 → 時刻 seed
    @test t isa Int && t > 0
    @test resolve_rnd_seed(ctx, 123, 777) == 777           # 明示 seed kwarg が優先
    # group1 offset (C vmcmain.c:257)
    ctx2 = MVMCOptimizers.ParallelContext(false, nothing, nothing, nothing,
                                          3, 4, 0, 1, 0, 1, 3)   # group1=3 の擬似 ctx
    @test resolve_rnd_seed(ctx2, 100, nothing) == 103
end
