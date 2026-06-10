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
