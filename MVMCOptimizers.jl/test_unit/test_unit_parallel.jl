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

using MVMCOptimizers: count_total_parameters, get_parameter_value,
                      set_parameter_value!, pack_parameters, unpack_parameters!
using MVMCExpertModeParsers

@testset "ModPara rnd_seed default matches C readdef.c:1967 (review F2)" begin
    # RndSeed 行が欠落した modpara.def は C では 11272 で走る。parser の kwdef
    # default が 12345 だと C parity が missing の行で崩れる（review 2026-06-11 F2）。
    @test MVMCExpertModeParsers.ModParaParameters().rnd_seed == 11272
end

@testset "PhysCal runner fails fast under MPI (review F7)" begin
    # R0 では PhysCal は MPI 未対応。mpiexec 検出下で走らせると全 rank が同一
    # output へ書くため、entry で error する（guard は parse より前なので
    # ダミーパスで検証できる）。
    withenv("JULIA_MVMC_MPI" => "1") do
        @test_throws ErrorException MVMCOptimizers.run_phys_cal_from_namelist(
            "nonexistent_namelist.def"; opt_para = "nonexistent.dat", mode = :real)
    end
end

@testset "parameter pack/unpack roundtrip (spec §5-2, F4)" begin
    fixture = joinpath(@__DIR__, "..", "..", "test", "integration", "reference",
                       "heisenberg_chain_real", "inputs", "namelist.def")
    data = MVMCExpertModeParsers.parse_expert_mode_files(fixture)
    n = count_total_parameters(data)
    @test n > 0

    para = pack_parameters(data)
    @test length(para) == n

    # update_parameter_value との整合: i 番目に δ を足すと pack[i] が δ 増える
    for i in (1, n)
        before = pack_parameters(data)[i]
        MVMCOptimizers.update_parameter_value(data, i, 0.25, -0.5)
        @test pack_parameters(data)[i] ≈ before + ComplexF64(0.25, -0.5)
    end

    # roundtrip: 摂動 → unpack で元に戻る
    original = pack_parameters(data)
    perturbed = original .+ ComplexF64(0.01, 0.02)
    unpack_parameters!(data, perturbed)
    @test pack_parameters(data) ≈ perturbed
    unpack_parameters!(data, original)
    @test pack_parameters(data) ≈ original

    @test_throws ArgumentError unpack_parameters!(data, original[1:max(n - 1, 0)])
end

@testset "duplicate idx invariant (plan review F7 / addendum C1)" begin
    fixture = joinpath(@__DIR__, "..", "..", "test", "integration", "reference",
                       "heisenberg_chain_real", "inputs", "namelist.def")
    data = MVMCExpertModeParsers.parse_expert_mode_files(fixture)

    # unpack 後、同一 idx の duplicate term は全て同値であること。
    para = pack_parameters(data) .+ ComplexF64(0.1, -0.1)
    unpack_parameters!(data, para)
    for idx in unique(t.idx for t in data.orbital_terms)
        vals = [t.value for t in data.orbital_terms if t.idx == idx]
        @test all(v -> v == vals[1], vals)
    end

    # orbital duplicate を人為的に不一致にすると fail-fast すること。
    dup_idx = findfirst(i -> count(t -> t.idx == data.orbital_terms[i].idx,
                                   data.orbital_terms) > 1,
                        eachindex(data.orbital_terms))
    if dup_idx !== nothing
        data.orbital_terms[dup_idx].value += ComplexF64(1.0, 0.0)
        @test_throws ErrorException unpack_parameters!(data, para)
    end

    # RBM duplicate も検査対象であること（addendum C1）。fixture に RBM がないため
    # synthetic に同一 idx・異値の term を 2 つ作る（fresh parse で orbital 側の
    # 不一致と混ざらないようにする）。
    data_rbm = MVMCExpertModeParsers.parse_expert_mode_files(fixture)
    push!(data_rbm.charge_rbm_phys_layer_terms,
          MVMCExpertModeParsers.ChargeRBMPhysLayerTerm(0, ComplexF64(0.1, 0.0), true, 0),
          MVMCExpertModeParsers.ChargeRBMPhysLayerTerm(1, ComplexF64(0.2, 0.0), true, 0))
    @test_throws ErrorException MVMCOptimizers.check_duplicate_consistency(data_rbm)
end

@testset "sync_modified_parameter!(ctx, data) serial == legacy" begin
    fixture = joinpath(@__DIR__, "..", "..", "test", "integration", "reference",
                       "heisenberg_chain_real", "inputs", "namelist.def")
    data_a = MVMCExpertModeParsers.parse_expert_mode_files(fixture)
    data_b = MVMCExpertModeParsers.parse_expert_mode_files(fixture)
    MVMCOptimizers.sync_modified_parameter!(data_a)                       # legacy
    MVMCOptimizers.sync_modified_parameter!(serial_context(), data_b)     # ctx 版
    @test pack_parameters(data_a) ≈ pack_parameters(data_b)
end
