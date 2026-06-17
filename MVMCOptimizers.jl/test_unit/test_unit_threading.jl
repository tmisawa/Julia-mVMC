using Test
using MVMCOptimizers

const MO = MVMCOptimizers

@testset "unit/threading: JULIA_MVMC_INNER_THREADS override" begin
    old = get(ENV, "JULIA_MVMC_INNER_THREADS", nothing)
    try
        delete!(ENV, "JULIA_MVMC_INNER_THREADS")
        @test !MO.vmc_inner_threading_requested(true)
        @test !MO.vmc_inner_threading_enabled(1024, true)

        ENV["JULIA_MVMC_INNER_THREADS"] = "0"
        @test !MO.vmc_inner_threading_requested(true)
        @test !MO.vmc_inner_threading_enabled(1024, true)

        ENV["JULIA_MVMC_INNER_THREADS"] = "1"
        expected = Threads.nthreads() > 1
        @test MO.vmc_inner_threading_requested(true) == expected
        @test MO.vmc_inner_threading_enabled(1024, true) == expected

        ENV["JULIA_MVMC_INNER_THREADS"] = "true"
        @test_throws ErrorException MO.vmc_inner_threading_requested(true)
    finally
        if old === nothing
            delete!(ENV, "JULIA_MVMC_INNER_THREADS")
        else
            ENV["JULIA_MVMC_INNER_THREADS"] = old
        end
    end
end

@testset "unit/threading: JULIA_MVMC_PFAPACK_THREADS override" begin
    old_inner = get(ENV, "JULIA_MVMC_INNER_THREADS", nothing)
    old_pfapack = get(ENV, "JULIA_MVMC_PFAPACK_THREADS", nothing)
    try
        ENV["JULIA_MVMC_INNER_THREADS"] = "1"
        delete!(ENV, "JULIA_MVMC_PFAPACK_THREADS")
        @test !MO.vmc_pfapack_threading_requested(true)

        ENV["JULIA_MVMC_PFAPACK_THREADS"] = "0"
        @test !MO.vmc_pfapack_threading_requested(true)

        ENV["JULIA_MVMC_PFAPACK_THREADS"] = "1"
        @test MO.vmc_pfapack_threading_requested(true) == (Threads.nthreads() > 1)

        ENV["JULIA_MVMC_PFAPACK_THREADS"] = "auto"
        @test_throws ErrorException MO.vmc_pfapack_threading_requested(true)
    finally
        if old_inner === nothing
            delete!(ENV, "JULIA_MVMC_INNER_THREADS")
        else
            ENV["JULIA_MVMC_INNER_THREADS"] = old_inner
        end
        if old_pfapack === nothing
            delete!(ENV, "JULIA_MVMC_PFAPACK_THREADS")
        else
            ENV["JULIA_MVMC_PFAPACK_THREADS"] = old_pfapack
        end
    end
end

@testset "unit/threading: inner copy helpers" begin
    n = 128
    real_src = [sin(0.03 * i) for i = 1:n]
    complex_dst_seq = fill(99.0 + 99.0im, n)
    complex_dst_thr = copy(complex_dst_seq)

    MO.copy_real_to_complex!(complex_dst_seq, real_src, n; threaded = false)
    MO.copy_real_to_complex!(complex_dst_thr, real_src, n; threaded = true)
    @test complex_dst_thr == complex_dst_seq
    @test complex_dst_thr == ComplexF64.(real_src, 0.0)

    complex_src = [
        ComplexF64(cos(0.05 * i), sin(0.07 * i))
        for i = 1:n
    ]
    real_dst_seq = fill(99.0, n)
    real_dst_thr = copy(real_dst_seq)

    MO.copy_complex_realpart!(real_dst_seq, complex_src, n; threaded = false)
    MO.copy_complex_realpart!(real_dst_thr, complex_src, n; threaded = true)
    @test real_dst_thr == real_dst_seq
    @test real_dst_thr == real.(complex_src)

    short_dst = fill(0.0 + 0.0im, 4)
    MO.copy_real_to_complex!(short_dst, real_src, n; threaded = true)
    @test short_dst == ComplexF64.(real_src[1:4], 0.0)
end

@testset "unit/threading: energy accumulator reduction" begin
    a = MO.VMCEnergyAccumulator()
    b = MO.VMCEnergyAccumulator()
    MO.accumulate_energy!(a, 2.0, 1.0 + 2.0im; sz = 3.0)
    MO.accumulate_energy!(b, 0.5, -2.0 + 1.0im; sz = -1.0)

    combined = MO.VMCEnergyAccumulator()
    MO.merge_energy_accumulator!(combined, a)
    MO.merge_energy_accumulator!(combined, b)

    @test combined.wc == 2.5 + 0.0im
    @test combined.etot == 1.0 + 4.5im
    @test combined.etot2 == 12.5 + 0.0im
    @test combined.sztot == 5.5 + 0.0im
    @test combined.sztot2 == 18.5 + 0.0im

    energy = MO.EnergyData()
    MO.merge_energy_accumulators!(energy, (a, b))
    @test energy.wc == combined.wc
    @test energy.etot == combined.etot
    @test energy.etot2 == combined.etot2
    @test energy.sztot == combined.sztot
    @test energy.sztot2 == combined.sztot2

    MO.clear_energy_accumulator!(combined)
    @test combined.wc == 0.0 + 0.0im
    @test combined.etot == 0.0 + 0.0im
end

@testset "unit/threading: SR accumulator reduction" begin
    sr = MO.SROptData(2, 3, false)
    a = MO.VMCSROptAccumulator(sr)
    b = MO.VMCSROptAccumulator(sr)

    a.sr_opt_oo[1:3] .= [1.0 + 1.0im, 2.0 + 0.0im, -1.0 + 2.0im]
    b.sr_opt_oo[1:3] .= [0.5 + 0.0im, 3.0 - 1.0im, 1.0 + 0.0im]
    a.sr_opt_ho[1:2] .= [2.0 + 0.0im, 1.0 - 1.0im]
    b.sr_opt_ho[1:2] .= [-1.0 + 1.0im, 4.0 + 0.0im]
    a.sr_opt_o_store[1:2] .= [3.0 + 0.0im, 4.0 + 0.0im]
    b.sr_opt_o_store[1:2] .= [5.0 + 0.0im, 6.0 + 0.0im]

    a.sr_opt_oo_real[1:3] .= [1.0, 2.0, 3.0]
    b.sr_opt_oo_real[1:3] .= [4.0, 5.0, 6.0]
    a.sr_opt_ho_real[1:2] .= [7.0, 8.0]
    b.sr_opt_ho_real[1:2] .= [9.0, 10.0]
    a.sr_opt_o_store_real[1:2] .= [11.0, 12.0]
    b.sr_opt_o_store_real[1:2] .= [13.0, 14.0]

    MO.merge_sropt_accumulators!(sr, (a, b))
    @test sr.sr_opt_oo[1:3] == [1.5 + 1.0im, 5.0 - 1.0im, 0.0 + 2.0im]
    @test sr.sr_opt_ho[1:2] == [1.0 + 1.0im, 5.0 - 1.0im]
    @test sr.sr_opt_o_store[1:2] == [8.0 + 0.0im, 10.0 + 0.0im]
    @test sr.sr_opt_oo_real[1:3] == [5.0, 7.0, 9.0]
    @test sr.sr_opt_ho_real[1:2] == [16.0, 18.0]
    @test sr.sr_opt_o_store_real[1:2] == [24.0, 26.0]

    MO.clear_sropt_accumulator!(a)
    @test all(iszero, a.sr_opt_oo)
    @test all(iszero, a.sr_opt_o)
    @test all(iszero, a.sr_opt_oo_real)
    @test all(iszero, a.sr_opt_o_real)

    sr.sr_opt_o_store[1:2] .= [1.0 + 2.0im, 3.0 + 4.0im]
    sr.sr_opt_o_store_real[1:2] .= [5.0, 6.0]
    MO.clear_sropt_store!(sr)
    @test all(iszero, sr.sr_opt_o_store)
    @test all(iszero, sr.sr_opt_o_store_real)
end

@testset "unit/threading: PhysCal accumulator reduction" begin
    phys = MO.PhysicalQuantities(2, 1, 2)
    phys.local_cis_ajs .= [99.0 + 0.0im, 100.0 + 0.0im]

    a = MO.VMCPhysAccumulator(phys)
    b = MO.VMCPhysAccumulator(phys)
    a.local_cis_ajs .= [101.0 + 0.0im, 102.0 + 0.0im]
    b.local_cis_ajs .= [201.0 + 0.0im, 202.0 + 0.0im]
    a.phys_cis_ajs .= [1.0 + 1.0im, 2.0 + 0.0im]
    b.phys_cis_ajs .= [3.0 + 0.0im, -1.0 + 2.0im]
    a.phys_cis_ajs_ckt_alt .= [4.0 + 0.0im]
    b.phys_cis_ajs_ckt_alt .= [5.0 + 1.0im]
    a.local_cis_ajs_ckt_alt_dc .= [301.0 + 0.0im, 302.0 + 0.0im]
    b.local_cis_ajs_ckt_alt_dc .= [401.0 + 0.0im, 402.0 + 0.0im]
    a.phys_cis_ajs_ckt_alt_dc .= [6.0 + 0.0im, 7.0 + 0.0im]
    b.phys_cis_ajs_ckt_alt_dc .= [-1.0 + 0.0im, 2.0 + 3.0im]

    MO.merge_phys_accumulators!(phys, (a, b))
    @test phys.phys_cis_ajs == [4.0 + 1.0im, 1.0 + 2.0im]
    @test phys.phys_cis_ajs_ckt_alt == [9.0 + 1.0im]
    @test phys.phys_cis_ajs_ckt_alt_dc == [5.0 + 0.0im, 9.0 + 3.0im]
    @test phys.local_cis_ajs == [99.0 + 0.0im, 100.0 + 0.0im]
    @test phys.local_cis_ajs_ckt_alt_dc == zeros(ComplexF64, 2)

    none_acc = MO.VMCPhysAccumulator(nothing)
    @test isempty(none_acc.local_cis_ajs)
    @test isempty(none_acc.phys_cis_ajs)
    @test MO.merge_phys_accumulators!(nothing, (none_acc,)) === nothing

    phys.cis_ajs_ckt_alt_idx = [(1, 2)]
    phys.local_cis_ajs .= [2.0 + 1.0im, 3.0 - 4.0im]
    factored = MO.VMCPhysAccumulator(phys)
    factored.local_cis_ajs .= [4.0 + 2.0im, -1.0 + 3.0im]
    MO.accumulate_factored_green!(factored, phys, 0.5)
    @test factored.phys_cis_ajs_ckt_alt[1] ≈
          0.5 * (factored.local_cis_ajs[1] * conj(factored.local_cis_ajs[2]))
end

@testset "unit/threading: counter and timer reductions" begin
    a = MO.VMCCounterAccumulator(3)
    b = MO.VMCCounterAccumulator(3)
    MO.record_counter!(a, 1)
    MO.record_counter!(a, 3, 4)
    MO.record_counter!(b, 2, 7)
    @test_throws ErrorException MO.record_counter!(a, 0)

    counters = zeros(Int, 3)
    MO.merge_counter_accumulators!(counters, (a, b))
    @test counters == [1, 7, 4]

    short_counters = zeros(Int, 1)
    MO.merge_counter_accumulator!(short_counters, b)
    @test short_counters == [0, 7, 0]

    parent = MO.CTimer(true)
    local1 = MO.CTimer(true)
    local2 = MO.CTimer(true)
    MO.ctimer_add_elapsed!(local1, 4, UInt64(10))
    MO.ctimer_add_elapsed!(local2, 4, UInt64(15))
    MO.ctimer_add_elapsed!(local2, 7, UInt64(20))
    MO.ctimer_merge_all!(parent, (local1, local2))
    @test parent.elapsed_ns[5] == UInt64(25)
    @test parent.elapsed_ns[8] == UInt64(20)

    disabled = MO.CTimer(false)
    MO.ctimer_add_elapsed!(disabled, 4, UInt64(999))
    MO.ctimer_merge!(disabled, parent)
    @test all(iszero, disabled.elapsed_ns)
    @test_throws ErrorException MO.ctimer_add_elapsed!(parent, MO.CTIMER_N, UInt64(1))
end

@testset "unit/threading: local accumulator construction and merge" begin
    state = MO.VMCOptimizationState(2, 1, 1, 2, 1, 2, true, false)
    state.phys_quantities = MO.PhysicalQuantities(1, 1, 1)
    parent_timer = MO.CTimer(true)
    local_acc = MO.VMCThreadAccumulator(state, parent_timer)

    @test MO.ctimer_enabled(local_acc.timer)
    @test length(local_acc.sr_opt.sr_opt_oo) == length(state.sr_opt.sr_opt_oo)
    @test length(local_acc.sr_opt.sr_opt_o) == length(state.sr_opt.sr_opt_o)
    @test length(local_acc.phys.local_cis_ajs) == 1
    @test length(local_acc.phys.phys_cis_ajs) == 1
    @test length(local_acc.main_cal_scratch.proj_cnt_new) ==
          length(state.electron_config.tmp_ele_proj_cnt)
    @test length(local_acc.main_cal_scratch.pf_m_new_real) ==
          length(state.slater_matrix.pf_m_real)

    MO.accumulate_energy!(local_acc.energy, 1.0, 3.0 + 0.0im)
    MO.record_counter!(local_acc.counter, 1, 5)
    local_acc.phys.phys_cis_ajs[1] = 6.0 + 0.0im
    MO.ctimer_add_elapsed!(local_acc.timer, 4, UInt64(13))

    MO.merge_thread_accumulator!(state, parent_timer, local_acc)
    @test state.energy.wc == 1.0 + 0.0im
    @test state.energy.etot == 3.0 + 0.0im
    @test state.electron_config.counter[1] == 5
    @test state.phys_quantities.phys_cis_ajs[1] == 6.0 + 0.0im
    @test parent_timer.elapsed_ns[5] == UInt64(13)

    cached = MO.main_cal_accumulator!(
        state,
        parent_timer;
        all_complex = true,
        use_sr_store = false,
        nsrcg = false,
        use_sr_opt = true,
    )
    @test cached === state.workspace.main_cal_accumulator
    MO.accumulate_energy!(cached.energy, 2.0, 7.0 + 0.0im)
    MO.record_counter!(cached.counter, 1, 9)
    cached.sr_opt.sr_opt_oo[1] = 11.0 + 0.0im
    cached.sr_opt.sr_opt_ho[1] = 12.0 + 0.0im
    MO.ctimer_add_elapsed!(cached.timer, 4, UInt64(99))

    reused = MO.main_cal_accumulator!(
        state,
        parent_timer;
        all_complex = true,
        use_sr_store = false,
        nsrcg = false,
        use_sr_opt = true,
    )
    @test reused === cached
    @test reused.energy.wc == 0.0 + 0.0im
    @test all(iszero, reused.counter.counter)
    @test iszero(reused.sr_opt.sr_opt_oo[1])
    @test iszero(reused.sr_opt.sr_opt_ho[1])
    @test reused.timer.elapsed_ns[5] == UInt64(0)

    real_state = MO.VMCOptimizationState(2, 1, 1, 2, 1, 2, false, false)
    real_acc = MO.main_cal_accumulator!(
        real_state,
        parent_timer;
        all_complex = true,
        use_sr_store = false,
        nsrcg = false,
        use_sr_opt = true,
    )
    real_acc.sr_opt.sr_opt_oo[1] = 21.0 + 0.0im
    real_acc = MO.main_cal_accumulator!(
        real_state,
        parent_timer;
        all_complex = false,
        use_sr_store = true,
        nsrcg = false,
        use_sr_opt = true,
    )
    @test iszero(real_acc.sr_opt.sr_opt_oo[1])

    sr_opt_size = real_state.sr_opt.sr_opt_size
    active_real = sr_opt_size * sr_opt_size
    real_acc.sr_opt.sr_opt_oo_real[end] = 31.0
    real_acc.sr_opt.sr_opt_ho_real[1] = 32.0
    real_acc.sr_opt.sr_opt_o_store_real[1] = 33.0
    real_acc = MO.main_cal_accumulator!(
        real_state,
        parent_timer;
        all_complex = false,
        use_sr_store = true,
        nsrcg = false,
        use_sr_opt = true,
    )
    @test iszero(real_acc.sr_opt.sr_opt_ho_real[1])
    @test iszero(real_acc.sr_opt.sr_opt_o_store_real[1])
    @test all(iszero, @view(real_acc.sr_opt.sr_opt_oo_real[(active_real+1):end]))
end
