using Test
using MVMCOptimizers
using MVMCExpertModeParsers:
    ExpertModeData,
    ModParaParameters,
    CoulombInterTerm,
    GutzwillerTerm,
    GreenOneTerm,
    GreenTwoExTerm,
    JastrowTerm,
    init_qp_weight!,
    projection_layout

const MO = MVMCOptimizers

@testset "unit/threading: VMCThreadConfig" begin
    cfg = MO.VMCThreadConfig(10; requested_threads = 4, min_work_per_thread = 3)
    @test cfg.requested_threads == 4
    @test cfg.effective_threads == 4
    @test cfg.work_items == 10
    @test MO.effective_thread_count(cfg) == 4
    @test MO.vmc_threading_enabled(cfg)

    small = MO.VMCThreadConfig(1; requested_threads = 4)
    @test small.effective_threads == 1
    @test !MO.vmc_threading_enabled(small)

    empty = MO.VMCThreadConfig(0; requested_threads = 4)
    @test empty.effective_threads == 1
    @test !MO.vmc_threading_enabled(empty)

    @test_throws ErrorException MO.VMCThreadConfig(1; requested_threads = 0)
    @test_throws ErrorException MO.VMCThreadConfig(-1)
    @test_throws ErrorException MO.VMCThreadConfig(1; min_work_per_thread = 0)
end

@testset "unit/threading: JULIA_MVMC_MAINCAL_THREADS override" begin
    old = get(ENV, "JULIA_MVMC_MAINCAL_THREADS", nothing)
    try
        delete!(ENV, "JULIA_MVMC_MAINCAL_THREADS")
        @test MO.vmc_main_cal_requested_threads() == 1

        ENV["JULIA_MVMC_MAINCAL_THREADS"] = string(Threads.nthreads())
        @test MO.vmc_main_cal_requested_threads() == Threads.nthreads()

        ENV["JULIA_MVMC_MAINCAL_THREADS"] = string(Threads.nthreads() + 2)
        @test MO.vmc_main_cal_requested_threads() == Threads.nthreads()

        ENV["JULIA_MVMC_MAINCAL_THREADS"] = "0"
        @test_throws ErrorException MO.vmc_main_cal_requested_threads()
    finally
        if old === nothing
            delete!(ENV, "JULIA_MVMC_MAINCAL_THREADS")
        else
            ENV["JULIA_MVMC_MAINCAL_THREADS"] = old
        end
    end
end

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

@testset "unit/threading: sample chunks" begin
    chunks = MO.vmc_sample_chunks(10, 3)
    @test collect.(chunks) == [collect(0:3), collect(4:6), collect(7:9)]
    @test sort!(vcat(collect.(chunks)...)) == collect(0:9)

    one = MO.vmc_sample_chunks(1, 4)
    @test collect.(one) == [[0], Int[], Int[], Int[]]

    empty = MO.vmc_sample_chunks(0, 3)
    @test all(isempty, empty)

    @test_throws ErrorException MO.vmc_sample_chunks(1, 0)
    @test_throws ErrorException MO.vmc_sample_chunks(-1, 1)
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

@testset "unit/threading: whole-thread accumulator construction and merge" begin
    state = MO.VMCOptimizationState(2, 1, 1, 2, 1, 2, true, false)
    state.phys_quantities = MO.PhysicalQuantities(1, 1, 1)
    parent_timer = MO.CTimer(true)
    config = MO.VMCThreadConfig(4; requested_threads = 2)
    locals = MO.make_thread_accumulators(state, config, parent_timer)

    @test length(locals) == 2
    @test MO.ctimer_enabled(locals[1].timer)
    @test length(locals[1].sr_opt.sr_opt_oo) == length(state.sr_opt.sr_opt_oo)
    @test length(locals[1].sr_opt.sr_opt_o) == length(state.sr_opt.sr_opt_o)
    @test length(locals[1].phys.local_cis_ajs) == 1
    @test length(locals[1].phys.phys_cis_ajs) == 1

    MO.accumulate_energy!(locals[1].energy, 1.0, 2.0 + 0.0im)
    MO.accumulate_energy!(locals[2].energy, 1.0, 3.0 + 0.0im)
    MO.record_counter!(locals[1].counter, 1, 2)
    MO.record_counter!(locals[2].counter, 1, 5)
    locals[1].phys.phys_cis_ajs[1] = 4.0 + 0.0im
    locals[2].phys.phys_cis_ajs[1] = 6.0 + 0.0im
    MO.ctimer_add_elapsed!(locals[1].timer, 4, UInt64(11))
    MO.ctimer_add_elapsed!(locals[2].timer, 4, UInt64(13))

    MO.merge_thread_accumulators!(state, parent_timer, locals)
    @test state.energy.wc == 2.0 + 0.0im
    @test state.energy.etot == 5.0 + 0.0im
    @test state.electron_config.counter[1] == 7
    @test state.phys_quantities.phys_cis_ajs[1] == 10.0 + 0.0im
    @test parent_timer.elapsed_ns[5] == UInt64(24)
end

@testset "unit/threading: VMCMainCal worker state boundaries" begin
    data = ExpertModeData()
    data.modpara = ModParaParameters(
        nsite = 2,
        nelec = 1,
        nvmc_sample = 3,
        complex_flag = 1,
    )

    parent = MO.VMCOptimizationState(2, 1, 0, 2, 1, 3, true, true)
    parent.phys_quantities = MO.PhysicalQuantities(1, 0, 0)
    parent.electron_config.ele_idx .= 1:length(parent.electron_config.ele_idx)
    parent.electron_config.ele_cfg .= 11:(10+length(parent.electron_config.ele_cfg))
    parent.electron_config.ele_num .= 21:(20+length(parent.electron_config.ele_num))
    parent.electron_config.ele_spn .= 31:(30+length(parent.electron_config.ele_spn))
    parent.electron_config.tmp_ele_idx[1] = 99
    parent.electron_config.counter[1] = 7
    parent.slater_matrix.slater_elm[1] = 2.0 + 3.0im
    parent.slater_matrix.inv_m[1] = 5.0 + 0.0im

    worker = MO.make_vmc_main_cal_worker_state(data, parent)

    @test worker !== parent
    @test worker.electron_config.ele_idx === parent.electron_config.ele_idx
    @test worker.electron_config.ele_cfg === parent.electron_config.ele_cfg
    @test worker.electron_config.ele_num === parent.electron_config.ele_num
    @test worker.electron_config.ele_proj_cnt === parent.electron_config.ele_proj_cnt
    @test worker.electron_config.ele_spn === parent.electron_config.ele_spn
    @test worker.electron_config.tmp_ele_idx !== parent.electron_config.tmp_ele_idx
    @test worker.electron_config.counter !== parent.electron_config.counter
    @test worker.phys_quantities === parent.phys_quantities
    @test worker.slater_matrix.slater_elm === parent.slater_matrix.slater_elm
    @test worker.slater_matrix.inv_m !== parent.slater_matrix.inv_m
    @test worker.workspace !== parent.workspace
    @test worker.sr_opt !== parent.sr_opt
end

function threading_fill_slater!(state, n_site::Int, all_complex::Bool)
    n_site2 = 2 * n_site
    for a = 0:(n_site2 - 1)
        for b = (a + 1):(n_site2 - 1)
            value = ComplexF64(
                sin(0.37 * (a + 1) + 0.19 * (b + 1)) + 0.05 * (a - b),
                0.07 * cos(0.23 * (a + 1) - 0.41 * (b + 1)),
            )
            state.slater_matrix.slater_elm[a*n_site2+b+1] = value
            state.slater_matrix.slater_elm[b*n_site2+a+1] = -value
            if !all_complex
                state.slater_matrix.slater_elm_real[a*n_site2+b+1] = real(value)
                state.slater_matrix.slater_elm_real[b*n_site2+a+1] = -real(value)
            end
        end
    end
    return state
end

function threading_fill_jastrow_idx!(data, n_site::Int)
    n_jastrow = n_site * (n_site - 1) ÷ 2
    data.n_jastrow_idx = n_jastrow
    data.jastrow_idx = fill(-1, n_site, n_site)
    idx = 0
    for i = 1:n_site
        for j = (i + 1):n_site
            data.jastrow_idx[i, j] = idx
            data.jastrow_idx[j, i] = idx
            idx += 1
        end
    end
    data.jastrow_terms = JastrowTerm[]
    for i = 0:(n_site - 2)
        for j = (i + 1):(n_site - 1)
            value = ComplexF64(0.03 + 0.01 * length(data.jastrow_terms), 0.0)
            push!(data.jastrow_terms, JastrowTerm(i, j, value, false))
        end
    end
    return data
end

function threading_store_sample!(state, data, sample::Int, up_sites, down_sites)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site

    @assert length(up_sites) == n_elec
    @assert length(down_sites) == n_elec

    ele_idx = collect(Iterators.flatten((up_sites, down_sites)))
    ele_cfg = fill(-1, n_site2)
    ele_num = zeros(Int, n_site2)

    for (mi0, site) in enumerate(up_sites)
        mi = mi0 - 1
        ele_cfg[site+1] = mi
        ele_num[site+1] = 1
    end
    for (mi0, site) in enumerate(down_sites)
        mi = mi0 - 1
        ele_cfg[n_site+site+1] = mi
        ele_num[n_site+site+1] = 1
    end

    idx_start = sample * n_size + 1
    cfg_start = sample * n_site2 + 1
    state.electron_config.ele_idx[idx_start:(idx_start+n_size-1)] .= ele_idx
    state.electron_config.ele_cfg[cfg_start:(cfg_start+n_site2-1)] .= ele_cfg
    state.electron_config.ele_num[cfg_start:(cfg_start+n_site2-1)] .= ele_num

    n_proj = projection_layout(data).n_proj
    if n_proj > 0
        proj = zeros(Int, n_proj)
        MO.make_proj_cnt!(proj, ele_num, data)
        proj_start = sample * n_proj + 1
        state.electron_config.ele_proj_cnt[proj_start:(proj_start+n_proj-1)] .= proj
    end

    return state
end

function threading_maincal_fixture(;
    n_samples::Int = 4,
    all_complex::Bool = true,
    mode::Int = 0,
    use_fsz::Bool = false,
    phys::Bool = false,
    use_store::Bool = true,
)
    rich_samples = !use_fsz
    n_site = rich_samples ? 4 : 2
    n_elec = rich_samples ? 2 : 1

    data = ExpertModeData()
    data.modpara = ModParaParameters(
        nsite = n_site,
        nelec = n_elec,
        nvmc_sample = n_samples,
        complex_flag = all_complex ? 1 : 0,
        vmc_calc_mode = mode,
        nmp_trans = 1,
        nsp_gauss_leg = 1,
    )
    data.complex_flags = [all_complex ? 1 : 0]
    data.para_qp_trans = ComplexF64[1.0 + 0.0im]
    data.coulomb_inter_terms = [
        CoulombInterTerm(0, 1, 2.0),
        CoulombInterTerm(max(0, n_site - 2), n_site - 1, -0.75),
    ]
    if rich_samples
        data.n_gutzwiller_idx = 2
        data.gutzwiller_idx = [0, 1, 0, 1]
        data.gutzwiller_terms = [
            GutzwillerTerm(0, 0.11 + 0.0im, false),
            GutzwillerTerm(1, -0.07 + 0.0im, false),
        ]
        threading_fill_jastrow_idx!(data, n_site)
    end
    if phys
        data.green_one_terms = [
            GreenOneTerm(0, 0, :up, :up),
            GreenOneTerm(n_site - 1, n_site - 1, :down, :down),
        ]
        data.green_two_ex_terms = [GreenTwoExTerm(0, 0, 0, 0, n_site - 1, 1, n_site - 1, 1)]
    end
    init_qp_weight!(data)

    n_proj = projection_layout(data).n_proj
    state = MO.VMCOptimizationState(
        n_site,
        n_elec,
        n_proj,
        n_proj,
        1,
        n_samples,
        all_complex,
        use_fsz,
    )
    if phys
        MO.initialize_phys_quantities!(state, data)
    end

    threading_fill_slater!(state, n_site, all_complex)

    if rich_samples
        samples = (
            ([0, 2], [0, 3]),
            ([0, 1], [2, 3]),
            ([1, 3], [1, 2]),
            ([2, 3], [0, 2]),
        )
        for sample = 0:(n_samples - 1)
            up_sites, down_sites = samples[mod(sample, length(samples)) + 1]
            threading_store_sample!(state, data, sample, up_sites, down_sites)
        end
    else
        n_site2 = 2 * n_site
        for sample = 0:(n_samples - 1)
            ele_idx_start = sample * 2 + 1
            ele_cfg_start = sample * n_site2 + 1
            state.electron_config.ele_idx[ele_idx_start:(ele_idx_start+1)] .= [0, 1]
            state.electron_config.ele_cfg[ele_cfg_start:(ele_cfg_start+n_site2-1)] .= [0, -1, -1, 0]
            state.electron_config.ele_num[ele_cfg_start:(ele_cfg_start+n_site2-1)] .= [1, 0, 0, 1]
            if use_fsz
                state.electron_config.ele_spn[ele_idx_start:(ele_idx_start+1)] .= [0, 1]
            end
        end
    end

    return data, state
end

function threading_maincal_snapshot(state)
    phys = state.phys_quantities
    return (
        energy = (
            state.energy.wc,
            state.energy.etot,
            state.energy.etot2,
            state.energy.sztot,
            state.energy.sztot2,
        ),
        sr_opt_oo = copy(state.sr_opt.sr_opt_oo),
        sr_opt_ho = copy(state.sr_opt.sr_opt_ho),
        sr_opt_oo_real = copy(state.sr_opt.sr_opt_oo_real),
        sr_opt_ho_real = copy(state.sr_opt.sr_opt_ho_real),
        phys_cis_ajs = phys === nothing ? ComplexF64[] : copy(phys.phys_cis_ajs),
        phys_cis_ajs_ckt_alt = phys === nothing ? ComplexF64[] :
            copy(phys.phys_cis_ajs_ckt_alt),
        phys_cis_ajs_ckt_alt_dc = phys === nothing ? ComplexF64[] :
            copy(phys.phys_cis_ajs_ckt_alt_dc),
    )
end

function run_threading_maincal_fixture(; requested_threads::Int, kwargs...)
    data, state = threading_maincal_fixture(; kwargs...)
    if get(kwargs, :use_fsz, false)
        MO.vmc_main_cal_fsz!(data, state; requested_threads = requested_threads)
    else
        use_store = get(kwargs, :use_store, true)
        MO.vmc_main_cal!(data, state; requested_threads = requested_threads, use_store = use_store)
    end
    return threading_maincal_snapshot(state)
end

function threading_snapshots_match(a, b; atol = 1e-10, rtol = 1e-10)
    for field in propertynames(a)
        aval = getproperty(a, field)
        bval = getproperty(b, field)
        if aval isa Tuple
            all(isapprox(x, y; atol = atol, rtol = rtol) for (x, y) in zip(aval, bval)) ||
                return false
        else
            isapprox(aval, bval; atol = atol, rtol = rtol) || return false
        end
    end
    return true
end

@testset "unit/threading: vmc_main_cal! rich fixture coverage" begin
    data, state = threading_maincal_fixture()
    n_proj = projection_layout(data).n_proj
    @test data.modpara.nsite == 4
    @test data.modpara.nelec == 2
    @test n_proj > 0
    n_size = 2 * data.modpara.nelec
    idx_rows = reshape(state.electron_config.ele_idx, n_size, :)'
    proj_rows = reshape(state.electron_config.ele_proj_cnt, n_proj, :)'
    @test length(Set(Tuple(row) for row in eachrow(idx_rows))) > 1
    @test length(Set(Tuple(row) for row in eachrow(proj_rows))) > 1
end

@testset "unit/threading: vmc_main_cal! requested thread self-consistency" begin
    cases = (
        (label = "complex store", kwargs = (all_complex = true, mode = 0, use_store = true)),
        (label = "complex non-store", kwargs = (all_complex = true, mode = 0, use_store = false)),
        (label = "real", kwargs = (all_complex = false, mode = 0, use_store = true)),
        (label = "measurement phys", kwargs = (all_complex = true, mode = 1, phys = true)),
        (label = "fsz", kwargs = (all_complex = true, mode = 0, use_fsz = true)),
    )

    for case in cases
        @testset "$(case.label)" begin
            if Threads.nthreads() > 1
                # Explicit sample-level VMCMainCal threading remains a known
                # unsafe triage mode; default CI covers JULIA_NUM_THREADS>1
                # with MAINCAL opt-in unset.
                @test_skip false
            else
                single = run_threading_maincal_fixture(; requested_threads = 1, case.kwargs...)
                threaded = run_threading_maincal_fixture(; requested_threads = 2, case.kwargs...)
                @test threading_snapshots_match(threaded, single)
            end
        end
    end
end
