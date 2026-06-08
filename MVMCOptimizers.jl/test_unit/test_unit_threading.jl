using Test
using MVMCOptimizers
using MVMCExpertModeParsers:
    ExpertModeData,
    ModParaParameters,
    CoulombInterTerm,
    GreenOneTerm,
    GreenTwoExTerm,
    init_qp_weight!

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

function threading_maincal_fixture(;
    n_samples::Int = 4,
    all_complex::Bool = true,
    mode::Int = 0,
    use_fsz::Bool = false,
    phys::Bool = false,
    use_store::Bool = true,
)
    data = ExpertModeData()
    data.modpara = ModParaParameters(
        nsite = 2,
        nelec = 1,
        nvmc_sample = n_samples,
        complex_flag = all_complex ? 1 : 0,
        vmc_calc_mode = mode,
        nmp_trans = 1,
        nsp_gauss_leg = 1,
    )
    data.complex_flags = [all_complex ? 1 : 0]
    data.para_qp_trans = ComplexF64[1.0 + 0.0im]
    data.coulomb_inter_terms = [CoulombInterTerm(0, 1, 2.0)]
    if phys
        data.green_one_terms = [
            GreenOneTerm(0, 0, :up, :up),
            GreenOneTerm(1, 1, :down, :down),
        ]
        data.green_two_ex_terms = [GreenTwoExTerm(0, 0, 0, 0, 1, 1, 1, 1)]
    end
    init_qp_weight!(data)

    state = MO.VMCOptimizationState(2, 1, 0, 0, 1, n_samples, all_complex, use_fsz)
    if phys
        MO.initialize_phys_quantities!(state, data)
    end

    n_site2 = 4
    state.slater_matrix.slater_elm[0*n_site2+3+1] = 1.0 + 0.0im
    state.slater_matrix.slater_elm[3*n_site2+0+1] = -1.0 + 0.0im
    if !all_complex
        state.slater_matrix.slater_elm_real[0*n_site2+3+1] = 1.0
        state.slater_matrix.slater_elm_real[3*n_site2+0+1] = -1.0
    end

    for sample = 0:(n_samples - 1)
        ele_idx_start = sample * 2 + 1
        ele_cfg_start = sample * 4 + 1
        state.electron_config.ele_idx[ele_idx_start:(ele_idx_start+1)] .= [0, 1]
        state.electron_config.ele_cfg[ele_cfg_start:(ele_cfg_start+3)] .= [0, -1, -1, 0]
        state.electron_config.ele_num[ele_cfg_start:(ele_cfg_start+3)] .= [1, 0, 0, 1]
        if use_fsz
            state.electron_config.ele_spn[ele_idx_start:(ele_idx_start+1)] .= [0, 1]
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
            single = run_threading_maincal_fixture(; requested_threads = 1, case.kwargs...)
            threaded = run_threading_maincal_fixture(; requested_threads = 2, case.kwargs...)
            @test threaded == single
        end
    end
end
