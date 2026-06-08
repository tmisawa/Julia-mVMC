using Test
using MVMCOptimizers
using MVMCExpertModeParsers
using MVMCExpertModeParsers:
    ExpertModeData,
    ModParaParameters,
    ChargeRBMPhysLayerTerm,
    ChargeRBMHiddenLayerTerm,
    ChargeRBMPhysHiddenTerm

if !isdefined(@__MODULE__, :make_minimal_data_for_rbm_diff_tests)
    include(joinpath(@__DIR__, "helpers", "mock_data.jl"))
end

@testset "unit/vmc_main_cal: set_projection_diff!" begin
    n_proj = 3
    ele_proj_cnt = [5, 0, -2]
    sr_opt_o = fill(99.0 + 99.0im, 2 * (n_proj + 1))

    MVMCOptimizers.set_projection_diff!(sr_opt_o, ele_proj_cnt, n_proj)

    @test sr_opt_o[1] == 1.0 + 0.0im
    @test sr_opt_o[2] == 0.0 + 0.0im
    for i = 0:(n_proj - 1)
        @test sr_opt_o[(i + 1) * 2 + 1] == ComplexF64(ele_proj_cnt[i + 1])
        @test sr_opt_o[(i + 1) * 2 + 2] == 0.0 + 0.0im
    end
end

@testset "unit/vmc_main_cal: set_rbm_diff!" begin
    @testset "physical layer only" begin
        data = make_minimal_data_for_rbm_diff_tests(nsite = 1, nneuron_charge = 0)
        data.charge_rbm_phys_layer_terms = [ChargeRBMPhysLayerTerm(0, 0.0 + 0.0im, false, 0)]

        rbm_cnt = ComplexF64[2.0 + 3.0im]
        ele_num = Int[1, 0]  # length = 2*nsite
        sr_opt_o = zeros(ComplexF64, 2)  # 2*n_rbm (n_rbm=1)

        MVMCOptimizers.set_rbm_diff!(sr_opt_o, rbm_cnt, ele_num, data)

        @test sr_opt_o[1] == rbm_cnt[1]
        @test sr_opt_o[2] == im * rbm_cnt[1]
    end

    @testset "physical + hidden + phys-hidden (charge)" begin
        data = make_minimal_data_for_rbm_diff_tests(nsite = 1, nneuron_charge = 1)
        data.charge_rbm_phys_layer_terms = [ChargeRBMPhysLayerTerm(0, 0.0 + 0.0im, false, 0)]
        data.charge_rbm_hidden_layer_terms = [ChargeRBMHiddenLayerTerm(0, 0.0 + 0.0im, false, 0)]
        data.charge_rbm_phys_hidden_terms = [ChargeRBMPhysHiddenTerm(0, 0, 0.0 + 0.0im, false, 0)]

        # rbm_cnt layout (see make_rbm_cnt):
        # [phys params... ; charge hidden neuron(s) ...]
        rbm_cnt_phys = 0.1 + 0.0im
        rbm_cnt_hidden = 0.7 + 0.0im
        rbm_cnt = ComplexF64[rbm_cnt_phys, rbm_cnt_hidden]

        ele_num = Int[1, 1]  # up=1, down=1 at site 0 => xi = 1
        sr_opt_o = zeros(ComplexF64, 6)  # 2*n_rbm with n_rbm = n_phys(1)+n_hidden(1)+n_ph(1) = 3

        MVMCOptimizers.set_rbm_diff!(sr_opt_o, rbm_cnt, ele_num, data)

        # Physical-layer derivative for idx0=0
        @test sr_opt_o[1] == rbm_cnt_phys
        @test sr_opt_o[2] == im * rbm_cnt_phys

        # Hidden-layer parameter derivative (adds tanh(hidden counter))
        t = tanh(rbm_cnt_hidden)
        @test sr_opt_o[3] == t
        @test sr_opt_o[4] == im * t

        # Phys-hidden coupling derivative: xi * tanh(hidden counter)
        xi = (ele_num[1] + ele_num[2] - 1)
        @test sr_opt_o[5] == xi * t
        @test sr_opt_o[6] == im * (xi * t)
    end
end

@testset "unit/vmc_main_cal: opt_trans_diff!" begin
    data = ExpertModeData()
    data.modpara = ModParaParameters(complex_flag = 1)
    data.opt_trans = [1.0 + 0.0im, 2.0 + 0.0im]
    data.qp_weights = MVMCExpertModeParsers.QuantumProjectionWeights()
    data.qp_weights.qp_fix_weight = ComplexF64[
        2.0 + 0.0im,
        3.0 + 0.0im,
    ]

    state = MVMCOptimizers.VMCOptimizationState(1, 1, 0, 2, 4, 1, true, false)
    state.slater_matrix.pf_m .= ComplexF64[
        1.0 + 1.0im,
        2.0 - 1.0im,
        -1.0 + 0.5im,
        0.25 - 2.0im,
    ]

    ip = 5.0 - 2.0im
    sr_opt_o = fill(99.0 + 99.0im, 2 * length(data.opt_trans))

    MVMCOptimizers.opt_trans_diff!(sr_opt_o, ip, data, state)

    expected1 =
        (data.qp_weights.qp_fix_weight[1] * state.slater_matrix.pf_m[1] +
         data.qp_weights.qp_fix_weight[2] * state.slater_matrix.pf_m[2]) / ip
    expected2 =
        (data.qp_weights.qp_fix_weight[1] * state.slater_matrix.pf_m[3] +
         data.qp_weights.qp_fix_weight[2] * state.slater_matrix.pf_m[4]) / ip

    @test sr_opt_o[1] == expected1
    @test sr_opt_o[2] == im * expected1
    @test sr_opt_o[3] == expected2
    @test sr_opt_o[4] == im * expected2
end

@testset "unit/vmc_main_cal: calculate_oo!" begin
    sr_opt_size = 3
    size_2 = 2 * sr_opt_size

    sr_opt_o = ComplexF64[
        1.0 + 0.0im,
        2.0 + 0.0im,
        3.0 + 4.0im,
        -1.0 + 2.0im,
        0.5 - 0.3im,
        0.0 + 1.0im,
    ]
    @test length(sr_opt_o) == size_2

    w = 0.25
    e = 1.2 - 0.5im

    sr_opt_oo = zeros(ComplexF64, size_2 * size_2)
    sr_opt_ho = zeros(ComplexF64, size_2)

    MVMCOptimizers.calculate_oo!(sr_opt_oo, sr_opt_ho, sr_opt_o, w, e, sr_opt_size)

    # First row: <O>
    for j = 0:(size_2 - 1)
        @test sr_opt_oo[j + 1] == w * sr_opt_o[j + 1]
        @test sr_opt_ho[j + 1] == e * (w * sr_opt_o[j + 1])
    end

    # Row i=1 (0-based) is intentionally not updated by calculate_oo!
    for j = 0:(size_2 - 1)
        @test sr_opt_oo[1 * size_2 + j + 1] == 0.0 + 0.0im
    end

    # Rows i=2..size_2-1: <O†O>
    for i = 2:(size_2 - 1)
        for j = 0:(size_2 - 1)
            idx = i * size_2 + j + 1
            @test sr_opt_oo[idx] == w * sr_opt_o[j + 1] * conj(sr_opt_o[i + 1])
        end
    end
end

@testset "unit/vmc_main_cal: calculate_oo_real!" begin
    sr_opt_size = 4
    sr_opt_o = Float64[1.0, -2.0, 0.5, 3.0]
    w = 0.25
    e = -1.5

    # Use non-zero initial values to ensure in-place += semantics are preserved.
    sr_opt_oo = fill(0.1, sr_opt_size * sr_opt_size)
    sr_opt_ho = fill(-0.2, sr_opt_size)

    MVMCOptimizers.calculate_oo_real!(sr_opt_oo, sr_opt_ho, sr_opt_o, w, e, sr_opt_size)

    expected_oo = fill(0.1, sr_opt_size * sr_opt_size)
    expected_oo .+= vec(w .* (sr_opt_o * transpose(sr_opt_o)))
    @test sr_opt_oo ≈ expected_oo

    expected_ho = fill(-0.2, sr_opt_size)
    expected_ho .+= (w * e) .* sr_opt_o
    @test sr_opt_ho ≈ expected_ho
end

@testset "unit/vmc_main_cal: threaded OO/store inner loops match sequential" begin
    sr_opt_size = 40
    size_2 = 2 * sr_opt_size
    w = 0.375
    e = 1.1 - 0.25im

    sr_opt_o = ComplexF64[
        ComplexF64(sin(0.17 * i), cos(0.11 * i))
        for i = 1:size_2
    ]

    oo_seq = fill(0.2 + 0.1im, size_2 * size_2)
    ho_seq = fill(-0.3 + 0.4im, size_2)
    oo_thr = copy(oo_seq)
    ho_thr = copy(ho_seq)

    MVMCOptimizers.calculate_oo!(
        oo_seq,
        ho_seq,
        sr_opt_o,
        w,
        e,
        sr_opt_size;
        threaded = false,
    )
    MVMCOptimizers.calculate_oo!(
        oo_thr,
        ho_thr,
        sr_opt_o,
        w,
        e,
        sr_opt_size;
        threaded = true,
    )
    @test oo_thr == oo_seq
    @test ho_thr == ho_seq

    sr_opt_o_real = [sin(0.13 * i) for i = 1:sr_opt_size]
    oo_real_seq = fill(0.2, sr_opt_size * sr_opt_size)
    ho_real_seq = fill(-0.3, sr_opt_size)
    oo_real_thr = copy(oo_real_seq)
    ho_real_thr = copy(ho_real_seq)
    MVMCOptimizers.calculate_oo_real!(
        oo_real_seq,
        ho_real_seq,
        sr_opt_o_real,
        w,
        real(e),
        sr_opt_size;
        threaded = false,
    )
    MVMCOptimizers.calculate_oo_real!(
        oo_real_thr,
        ho_real_thr,
        sr_opt_o_real,
        w,
        real(e),
        sr_opt_size;
        threaded = true,
    )
    @test oo_real_thr == oo_real_seq
    @test ho_real_thr == ho_real_seq

    n_samples = 5
    sample = 2
    store_seq = fill(0.0 + 0.0im, size_2 * n_samples)
    store_thr = copy(store_seq)
    ho_store_seq = fill(0.1 - 0.2im, size_2)
    ho_store_thr = copy(ho_store_seq)
    dummy_oo = ComplexF64[]
    MVMCOptimizers.calculate_oo_store!(
        dummy_oo,
        ho_store_seq,
        store_seq,
        sr_opt_o,
        w,
        e,
        sample,
        sr_opt_size;
        threaded = false,
    )
    MVMCOptimizers.calculate_oo_store!(
        dummy_oo,
        ho_store_thr,
        store_thr,
        sr_opt_o,
        w,
        e,
        sample,
        sr_opt_size;
        threaded = true,
    )
    @test store_thr == store_seq
    @test ho_store_thr == ho_store_seq

    final_seq = fill(0.0 + 0.0im, size_2 * size_2)
    final_thr = copy(final_seq)
    full_store = ComplexF64[
        ComplexF64(sin(0.07 * i), cos(0.05 * i))
        for i = 1:(size_2 * n_samples)
    ]
    MVMCOptimizers.finalize_oo_store!(
        final_seq,
        full_store,
        sr_opt_size,
        n_samples;
        threaded = false,
    )
    MVMCOptimizers.finalize_oo_store!(
        final_thr,
        full_store,
        sr_opt_size,
        n_samples;
        threaded = true,
    )
    @test final_thr == final_seq
end
