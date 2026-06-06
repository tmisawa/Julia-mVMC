using Test
using MVMCOptimizers
using MVMCExpertModeParsers
using MVMCExpertModeParsers:
    ExpertModeData,
    ModParaParameters,
    GutzwillerTerm,
    JastrowTerm,
    OrbitalTerm,
    DoublonHolon2SiteIndex,
    DoublonHolon4SiteIndex,
    ChargeRBMPhysLayerTerm,
    SpinRBMPhysLayerTerm,
    GeneralRBMPhysLayerTerm,
    ChargeRBMHiddenLayerTerm,
    SpinRBMHiddenLayerTerm,
    GeneralRBMHiddenLayerTerm,
    ChargeRBMPhysHiddenTerm,
    SpinRBMPhysHiddenTerm,
    GeneralRBMPhysHiddenTerm

function make_mock_data_for_stochastic_opt_tests()
    data = ExpertModeData()

    # Proj terms
    data.gutzwiller_terms = [GutzwillerTerm(0, 1.0 + 0.0im, false)]
    data.jastrow_terms = [JastrowTerm(0, 1, 2.0 + 0.0im, false)]

    # RBM terms (9 sections, idx=0 each)
    # Keep duplicate idx in some sections to verify "update all terms with same idx".
    data.charge_rbm_phys_layer_terms = [
        ChargeRBMPhysLayerTerm(0, 10.0 + 0.0im, false, 0),
        ChargeRBMPhysLayerTerm(1, 11.0 + 0.0im, false, 0),
    ]
    data.spin_rbm_phys_layer_terms = [SpinRBMPhysLayerTerm(0, 20.0 + 0.0im, false, 0)]
    data.general_rbm_phys_layer_terms = [GeneralRBMPhysLayerTerm(0, 1, 30.0 + 0.0im, false, 0)]

    data.charge_rbm_hidden_layer_terms = [ChargeRBMHiddenLayerTerm(0, 40.0 + 0.0im, false, 0)]
    data.spin_rbm_hidden_layer_terms = [SpinRBMHiddenLayerTerm(0, 50.0 + 0.0im, false, 0)]
    data.general_rbm_hidden_layer_terms = [GeneralRBMHiddenLayerTerm(0, 60.0 + 0.0im, false, 0)]

    data.charge_rbm_phys_hidden_terms = [ChargeRBMPhysHiddenTerm(0, 0, 70.0 + 0.0im, false, 0)]
    data.spin_rbm_phys_hidden_terms = [SpinRBMPhysHiddenTerm(0, 0, 80.0 + 0.0im, false, 0)]
    data.general_rbm_phys_hidden_terms = [
        GeneralRBMPhysHiddenTerm(0, 1, 0, 90.0 + 0.0im, false, 0),
        GeneralRBMPhysHiddenTerm(1, 0, 0, 91.0 + 0.0im, false, 0),
    ]

    # Slater terms (idx=0 duplicated, idx=1 single)
    data.orbital_terms = [
        OrbitalTerm(0, 1, 0, 100.0 + 0.0im, false, 1),
        OrbitalTerm(1, 2, 0, 101.0 + 0.0im, false, 1),
        OrbitalTerm(0, 2, 1, 110.0 + 0.0im, false, 1),
    ]

    # Needed only for test readability (not required by update_parameter_value itself)
    data.modpara = ModParaParameters(n_orbital_idx = 2)

    return data
end

@testset "unit/stochastic_opt: update_parameter_value DH write-back and shifted offsets" begin
    data = make_mock_data_for_stochastic_opt_tests()
    data.doublon_holon_2site_indices = [
        DoublonHolon2SiteIndex([
            1 0
            0 1
        ]),
    ]
    data.doublon_holon_2site_params = [ComplexF64(1000 + i, 0) for i = 1:6]
    data.doublon_holon_2site_opt_flags = fill(true, 6)

    delta = -0.125 + 0.375im
    orig_charge_phys = [t.value for t in data.charge_rbm_phys_layer_terms]
    orig_orbital = [t.value for t in data.orbital_terms]

    # NProj = Gutz(1) + Jastrow(1) + DH2(6) = 8. DH occupies para_idx 3..8.
    MVMCOptimizers.update_parameter_value(data, 3, real(delta), imag(delta))
    @test data.doublon_holon_2site_params[1] == 1001.0 + 0.0im + delta
    @test [t.value for t in data.charge_rbm_phys_layer_terms] == orig_charge_phys
    @test [t.value for t in data.orbital_terms] == orig_orbital

    MVMCOptimizers.update_parameter_value(data, 8, real(delta), imag(delta))
    @test data.doublon_holon_2site_params[6] == 1006.0 + 0.0im + delta

    # RBM and Slater offsets move after the DH projection slice.
    MVMCOptimizers.update_parameter_value(data, 9, real(delta), imag(delta))
    @test data.charge_rbm_phys_layer_terms[1].value == orig_charge_phys[1] + delta
    @test data.charge_rbm_phys_layer_terms[2].value == orig_charge_phys[2] + delta

    MVMCOptimizers.update_parameter_value(data, 18, real(delta), imag(delta))
    @test data.orbital_terms[1].value == orig_orbital[1] + delta
    @test data.orbital_terms[2].value == orig_orbital[2] + delta
    @test data.orbital_terms[3].value == orig_orbital[3]
end

@testset "unit/stochastic_opt: update_parameter_value DH4 write-back" begin
    data = make_mock_data_for_stochastic_opt_tests()
    data.doublon_holon_4site_indices = [
        DoublonHolon4SiteIndex([
            1 0 1 0
            0 1 0 1
        ]),
    ]
    data.doublon_holon_4site_params = [ComplexF64(2000 + i, 0) for i = 1:10]
    data.doublon_holon_4site_opt_flags = fill(true, 10)

    delta = 0.625 - 0.25im
    orig_charge_phys = [t.value for t in data.charge_rbm_phys_layer_terms]
    orig_orbital = [t.value for t in data.orbital_terms]

    # NProj = Gutz(1) + Jastrow(1) + DH4(10) = 12. DH4 occupies para_idx 3..12.
    MVMCOptimizers.update_parameter_value(data, 3, real(delta), imag(delta))
    @test data.doublon_holon_4site_params[1] == 2001.0 + 0.0im + delta
    @test [t.value for t in data.charge_rbm_phys_layer_terms] == orig_charge_phys
    @test [t.value for t in data.orbital_terms] == orig_orbital

    MVMCOptimizers.update_parameter_value(data, 12, real(delta), imag(delta))
    @test data.doublon_holon_4site_params[10] == 2010.0 + 0.0im + delta

    # RBM and Slater offsets move after the DH4 projection slice.
    MVMCOptimizers.update_parameter_value(data, 13, real(delta), imag(delta))
    @test data.charge_rbm_phys_layer_terms[1].value == orig_charge_phys[1] + delta
    @test data.charge_rbm_phys_layer_terms[2].value == orig_charge_phys[2] + delta

    MVMCOptimizers.update_parameter_value(data, 22, real(delta), imag(delta))
    @test data.orbital_terms[1].value == orig_orbital[1] + delta
    @test data.orbital_terms[2].value == orig_orbital[2] + delta
    @test data.orbital_terms[3].value == orig_orbital[3]
end

@testset "unit/stochastic_opt: SR enumeration writes DH without touching RBM or Slater" begin
    data = make_mock_data_for_stochastic_opt_tests()
    data.modpara = ModParaParameters(
        nsite = 2,
        nelec = 1,
        nvmc_sample = 1,
        n_orbital_idx = 2,
        dsr_opt_red_cut = 0.0,
        dsr_opt_sta_del = 0.0,
        dsr_opt_step_dt = 0.25,
    )
    data.complex_flags = [1]
    data.doublon_holon_2site_indices = [
        DoublonHolon2SiteIndex([
            1 0
            0 1
        ]),
    ]
    data.doublon_holon_2site_params = [ComplexF64(3000 + i, 0) for i = 1:6]
    data.doublon_holon_2site_opt_flags = fill(true, 6)

    layout = MVMCExpertModeParsers.projection_layout(data)
    n_rbm = MVMCExpertModeParsers.count_rbm_parameters(data)
    n_para = layout.n_proj + n_rbm + data.modpara.n_orbital_idx
    state = MVMCOptimizers.VMCOptimizationState(
        data.modpara.nsite,
        data.modpara.nelec,
        layout.n_proj,
        n_para,
        1,
        data.modpara.nvmc_sample,
        true,
        false,
    )

    target_para_idx = layout.dh2_offset + 1  # 1-based parameter index
    target_pi = 2 * (target_para_idx - 1)    # C/Julia SR real component index, 0-based
    data.optimization_flags = falses(2 * n_para)
    data.optimization_flags[target_pi + 1] = true

    sr_opt_size = state.sr_opt.sr_opt_size
    diag_idx = (target_pi + 2) * (2 * sr_opt_size) + (target_pi + 2) + 1
    state.sr_opt.sr_opt_oo[diag_idx] = 2.0 + 0.0im
    state.sr_opt.sr_opt_ho[target_pi + 3] = 3.0 + 0.0im

    orig_dh = copy(data.doublon_holon_2site_params)
    orig_rbm = (
        [t.value for t in data.charge_rbm_phys_layer_terms],
        [t.value for t in data.spin_rbm_phys_layer_terms],
        [t.value for t in data.general_rbm_phys_layer_terms],
        [t.value for t in data.charge_rbm_hidden_layer_terms],
        [t.value for t in data.spin_rbm_hidden_layer_terms],
        [t.value for t in data.general_rbm_hidden_layer_terms],
        [t.value for t in data.charge_rbm_phys_hidden_terms],
        [t.value for t in data.spin_rbm_phys_hidden_terms],
        [t.value for t in data.general_rbm_phys_hidden_terms],
    )
    orig_orbital = [t.value for t in data.orbital_terms]

    @test MVMCOptimizers.stochastic_opt!(data, state) == 0

    expected_delta = -0.75 + 0.0im
    @test data.doublon_holon_2site_params[1] ≈ orig_dh[1] + expected_delta atol = 1e-14
    @test data.doublon_holon_2site_params[2:end] == orig_dh[2:end]
    @test (
        [t.value for t in data.charge_rbm_phys_layer_terms],
        [t.value for t in data.spin_rbm_phys_layer_terms],
        [t.value for t in data.general_rbm_phys_layer_terms],
        [t.value for t in data.charge_rbm_hidden_layer_terms],
        [t.value for t in data.spin_rbm_hidden_layer_terms],
        [t.value for t in data.general_rbm_hidden_layer_terms],
        [t.value for t in data.charge_rbm_phys_hidden_terms],
        [t.value for t in data.spin_rbm_phys_hidden_terms],
        [t.value for t in data.general_rbm_phys_hidden_terms],
    ) == orig_rbm
    @test [t.value for t in data.orbital_terms] == orig_orbital
end

@testset "unit/stochastic_opt: SR enumeration writes OptTrans after Slater" begin
    data = ExpertModeData()
    data.modpara = ModParaParameters(
        nsite = 1,
        nelec = 1,
        nvmc_sample = 1,
        n_orbital_idx = 1,
        complex_flag = 1,
        dsr_opt_red_cut = 0.0,
        dsr_opt_sta_del = 0.0,
        dsr_opt_step_dt = 0.25,
    )
    data.orbital_terms = [
        OrbitalTerm(0, 0, 0, 10.0 + 0.0im, false, 1),
    ]
    data.opt_trans = ComplexF64[
        1.0 + 0.0im,
        2.0 + 0.0im,
    ]
    data.qp_weights = MVMCExpertModeParsers.QuantumProjectionWeights()
    data.qp_weights.qp_fix_weight = ComplexF64[2.0 + 0.0im]

    layout = MVMCExpertModeParsers.projection_layout(data)
    n_para =
        layout.n_proj +
        data.modpara.n_orbital_idx +
        MVMCExpertModeParsers.count_opt_trans_parameters(data)
    state = MVMCOptimizers.VMCOptimizationState(
        data.modpara.nsite,
        data.modpara.nelec,
        layout.n_proj,
        n_para,
        1,
        data.modpara.nvmc_sample,
        true,
        false,
    )

    # Target OptTrans[2], after the Slater parameter:
    # para_idx = NProj(0) + NSlater(1) + optidx(2) = 3.
    target_para_idx = 3
    target_pi = 2 * (target_para_idx - 1)
    data.optimization_flags = falses(2 * n_para)
    data.optimization_flags[target_pi + 1] = true

    sr_opt_size = state.sr_opt.sr_opt_size
    diag_idx = (target_pi + 2) * (2 * sr_opt_size) + (target_pi + 2) + 1
    state.sr_opt.sr_opt_oo[diag_idx] = 2.0 + 0.0im
    state.sr_opt.sr_opt_ho[target_pi + 3] = 3.0 + 0.0im

    orig_orbital = [t.value for t in data.orbital_terms]
    orig_opt_trans = copy(data.opt_trans)

    info = MVMCOptimizers.stochastic_opt!(data, state)

    @test info == 0
    # g = -0.25 * 2 * 3 = -1.5; S = 2 => delta = -0.75.
    @test data.opt_trans[1] == orig_opt_trans[1]
    @test data.opt_trans[2] == orig_opt_trans[2] - 0.75
    @test [t.value for t in data.orbital_terms] == orig_orbital
    @test data.qp_weights.qp_full_weight == ComplexF64[
        2.0 + 0.0im,
        2.5 + 0.0im,
    ]
end

@testset "unit/stochastic_opt: get_opt_flag_for_parameter" begin
    data = ExpertModeData()
    @test MVMCOptimizers.get_opt_flag_for_parameter(data, 0) == 0
    @test MVMCOptimizers.get_opt_flag_for_parameter(data, 10) == 0

    data.optimization_flags = [true, false, true]
    @test MVMCOptimizers.get_opt_flag_for_parameter(data, 0) == 1
    @test MVMCOptimizers.get_opt_flag_for_parameter(data, 1) == 0
    @test MVMCOptimizers.get_opt_flag_for_parameter(data, 2) == 1
    @test MVMCOptimizers.get_opt_flag_for_parameter(data, 3) == 0
end

@testset "unit/stochastic_opt: update_parameter_value Proj/RBM/Slater mapping" begin
    data = make_mock_data_for_stochastic_opt_tests()
    delta = 0.25 - 0.5im

    # NProj = 2
    MVMCOptimizers.update_parameter_value(data, 1, real(delta), imag(delta))
    @test data.gutzwiller_terms[1].value == 1.0 + 0.0im + delta
    @test data.jastrow_terms[1].value == 2.0 + 0.0im

    MVMCOptimizers.update_parameter_value(data, 2, real(delta), imag(delta))
    @test data.jastrow_terms[1].value == 2.0 + 0.0im + delta

    # NRBM = 9 (one idx per section); para_idx 3..11 map to RBM block.
    # First RBM parameter -> charge phys idx=0 (all matching terms updated).
    MVMCOptimizers.update_parameter_value(data, 3, real(delta), imag(delta))
    @test data.charge_rbm_phys_layer_terms[1].value == 10.0 + 0.0im + delta
    @test data.charge_rbm_phys_layer_terms[2].value == 11.0 + 0.0im + delta
    @test data.spin_rbm_phys_layer_terms[1].value == 20.0 + 0.0im

    # Last RBM parameter -> general phys-hidden idx=0 (all matching terms updated).
    MVMCOptimizers.update_parameter_value(data, 11, real(delta), imag(delta))
    @test data.general_rbm_phys_hidden_terms[1].value == 90.0 + 0.0im + delta
    @test data.general_rbm_phys_hidden_terms[2].value == 91.0 + 0.0im + delta
    @test data.charge_rbm_phys_hidden_terms[1].value == 70.0 + 0.0im

    # Slater starts at para_idx = NProj + NRBM + 1 = 12
    MVMCOptimizers.update_parameter_value(data, 12, real(delta), imag(delta))
    @test data.orbital_terms[1].value == 100.0 + 0.0im + delta
    @test data.orbital_terms[2].value == 101.0 + 0.0im + delta
    @test data.orbital_terms[3].value == 110.0 + 0.0im

    MVMCOptimizers.update_parameter_value(data, 13, real(delta), imag(delta))
    @test data.orbital_terms[3].value == 110.0 + 0.0im + delta

    data.opt_trans = [1.0 + 0.0im, 2.0 + 0.0im]
    MVMCOptimizers.update_parameter_value(data, 14, real(delta), imag(delta))
    @test data.opt_trans[1] == 1.0 + 0.0im + delta
    @test data.opt_trans[2] == 2.0 + 0.0im

    MVMCOptimizers.update_parameter_value(data, 15, real(delta), imag(delta))
    @test data.opt_trans[2] == 2.0 + 0.0im + delta
end

@testset "unit/stochastic_opt: build_s_matrix_and_g_vector!" begin
    # Small hand-checkable layout.
    # sr_opt_size=2 -> lda_oo = 4, and function reads:
    #   OO[0][pi+1] at sr_opt_oo[pi+3]
    #   OO[pi+1][pj+1] at sr_opt_oo[(pi+2)*4 + (pj+2) + 1]
    sr_opt_size = 2
    sr_opt_oo = zeros(ComplexF64, 16)
    sr_opt_ho = zeros(ComplexF64, 8)

    # OO first row
    sr_opt_oo[3] = 0.5 + 0im    # OO[0][1]
    sr_opt_oo[4] = -0.2 + 0im   # OO[0][2]

    # OO block used by S
    sr_opt_oo[11] = 2.0 + 0im   # OO[1][1]
    sr_opt_oo[12] = 0.3 + 0im   # OO[1][2]
    sr_opt_oo[15] = 0.4 + 0im   # OO[2][1]
    sr_opt_oo[16] = 1.5 + 0im   # OO[2][2]

    # HO
    sr_opt_ho[1] = 1.2 + 0im    # HO[0]
    sr_opt_ho[3] = 2.3 + 0im    # HO[1]
    sr_opt_ho[4] = -0.7 + 0im   # HO[2]

    smat_to_para_idx = [0, 1]
    S = zeros(Float64, 2, 2)
    g = zeros(Float64, 2)

    dsr_opt_sta_del = 0.1
    dsr_opt_step_dt = 0.05

    MVMCOptimizers.build_s_matrix_and_g_vector!(
        S,
        g,
        smat_to_para_idx,
        sr_opt_oo,
        sr_opt_ho,
        sr_opt_size,
        dsr_opt_sta_del,
        dsr_opt_step_dt,
    )

    @test S[1, 1] ≈ 1.925 atol = 1e-12  # (2.0 - 0.5^2) * 1.1
    @test S[1, 2] ≈ 0.4 atol = 1e-12    # 0.3 - 0.5*(-0.2)
    @test S[2, 1] ≈ 0.5 atol = 1e-12    # 0.4 - (-0.2)*0.5
    @test S[2, 2] ≈ 1.606 atol = 1e-12  # (1.5 - (-0.2)^2) * 1.1

    @test g[1] ≈ -0.17 atol = 1e-12     # -0.1*(2.3 - 1.2*0.5)
    @test g[2] ≈ 0.046 atol = 1e-12     # -0.1*(-0.7 - 1.2*(-0.2))
end
