using Test
using MVMCOptimizers
using MVMCExpertModeParsers
using MVMCExpertModeParsers:
    ExpertModeData,
    ModParaParameters,
    ChargeRBMPhysLayerTerm,
    SpinRBMPhysLayerTerm,
    GeneralRBMPhysLayerTerm,
    ChargeRBMHiddenLayerTerm,
    SpinRBMHiddenLayerTerm,
    GeneralRBMHiddenLayerTerm,
    ChargeRBMPhysHiddenTerm,
    SpinRBMPhysHiddenTerm,
    GeneralRBMPhysHiddenTerm

if !isdefined(@__MODULE__, :make_ele_num)
    include(joinpath(@__DIR__, "helpers", "mock_state.jl"))
end
if !isdefined(@__MODULE__, :make_mock_data_for_rbm_tests)
    include(joinpath(@__DIR__, "helpers", "mock_data.jl"))
end

@testset "unit/vmc_sampling: has_rbm_terms" begin
    data = ExpertModeData()
    data.modpara = ModParaParameters(nsite = 2)
    @test !MVMCOptimizers.has_rbm_terms(data)

    data.charge_rbm_hidden_layer_terms = [ChargeRBMHiddenLayerTerm(0, 0.1 + 0im, false, 0)]
    @test MVMCOptimizers.has_rbm_terms(data)
end

@testset "unit/vmc_sampling: make_rbm_cnt / update_rbm_cnt_hopping! / log_rbm_* consistency" begin
    data = make_mock_data_for_rbm_tests(nsite = 4)
    nsite = data.modpara.nsite

    # Electron occupation (length = 2*nsite): up[1:nsite], down[nsite+1:2*nsite]
    # up: sites 0,2 ; down: sites 1,3
    ele_num0 = make_ele_num(nsite; up_sites = [0, 2], down_sites = [1, 3])

    @test MVMCOptimizers.log_rbm_val(ele_num0, data) isa ComplexF64

    rbm_cnt0 = MVMCOptimizers.make_rbm_cnt(ele_num0, data)
    @test !isempty(rbm_cnt0)

    # One hop: up electron 0 -> 1 (s=0)
    ri, rj, s = 0, 1, 0
    @test ele_num0[ri + 1] == 1
    @test ele_num0[rj + 1] == 0
    ele_num1 = apply_hop(ele_num0, ri, rj, s, nsite)

    rbm_cnt1_full = MVMCOptimizers.make_rbm_cnt(ele_num1, data)
    rbm_cnt1_inc = similar(rbm_cnt0)
    MVMCOptimizers.update_rbm_cnt_hopping!(rbm_cnt1_inc, rbm_cnt0, ri, rj, s, data)
    @test rbm_cnt1_inc ≈ rbm_cnt1_full rtol = 0 atol = 1e-12

    # log ratio should match log-value difference
    d1_ratio = MVMCOptimizers.log_rbm_ratio(rbm_cnt1_full, rbm_cnt0, data)
    d1_val = MVMCOptimizers.log_rbm_val(ele_num1, data) - MVMCOptimizers.log_rbm_val(ele_num0, data)
    @test d1_ratio ≈ d1_val rtol = 0 atol = 1e-12

    # No-op hop (ri==rj) keeps counters unchanged
    rbm_cnt_noop = similar(rbm_cnt0)
    MVMCOptimizers.update_rbm_cnt_hopping!(rbm_cnt_noop, rbm_cnt0, 2, 2, 0, data)
    @test rbm_cnt_noop == rbm_cnt0

    # Two hops (exchange-like): apply two incremental updates, compare to full recompute.
    # Hop 2: down electron 3 -> 0 (s=1)
    ri2, rj2, s2 = 3, 0, 1
    @test ele_num1[ri2 + s2 * nsite + 1] == 1
    @test ele_num1[rj2 + s2 * nsite + 1] == 0
    ele_num2 = apply_hop(ele_num1, ri2, rj2, s2, nsite)

    rbm_cnt2_full = MVMCOptimizers.make_rbm_cnt(ele_num2, data)
    rbm_cnt2_inc = similar(rbm_cnt0)
    MVMCOptimizers.update_rbm_cnt_hopping!(rbm_cnt2_inc, rbm_cnt0, ri, rj, s, data)
    MVMCOptimizers.update_rbm_cnt_hopping!(rbm_cnt2_inc, rbm_cnt2_inc, ri2, rj2, s2, data)
    @test rbm_cnt2_inc ≈ rbm_cnt2_full rtol = 0 atol = 1e-12

    d2_ratio = MVMCOptimizers.log_rbm_ratio(rbm_cnt2_full, rbm_cnt0, data)
    d2_val = MVMCOptimizers.log_rbm_val(ele_num2, data) - MVMCOptimizers.log_rbm_val(ele_num0, data)
    @test d2_ratio ≈ d2_val rtol = 0 atol = 1e-12
end
