using Test
using MVMCExpertModeParsers

@testset "contract/read_input_parameters: RBM layout (count/offsets)" begin
    data = MVMCExpertModeParsers.ExpertModeData()

    # Build sparse idx patterns to ensure we count `maximum(idx)+1` (not just number of terms).
    # Section order must remain:
    #   charge_phys, spin_phys, general_phys,
    #   charge_hidden, spin_hidden, general_hidden,
    #   charge_phys_hidden, spin_phys_hidden, general_phys_hidden
    data.charge_rbm_phys_layer_terms = [
        MVMCExpertModeParsers.ChargeRBMPhysLayerTerm(0, 0.0 + 0.0im, false, 0),
        MVMCExpertModeParsers.ChargeRBMPhysLayerTerm(1, 0.0 + 0.0im, false, 2),
    ] # => nparam = 3
    data.spin_rbm_phys_layer_terms = [
        MVMCExpertModeParsers.SpinRBMPhysLayerTerm(0, 0.0 + 0.0im, false, 1),
    ] # => nparam = 2
    data.general_rbm_phys_layer_terms = MVMCExpertModeParsers.GeneralRBMPhysLayerTerm[]

    data.charge_rbm_hidden_layer_terms = [
        MVMCExpertModeParsers.ChargeRBMHiddenLayerTerm(0, 0.0 + 0.0im, false, 0),
    ] # => nparam = 1
    data.spin_rbm_hidden_layer_terms = MVMCExpertModeParsers.SpinRBMHiddenLayerTerm[]
    data.general_rbm_hidden_layer_terms = [
        MVMCExpertModeParsers.GeneralRBMHiddenLayerTerm(0, 0.0 + 0.0im, false, 0),
        MVMCExpertModeParsers.GeneralRBMHiddenLayerTerm(1, 0.0 + 0.0im, false, 3),
    ] # => nparam = 4

    data.charge_rbm_phys_hidden_terms = MVMCExpertModeParsers.ChargeRBMPhysHiddenTerm[]
    data.spin_rbm_phys_hidden_terms = [
        MVMCExpertModeParsers.SpinRBMPhysHiddenTerm(0, 0, 0.0 + 0.0im, false, 0),
    ] # => nparam = 1
    data.general_rbm_phys_hidden_terms = MVMCExpertModeParsers.GeneralRBMPhysHiddenTerm[]

    @test MVMCExpertModeParsers.count_rbm_parameters(data) == 11

    # Verify `set_rbm_opt_flags!` uses the same global indexing we expect from the section offsets.
    # Offsets (nparam): 3,2,0,1,0,4,0,1,0 => cumulative before spin_phys_hidden = 10.
    opt_flags = Dict{Int,Int}(
        0 => 1,   # charge_phys idx=0
        2 => 0,   # charge_phys idx=2
        10 => 1,  # spin_phys_hidden local idx=0
    )

    # No Proj terms in this fixture, so RBM block starts at fidx_offset = 0.
    MVMCExpertModeParsers.set_rbm_opt_flags!(data, opt_flags, 0; is_complex=false)

    @test data.optimization_flags[2 * 0 + 1] == true
    @test data.optimization_flags[2 * 0 + 2] == false
    @test data.optimization_flags[2 * 2 + 1] == false
    @test data.optimization_flags[2 * 2 + 2] == false
    @test data.optimization_flags[2 * 10 + 1] == true
    @test data.optimization_flags[2 * 10 + 2] == false
end
