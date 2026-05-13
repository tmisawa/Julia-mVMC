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

"""
Helpers for unit tests: minimal `ExpertModeData` fixtures.

Keep these fixtures small and deterministic: they are meant for unit/contract
tests (fast, no file I/O).
"""

function make_mock_data_for_rbm_tests(; nsite::Int = 4)
    data = ExpertModeData()
    data.modpara = ModParaParameters(
        nsite = nsite,
        nneuron_charge = 2,
        nneuron_spin = 2,
        nneuron_general = 2,
        nblock_size_rbm_ratio = 2,
    )

    # Physical-layer RBM parameters (idx layout is tested; values used by log_rbm_val/log_rbm_ratio).
    # NOTE: Terms that share the same `idx` must share the same parameter value.
    charge_phys_v = 0.31 + 0.00im
    spin_phys_v = 0.29 + 0.00im
    general_phys_v = 0.07 + 0.00im
    data.charge_rbm_phys_layer_terms = [
        ChargeRBMPhysLayerTerm(0, charge_phys_v, false, 0),
        ChargeRBMPhysLayerTerm(1, charge_phys_v, false, 0),
    ]
    data.spin_rbm_phys_layer_terms = [
        SpinRBMPhysLayerTerm(0, spin_phys_v, false, 0),
        SpinRBMPhysLayerTerm(2, spin_phys_v, false, 0),
    ]
    data.general_rbm_phys_layer_terms = [
        GeneralRBMPhysLayerTerm(0, 0, general_phys_v, false, 0),
        GeneralRBMPhysLayerTerm(3, 1, general_phys_v, false, 0),
    ]

    # Hidden-layer biases (site = neuron index).
    data.charge_rbm_hidden_layer_terms = [
        ChargeRBMHiddenLayerTerm(0, -0.20 + 0.00im, false, 0),
        ChargeRBMHiddenLayerTerm(1, 0.10 + 0.00im, false, 0),
    ]
    data.spin_rbm_hidden_layer_terms = [
        SpinRBMHiddenLayerTerm(0, 0.30 + 0.00im, false, 0),
        SpinRBMHiddenLayerTerm(1, -0.15 + 0.00im, false, 0),
    ]
    data.general_rbm_hidden_layer_terms = [
        GeneralRBMHiddenLayerTerm(0, -0.25 + 0.00im, false, 0),
        GeneralRBMHiddenLayerTerm(1, 0.05 + 0.00im, false, 0),
    ]

    # Phys-hidden couplings (site1 = physical site, site2 = neuron index).
    data.charge_rbm_phys_hidden_terms = [
        ChargeRBMPhysHiddenTerm(0, 0, 0.08 + 0.00im, false, 0),
        ChargeRBMPhysHiddenTerm(1, 1, -0.06 + 0.00im, false, 0),
    ]
    data.spin_rbm_phys_hidden_terms = [
        SpinRBMPhysHiddenTerm(2, 0, 0.09 + 0.00im, false, 0),
        SpinRBMPhysHiddenTerm(3, 1, -0.04 + 0.00im, false, 0),
    ]
    data.general_rbm_phys_hidden_terms = [
        GeneralRBMPhysHiddenTerm(0, 0, 0, 0.03 + 0.00im, false, 0),
        GeneralRBMPhysHiddenTerm(3, 1, 1, -0.02 + 0.00im, false, 0),
    ]

    return data
end

function make_minimal_data_for_rbm_diff_tests(; nsite::Int = 1, nneuron_charge::Int = 1)
    data = ExpertModeData()
    data.modpara = ModParaParameters(nsite = nsite, nneuron_charge = nneuron_charge)
    return data
end

function make_mock_data_for_proj_tests(; nsite::Int = 4)
    data = ExpertModeData()
    data.modpara = ModParaParameters(nsite = nsite)

    # Gutzwiller index table: site -> idx (0-based)
    data.n_gutzwiller_idx = 2
    data.gutzwiller_idx = [0, 1, 0, 1][1:nsite]

    # Jastrow index table: site-pair -> idx (0-based)
    # For nsite=4, pairs = 6 => indices 0..5.
    n_j = nsite * (nsite - 1) ÷ 2
    data.n_jastrow_idx = n_j
    jmat = fill(-1, nsite, nsite)
    idx = 0
    for i = 1:nsite
        for j = (i + 1):nsite
            jmat[i, j] = idx
            jmat[j, i] = idx
            idx += 1
        end
    end
    data.jastrow_idx = jmat

    return data
end

