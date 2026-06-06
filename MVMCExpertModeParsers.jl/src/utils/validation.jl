"""
Validation Utilities for Expert Mode Parsing

Parameter validation and consistency checking functions.
"""

"""
    validate_modpara_params(params::ModParaParameters) -> ValidationResult

Validate ModPara parameters for consistency.
"""
function validate_modpara_params(params::ModParaParameters)::ValidationResult
    errors = String[]
    warnings = String[]

    # Basic parameter validation
    if params.nsite <= 0
        push!(errors, "NSite must be positive")
    end

    if params.nelec < 0
        push!(errors, "NElec must be non-negative")
    end

    if params.nlocspin < 0
        push!(errors, "NLocSpin must be non-negative")
    end

    if params.ncond != -1 && params.ncond < 0
        push!(errors, "NCond must be non-negative or -1")
    end

    # Electron number consistency
    if params.ncond != -1
        if params.ncond % 2 != 0
            push!(errors, "NCond must be even")
        end
        expected_nelec = (params.nlocspin + params.ncond) ÷ 2
        if params.nelec != expected_nelec
            push!(
                warnings,
                "NElec ($(params.nelec)) differs from expected value ($expected_nelec) based on NCond",
            )
        end
    end

    # Spin consistency
    if params.two_sz != 0 && params.two_sz % 2 != 0 && params.two_sz != -1
        push!(errors, "2Sz must be even or -1")
    end

    # Local spin consistency
    if params.nlocspin > 0
        if params.nlocspin == 2 * params.nelec && params.nex_update_path != 2
            push!(errors, "NExUpdatePath must be 2 when 2*Ne = NLocalSpin (spin system)")
        elseif params.nex_update_path == 0
            push!(errors, "NExUpdatePath must be 1")
        elseif params.nlocspin > 2 * params.nelec
            push!(errors, "2*Ne must satisfy 2*Ne >= NLocalSpin")
        end
    end

    # VMC parameters
    if params.nvmc_sample <= 0
        push!(errors, "NVMCSample must be positive")
    end

    if params.nvmc_interval <= 0
        push!(errors, "NVMCInterval must be positive")
    end

    if params.nvmc_warmup < 0
        push!(errors, "NVMCWarmUp must be non-negative")
    end

    # SR parameters
    if params.nsr_opt_itr_step <= 0
        push!(errors, "NSROptItrStep must be positive")
    end

    if params.nsr_opt_itr_smp <= 0
        push!(errors, "NSROptItrSmp must be positive")
    end

    if params.dsr_opt_red_cut < 0
        push!(errors, "DSROptRedCut must be non-negative")
    end

    if params.dsr_opt_step_dt <= 0
        push!(errors, "DSROptStepDt must be positive")
    end

    # RBM parameters
    if params.nblock_size_rbm_ratio > 0 && params.nblock_size_rbm_ratio % 8 != 0
        push!(warnings, "NBlockSize_RBMRatio should be multiple of 8")
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_transfer_terms(terms::Vector{TransferTerm}, nsite::Int) -> ValidationResult

Validate transfer terms for consistency.
"""
function validate_transfer_terms(terms::Vector{TransferTerm}, nsite::Int)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site1 < 0 || term.site1 >= nsite
            push!(
                errors,
                "Transfer term $i: site1 ($(term.site1)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site2 < 0 || term.site2 >= nsite
            push!(
                errors,
                "Transfer term $i: site2 ($(term.site2)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site1 == term.site2
            push!(warnings, "Transfer term $i: site1 == site2 (diagonal term)")
        end

        if term.spin1 < 0 || term.spin1 > 1 || term.spin2 < 0 || term.spin2 > 1
            push!(errors, "Transfer term $i: invalid spin indices $(term.spin1), $(term.spin2) (must be 0 or 1)")
        elseif term.spin1 != term.spin2
            push!(warnings, "Transfer term $i: spin indices differ ($(term.spin1) != $(term.spin2))")
        end

        if !(term.spin in [:up, :down, :both])
            push!(errors, "Transfer term $i: invalid spin symbol $(term.spin)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_coulomb_intra_terms(terms::Vector{CoulombIntraTerm}, nsite::Int) -> ValidationResult

Validate Coulomb intra terms for consistency.
"""
function validate_coulomb_intra_terms(
    terms::Vector{CoulombIntraTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site < 0 || term.site >= nsite
            push!(
                errors,
                "CoulombIntra term $i: site ($(term.site)) out of range [0, $(nsite-1)]",
            )
        end

        if term.value < 0
            push!(warnings, "CoulombIntra term $i: negative value $(term.value)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_coulomb_inter_terms(terms::Vector{CoulombInterTerm}, nsite::Int) -> ValidationResult

Validate Coulomb inter terms for consistency.
"""
function validate_coulomb_inter_terms(
    terms::Vector{CoulombInterTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site1 < 0 || term.site1 >= nsite
            push!(
                errors,
                "CoulombInter term $i: site1 ($(term.site1)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site2 < 0 || term.site2 >= nsite
            push!(
                errors,
                "CoulombInter term $i: site2 ($(term.site2)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site1 == term.site2
            push!(warnings, "CoulombInter term $i: site1 == site2 (on-site term)")
        end

        if term.value < 0
            push!(warnings, "CoulombInter term $i: negative value $(term.value)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_gutzwiller_terms(terms::Vector{GutzwillerTerm}, nsite::Int) -> ValidationResult

Validate Gutzwiller terms for consistency.
"""
function validate_gutzwiller_terms(
    terms::Vector{GutzwillerTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site < 0 || term.site >= nsite
            push!(
                errors,
                "Gutzwiller term $i: site ($(term.site)) out of range [0, $(nsite-1)]",
            )
        end

        if real(term.value) < 0
            push!(warnings, "Gutzwiller term $i: negative real part $(real(term.value))")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_jastrow_terms(terms::Vector{JastrowTerm}, nsite::Int) -> ValidationResult

Validate Jastrow terms for consistency.
"""
function validate_jastrow_terms(terms::Vector{JastrowTerm}, nsite::Int)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site1 < 0 || term.site1 >= nsite
            push!(
                errors,
                "Jastrow term $i: site1 ($(term.site1)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site2 < 0 || term.site2 >= nsite
            push!(
                errors,
                "Jastrow term $i: site2 ($(term.site2)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site1 == term.site2
            push!(warnings, "Jastrow term $i: site1 == site2 (diagonal term)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_orbital_terms(terms::Vector{OrbitalTerm}, nsite::Int) -> ValidationResult

Validate orbital terms for consistency.
"""
function validate_orbital_terms(terms::Vector{OrbitalTerm}, nsite::Int)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site1 < 0 || term.site1 >= nsite
            push!(
                errors,
                "Orbital term $i: site1 ($(term.site1)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site2 < 0 || term.site2 >= nsite
            push!(
                errors,
                "Orbital term $i: site2 ($(term.site2)) out of range [0, $(nsite-1)]",
            )
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_charge_rbm_phys_layer_terms(terms::Vector{ChargeRBMPhysLayerTerm}, nsite::Int) -> ValidationResult

Validate ChargeRBM_PhysLayer terms for consistency.
"""
function validate_charge_rbm_phys_layer_terms(
    terms::Vector{ChargeRBMPhysLayerTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site < 0 || term.site >= nsite
            push!(
                errors,
                "ChargeRBM_PhysLayer term $i: site ($(term.site)) out of range [0, $(nsite-1)]",
            )
        end

        if abs(term.value) > 1e10
            push!(warnings, "ChargeRBM_PhysLayer term $i: very large value $(term.value)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_spin_rbm_phys_layer_terms(terms::Vector{SpinRBMPhysLayerTerm}, nsite::Int) -> ValidationResult

Validate SpinRBM_PhysLayer terms for consistency.
"""
function validate_spin_rbm_phys_layer_terms(
    terms::Vector{SpinRBMPhysLayerTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site < 0 || term.site >= nsite
            push!(
                errors,
                "SpinRBM_PhysLayer term $i: site ($(term.site)) out of range [0, $(nsite-1)]",
            )
        end

        if abs(term.value) > 1e10
            push!(warnings, "SpinRBM_PhysLayer term $i: very large value $(term.value)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_general_rbm_phys_layer_terms(terms::Vector{GeneralRBMPhysLayerTerm}, nsite::Int) -> ValidationResult

Validate GeneralRBM_PhysLayer terms for consistency.
"""
function validate_general_rbm_phys_layer_terms(
    terms::Vector{GeneralRBMPhysLayerTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site < 0 || term.site >= nsite
            push!(
                errors,
                "GeneralRBM_PhysLayer term $i: site ($(term.site)) out of range [0, $(nsite-1)]",
            )
        end

        if abs(term.value) > 1e10
            push!(warnings, "GeneralRBM_PhysLayer term $i: very large value $(term.value)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_charge_rbm_hidden_layer_terms(terms::Vector{ChargeRBMHiddenLayerTerm}, nsite::Int) -> ValidationResult

Validate ChargeRBM_HiddenLayer terms for consistency.
"""
function validate_charge_rbm_hidden_layer_terms(
    terms::Vector{ChargeRBMHiddenLayerTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site < 0 || term.site >= nsite
            push!(
                errors,
                "ChargeRBM_HiddenLayer term $i: site ($(term.site)) out of range [0, $(nsite-1)]",
            )
        end

        if abs(term.value) > 1e10
            push!(warnings, "ChargeRBM_HiddenLayer term $i: very large value $(term.value)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_spin_rbm_hidden_layer_terms(terms::Vector{SpinRBMHiddenLayerTerm}, nsite::Int) -> ValidationResult

Validate SpinRBM_HiddenLayer terms for consistency.
"""
function validate_spin_rbm_hidden_layer_terms(
    terms::Vector{SpinRBMHiddenLayerTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site < 0 || term.site >= nsite
            push!(
                errors,
                "SpinRBM_HiddenLayer term $i: site ($(term.site)) out of range [0, $(nsite-1)]",
            )
        end

        if abs(term.value) > 1e10
            push!(warnings, "SpinRBM_HiddenLayer term $i: very large value $(term.value)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_general_rbm_hidden_layer_terms(terms::Vector{GeneralRBMHiddenLayerTerm}, nsite::Int) -> ValidationResult

Validate GeneralRBM_HiddenLayer terms for consistency.
"""
function validate_general_rbm_hidden_layer_terms(
    terms::Vector{GeneralRBMHiddenLayerTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site < 0 || term.site >= nsite
            push!(
                errors,
                "GeneralRBM_HiddenLayer term $i: site ($(term.site)) out of range [0, $(nsite-1)]",
            )
        end

        if abs(term.value) > 1e10
            push!(
                warnings,
                "GeneralRBM_HiddenLayer term $i: very large value $(term.value)",
            )
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_charge_rbm_phys_hidden_terms(terms::Vector{ChargeRBMPhysHiddenTerm}, nsite::Int) -> ValidationResult

Validate ChargeRBM_PhysHidden terms for consistency.
"""
function validate_charge_rbm_phys_hidden_terms(
    terms::Vector{ChargeRBMPhysHiddenTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site1 < 0 || term.site1 >= nsite
            push!(
                errors,
                "ChargeRBM_PhysHidden term $i: site1 ($(term.site1)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site2 < 0 || term.site2 >= nsite
            push!(
                errors,
                "ChargeRBM_PhysHidden term $i: site2 ($(term.site2)) out of range [0, $(nsite-1)]",
            )
        end

        if abs(term.value) > 1e10
            push!(warnings, "ChargeRBM_PhysHidden term $i: very large value $(term.value)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_spin_rbm_phys_hidden_terms(terms::Vector{SpinRBMPhysHiddenTerm}, nsite::Int) -> ValidationResult

Validate SpinRBM_PhysHidden terms for consistency.
"""
function validate_spin_rbm_phys_hidden_terms(
    terms::Vector{SpinRBMPhysHiddenTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site1 < 0 || term.site1 >= nsite
            push!(
                errors,
                "SpinRBM_PhysHidden term $i: site1 ($(term.site1)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site2 < 0 || term.site2 >= nsite
            push!(
                errors,
                "SpinRBM_PhysHidden term $i: site2 ($(term.site2)) out of range [0, $(nsite-1)]",
            )
        end

        if abs(term.value) > 1e10
            push!(warnings, "SpinRBM_PhysHidden term $i: very large value $(term.value)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_general_rbm_phys_hidden_terms(terms::Vector{GeneralRBMPhysHiddenTerm}, nsite::Int) -> ValidationResult

Validate GeneralRBM_PhysHidden terms for consistency.
"""
function validate_general_rbm_phys_hidden_terms(
    terms::Vector{GeneralRBMPhysHiddenTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site1 < 0 || term.site1 >= nsite
            push!(
                errors,
                "GeneralRBM_PhysHidden term $i: site1 ($(term.site1)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site2 < 0 || term.site2 >= nsite
            push!(
                errors,
                "GeneralRBM_PhysHidden term $i: site2 ($(term.site2)) out of range [0, $(nsite-1)]",
            )
        end

        if abs(term.value) > 1e10
            push!(warnings, "GeneralRBM_PhysHidden term $i: very large value $(term.value)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_doublon_holon_2site_terms(terms::Vector{DoublonHolon2SiteTerm}, nsite::Int) -> ValidationResult

Deprecated compatibility shim for the pre-DH-1 value-bearing DH2 term API.
Runtime data validates `DoublonHolon2SiteIndex` tables instead.
"""
function validate_doublon_holon_2site_terms(
    terms::Vector{DoublonHolon2SiteTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site1 < 0 || term.site1 >= nsite
            push!(
                errors,
                "DoublonHolon2Site term $i: site1 ($(term.site1)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site2 < 0 || term.site2 >= nsite
            push!(
                errors,
                "DoublonHolon2Site term $i: site2 ($(term.site2)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site1 == term.site2
            push!(warnings, "DoublonHolon2Site term $i: site1 == site2 (diagonal term)")
        end

        if abs(term.value) > 1e10
            push!(warnings, "DoublonHolon2Site term $i: very large value $(term.value)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_doublon_holon_4site_terms(terms::Vector{DoublonHolon4SiteTerm}, nsite::Int) -> ValidationResult

Deprecated compatibility shim for the pre-DH-1 value-bearing DH4 term API.
Runtime data validates `DoublonHolon4SiteIndex` tables instead.
"""
function validate_doublon_holon_4site_terms(
    terms::Vector{DoublonHolon4SiteTerm},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (i, term) in enumerate(terms)
        if term.site1 < 0 || term.site1 >= nsite
            push!(
                errors,
                "DoublonHolon4Site term $i: site1 ($(term.site1)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site2 < 0 || term.site2 >= nsite
            push!(
                errors,
                "DoublonHolon4Site term $i: site2 ($(term.site2)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site3 < 0 || term.site3 >= nsite
            push!(
                errors,
                "DoublonHolon4Site term $i: site3 ($(term.site3)) out of range [0, $(nsite-1)]",
            )
        end

        if term.site4 < 0 || term.site4 >= nsite
            push!(
                errors,
                "DoublonHolon4Site term $i: site4 ($(term.site4)) out of range [0, $(nsite-1)]",
            )
        end

        if length(unique([term.site1, term.site2, term.site3, term.site4])) < 4
            push!(warnings, "DoublonHolon4Site term $i: duplicate sites in 4-site term")
        end

        if abs(term.value) > 1e10
            push!(warnings, "DoublonHolon4Site term $i: very large value $(term.value)")
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

function validate_doublon_holon_2site_indices(
    indices::Vector{DoublonHolon2SiteIndex},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (idx, table) in enumerate(indices)
        if size(table.neighbors, 1) != nsite || size(table.neighbors, 2) != 2
            push!(errors, "DH2 index $(idx - 1): neighbors must be $nsite x 2")
            continue
        end
        for site = 1:nsite, col = 1:2
            neighbor = table.neighbors[site, col]
            if neighbor < 0 || neighbor >= nsite
                push!(
                    errors,
                    "DH2 index $(idx - 1) site $(site - 1) neighbor $(col - 1)=$neighbor out of range [0, $(nsite - 1)]",
                )
            end
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

function validate_doublon_holon_4site_indices(
    indices::Vector{DoublonHolon4SiteIndex},
    nsite::Int,
)::ValidationResult
    errors = String[]
    warnings = String[]

    for (idx, table) in enumerate(indices)
        if size(table.neighbors, 1) != nsite || size(table.neighbors, 2) != 4
            push!(errors, "DH4 index $(idx - 1): neighbors must be $nsite x 4")
            continue
        end
        for site = 1:nsite, col = 1:4
            neighbor = table.neighbors[site, col]
            if neighbor < 0 || neighbor >= nsite
                push!(
                    errors,
                    "DH4 index $(idx - 1) site $(site - 1) neighbor $(col - 1)=$neighbor out of range [0, $(nsite - 1)]",
                )
            end
        end
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end

"""
    validate_expert_mode_data(data::ExpertModeData) -> ValidationResult

Validate all Expert Mode data for consistency.
"""
function validate_expert_mode_data(data::ExpertModeData)::ValidationResult
    errors = String[]
    warnings = String[]

    # Validate ModPara parameters
    modpara_result = validate_modpara_params(data.modpara)
    append!(errors, modpara_result.errors)
    append!(warnings, modpara_result.warnings)

    nsite = data.modpara.nsite

    # Validate transfer terms
    if !isempty(data.transfer_terms)
        transfer_result = validate_transfer_terms(data.transfer_terms, nsite)
        append!(errors, transfer_result.errors)
        append!(warnings, transfer_result.warnings)
    end

    # Validate Coulomb terms
    if !isempty(data.coulomb_intra_terms)
        coulomb_intra_result = validate_coulomb_intra_terms(data.coulomb_intra_terms, nsite)
        append!(errors, coulomb_intra_result.errors)
        append!(warnings, coulomb_intra_result.warnings)
    end

    if !isempty(data.coulomb_inter_terms)
        coulomb_inter_result = validate_coulomb_inter_terms(data.coulomb_inter_terms, nsite)
        append!(errors, coulomb_inter_result.errors)
        append!(warnings, coulomb_inter_result.warnings)
    end

    # Validate variational terms
    if !isempty(data.gutzwiller_terms)
        gutzwiller_result = validate_gutzwiller_terms(data.gutzwiller_terms, nsite)
        append!(errors, gutzwiller_result.errors)
        append!(warnings, gutzwiller_result.warnings)
    end

    if !isempty(data.jastrow_terms)
        jastrow_result = validate_jastrow_terms(data.jastrow_terms, nsite)
        append!(errors, jastrow_result.errors)
        append!(warnings, jastrow_result.warnings)
    end

    if !isempty(data.orbital_terms)
        orbital_result = validate_orbital_terms(data.orbital_terms, nsite)
        append!(errors, orbital_result.errors)
        append!(warnings, orbital_result.warnings)
    end

    # Validate RBM terms
    if !isempty(data.charge_rbm_phys_layer_terms)
        charge_rbm_phys_layer_result =
            validate_charge_rbm_phys_layer_terms(data.charge_rbm_phys_layer_terms, nsite)
        append!(errors, charge_rbm_phys_layer_result.errors)
        append!(warnings, charge_rbm_phys_layer_result.warnings)
    end

    if !isempty(data.spin_rbm_phys_layer_terms)
        spin_rbm_phys_layer_result =
            validate_spin_rbm_phys_layer_terms(data.spin_rbm_phys_layer_terms, nsite)
        append!(errors, spin_rbm_phys_layer_result.errors)
        append!(warnings, spin_rbm_phys_layer_result.warnings)
    end

    if !isempty(data.general_rbm_phys_layer_terms)
        general_rbm_phys_layer_result =
            validate_general_rbm_phys_layer_terms(data.general_rbm_phys_layer_terms, nsite)
        append!(errors, general_rbm_phys_layer_result.errors)
        append!(warnings, general_rbm_phys_layer_result.warnings)
    end

    if !isempty(data.charge_rbm_hidden_layer_terms)
        charge_rbm_hidden_layer_result = validate_charge_rbm_hidden_layer_terms(
            data.charge_rbm_hidden_layer_terms,
            nsite,
        )
        append!(errors, charge_rbm_hidden_layer_result.errors)
        append!(warnings, charge_rbm_hidden_layer_result.warnings)
    end

    if !isempty(data.spin_rbm_hidden_layer_terms)
        spin_rbm_hidden_layer_result =
            validate_spin_rbm_hidden_layer_terms(data.spin_rbm_hidden_layer_terms, nsite)
        append!(errors, spin_rbm_hidden_layer_result.errors)
        append!(warnings, spin_rbm_hidden_layer_result.warnings)
    end

    if !isempty(data.general_rbm_hidden_layer_terms)
        general_rbm_hidden_layer_result = validate_general_rbm_hidden_layer_terms(
            data.general_rbm_hidden_layer_terms,
            nsite,
        )
        append!(errors, general_rbm_hidden_layer_result.errors)
        append!(warnings, general_rbm_hidden_layer_result.warnings)
    end

    if !isempty(data.charge_rbm_phys_hidden_terms)
        charge_rbm_phys_hidden_result =
            validate_charge_rbm_phys_hidden_terms(data.charge_rbm_phys_hidden_terms, nsite)
        append!(errors, charge_rbm_phys_hidden_result.errors)
        append!(warnings, charge_rbm_phys_hidden_result.warnings)
    end

    if !isempty(data.spin_rbm_phys_hidden_terms)
        spin_rbm_phys_hidden_result =
            validate_spin_rbm_phys_hidden_terms(data.spin_rbm_phys_hidden_terms, nsite)
        append!(errors, spin_rbm_phys_hidden_result.errors)
        append!(warnings, spin_rbm_phys_hidden_result.warnings)
    end

    if !isempty(data.general_rbm_phys_hidden_terms)
        general_rbm_phys_hidden_result = validate_general_rbm_phys_hidden_terms(
            data.general_rbm_phys_hidden_terms,
            nsite,
        )
        append!(errors, general_rbm_phys_hidden_result.errors)
        append!(warnings, general_rbm_phys_hidden_result.warnings)
    end

    # Validate Doublon-Holon index tables
    if !isempty(data.doublon_holon_2site_indices)
        doublon_holon_2site_result =
            validate_doublon_holon_2site_indices(data.doublon_holon_2site_indices, nsite)
        append!(errors, doublon_holon_2site_result.errors)
        append!(warnings, doublon_holon_2site_result.warnings)
    end

    if !isempty(data.doublon_holon_4site_indices)
        doublon_holon_4site_result =
            validate_doublon_holon_4site_indices(data.doublon_holon_4site_indices, nsite)
        append!(errors, doublon_holon_4site_result.errors)
        append!(warnings, doublon_holon_4site_result.warnings)
    end

    return ValidationResult(length(errors) == 0, errors, warnings)
end
