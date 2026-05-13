"""
Green Function Calculation for VMCPhysCal

Calculate and accumulate Green's functions for physical quantity measurement.
"""

"""
    initialize_phys_quantities!(state::VMCOptimizationState, data::ExpertModeData)

Initialize PhysicalQuantities in VMCOptimizationState based on ExpertModeData.
"""
function initialize_phys_quantities!(state::VMCOptimizationState, data::ExpertModeData)
    n_cis_ajs = length(data.green_one_terms)

    # Calculate NCisAjsCktAlt from green_two_terms
    # This is the number of product-type 2-body Green's functions
    # For now, we'll use green_two_terms length as a proxy
    # TODO: Parse cisajscktalt.def to get exact count
    n_cis_ajs_ckt_alt = 0  # Will be set from data if available

    n_cis_ajs_ckt_alt_dc = length(data.green_two_terms)

    state.phys_quantities =
        PhysicalQuantities(n_cis_ajs, n_cis_ajs_ckt_alt, n_cis_ajs_ckt_alt_dc)
end

"""
    reset_phys_quantities!(state::VMCOptimizationState)

Reset physical quantities for a new sampling run.
"""
function reset_phys_quantities!(state::VMCOptimizationState)
    if state.phys_quantities === nothing
        return
    end

    phys = state.phys_quantities
    fill!(phys.local_cis_ajs, 0.0 + 0.0im)
    fill!(phys.phys_cis_ajs, 0.0 + 0.0im)
    fill!(phys.phys_cis_ajs_ckt_alt, 0.0 + 0.0im)
    fill!(phys.local_cis_ajs_ckt_alt_dc, 0.0 + 0.0im)
    fill!(phys.phys_cis_ajs_ckt_alt_dc, 0.0 + 0.0im)
end

"""
    calculate_green_func!(data::ExpertModeData, state::VMCOptimizationState,
                         w::Float64, ip::ComplexF64,
                         ele_idx::Vector{Int}, ele_cfg::Vector{Int},
                         ele_num::Vector{Int}, ele_proj_cnt::Vector{Int})

Calculate Green's functions and accumulate with weight w.
Equivalent to C's `CalculateGreenFunc()` in calgrn.c.

# Arguments
- `w`: Weight for this sample
- `ip`: Inner product <ψ|x>
- `ele_idx`, `ele_cfg`, `ele_num`, `ele_proj_cnt`: Electron configuration
"""
function calculate_green_func!(
    data::ExpertModeData,
    state::VMCOptimizationState,
    w::Float64,
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
)
    if state.phys_quantities === nothing
        error("PhysicalQuantities not initialized. Call initialize_phys_quantities! first.")
    end

    phys = state.phys_quantities
    all_complex = get_all_complex_flag(data)

    # 1-body Green's function: <c†_{ri,s} c_{rj,s}>
    # C implementation (sz-conserved): uses CisAjsIdx[idx][3] = sj (spin2)
    # For sz-conserved case, spin1 == spin2, so we use spin2 to match C
    for (idx, term) in enumerate(data.green_one_terms)
        ri = term.site1
        rj = term.site2
        # Use spin2 (sj) to match C implementation
        # C: s = CisAjsIdx[idx][3] where format is (ri, si, rj, sj)
        s = term.spin2 == :up ? 0 : 1

        # Calculate GreenFunc1
        local_val = green_func1(
            ri,
            rj,
            s,
            ip,
            ele_idx,
            ele_cfg,
            ele_num,
            ele_proj_cnt,
            data,
            state;
            all_complex = all_complex,
        )

        phys.local_cis_ajs[idx] = local_val
        phys.phys_cis_ajs[idx] += w * local_val
    end

    # 2-body Green's function (direct): <c†_i c_j c†_k c_l>
    for (idx, term) in enumerate(data.green_two_terms)
        ri = term.site1
        rj = term.site2
        rk = term.site3
        rl = term.site4
        s = term.spin1 == :up ? 0 : 1
        t = term.spin3 == :up ? 0 : 1

        # Calculate GreenFunc2
        local_val = green_func2(
            ri,
            rj,
            rk,
            rl,
            s,
            t,
            ip,
            ele_idx,
            ele_cfg,
            ele_num,
            ele_proj_cnt,
            data,
            state;
            all_complex = all_complex,
        )

        phys.local_cis_ajs_ckt_alt_dc[idx] = local_val
        phys.phys_cis_ajs_ckt_alt_dc[idx] += w * local_val
    end

    # 2-body Green's function (product): <c†_i c_j> × <c†_k c_l>
    # This requires CisAjsCktAltIdx which maps to pairs of 1-body Green's function indices
    # For now, we'll skip this if not available
    # TODO: Parse cisajscktalt.def to get CisAjsCktAltIdx
    if length(phys.phys_cis_ajs_ckt_alt) > 0
        # This will be implemented when we have CisAjsCktAltIdx
        # For now, leave empty
    end
end
