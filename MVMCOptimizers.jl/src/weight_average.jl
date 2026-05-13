"""
Weighted Average Functions

Calculate weighted averages for energy and SR optimization data.
"""

"""
    weight_average_we!(state::VMCOptimizationState)

Calculate weighted average of energy.
Equivalent to C's `WeightAverageWE()`.

Normalizes Wc, Etot, Etot2, Sztot, Sztot2 by Wc.
"""
function weight_average_we!(state::VMCOptimizationState)
    if abs(state.energy.wc) < 1e-15
        @warn "Weight Wc is too small: $(state.energy.wc), etot=$(state.energy.etot)"
        return
    end

    # Debug: Check if etot is zero before normalization
    if abs(state.energy.etot) < 1e-15
        @warn "Etot is zero before normalization: wc=$(state.energy.wc), etot=$(state.energy.etot)" maxlog=3
    end

    inv_w = 1.0 / state.energy.wc
    state.energy.etot *= inv_w
    state.energy.etot2 *= inv_w
    state.energy.sztot *= inv_w
    state.energy.sztot2 *= inv_w
end

"""
    weight_average_sr_opt!(state::VMCOptimizationState)

Calculate weighted average of SR optimization data (complex version).
Equivalent to C's `WeightAverageSROpt()`.

Normalizes SROptOO and SROptHO by Wc.
Note: C implementation explicitly excludes SROptO from normalization.
"""
function weight_average_sr_opt!(state::VMCOptimizationState)
    if abs(state.energy.wc) < 1e-15
        @warn "Weight Wc is too small: $(state.energy.wc)"
        return
    end

    inv_w = 1.0 / state.energy.wc

    # Normalize SROptOO and SROptHO (not SROptO - per C implementation)
    state.sr_opt.sr_opt_oo .*= inv_w
    state.sr_opt.sr_opt_ho .*= inv_w
end

"""
    weight_average_sr_opt_real!(state::VMCOptimizationState)

Calculate weighted average of SR optimization data (real version).
Equivalent to C's `WeightAverageSROpt_real()`.

Normalizes SROptOO_real and SROptHO_real by Wc.
Note: C implementation explicitly excludes SROptO_real from normalization.
"""
function weight_average_sr_opt_real!(state::VMCOptimizationState)
    if abs(state.energy.wc) < 1e-15
        @warn "Weight Wc is too small: $(state.energy.wc)"
        return
    end

    inv_w = 1.0 / state.energy.wc

    # Normalize SROptOO_real and SROptHO_real (not SROptO_real - per C implementation)
    state.sr_opt.sr_opt_oo_real .*= inv_w
    state.sr_opt.sr_opt_ho_real .*= inv_w
end

"""
    weight_average_green_func!(state::VMCOptimizationState)

Calculate weighted average of Green's functions.
Equivalent to C's `WeightAverageGreenFunc()` in average.c.

Normalizes all Green's function arrays by Wc.
"""
function weight_average_green_func!(state::VMCOptimizationState)
    if state.phys_quantities === nothing
        return  # No physical quantities to average
    end

    if abs(state.energy.wc) < 1e-15
        @warn "Weight Wc is too small: $(state.energy.wc), cannot normalize Green's functions"
        return
    end

    phys = state.phys_quantities
    inv_w = 1.0 / state.energy.wc

    # Normalize all Green's function arrays
    phys.phys_cis_ajs .*= inv_w
    phys.phys_cis_ajs_ckt_alt .*= inv_w
    phys.phys_cis_ajs_ckt_alt_dc .*= inv_w
end
