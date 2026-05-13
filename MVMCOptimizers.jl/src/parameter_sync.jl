"""
Parameter Synchronization Functions

Synchronize modified parameters (MPI communication not needed for single process).
"""

using MVMCExpertModeParsers: is_gutzwiller_optimized, is_jastrow_optimized

const D_AMP_MAX = 4.0  # Maximum amplitude for Slater parameters (D_AmpMax)

"""
    sync_modified_parameter!(data::ExpertModeData)

Sync modified parameters.
Equivalent to C's `SyncModifiedParameter()`.

This function:
- Rescales Slater parameters to ensure maximum amplitude <= D_AMP_MAX (4.0)
- Similar to C implementation's SyncModifiedParameter() after parameter updates

# Arguments
- `data::ExpertModeData`: Expert Mode data structure to modify

# Notes
- Only rescales Slater (Orbital) parameters
- Finds maximum absolute value and scales all parameters proportionally
- Ensures max(|Slater[i]|) <= D_AMP_MAX
"""
function sync_modified_parameter!(data::ExpertModeData)
    # Shift Gutzwiller/Jastrow if all are optimized (C: FlagShiftGJ / shiftGJ)
    n_gutz = length(data.gutzwiller_terms)
    n_jast = length(data.jastrow_terms)
    flag_shift_gj = false
    if n_gutz > 0 && n_jast > 0
        if isempty(data.optimization_flags)
            # Default: treat all parameters as optimized
            flag_shift_gj = true
        else
            all_gutz = all(i -> is_gutzwiller_optimized(data, i), 1:n_gutz)
            all_jast = all(i -> is_jastrow_optimized(data, i), 1:n_jast)
            flag_shift_gj = all_gutz && all_jast
        end
    end
    if flag_shift_gj
        shift = 0.0
        n_all = n_gutz + n_jast
        for term in data.gutzwiller_terms
            shift += real(term.value)
        end
        for term in data.jastrow_terms
            shift += real(term.value)
        end
        if n_all > 0
            shift /= n_all
            for term in data.gutzwiller_terms
                term.value -= shift
            end
            for term in data.jastrow_terms
                term.value -= shift
            end
        end
    end

    # Rescale Slater parameters (C: D_AmpMax/xmax)
    xmax = 0.0
    for term in data.orbital_terms
        abs_val = abs(term.value)
        if abs_val > xmax
            xmax = abs_val
        end
    end
    if xmax > 0.0
        ratio = D_AMP_MAX / xmax
        for term in data.orbital_terms
            term.value *= ratio
        end
    end

    # OptTrans normalization (C parameter.c:163-175, gated by FlagOptTrans>0)
    # is intentionally NOT performed here. The C code rescales the dedicated
    # OptTrans[] array, which has no Julia counterpart yet. The previous
    # implementation rescaled `data.para_qp_trans` instead, but that array
    # holds the QPTrans phase factors consumed by qp_weight_update.jl — a
    # different quantity. Rescaling it would silently corrupt the quantum
    # projection weights. When OptTrans is added to Julia, normalize that
    # field here (gated by an analogue of FlagOptTrans), not para_qp_trans.

    return 0
end
