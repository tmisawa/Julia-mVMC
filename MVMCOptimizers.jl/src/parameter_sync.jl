"""
Parameter Synchronization Functions

Synchronize modified parameters (MPI communication not needed for single process).
"""

using MVMCExpertModeParsers: is_gutzwiller_optimized, is_jastrow_optimized

const D_AMP_MAX = 4.0  # Maximum amplitude for Slater parameters (D_AmpMax)

function _real_opt_flag(data::ExpertModeData, fidx0::Int)::Bool
    flag_idx = 2 * fidx0 + 1
    flag_idx <= length(data.optimization_flags) || return false
    return data.optimization_flags[flag_idx]
end

function _all_gutzwiller_real_optimized(data::ExpertModeData, layout)::Bool
    layout.n_gutzwiller > 0 || return false
    return all(i -> _real_opt_flag(data, layout.gutzwiller_offset + i - 1), 1:layout.n_gutzwiller)
end

function flag_shift_dh2(data::ExpertModeData, layout = MVMCExpertModeParsers.projection_layout(data))::Bool
    layout.n_dh2 > 0 || return false
    _all_gutzwiller_real_optimized(data, layout) || return false
    return all(i -> _real_opt_flag(data, layout.dh2_offset + i - 1), 1:(6 * layout.n_dh2))
end

function flag_shift_dh4(data::ExpertModeData, layout = MVMCExpertModeParsers.projection_layout(data))::Bool
    layout.n_dh4 > 0 || return false
    _all_gutzwiller_real_optimized(data, layout) || return false
    return all(i -> _real_opt_flag(data, layout.dh4_offset + i - 1), 1:(10 * layout.n_dh4))
end

function shift_dh2!(data::ExpertModeData, layout = MVMCExpertModeParsers.projection_layout(data))::Float64
    layout.n_dh2 == 0 && return 0.0
    params = data.doublon_holon_2site_params
    g_shift = 0.0
    for group0 = 0:(2 * layout.n_dh2 - 1)
        i0 = group0 + 1
        i1 = group0 + 2 * layout.n_dh2 + 1
        i2 = group0 + 4 * layout.n_dh2 + 1
        if i2 <= length(params)
            shift = (real(params[i0]) + real(params[i1]) + real(params[i2])) / 3.0
            params[i0] -= shift
            params[i1] -= shift
            params[i2] -= shift
            g_shift += shift
        end
    end
    return g_shift
end

function shift_dh4!(data::ExpertModeData, layout = MVMCExpertModeParsers.projection_layout(data))::Float64
    layout.n_dh4 == 0 && return 0.0
    params = data.doublon_holon_4site_params
    g_shift = 0.0
    for group0 = 0:(2 * layout.n_dh4 - 1)
        indices = ntuple(k -> group0 + 2 * (k - 1) * layout.n_dh4 + 1, 5)
        if indices[end] <= length(params)
            shift = sum(real(params[i]) for i in indices) / 5.0
            for i in indices
                params[i] -= shift
            end
            g_shift += shift
        end
    end
    return g_shift
end

"""
    sync_modified_parameter!(data::ExpertModeData; shift_correlations::Bool = true)

Sync modified parameters.
Equivalent to C's `SyncModifiedParameter()`.

This function:
- Shifts DH/Gutzwiller/Jastrow parameters when `shift_correlations` is true
- Rescales Slater parameters to ensure maximum amplitude <= D_AMP_MAX (4.0)
- Similar to C implementation's SyncModifiedParameter() after parameter updates

# Arguments
- `data::ExpertModeData`: Expert Mode data structure to modify
- `shift_correlations::Bool`: apply DH/Gutzwiller/Jastrow gauge shifts

# Notes
- C enables DH/Gutzwiller/Jastrow shift flags only in optimization mode
- Finds maximum absolute value and scales all parameters proportionally
- Ensures max(|Slater[i]|) <= D_AMP_MAX
"""
function sync_modified_parameter!(data::ExpertModeData; shift_correlations::Bool = true)
    layout = MVMCExpertModeParsers.projection_layout(data)

    if shift_correlations
        g_shift = 0.0
        flag_shift_dh2(data, layout) && (g_shift += shift_dh2!(data, layout))
        flag_shift_dh4(data, layout) && (g_shift += shift_dh4!(data, layout))
        if g_shift != 0.0
            for term in data.gutzwiller_terms
                term.value += g_shift
            end
        end

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

    # Normalize OptTrans (C parameter.c:163-175, gated by FlagOptTrans>0).
    if !isempty(data.opt_trans)
        xmax = maximum(abs, data.opt_trans)
        if xmax > 0.0
            ratio = 1.0 / xmax
            for i in eachindex(data.opt_trans)
                data.opt_trans[i] *= ratio
            end
        end
    end

    return 0
end
