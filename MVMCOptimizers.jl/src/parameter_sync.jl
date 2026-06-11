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

"""
    count_total_parameters(data) -> Int

C の NPara 相当（Proj + RBM + Slater idx + OptTrans）。
"""
function count_total_parameters(data::ExpertModeData)::Int
    return MVMCExpertModeParsers.count_variational_parameters(data)
end

"""
    set_parameter_value!(data, para_idx, value)

flat index の値を `value` に直接代入する。RBM / orbital の duplicate idx は同じ
`value` へまとめて代入するため、C の contiguous `Para` broadcast と同じく rank-local
divergence を root 値へ修復できる。`update_parameter_value` は別の delta 加算 helper
を使うが、同じ flat-index layout に従う。
"""
function set_parameter_value!(data::ExpertModeData, para_idx::Int, value::ComplexF64)
    _set_parameter_value_direct!(data, para_idx, value)
    return data
end

"""
duplicate idx invariant: 同一 idx の duplicate term が全て同値であることを確認する。
R1 の direct-assignment unpack は不一致 duplicate を修復できるため、この check は
手動検査用に残す。RBM を使わない入力では各 section は empty vector なので no-op。
"""
function _duplicate_checked_sections(data::ExpertModeData)
    return (
        (:orbital_terms, data.orbital_terms),
        (:charge_rbm_phys_layer_terms, _rbm_parameter_sections(data)[1]),
        (:spin_rbm_phys_layer_terms, _rbm_parameter_sections(data)[2]),
        (:general_rbm_phys_layer_terms, _rbm_parameter_sections(data)[3]),
        (:charge_rbm_hidden_layer_terms, _rbm_parameter_sections(data)[4]),
        (:spin_rbm_hidden_layer_terms, _rbm_parameter_sections(data)[5]),
        (:general_rbm_hidden_layer_terms, _rbm_parameter_sections(data)[6]),
        (:charge_rbm_phys_hidden_terms, _rbm_parameter_sections(data)[7]),
        (:spin_rbm_phys_hidden_terms, _rbm_parameter_sections(data)[8]),
        (:general_rbm_phys_hidden_terms, _rbm_parameter_sections(data)[9]),
    )
end

function check_duplicate_consistency(data::ExpertModeData)
    for (name, terms) in _duplicate_checked_sections(data)
        by_idx = Dict{Int,ComplexF64}()
        for term in terms
            v = get(by_idx, term.idx, nothing)
            if v === nothing
                by_idx[term.idx] = term.value
            elseif v != term.value
                error("$name idx=$(term.idx) has inconsistent duplicate values " *
                      "($v vs $(term.value)).")
            end
        end
    end
    return nothing
end

"C の contiguous `Para` 相当へ pack（spec §5-2 の pack/unpack helper）。"
function pack_parameters(data::ExpertModeData)
    counts = _parameter_count_breakdown(data)
    para = zeros(ComplexF64, counts.n_para)
    layout = counts.layout
    n_proj = counts.n_proj
    n_rbm = counts.n_rbm
    n_orbital_idx = counts.n_orbital_idx
    n_opt_trans = counts.n_opt_trans

    for i = 1:min(layout.n_gutzwiller, length(data.gutzwiller_terms))
        para[layout.gutzwiller_offset+i] = ComplexF64(data.gutzwiller_terms[i].value)
    end
    for i = 1:min(layout.n_jastrow, length(data.jastrow_terms))
        para[layout.jastrow_offset+i] = ComplexF64(data.jastrow_terms[i].value)
    end
    for i = 1:min(6 * layout.n_dh2, length(data.doublon_holon_2site_params))
        para[layout.dh2_offset+i] = ComplexF64(data.doublon_holon_2site_params[i])
    end
    for i = 1:min(10 * layout.n_dh4, length(data.doublon_holon_4site_params))
        para[layout.dh4_offset+i] = ComplexF64(data.doublon_holon_4site_params[i])
    end

    rbm_offset = n_proj
    for terms in _rbm_parameter_sections(data)
        n_section = _parameter_section_width(terms)
        for term in terms
            idx = rbm_offset + term.idx + 1
            if 1 <= idx <= length(para)
                para[idx] = ComplexF64(term.value)
            end
        end
        rbm_offset += n_section
    end

    slater_start = n_proj + n_rbm + 1
    for term in data.orbital_terms
        idx = slater_start + term.idx
        if slater_start <= idx < slater_start + n_orbital_idx
            para[idx] = ComplexF64(term.value)
        end
    end

    opt_trans_start = n_proj + n_rbm + n_orbital_idx + 1
    for i = 1:min(n_opt_trans, length(data.opt_trans))
        para[opt_trans_start+i-1] = ComplexF64(data.opt_trans[i])
    end
    return para
end

function unpack_parameters!(data::ExpertModeData, para::AbstractVector{ComplexF64})
    counts = _parameter_count_breakdown(data)
    n = counts.n_para
    length(para) == n || throw(ArgumentError(
        "parameter vector length $(length(para)) != NPara $n"))
    layout = counts.layout
    n_proj = counts.n_proj
    n_rbm = counts.n_rbm
    n_orbital_idx = counts.n_orbital_idx
    n_opt_trans = counts.n_opt_trans

    for i = 1:min(layout.n_gutzwiller, length(data.gutzwiller_terms))
        data.gutzwiller_terms[i].value = para[layout.gutzwiller_offset+i]
    end
    for i = 1:min(layout.n_jastrow, length(data.jastrow_terms))
        data.jastrow_terms[i].value = para[layout.jastrow_offset+i]
    end
    for i = 1:min(6 * layout.n_dh2, length(data.doublon_holon_2site_params))
        data.doublon_holon_2site_params[i] = para[layout.dh2_offset+i]
    end
    for i = 1:min(10 * layout.n_dh4, length(data.doublon_holon_4site_params))
        data.doublon_holon_4site_params[i] = para[layout.dh4_offset+i]
    end

    rbm_offset = n_proj
    for terms in _rbm_parameter_sections(data)
        n_section = _parameter_section_width(terms)
        for term in terms
            idx = rbm_offset + term.idx + 1
            if 1 <= idx <= length(para)
                term.value = para[idx]
            end
        end
        rbm_offset += n_section
    end

    slater_start = n_proj + n_rbm + 1
    for term in data.orbital_terms
        idx = slater_start + term.idx
        if slater_start <= idx < slater_start + n_orbital_idx
            term.value = para[idx]
        end
    end

    opt_trans_start = n_proj + n_rbm + n_orbital_idx + 1
    for i = 1:min(n_opt_trans, length(data.opt_trans))
        data.opt_trans[i] = para[opt_trans_start+i-1]
    end
    if n_opt_trans > 0 && data.qp_weights !== nothing
        MVMCExpertModeParsers.update_qp_weight!(data.qp_weights, data.opt_trans)
    end
    return data
end

"""
    sync_modified_parameter!(ctx, data; shift_correlations=true)

C `SyncModifiedParameter(comm0)`（mVMC/src/mVMC/parameter.c:134-175）の MPI-aware 版:
rank0 の parameter vector を comm0 で bcast（parameter.c:142）してから、全 rank が
同じ値に対して既存の local sync（DH/GJ shift → Slater rescale → OptTrans normalize）
を実行する。serial（`ctx.is_mpi == false`）では bcast を完全に skip し、既存実装と
bit 同一の経路になる。
"""
function sync_modified_parameter!(ctx::ParallelContext, data::ExpertModeData;
                                  shift_correlations::Bool = true)
    if ctx.is_mpi
        para = pack_parameters(data)
        bcast!(ctx, para; root = 0, which = :comm0)
        unpack_parameters!(data, para)
    end
    return sync_modified_parameter!(data; shift_correlations = shift_correlations)
end
