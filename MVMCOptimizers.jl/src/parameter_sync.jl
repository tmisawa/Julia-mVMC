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
    layout = MVMCExpertModeParsers.projection_layout(data)
    n_rbm = has_rbm_terms(data) ? MVMCExpertModeParsers.count_rbm_parameters(data) : 0
    return layout.n_proj + n_rbm +
           MVMCExpertModeParsers.count_orbital_parameters(data) +
           MVMCExpertModeParsers.count_opt_trans_parameters(data)
end

"""
    set_parameter_value!(data, para_idx, value)

flat index の値を `value` にする。書き込みは既存 `update_parameter_value`
（optimizer と同じ経路）への差分適用で行い、対応の二重実装を避ける。

**制限（plan review F7）:** `update_parameter_value` は同一 idx を持つ全
duplicate term（RBM 各 section / `orbital_terms`）へ同じ delta を加算するため、
duplicate 同士が何らかの理由で不一致になっている場合、この差分適用では root 値へ
「代入」修復できない（C の contiguous `Para` bcast より修復力が弱い）。R0 では
rank-local mutation が存在せず divergence は起きないため、delta 方式 + 下の
invariant check（fail-fast）で対応する。R1 以降で rank-local mutation が増える
場合は direct-assignment setter への refactor
（`update_parameter_value` の layout walker 共通化）を再検討する。
"""
function set_parameter_value!(data::ExpertModeData, para_idx::Int, value::ComplexF64)
    delta = value - get_parameter_value(data, para_idx)
    update_parameter_value(data, para_idx, real(delta), imag(delta))
    return data
end

"""
duplicate idx invariant（plan review F7 / addendum C1）: 同一 idx の duplicate term が
全て同値であることを確認する。不一致なら error（silent な部分修復より fail-fast を
選ぶ）。`unpack_parameters!` の冒頭で呼ぶ。

検査対象は `update_parameter_value` が「同一 idx の全 term へ同じ delta」を適用する
全 section: `orbital_terms` + RBM 9 sections（stochastic_opt.jl の `rbm_sections`
tuple と同一順）。RBM を使わない入力では各 section は empty vector なので no-op。
"""
function _duplicate_checked_sections(data::ExpertModeData)
    return (
        (:orbital_terms, data.orbital_terms),
        (:charge_rbm_phys_layer_terms, data.charge_rbm_phys_layer_terms),
        (:spin_rbm_phys_layer_terms, data.spin_rbm_phys_layer_terms),
        (:general_rbm_phys_layer_terms, data.general_rbm_phys_layer_terms),
        (:charge_rbm_hidden_layer_terms, data.charge_rbm_hidden_layer_terms),
        (:spin_rbm_hidden_layer_terms, data.spin_rbm_hidden_layer_terms),
        (:general_rbm_hidden_layer_terms, data.general_rbm_hidden_layer_terms),
        (:charge_rbm_phys_hidden_terms, data.charge_rbm_phys_hidden_terms),
        (:spin_rbm_phys_hidden_terms, data.spin_rbm_phys_hidden_terms),
        (:general_rbm_phys_hidden_terms, data.general_rbm_phys_hidden_terms),
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
                      "($v vs $(term.value)); delta-based unpack cannot repair this. " *
                      "See plan review F7 / addendum C1.")
            end
        end
    end
    return nothing
end

"C の contiguous `Para` 相当へ pack（spec §5-2 の pack/unpack helper）。"
pack_parameters(data::ExpertModeData) =
    ComplexF64[get_parameter_value(data, i) for i in 1:count_total_parameters(data)]

function unpack_parameters!(data::ExpertModeData, para::AbstractVector{ComplexF64})
    n = count_total_parameters(data)
    length(para) == n || throw(ArgumentError(
        "parameter vector length $(length(para)) != NPara $n"))
    check_duplicate_consistency(data)   # F7: 差分適用の前提（duplicate 同値）を保証
    for i in 1:n
        set_parameter_value!(data, i, para[i])
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
