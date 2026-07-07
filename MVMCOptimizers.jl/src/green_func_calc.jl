"""
Green Function Calculation for VMCPhysCal

Calculate and accumulate Green's functions for physical quantity measurement.
"""

# Symbol spin (:up/:down) -> integer (0/1), matching C's spin encoding.
# Reject anything else (e.g. :both or a typo) instead of silently mapping it to
# down spin, so malformed data surfaces.
function _spin_int(s::Symbol)
    if s == :up
        return 0
    elseif s == :down
        return 1
    end
    error("one-body Green spin must be :up or :down, got :$s")
end

"""
    build_canonical_cis_ajs_idx(green_one_terms, green_two_ex_terms, n_site)
        -> Vector{NTuple{4,Int}}

Build the C-compatible canonical one-body Green list of `(ri, si, rj, sj)`.

- When `green_two_ex_terms` is empty, this is exactly `greenone.def` order with
  no de-duplication (C reads cisajs.def directly, `IndirectGFOn = false`).
- When `green_two_ex_terms` is non-empty, the list is de-duplicated by
  `(ri, si, rj, sj)`: explicit `greenone.def` terms first, then for each
  factored term its first constituent `(i1,j1)` and second constituent
  `(i2,j2)` appended if not already present (C's `CountOneBodyGForLanczos` /
  `iOneBodyGIdx`).

Sites are validated against `[0, n_site)`.
"""
function build_canonical_cis_ajs_idx(green_one_terms, green_two_ex_terms, n_site::Int)
    canonical = NTuple{4,Int}[]
    if isempty(green_two_ex_terms)
        for t in green_one_terms
            push!(canonical, (t.site1, _spin_int(t.spin1), t.site2, _spin_int(t.spin2)))
        end
    else
        seen = Dict{NTuple{4,Int},Int}()
        function add!(key::NTuple{4,Int})
            if !haskey(seen, key)
                push!(canonical, key)
                seen[key] = length(canonical)  # 1-based
            end
        end
        for t in green_one_terms
            add!((t.site1, _spin_int(t.spin1), t.site2, _spin_int(t.spin2)))
        end
        for t in green_two_ex_terms
            add!((t.site_i1, t.spin_i1, t.site_j1, t.spin_j1))
            add!((t.site_i2, t.spin_i2, t.site_j2, t.spin_j2))
        end
    end
    for (ri, si, rj, sj) in canonical
        if ri < 0 || ri >= n_site || rj < 0 || rj >= n_site
            error("one-body Green site out of range [0,$n_site): (ri=$ri, rj=$rj)")
        end
        if si < 0 || si > 1 || sj < 0 || sj > 1
            error("one-body Green spin must be 0 or 1: (si=$si, sj=$sj)")
        end
    end
    return canonical
end

"""
    resolve_cis_ajs_ckt_alt_idx(cis_ajs_idx, green_two_ex_terms)
        -> Vector{Tuple{Int,Int}}

Resolve each factored term to a 1-based `(idx0, idx1)` pair into `cis_ajs_idx`.
`idx0` is the first one-body constituent `(i1,j1)`, `idx1` the second `(i2,j2)`.
After `build_canonical_cis_ajs_idx` every constituent is present, so a miss is
an internal invariant error.
"""
function resolve_cis_ajs_ckt_alt_idx(cis_ajs_idx, green_two_ex_terms)
    lookup = Dict{NTuple{4,Int},Int}()
    for (i, key) in enumerate(cis_ajs_idx)
        # First occurrence wins (canonical is already de-duplicated).
        get!(lookup, key, i)
    end
    pairs = Tuple{Int,Int}[]
    for t in green_two_ex_terms
        k0 = (t.site_i1, t.spin_i1, t.site_j1, t.spin_j1)
        k1 = (t.site_i2, t.spin_i2, t.site_j2, t.spin_j2)
        idx0 = get(lookup, k0, 0)
        idx1 = get(lookup, k1, 0)
        if idx0 == 0 || idx1 == 0
            error("internal: factored constituent missing from canonical one-body list")
        end
        push!(pairs, (idx0, idx1))
    end
    return pairs
end

function validate_lanczos_mode2_cis_ajs_unique(cis_ajs_idx)
    seen = Set{NTuple{4,Int}}()
    for key in cis_ajs_idx
        if key in seen
            error(
                "NLanczosMode = 2 does not support duplicate OneBodyG entries " *
                "without TwoBodyGEx in Julia-mVMC; duplicate entry = $key. " *
                "Remove duplicate greenone.def rows before running mode2.",
            )
        end
        push!(seen, key)
    end
    return nothing
end

"""
    validate_factored_green_supported(data::ExpertModeData)

Reject the factored two-body Green (`TwoBodyGEx`) in FSZ / general-orbital mode.
The FSZ measurement path is not yet wired for Green functions (a separate spec),
so producing factored output there would be silently wrong. Non-mutating;
returns `nothing` when supported.
"""
function validate_factored_green_supported(data::ExpertModeData)
    if !isempty(data.green_two_ex_terms) && data.i_flg_orbital_general != 0
        error(
            "TwoBodyGEx (factored two-body Green) is not supported in FSZ / " *
            "general-orbital mode (i_flg_orbital_general = $(data.i_flg_orbital_general)). " *
            "The FSZ Green measurement path is not yet implemented; use sz-conserved " *
            "mode or remove the TwoBodyGEx input.",
        )
    end
    return nothing
end

"""
    initialize_phys_quantities!(state::VMCOptimizationState, data::ExpertModeData)

Initialize PhysicalQuantities, including the canonical one-body list and the
resolved factored index pairs.
"""
function initialize_phys_quantities!(state::VMCOptimizationState, data::ExpertModeData)
    n_site = data.modpara.nsite
    cis_ajs_idx =
        build_canonical_cis_ajs_idx(data.green_one_terms, data.green_two_ex_terms, n_site)
    cis_ajs_ckt_alt_idx = resolve_cis_ajs_ckt_alt_idx(cis_ajs_idx, data.green_two_ex_terms)
    if data.modpara.lanczos_mode > 1 && isempty(data.green_two_ex_terms)
        validate_lanczos_mode2_cis_ajs_unique(cis_ajs_idx)
    end

    n_cis_ajs = length(cis_ajs_idx)
    n_cis_ajs_ckt_alt = length(cis_ajs_ckt_alt_idx)  # == length(green_two_ex_terms)
    n_cis_ajs_ckt_alt_dc = length(data.green_two_terms)

    phys = PhysicalQuantities(n_cis_ajs, n_cis_ajs_ckt_alt, n_cis_ajs_ckt_alt_dc)
    phys.cis_ajs_idx = cis_ajs_idx
    phys.cis_ajs_ckt_alt_idx = cis_ajs_ckt_alt_idx
    state.phys_quantities = phys
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
    fill!(phys.phys_lanczos_qqqq, 0.0 + 0.0im)
    fill!(phys.phys_lanczos_qcisajsq, 0.0 + 0.0im)
    fill!(phys.phys_lanczos_qcisajscktaltq, 0.0 + 0.0im)
    fill!(phys.phys_lanczos_qcisajscktaltq_dc, 0.0 + 0.0im)
end

"""
    accumulate_factored_green!(phys::PhysicalQuantities, w::Float64)

Accumulate the factored two-body Green:
`phys_cis_ajs_ckt_alt[idx] += w * local_cis_ajs[idx0] * conj(local_cis_ajs[idx1])`,
faithful to C's `calgrn.c` `PhysCisAjsCktAlt` loop. Indices are 1-based.
"""
function accumulate_factored_green!(phys::PhysicalQuantities, w::Float64)
    @inbounds for (idx, (idx0, idx1)) in enumerate(phys.cis_ajs_ckt_alt_idx)
        phys.phys_cis_ajs_ckt_alt[idx] +=
            w * phys.local_cis_ajs[idx0] * conj(phys.local_cis_ajs[idx1])
    end
    return nothing
end

function accumulate_factored_green!(
    dst::PhysicalQuantities,
    phys::PhysicalQuantities,
    w::Float64,
)
    @inbounds for (idx, (idx0, idx1)) in enumerate(phys.cis_ajs_ckt_alt_idx)
        dst.phys_cis_ajs_ckt_alt[idx] +=
            w * phys.local_cis_ajs[idx0] * conj(phys.local_cis_ajs[idx1])
    end
    return nothing
end

function accumulate_factored_green!(
    acc::VMCPhysAccumulator,
    phys::PhysicalQuantities,
    w::Float64,
)
    @inbounds for (idx, (idx0, idx1)) in enumerate(phys.cis_ajs_ckt_alt_idx)
        acc.phys_cis_ajs_ckt_alt[idx] +=
            w * acc.local_cis_ajs[idx0] * conj(acc.local_cis_ajs[idx1])
    end
    return nothing
end

@inline function accumulate_phys_cis_ajs!(
    phys::PhysicalQuantities,
    idx::Integer,
    w::Float64,
    local_val::ComplexF64,
)
    phys.local_cis_ajs[idx] = local_val
    phys.phys_cis_ajs[idx] += w * local_val
    return nothing
end

@inline function accumulate_phys_cis_ajs!(
    dst::PhysicalQuantities,
    phys::PhysicalQuantities,
    idx::Integer,
    w::Float64,
    local_val::ComplexF64,
)
    phys.local_cis_ajs[idx] = local_val
    dst.phys_cis_ajs[idx] += w * local_val
    return nothing
end

@inline function accumulate_phys_cis_ajs!(
    acc::VMCPhysAccumulator,
    phys::PhysicalQuantities,
    idx::Integer,
    w::Float64,
    local_val::ComplexF64,
)
    acc.local_cis_ajs[idx] = local_val
    acc.phys_cis_ajs[idx] += w * local_val
    return nothing
end

@inline function accumulate_phys_cis_ajs_ckt_alt_dc!(
    phys::PhysicalQuantities,
    idx::Integer,
    w::Float64,
    local_val::ComplexF64,
)
    phys.local_cis_ajs_ckt_alt_dc[idx] = local_val
    phys.phys_cis_ajs_ckt_alt_dc[idx] += w * local_val
    return nothing
end

@inline function accumulate_phys_cis_ajs_ckt_alt_dc!(
    dst::PhysicalQuantities,
    phys::PhysicalQuantities,
    idx::Integer,
    w::Float64,
    local_val::ComplexF64,
)
    phys.local_cis_ajs_ckt_alt_dc[idx] = local_val
    dst.phys_cis_ajs_ckt_alt_dc[idx] += w * local_val
    return nothing
end

@inline function accumulate_phys_cis_ajs_ckt_alt_dc!(
    acc::VMCPhysAccumulator,
    phys::PhysicalQuantities,
    idx::Integer,
    w::Float64,
    local_val::ComplexF64,
)
    acc.local_cis_ajs_ckt_alt_dc[idx] = local_val
    acc.phys_cis_ajs_ckt_alt_dc[idx] += w * local_val
    return nothing
end

const LANCZOS_N_LS_HAM = 2
const LANCZOS_QQQQ_LENGTH = LANCZOS_N_LS_HAM^4

struct LanczosElectronConfig
    ele_idx::Vector{Int}
    ele_cfg::Vector{Int}
    ele_num::Vector{Int}
    ele_proj_cnt::Vector{Int}
end

function _lanczos_config(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
)
    return LanczosElectronConfig(copy(ele_idx), copy(ele_cfg), copy(ele_num), copy(ele_proj_cnt))
end

_copy_lanczos_config(cfg::LanczosElectronConfig) =
    LanczosElectronConfig(copy(cfg.ele_idx), copy(cfg.ele_cfg), copy(cfg.ele_num), copy(cfg.ele_proj_cnt))

function _lanczos_same_config(a::LanczosElectronConfig, b::LanczosElectronConfig)
    return a.ele_idx == b.ele_idx &&
           a.ele_cfg == b.ele_cfg &&
           a.ele_num == b.ele_num &&
           a.ele_proj_cnt == b.ele_proj_cnt
end

function _lanczos_recompute_proj_cnt!(cfg::LanczosElectronConfig, data::ExpertModeData)
    make_proj_cnt!(cfg.ele_proj_cnt, cfg.ele_num, data)
    return cfg
end

@inline _lanczos_value(value, all_complex::Bool) =
    all_complex ? ComplexF64(value) : ComplexF64(real(value), 0.0)

function _lanczos_hamiltonian0(data::ExpertModeData, cfg::LanczosElectronConfig)::ComplexF64
    n_site = data.modpara.nsite
    n0 = @view cfg.ele_num[1:n_site]
    n1 = @view cfg.ele_num[(n_site+1):(2*n_site)]

    e = 0.0 + 0.0im
    for term in data.coulomb_intra_terms
        ri = term.site
        if 0 <= ri < n_site
            e += term.value * n0[ri+1] * n1[ri+1]
        end
    end
    for term in data.coulomb_inter_terms
        ri = term.site1
        rj = term.site2
        if 0 <= ri < n_site && 0 <= rj < n_site
            e += term.value * (n0[ri+1] + n1[ri+1]) * (n0[rj+1] + n1[rj+1])
        end
    end
    for term in data.hund_terms
        ri = term.site1
        rj = term.site2
        if 0 <= ri < n_site && 0 <= rj < n_site
            e -= term.value * (n0[ri+1] * n0[rj+1] + n1[ri+1] * n1[rj+1])
        end
    end
    return e
end

function _lanczos_preserve_m_all!(f, state::VMCOptimizationState)
    sm = state.slater_matrix
    inv_m = copy(sm.inv_m)
    pf_m = copy(sm.pf_m)
    inv_m_real = copy(sm.inv_m_real)
    pf_m_real = copy(sm.pf_m_real)
    try
        return f()
    finally
        copyto!(sm.inv_m, inv_m)
        copyto!(sm.pf_m, pf_m)
        copyto!(sm.inv_m_real, inv_m_real)
        copyto!(sm.pf_m_real, pf_m_real)
    end
end

function _lanczos_overlap_ratio!(
    state::VMCOptimizationState,
    cfg::LanczosElectronConfig,
    original::LanczosElectronConfig,
    ip::ComplexF64,
    data::ExpertModeData;
    all_complex::Bool,
)::ComplexF64
    abs(ip) > 0.0 || error("Lanczos overlap ratio cannot be computed for zero inner product")
    n_qp_full = get_n_qp_full(data)
    z = proj_rbm_ratio(cfg.ele_proj_cnt, original.ele_proj_cnt, cfg.ele_num, original.ele_num, data)

    ip_new = if all_complex
        info = calculate_m_all_fcmp!(cfg.ele_idx, 1, n_qp_full + 1, data, state; threaded = false)
        info == 0 || return 0.0 + 0.0im
        calculate_ip_fcmp(state.slater_matrix.pf_m, 1, n_qp_full + 1, data; reduce = :none)
    else
        info = calculate_m_all_real!(cfg.ele_idx, 1, n_qp_full + 1, data, state; threaded = false)
        info == 0 || return 0.0 + 0.0im
        ComplexF64(
            calculate_ip_real(
                state.slater_matrix.pf_m_real,
                1,
                n_qp_full + 1,
                data;
                reduce = :none,
            ),
            0.0,
        )
    end

    return conj(z * ip_new / ip)
end

function _lanczos_scale_terms(terms, scale::ComplexF64)
    scaled = Vector{Tuple{ComplexF64,LanczosElectronConfig}}()
    for (coef, cfg) in terms
        new_coef = scale * coef
        if abs(new_coef) > 0
            push!(scaled, (new_coef, cfg))
        end
    end
    return scaled
end

function _lanczos_apply_one_body(
    cfg::LanczosElectronConfig,
    ri::Int,
    rj::Int,
    s::Int,
    data::ExpertModeData,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    rsi = ri + s * n_site
    rsj = rj + s * n_site

    if ri == rj
        if cfg.ele_num[rsi+1] == 1
            return [(1.0 + 0.0im, _copy_lanczos_config(cfg))]
        end
        return Tuple{ComplexF64,LanczosElectronConfig}[]
    end
    if cfg.ele_num[rsi+1] == 1 || cfg.ele_num[rsj+1] == 0
        return Tuple{ComplexF64,LanczosElectronConfig}[]
    end

    mj = cfg.ele_cfg[rsj+1]
    mj < 0 && return Tuple{ComplexF64,LanczosElectronConfig}[]
    moved = _copy_lanczos_config(cfg)
    moved.ele_idx[mj+s*n_elec+1] = ri
    moved.ele_cfg[rsj+1] = -1
    moved.ele_cfg[rsi+1] = mj
    moved.ele_num[rsj+1] = 0
    moved.ele_num[rsi+1] = 1
    _lanczos_recompute_proj_cnt!(moved, data)
    return [(1.0 + 0.0im, moved)]
end

function _lanczos_diagonal_term(coef::Integer, cfg::LanczosElectronConfig)
    coef == 0 && return Tuple{ComplexF64,LanczosElectronConfig}[]
    return [(ComplexF64(coef, 0.0), _copy_lanczos_config(cfg))]
end

function _lanczos_apply_two_body(
    cfg::LanczosElectronConfig,
    ri::Int,
    rj::Int,
    rk::Int,
    rl::Int,
    s::Int,
    t::Int,
    data::ExpertModeData,
)
    n_site = data.modpara.nsite
    rsi = ri + s * n_site
    rsj = rj + s * n_site
    rtk = rk + t * n_site
    rtl = rl + t * n_site

    if s == t
        if rk == rl
            cfg.ele_num[rtk+1] == 0 && return Tuple{ComplexF64,LanczosElectronConfig}[]
            return _lanczos_apply_one_body(cfg, ri, rj, s, data)
        elseif rj == rl
            return Tuple{ComplexF64,LanczosElectronConfig}[]
        elseif ri == rl
            if cfg.ele_num[rsi+1] == 0
                return Tuple{ComplexF64,LanczosElectronConfig}[]
            elseif rj == rk
                return _lanczos_diagonal_term(1 - cfg.ele_num[rsj+1], cfg)
            else
                return _lanczos_scale_terms(_lanczos_apply_one_body(cfg, rk, rj, s, data), -1.0 + 0.0im)
            end
        elseif rj == rk
            if cfg.ele_num[rsj+1] == 1
                return Tuple{ComplexF64,LanczosElectronConfig}[]
            elseif ri == rl
                return _lanczos_diagonal_term(cfg.ele_num[rsi+1], cfg)
            else
                return _lanczos_apply_one_body(cfg, ri, rl, s, data)
            end
        elseif ri == rk
            return Tuple{ComplexF64,LanczosElectronConfig}[]
        elseif ri == rj
            cfg.ele_num[rsi+1] == 0 && return Tuple{ComplexF64,LanczosElectronConfig}[]
            return _lanczos_apply_one_body(cfg, rk, rl, s, data)
        end
    else
        if rk == rl
            if cfg.ele_num[rtk+1] == 0
                return Tuple{ComplexF64,LanczosElectronConfig}[]
            elseif ri == rj
                return _lanczos_diagonal_term(cfg.ele_num[rsi+1], cfg)
            else
                return _lanczos_apply_one_body(cfg, ri, rj, s, data)
            end
        elseif ri == rj
            cfg.ele_num[rsi+1] == 0 && return Tuple{ComplexF64,LanczosElectronConfig}[]
            return _lanczos_apply_one_body(cfg, rk, rl, t, data)
        end
    end

    if cfg.ele_num[rsi+1] == 1 ||
       cfg.ele_num[rsj+1] == 0 ||
       cfg.ele_num[rtk+1] == 1 ||
       cfg.ele_num[rtl+1] == 0
        return Tuple{ComplexF64,LanczosElectronConfig}[]
    end

    results = Vector{Tuple{ComplexF64,LanczosElectronConfig}}()
    for (coef1, cfg1) in _lanczos_apply_one_body(cfg, rk, rl, t, data)
        for (coef2, cfg2) in _lanczos_apply_one_body(cfg1, ri, rj, s, data)
            push!(results, (coef1 * coef2, cfg2))
        end
    end
    return results
end

function _lanczos_transfer_spins(term)
    term.spin == :both && return (0, 1)
    return (term.spin1,)
end

function _lanczos_local_hamiltonian_from_config!(
    cfg::LanczosElectronConfig,
    original::LanczosElectronConfig,
    original_h1::ComplexF64,
    ip::ComplexF64,
    data::ExpertModeData,
    state::VMCOptimizationState;
    all_complex::Bool,
)::ComplexF64
    if _lanczos_same_config(cfg, original)
        return original_h1
    end

    value =
        _lanczos_hamiltonian0(data, cfg) *
        _lanczos_overlap_ratio!(state, cfg, original, ip, data; all_complex = all_complex)

    for term in data.transfer_terms
        coeff = -_lanczos_value(term.value, all_complex)
        for s in _lanczos_transfer_spins(term)
            for (op_coef, moved) in _lanczos_apply_one_body(cfg, term.site1, term.site2, s, data)
                value += coeff * op_coef *
                         _lanczos_overlap_ratio!(
                             state,
                             moved,
                             original,
                             ip,
                             data;
                             all_complex = all_complex,
                         )
            end
        end
    end

    for term in data.pair_hop_terms
        coeff = _lanczos_value(term.value, all_complex)
        for (op_coef, moved) in
            _lanczos_apply_two_body(cfg, term.site1, term.site2, term.site1, term.site2, 0, 1, data)
            value += coeff * op_coef *
                     _lanczos_overlap_ratio!(
                         state,
                         moved,
                         original,
                         ip,
                         data;
                         all_complex = all_complex,
                     )
        end
    end

    for term in data.exchange_terms
        coeff = _lanczos_value(term.value, all_complex)
        for (op_coef, moved) in
            _lanczos_apply_two_body(cfg, term.site1, term.site2, term.site2, term.site1, 0, 1, data)
            value += coeff * op_coef *
                     _lanczos_overlap_ratio!(
                         state,
                         moved,
                         original,
                         ip,
                         data;
                         all_complex = all_complex,
                     )
        end
        for (op_coef, moved) in
            _lanczos_apply_two_body(cfg, term.site1, term.site2, term.site2, term.site1, 1, 0, data)
            value += coeff * op_coef *
                     _lanczos_overlap_ratio!(
                         state,
                         moved,
                         original,
                         ip,
                         data;
                         all_complex = all_complex,
                     )
        end
    end

    for term in data.inter_all_terms
        coeff = _lanczos_value(term.value, all_complex)
        for (op_coef, moved) in _lanczos_apply_two_body(
            cfg,
            term.site0,
            term.site1,
            term.site2,
            term.site3,
            term.spin1,
            term.spin3,
            data,
        )
            value += coeff * op_coef *
                     _lanczos_overlap_ratio!(
                         state,
                         moved,
                         original,
                         ip,
                         data;
                         all_complex = all_complex,
                     )
        end
    end

    return value
end

function calculate_lanczos_h2!(
    h1::ComplexF64,
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState;
    all_complex::Bool,
)::ComplexF64
    original = _lanczos_config(ele_idx, ele_cfg, ele_num, ele_proj_cnt)

    return _lanczos_preserve_m_all!(state) do
        h2 = h1 * _lanczos_hamiltonian0(data, original)

        for term in data.transfer_terms
            coeff = -_lanczos_value(term.value, all_complex)
            for s in _lanczos_transfer_spins(term)
                for (op_coef, moved) in
                    _lanczos_apply_one_body(original, term.site1, term.site2, s, data)
                    h2 += coeff * op_coef *
                          _lanczos_local_hamiltonian_from_config!(
                              moved,
                              original,
                              h1,
                              ip,
                              data,
                              state;
                              all_complex = all_complex,
                          )
                end
            end
        end

        for term in data.pair_hop_terms
            coeff = _lanczos_value(term.value, all_complex)
            for (op_coef, moved) in _lanczos_apply_two_body(
                original,
                term.site1,
                term.site2,
                term.site1,
                term.site2,
                0,
                1,
                data,
            )
                h2 += coeff * op_coef *
                      _lanczos_local_hamiltonian_from_config!(
                          moved,
                          original,
                          h1,
                          ip,
                          data,
                          state;
                          all_complex = all_complex,
                      )
            end
        end

        for term in data.exchange_terms
            coeff = _lanczos_value(term.value, all_complex)
            for (op_coef, moved) in _lanczos_apply_two_body(
                original,
                term.site1,
                term.site2,
                term.site2,
                term.site1,
                0,
                1,
                data,
            )
                h2 += coeff * op_coef *
                      _lanczos_local_hamiltonian_from_config!(
                          moved,
                          original,
                          h1,
                          ip,
                          data,
                          state;
                          all_complex = all_complex,
                      )
            end
            for (op_coef, moved) in _lanczos_apply_two_body(
                original,
                term.site1,
                term.site2,
                term.site2,
                term.site1,
                1,
                0,
                data,
            )
                h2 += coeff * op_coef *
                      _lanczos_local_hamiltonian_from_config!(
                          moved,
                          original,
                          h1,
                          ip,
                          data,
                          state;
                          all_complex = all_complex,
                      )
            end
        end

        for term in data.inter_all_terms
            coeff = _lanczos_value(term.value, all_complex)
            for (op_coef, moved) in _lanczos_apply_two_body(
                original,
                term.site0,
                term.site1,
                term.site2,
                term.site3,
                term.spin1,
                term.spin3,
                data,
            )
                h2 += coeff * op_coef *
                      _lanczos_local_hamiltonian_from_config!(
                          moved,
                          original,
                          h1,
                          ip,
                          data,
                          state;
                          all_complex = all_complex,
                      )
            end
        end

        h2
    end
end

function accumulate_lanczos_qqqq!(
    dst::Union{PhysicalQuantities,VMCPhysAccumulator},
    w::Float64,
    h1::ComplexF64,
    h2::ComplexF64;
    all_complex::Bool,
)
    length(dst.phys_lanczos_qqqq) == LANCZOS_QQQQ_LENGTH ||
        throw(
            ArgumentError(
                "phys_lanczos_qqqq must have length $LANCZOS_QQQQ_LENGTH; " *
                "got $(length(dst.phys_lanczos_qqqq)).",
            ),
        )

    lslq = (1.0 + 0.0im, h1, h1, h2)
    @inbounds for i0 = 0:(LANCZOS_QQQQ_LENGTH-1)
        rj = i0 % LANCZOS_N_LS_HAM
        ri = (i0 ÷ LANCZOS_N_LS_HAM) % LANCZOS_N_LS_HAM
        rp = (i0 ÷ (LANCZOS_N_LS_HAM^2)) % LANCZOS_N_LS_HAM
        rq = (i0 ÷ (LANCZOS_N_LS_HAM^3)) % LANCZOS_N_LS_HAM
        left = lslq[rq*LANCZOS_N_LS_HAM+ri+1]
        right = lslq[rp*LANCZOS_N_LS_HAM+rj+1]
        dst.phys_lanczos_qqqq[i0+1] += w * (all_complex ? conj(left) : left) * right
    end
    return nothing
end

@inline function _lanczos_local_value(value::ComplexF64, all_complex::Bool)
    return all_complex ? value : ComplexF64(real(value), 0.0)
end

@inline _local_cis_ajs(dst::PhysicalQuantities, phys::PhysicalQuantities) = phys.local_cis_ajs
@inline _local_cis_ajs(dst::VMCPhysAccumulator, phys::PhysicalQuantities) = dst.local_cis_ajs

@inline _local_cis_ajs_ckt_alt_dc(dst::PhysicalQuantities, phys::PhysicalQuantities) =
    phys.local_cis_ajs_ckt_alt_dc
@inline _local_cis_ajs_ckt_alt_dc(dst::VMCPhysAccumulator, phys::PhysicalQuantities) =
    dst.local_cis_ajs_ckt_alt_dc

function _lanczos_hca_from_config!(
    ri::Int,
    rj::Int,
    s::Int,
    h1::ComplexF64,
    ip::ComplexF64,
    original::LanczosElectronConfig,
    data::ExpertModeData,
    state::VMCOptimizationState;
    all_complex::Bool,
)::ComplexF64
    value = 0.0 + 0.0im
    for (op_coef, moved) in _lanczos_apply_one_body(original, ri, rj, s, data)
        value += op_coef *
                 _lanczos_local_hamiltonian_from_config!(
                     moved,
                     original,
                     h1,
                     ip,
                     data,
                     state;
                     all_complex = all_complex,
                 )
    end
    return _lanczos_local_value(value, all_complex)
end

function _lanczos_hcaca_from_config!(
    ri::Int,
    rj::Int,
    rk::Int,
    rl::Int,
    s::Int,
    t::Int,
    h1::ComplexF64,
    ip::ComplexF64,
    original::LanczosElectronConfig,
    data::ExpertModeData,
    state::VMCOptimizationState;
    all_complex::Bool,
)::ComplexF64
    value = 0.0 + 0.0im
    for (op_coef, moved) in _lanczos_apply_two_body(original, ri, rj, rk, rl, s, t, data)
        value += op_coef *
                 _lanczos_local_hamiltonian_from_config!(
                     moved,
                     original,
                     h1,
                     ip,
                     data,
                     state;
                     all_complex = all_complex,
                 )
    end
    return _lanczos_local_value(value, all_complex)
end

function accumulate_lanczos_green!(
    dst::Union{PhysicalQuantities,VMCPhysAccumulator},
    phys::PhysicalQuantities,
    data::ExpertModeData,
    state::VMCOptimizationState,
    w::Float64,
    h1::ComplexF64,
    h2::ComplexF64,
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int};
    all_complex::Bool,
)
    n_cis_ajs = length(phys.cis_ajs_idx)
    n_cis_ajs_ckt_alt = length(phys.cis_ajs_ckt_alt_idx)
    n_cis_ajs_ckt_alt_dc = length(data.green_two_terms)

    length(dst.phys_lanczos_qcisajsq) == LANCZOS_N_LS_HAM^2 * n_cis_ajs ||
        throw(
            ArgumentError(
                "phys_lanczos_qcisajsq has length $(length(dst.phys_lanczos_qcisajsq)); " *
                "expected $(LANCZOS_N_LS_HAM^2 * n_cis_ajs).",
            ),
        )
    length(dst.phys_lanczos_qcisajscktaltq) == LANCZOS_N_LS_HAM^2 * n_cis_ajs_ckt_alt ||
        throw(
            ArgumentError(
                "phys_lanczos_qcisajscktaltq has length " *
                "$(length(dst.phys_lanczos_qcisajscktaltq)); expected " *
                "$(LANCZOS_N_LS_HAM^2 * n_cis_ajs_ckt_alt).",
            ),
        )
    length(dst.phys_lanczos_qcisajscktaltq_dc) ==
    LANCZOS_N_LS_HAM^2 * n_cis_ajs_ckt_alt_dc ||
        throw(
            ArgumentError(
                "phys_lanczos_qcisajscktaltq_dc has length " *
                "$(length(dst.phys_lanczos_qcisajscktaltq_dc)); expected " *
                "$(LANCZOS_N_LS_HAM^2 * n_cis_ajs_ckt_alt_dc).",
            ),
        )

    original = _lanczos_config(ele_idx, ele_cfg, ele_num, ele_proj_cnt)
    lslq = (1.0 + 0.0im, h1, h1, h2)
    local_cis_ajs = _local_cis_ajs(dst, phys)
    local_cis_ajs_ckt_alt_dc = _local_cis_ajs_ckt_alt_dc(dst, phys)

    return _lanczos_preserve_m_all!(state) do
        lslca = Vector{ComplexF64}(undef, LANCZOS_N_LS_HAM * n_cis_ajs)
        @inbounds for idx in 1:n_cis_ajs
            lslca[idx] = _lanczos_local_value(local_cis_ajs[idx], all_complex)
        end
        @inbounds for (idx, (ri, _si, rj, sj)) in enumerate(phys.cis_ajs_idx)
            lslca[n_cis_ajs+idx] =
                _lanczos_hca_from_config!(
                    ri,
                    rj,
                    sj,
                    h1,
                    ip,
                    original,
                    data,
                    state;
                    all_complex = all_complex,
                )
        end

        @inbounds for rq = 0:(LANCZOS_N_LS_HAM-1), rp = 0:(LANCZOS_N_LS_HAM-1)
            right_q = _lanczos_local_value(lslq[rp*LANCZOS_N_LS_HAM+1], all_complex)
            for idx in 1:n_cis_ajs
                out_idx = idx + n_cis_ajs * (rp + LANCZOS_N_LS_HAM * rq)
                left = lslca[rq*n_cis_ajs+idx]
                dst.phys_lanczos_qcisajsq[out_idx] +=
                    w * (all_complex ? conj(left) : left) * right_q
            end
        end

        @inbounds for rq = 0:(LANCZOS_N_LS_HAM-1), rp = 0:(LANCZOS_N_LS_HAM-1)
            for (idx, (idx0, idx1)) in enumerate(phys.cis_ajs_ckt_alt_idx)
                out_idx = idx + n_cis_ajs_ckt_alt * (rp + LANCZOS_N_LS_HAM * rq)
                left = lslca[rq*n_cis_ajs+idx0]
                right = lslca[rp*n_cis_ajs+idx1]
                dst.phys_lanczos_qcisajscktaltq[out_idx] +=
                    w * (all_complex ? conj(left) : left) * right
            end
        end

        @inbounds for rq = 0:(LANCZOS_N_LS_HAM-1), rp = 0:(LANCZOS_N_LS_HAM-1)
            right_q = _lanczos_local_value(lslq[rp*LANCZOS_N_LS_HAM+1], all_complex)
            for (idx, term) in enumerate(data.green_two_terms)
                ri = term.site1
                rj = term.site2
                rk = term.site3
                rl = term.site4
                s = term.spin1 == :up ? 0 : 1
                t = term.spin3 == :up ? 0 : 1
                phys_value =
                    if rq == 0
                        _lanczos_local_value(local_cis_ajs_ckt_alt_dc[idx], all_complex)
                    else
                        _lanczos_hcaca_from_config!(
                            ri,
                            rj,
                            rk,
                            rl,
                            s,
                            t,
                            h1,
                            ip,
                            original,
                            data,
                            state;
                            all_complex = all_complex,
                        )
                    end
                out_idx = idx + n_cis_ajs_ckt_alt_dc * (rp + LANCZOS_N_LS_HAM * rq)
                dst.phys_lanczos_qcisajscktaltq_dc[out_idx] += w * right_q * phys_value
            end
        end
        nothing
    end
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
    return calculate_green_func_into!(
        data,
        state,
        state.phys_quantities,
        w,
        ip,
        ele_idx,
        ele_cfg,
        ele_num,
        ele_proj_cnt,
    )
end

function calculate_green_func!(
    data::ExpertModeData,
    state::VMCOptimizationState,
    acc::VMCPhysAccumulator,
    w::Float64,
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
)
    return calculate_green_func_into!(
        data,
        state,
        acc,
        w,
        ip,
        ele_idx,
        ele_cfg,
        ele_num,
        ele_proj_cnt,
    )
end

function calculate_green_func_into!(
    data::ExpertModeData,
    state::VMCOptimizationState,
    dst::Union{PhysicalQuantities,VMCPhysAccumulator,Nothing},
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

    dst === nothing &&
        error("PhysicalQuantities not initialized. Call initialize_phys_quantities! first.")

    phys = state.phys_quantities
    all_complex = get_all_complex_flag(data)

    # 1-body Green's function: <c†_{ri,s} c_{rj,s}>, over the canonical list
    # (which, when TwoBodyGEx is present, includes appended factored
    # constituents in C order). C uses s = CisAjsIdx[idx][3] = sj (annihilation
    # spin); for sz-conserved inputs si == sj.
    for (idx, (ri, si, rj, sj)) in enumerate(phys.cis_ajs_idx)
        # Calculate GreenFunc1 with s = sj to match C
        local_val = green_func1(
            ri,
            rj,
            sj,
            ip,
            ele_idx,
            ele_cfg,
            ele_num,
            ele_proj_cnt,
            data,
            state;
            all_complex = all_complex,
        )

        accumulate_phys_cis_ajs!(dst, phys, idx, w, local_val)
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

        accumulate_phys_cis_ajs_ckt_alt_dc!(dst, phys, idx, w, local_val)
    end

    # 2-body Green's function (product): <c†_i c_j> × conj(<c†_k c_l>)
    # using the resolved 1-based index pairs into the one-body Greens above.
    accumulate_factored_green!(dst, phys, w)
end
