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
