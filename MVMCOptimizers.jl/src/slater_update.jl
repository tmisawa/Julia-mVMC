"""
Slater Matrix Update Functions

Update Slater matrix elements from variational parameters.
"""

"""
    build_orbital_idx_sgn_matrices(data::ExpertModeData) -> (Matrix{Int}, Matrix{Int}, Vector{ComplexF64})

Build OrbitalIdx and OrbitalSgn matrices from orbital_terms, and extract Slater parameter values.
Returns (OrbitalIdx, OrbitalSgn, Slater) where:
- OrbitalIdx[i+1, j+1]: orbital index for site pair (i, j) (1-based indexing)
- OrbitalSgn[i+1, j+1]: sign for site pair (i, j) (1-based indexing, +1 or -1)
- Slater: vector of orbital parameter values

# Note
This function constructs the matrices from orbital_terms and uses OrbitalSgn from ExpertModeData
if available. The orbital_terms should contain all site pairs (i, j) with their corresponding orbital indices.
"""
function build_orbital_idx_sgn_matrices(
    data::ExpertModeData,
)::Tuple{Matrix{Int},Matrix{Int},Vector{ComplexF64}}
    n_site = data.modpara.nsite
    n_orbital = length(data.orbital_terms)

    # Initialize matrices
    orbital_idx = fill(-1, n_site, n_site)

    # Use OrbitalSgn from ExpertModeData if available, otherwise build from orbital_terms
    if data.orbital_sgn !== nothing && size(data.orbital_sgn) == (n_site, n_site)
        orbital_sgn = data.orbital_sgn
    else
        # Fallback: build from orbital_terms
        orbital_sgn = zeros(Int, n_site, n_site)
        for term in data.orbital_terms
            site1 = term.site1
            site2 = term.site2
            if 0 <= site1 < n_site && 0 <= site2 < n_site
                orbital_sgn[site1+1, site2+1] = term.sign
            end
        end
    end

    # Find the maximum orbital index to size the slater array correctly
    max_orbital_idx = -1
    for term in data.orbital_terms
        if term.idx > max_orbital_idx
            max_orbital_idx = term.idx
        end
    end

    # Slater array size should be max_orbital_idx + 1 (since indices are 0-based)
    n_slater = max_orbital_idx >= 0 ? (max_orbital_idx + 1) : n_orbital
    slater = zeros(ComplexF64, n_slater)

    # Build mapping from (site1, site2) to orbital index
    # orbital_terms contains: site1, site2, value, idx
    # The orbital index (term.idx) is the unique parameter index from the input file
    for term in data.orbital_terms
        site1 = term.site1
        site2 = term.site2

        # Store orbital index mapping using term.idx (0-based)
        # C uses OrbitalIdx[tri][trj] which returns 0-based index
        if 0 <= site1 < n_site && 0 <= site2 < n_site
            orbital_idx[site1+1, site2+1] = term.idx  # Use the actual orbital index from file
        end

        # Store Slater parameter value at the correct index
        # term.idx is 0-based, so we use term.idx + 1 for Julia 1-based indexing
        # Only copy non-zero values to avoid overwriting with zeros from duplicate idx entries
        if term.idx >= 0 && term.idx < n_slater && abs(term.value) > 1e-14
            slater[term.idx+1] = term.value
        end
    end

    @debug "build_orbital_idx_sgn_matrices: n_orbital=$n_orbital, n_slater=$n_slater, max_orbital_idx=$max_orbital_idx"
    @debug "  slater[1:5]=$(slater[1:min(5, length(slater))])"
    @debug "  orbital_idx[1:3, 1:3]=$(orbital_idx[1:min(3, n_site), 1:min(3, n_site)])"

    return orbital_idx, orbital_sgn, slater
end

"""
    build_qp_trans_matrices(data::ExpertModeData) -> (Vector{Vector{Int}}, Vector{Vector{Int}}, Vector{Vector{Int}}, Vector{Vector{Int}}, Vector{Vector{Int}})

Build QPTrans, QPTransInv, QPTransSgn, QPOptTrans, and QPOptTransSgn matrices.
Equivalent to C's usage of QPTrans, QPTransInv, QPTransSgn, QPOptTrans, QPOptTransSgn.

# Behavior (matching C implementation)
- QPTrans, QPTransInv, QPTransSgn: **Required** (must be parsed from qptransidx.def).
  If not found, throws an error (no fallback).
- QPOptTrans, QPOptTransSgn: If parsed data exists, use it (FlagOptTrans > 0 equivalent).
  Otherwise, use identity mapping (FlagOptTrans <= 0 equivalent).

# Returns
Tuple of (QPTrans, QPTransInv, QPTransSgn, QPOptTrans, QPOptTransSgn)
"""
function build_qp_trans_matrices(
    data::ExpertModeData,
)::Tuple{
    Vector{Vector{Int}},
    Vector{Vector{Int}},
    Vector{Vector{Int}},
    Vector{Vector{Int}},
    Vector{Vector{Int}},
}
    n_site = data.modpara.nsite
    n_qp_trans = max(1, data.n_qp_trans)
    n_qp_opt_trans = max(1, data.n_qp_opt_trans)

    # QPTrans: Must exist (read from qptransidx.def) - no fallback
    # Equivalent to C's GetInfoTransSym() which is required
    if isempty(data.qp_trans) ||
       isempty(data.qp_trans_inv) ||
       isempty(data.qp_trans_sgn) ||
       length(data.qp_trans) < n_qp_trans ||
       length(data.qp_trans_inv) < n_qp_trans ||
       length(data.qp_trans_sgn) < n_qp_trans
        @error "QPTrans mappings not found. qptransidx.def must be parsed first via build_qp_trans_mappings!()."
        throw(
            ArgumentError(
                "QPTrans mappings are required. Please ensure qptransidx.def is parsed before calling update_slater_elm_fcmp!().",
            ),
        )
    end

    # Use parsed data from ExpertModeData
    qp_trans = data.qp_trans[1:n_qp_trans]
    qp_trans_inv = data.qp_trans_inv[1:n_qp_trans]
    qp_trans_sgn = data.qp_trans_sgn[1:n_qp_trans]

    # QPOptTrans: Use parsed data if available, otherwise use identity mapping
    # Equivalent to C's behavior: FlagOptTrans > 0 uses file, FlagOptTrans <= 0 uses identity
    if !isempty(data.qp_opt_trans) &&
       !isempty(data.qp_opt_trans_sgn) &&
       length(data.qp_opt_trans) >= n_qp_opt_trans &&
       length(data.qp_opt_trans_sgn) >= n_qp_opt_trans
        # Use parsed data from ExpertModeData (FlagOptTrans > 0 equivalent)
        qp_opt_trans = data.qp_opt_trans[1:n_qp_opt_trans]
        qp_opt_trans_sgn = data.qp_opt_trans_sgn[1:n_qp_opt_trans]
    else
        # Fallback: identity mapping (FlagOptTrans <= 0 equivalent)
        # This matches C's initialization in readdef.c:1138-1143
        qp_opt_trans = Vector{Vector{Int}}()
        qp_opt_trans_sgn = Vector{Vector{Int}}()

        for optidx = 1:n_qp_opt_trans
            opt_trans = collect(0:(n_site-1))  # Identity mapping: QPOptTrans[0][i] = i
            opt_trans_sgn = ones(Int, n_site)  # All signs are +1: QPOptTransSgn[0][i] = 1

            push!(qp_opt_trans, opt_trans)
            push!(qp_opt_trans_sgn, opt_trans_sgn)
        end
    end

    return qp_trans, qp_trans_inv, qp_trans_sgn, qp_opt_trans, qp_opt_trans_sgn
end

"""
    update_slater_elm_fcmp!(data::ExpertModeData, state::VMCOptimizationState)

Update Slater matrix elements from variational parameters (sz-conserved version).
Equivalent to C's `UpdateSlaterElm_fcmp()`.

This function implements the full Slater matrix update algorithm:
1. Loop over quantum projections (NQPFull)
2. For each projection, calculate Slater matrix elements
3. Apply translation and point group symmetries
4. Apply spin projection (cos/sin factors)

# Note
- OrbitalSgn is read from orbitalidx.def via ExpertModeData.orbital_sgn (if available).
- QPTrans mappings are read from qptransidx.def via ExpertModeData.qp_trans, etc. (if available).

# Algorithm
For each quantum projection index qpidx:
- Decompose qpidx into optidx, mpidx, spidx
- Get translation mappings (QPTrans, QPOptTrans) and signs
- Get spin projection factors (cos/sin)
- For each site pair (ri, rj):
  - Apply translation: ri -> ori -> tri, rj -> orj -> trj
  - Get orbital parameters: Slater[OrbitalIdx[tri][trj]]
  - Apply signs and spin projection
  - Store in SlaterElm matrix
"""
function update_slater_elm_fcmp!(data::ExpertModeData, state::VMCOptimizationState)
    # Get parameters
    n_site = data.modpara.nsite
    n_site2 = 2 * n_site
    n_sp_gauss_leg = max(1, data.modpara.nsp_gauss_leg)
    n_mp_trans = max(1, data.modpara.nmp_trans)
    n_qp_opt_trans = max(1, data.n_qp_opt_trans)

    # Counter for missing orbital indices (for debugging)
    missing_orbital_count = 0

    # Calculate NQPFix and NQPFull
    n_qp_fix = n_sp_gauss_leg * n_mp_trans
    n_qp_full = n_qp_fix * n_qp_opt_trans

    # Check quantum projection weights
    if data.qp_weights === nothing
        @error "Quantum projection weights not initialized. Call init_qp_weight!(data) first."
        return
    end

    weights = data.qp_weights
    if length(weights.spgl_cos) < n_sp_gauss_leg
        @error "Quantum projection weights have insufficient spin projection data."
        return
    end

    # Debug: Check orbital_terms values
    if !isempty(data.orbital_terms)
        nonzero_terms = [t for t in data.orbital_terms if abs(t.value) > 1e-14]
        @debug "update_slater_elm_fcmp!: orbital_terms count=$(length(data.orbital_terms)), nonzero_values=$(length(nonzero_terms))"
        if !isempty(data.orbital_terms)
            first_5 = [
                (t.site1, t.site2, t.idx, t.value) for
                t in data.orbital_terms[1:min(5, length(data.orbital_terms))]
            ]
            @debug "  first 5 orbital_terms: $first_5"
        end
    end

    # Build orbital index and sign matrices
    orbital_idx, orbital_sgn, slater = build_orbital_idx_sgn_matrices(data)
    @debug "update_slater_elm_fcmp!: slater length=$(length(slater)), slater[1:5]=$(slater[1:min(5, length(slater))])"

    # Build QPTrans matrices
    qp_trans, qp_trans_inv, qp_trans_sgn, qp_opt_trans, qp_opt_trans_sgn =
        build_qp_trans_matrices(data)

    # Main loop over quantum projections
    for qp_idx = 1:n_qp_full
        # Decompose qpidx: qpidx = optidx*NQPFix + NSPGaussLeg*mpidx + spidx
        # Julia uses 1-based indexing, so subtract 1 for calculations
        qpidx_0 = qp_idx - 1  # 0-based index

        optidx = qpidx_0 ÷ n_qp_fix
        remainder = qpidx_0 % n_qp_fix
        mpidx = remainder ÷ n_sp_gauss_leg
        spidx = remainder % n_sp_gauss_leg

        # Convert to 1-based for Julia arrays
        optidx_1 = optidx + 1
        mpidx_1 = mpidx + 1
        spidx_1 = spidx + 1

        # Get translation mappings
        if optidx_1 > length(qp_opt_trans) || mpidx_1 > length(qp_trans)
            @warn "QPTrans index out of bounds: optidx=$optidx_1, mpidx=$mpidx_1"
            continue
        end

        xqp_opt = qp_opt_trans[optidx_1]
        xqp_opt_sgn = qp_opt_trans_sgn[optidx_1]
        xqp = qp_trans[mpidx_1]
        xqp_sgn = qp_trans_sgn[mpidx_1]

        # Get spin projection factors
        cs = weights.spgl_cos_sin[spidx_1]
        cc = weights.spgl_cos_cos[spidx_1]
        ss = weights.spgl_sin_sin[spidx_1]

        # Debug: Check spin projection factors
        if qp_idx == 1
            @debug "update_slater_elm_fcmp!: qp_idx=1, spidx=$spidx, spidx_1=$spidx_1, cs=$cs, cc=$cc, ss=$ss"
        end

        # Get SlaterElm offset for this quantum projection
        slater_elm_offset = (qp_idx - 1) * n_site2 * n_site2

        # Loop over site pairs
        for ri = 0:(n_site-1)  # 0-based indexing to match C
            # Apply OptTrans and Trans: ri -> ori -> tri
            ori = xqp_opt[ri+1]  # Convert to 1-based for Julia array access
            tri = xqp[ori+1]
            sgni = xqp_sgn[ori+1] * xqp_opt_sgn[ri+1]

            # Spin indices
            rsi0 = ri  # up
            rsi1 = ri + n_site  # down

            # Get SlaterElm row pointers
            slater_elm_i0_offset = slater_elm_offset + rsi0 * n_site2
            slater_elm_i1_offset = slater_elm_offset + rsi1 * n_site2

            for rj = 0:(n_site-1)  # 0-based indexing
                # Apply OptTrans and Trans: rj -> orj -> trj
                orj = xqp_opt[rj+1]
                trj = xqp[orj+1]
                sgnj = xqp_sgn[orj+1] * xqp_opt_sgn[rj+1]

                # Spin indices
                rsj0 = rj  # up
                rsj1 = rj + n_site  # down

                # Get orbital indices (convert to 1-based for Julia)
                tri_1 = tri + 1
                trj_1 = trj + 1

                # Julia arrays are 1-based, so check both upper and lower bounds
                if tri_1 < 1 || tri_1 > n_site || trj_1 < 1 || trj_1 > n_site
                    @warn "Translated site index out of bounds: tri=$tri, trj=$trj (n_site=$n_site)" maxlog=1
                    continue
                end

                # Get orbital index and sign
                orb_idx_ij = orbital_idx[tri_1, trj_1]
                orb_sgn_ij = orbital_sgn[tri_1, trj_1]
                orb_idx_ji = orbital_idx[trj_1, tri_1]
                orb_sgn_ji = orbital_sgn[trj_1, tri_1]

                # Debug: Check orbital_idx for ri=0, rj=0 and ri=0, rj=1 (off-diagonal)
                if qp_idx == 1 && ri == 0 && rj == 0
                    @debug "update_slater_elm_fcmp!: qp_idx=1, ri=0, rj=0, tri=$tri, trj=$trj, orb_idx_ij=$orb_idx_ij, orb_idx_ji=$orb_idx_ji, orb_sgn_ij=$orb_sgn_ij, orb_sgn_ji=$orb_sgn_ji, sgni=$sgni, sgnj=$sgnj"
                end
                if qp_idx == 1 && ri == 0 && rj == 1
                    @debug "update_slater_elm_fcmp!: qp_idx=1, ri=0, rj=1, tri=$tri, trj=$trj, orb_idx_ij=$orb_idx_ij, orb_idx_ji=$orb_idx_ji, orb_sgn_ij=$orb_sgn_ij, orb_sgn_ji=$orb_sgn_ji"
                end
                if qp_idx == 1 && ri == 15 && rj == 9
                    @debug "update_slater_elm_fcmp!: qp_idx=1, ri=15, rj=9, tri=$tri, trj=$trj, orb_idx_ij=$orb_idx_ij, orb_idx_ji=$orb_idx_ji, orb_sgn_ij=$orb_sgn_ij, orb_sgn_ji=$orb_sgn_ji"
                end

                # Skip if orbital index is not defined (same as C implementation)
                # C code directly accesses OrbitalIdx[tri][trj] without checking,
                # but if the index is -1, the Slater parameter would be undefined.
                # In practice, undefined orbital pairs are skipped by setting the
                # Slater matrix element to 0 (which is the default).
                if orb_idx_ij < 0 || orb_idx_ji < 0
                    # Skip this site pair - Slater matrix element remains 0
                    missing_orbital_count += 1
                    continue
                end

                # Get Slater parameter values (convert to 1-based for Julia)
                slt_ij = slater[orb_idx_ij+1] * Float64(orb_sgn_ij * sgni * sgnj)
                slt_ji = slater[orb_idx_ji+1] * Float64(orb_sgn_ji * sgni * sgnj)

                # Debug: Check slater values for ri=0, rj=0 and off-diagonal
                if qp_idx == 1 && ri == 0 && rj == 0
                    @debug "update_slater_elm_fcmp!: qp_idx=1, ri=0, rj=0, slater[1]=$(slater[1]), slt_ij=$slt_ij, slt_ji=$slt_ji, cs=$cs"
                end
                if qp_idx == 1 && ri == 0 && rj == 1
                    @debug "update_slater_elm_fcmp!: qp_idx=1, ri=0, rj=1, orb_idx_ij=$orb_idx_ij, orb_idx_ji=$orb_idx_ji, slater[orb_idx_ij+1]=$(slater[orb_idx_ij + 1]), slater[orb_idx_ji+1]=$(slater[orb_idx_ji + 1]), slt_ij=$slt_ij, slt_ji=$slt_ji, (slt_ij-slt_ji)*cs=$((slt_ij - slt_ji) * cs)"
                end
                if qp_idx == 1 && ri == 15 && rj == 9
                    @debug "update_slater_elm_fcmp!: qp_idx=1, ri=15, rj=9, orb_idx_ij=$orb_idx_ij, orb_idx_ji=$orb_idx_ji, slater[orb_idx_ij+1]=$(slater[orb_idx_ij + 1]), slater[orb_idx_ji+1]=$(slater[orb_idx_ji + 1]), slt_ij=$slt_ij, slt_ji=$slt_ji, (slt_ij-slt_ji)*cs=$((slt_ij - slt_ji) * cs)"
                end

                # Calculate Slater matrix elements with spin projection
                # C code:
                # sltE_i0[rsj0] = -(slt_ij - slt_ji)*cs;   // up   - up
                # sltE_i0[rsj1] =   slt_ij*cc + slt_ji*ss; // up   - down
                # sltE_i1[rsj0] = -slt_ij*ss - slt_ji*cc;  // down - up
                # sltE_i1[rsj1] =  (slt_ij - slt_ji)*cs;   // down - down

                idx_i0_j0 = slater_elm_i0_offset + rsj0 + 1  # Convert to 1-based
                idx_i0_j1 = slater_elm_i0_offset + rsj1 + 1
                idx_i1_j0 = slater_elm_i1_offset + rsj0 + 1
                idx_i1_j1 = slater_elm_i1_offset + rsj1 + 1

                if idx_i0_j0 <= length(state.slater_matrix.slater_elm)
                    state.slater_matrix.slater_elm[idx_i0_j0] = -(slt_ij - slt_ji) * cs  # up - up
                end
                if idx_i0_j1 <= length(state.slater_matrix.slater_elm)
                    state.slater_matrix.slater_elm[idx_i0_j1] = slt_ij * cc + slt_ji * ss  # up - down
                end
                if idx_i1_j0 <= length(state.slater_matrix.slater_elm)
                    state.slater_matrix.slater_elm[idx_i1_j0] = -slt_ij * ss - slt_ji * cc  # down - up
                end
                if idx_i1_j1 <= length(state.slater_matrix.slater_elm)
                    state.slater_matrix.slater_elm[idx_i1_j1] = (slt_ij - slt_ji) * cs  # down - down
                end
            end
        end
    end

    # Log statistics about missing orbital indices (for debugging)
    if missing_orbital_count > 0
        total_pairs = n_qp_full * n_site * n_site
        @debug "update_slater_elm_fcmp!: $missing_orbital_count out of $total_pairs site pairs had undefined orbital indices (skipped)"
    end
end

"""
    build_orbital_idx_sgn_matrices_fsz(data::ExpertModeData) -> (Matrix{Int}, Matrix{Int}, Vector{ComplexF64})

Build OrbitalIdx and OrbitalSgn matrices for FSZ mode (2*Nsite × 2*Nsite).
Returns (OrbitalIdx, OrbitalSgn, Slater) where:
- OrbitalIdx[i+1, j+1]: orbital index for spin-site pair (i, j) (1-based indexing, i,j ∈ [0, 2*Nsite))
- OrbitalSgn[i+1, j+1]: sign for spin-site pair (i, j) (1-based indexing, +1 or -1)
- Slater: vector of orbital parameter values

For FSZ with OrbitalAntiParallel + OrbitalParallel files:
- Anti-parallel terms go into up-down and down-up blocks
- Parallel terms go into up-up and down-down blocks
"""
function build_orbital_idx_sgn_matrices_fsz(
    data::ExpertModeData,
)::Tuple{Matrix{Int},Matrix{Int},Vector{ComplexF64}}
    n_site = data.modpara.nsite
    n_site2 = 2 * n_site

    # Initialize 2*Nsite × 2*Nsite matrices
    orbital_idx = fill(-1, n_site2, n_site2)
    orbital_sgn = zeros(Int, n_site2, n_site2)

    # Default skew-symmetric signs: OrbitalSgn[i][j] = +1 for i < j, -1 for i > j
    for i = 1:n_site2
        for j = (i+1):n_site2
            orbital_sgn[i, j] = 1
            orbital_sgn[j, i] = -1
        end
    end

    # Determine whether we have both anti-parallel and parallel files
    has_anti_parallel = data.i_flg_orbital_anti_parallel == 1
    has_parallel = data.i_flg_orbital_parallel == 1

    # Find the maximum orbital index across all terms
    # For FSZ with anti-parallel + parallel, we need to account for the offset
    max_orbital_idx = -1
    n_anti_parallel_terms_count = n_site * n_site
    term_idx = 0
    max_anti_parallel_idx = -1
    max_parallel_idx = -1
    for term in data.orbital_terms
        term_idx += 1
        if has_anti_parallel && has_parallel && term_idx <= n_anti_parallel_terms_count
            if term.idx > max_anti_parallel_idx
                max_anti_parallel_idx = term.idx
            end
        else
            if term.idx > max_parallel_idx
                max_parallel_idx = term.idx
            end
        end
    end

    # Find max orbital index to determine Slater array size
    # With the pre-offset indices from MVMCExpertModeParsers.jl, we can just use the max idx
    max_idx = -1
    for term in data.orbital_terms
        if term.idx > max_idx
            max_idx = term.idx
        end
    end
    n_slater = max(max_idx + 1, 1)
    slater = zeros(ComplexF64, n_slater)

    # Determine the boundary between anti-parallel and parallel terms
    # Anti-parallel terms have indices 0..(n_anti_parallel-1)
    # Parallel terms have indices n_anti_parallel..(n_slater-1) with interleaved up/down
    n_anti_parallel_terms = n_site * n_site  # Number of anti-parallel term entries
    n_anti_parallel_idx = 0  # Maximum anti-parallel index + 1
    term_count = 0
    for term in data.orbital_terms
        term_count += 1
        if term_count <= n_anti_parallel_terms
            if term.idx >= n_anti_parallel_idx
                n_anti_parallel_idx = term.idx + 1
            end
        end
    end

    # Process orbital_terms
    # With pre-offset indices from MVMCExpertModeParsers.jl:
    # - Anti-parallel terms (first n_site*n_site): indices 0..143, go to up-down/down-up blocks
    # - Parallel up-up terms: even indices >= n_anti_parallel_idx, go to up-up block
    # - Parallel down-down terms: odd indices >= n_anti_parallel_idx, go to down-down block
    term_count = 0
    for term in data.orbital_terms
        term_count += 1
        site1 = term.site1
        site2 = term.site2
        actual_idx = term.idx

        # Skip invalid sites
        if site1 < 0 || site2 < 0
            continue
        end

        # Store Slater parameter value
        if actual_idx >= 0 && actual_idx < n_slater
            slater[actual_idx+1] = term.value
        end

        if has_anti_parallel && has_parallel
            # Both files exist
            if term_count <= n_anti_parallel_terms
                # Anti-parallel term: goes into up-down and down-up blocks
                if 0 <= site1 < n_site && 0 <= site2 < n_site
                    # up-down: row=site1 (up), col=site2+Nsite (down)
                    orbital_idx[site1+1, site2+n_site+1] = actual_idx
                    orbital_sgn[site1+1, site2+n_site+1] = term.sign

                    # down-up: row=site1+Nsite (down), col=site2 (up) - antisymmetric
                    orbital_idx[site1+n_site+1, site2+1] = actual_idx
                    orbital_sgn[site1+n_site+1, site2+1] = -term.sign
                end
            else
                # Parallel term: indices are already pre-offset and interleaved
                # Determine if this is up-up (even relative idx) or down-down (odd relative idx)
                rel_idx = actual_idx - n_anti_parallel_idx
                is_up_up = (rel_idx % 2 == 0)

                if 0 <= site1 < n_site && 0 <= site2 < n_site
                    if is_up_up
                        # up-up block
                        orbital_idx[site1+1, site2+1] = actual_idx
                        orbital_sgn[site1+1, site2+1] = term.sign
                        # Antisymmetric counterpart
                        orbital_idx[site2+1, site1+1] = actual_idx
                        orbital_sgn[site2+1, site1+1] = -term.sign
                    else
                        # down-down block
                        orbital_idx[site1+n_site+1, site2+n_site+1] = actual_idx
                        orbital_sgn[site1+n_site+1, site2+n_site+1] = term.sign
                        # Antisymmetric counterpart
                        orbital_idx[site2+n_site+1, site1+n_site+1] = actual_idx
                        orbital_sgn[site2+n_site+1, site1+n_site+1] = -term.sign
                    end
                end
            end
        elseif has_anti_parallel
            # Only anti-parallel file: up-down block only
            # Store Slater parameter value
            if term.idx >= 0 && term.idx < n_slater
                slater[term.idx+1] = term.value
            end
            if 0 <= site1 < n_site && 0 <= site2 < n_site
                orbital_idx[site1+1, site2+n_site+1] = term.idx
                orbital_sgn[site1+1, site2+n_site+1] = term.sign
                orbital_idx[site1+n_site+1, site2+1] = term.idx
                orbital_sgn[site1+n_site+1, site2+1] = -term.sign
            end
        else
            # General format: sites may already include spin info
            # Store Slater parameter value
            if term.idx >= 0 && term.idx < n_slater
                slater[term.idx+1] = term.value
            end
            if 0 <= site1 < n_site2 && 0 <= site2 < n_site2
                orbital_idx[site1+1, site2+1] = term.idx
                orbital_sgn[site1+1, site2+1] = term.sign
                orbital_idx[site2+1, site1+1] = term.idx
                orbital_sgn[site2+1, site1+1] = -term.sign
            end
        end
    end

    return orbital_idx, orbital_sgn, slater
end

"""
    update_slater_elm_fsz!(data::ExpertModeData, state::VMCOptimizationState)

Update Slater matrix elements from variational parameters (fsz version).
Equivalent to C's `UpdateSlaterElm_fsz()`.

# Note
FSZ mode does NOT use spin projection (SPGaussLeg must be 1).
The OrbitalIdx matrix is 2*Nsite × 2*Nsite, handling all spin combinations explicitly.
"""
function update_slater_elm_fsz!(data::ExpertModeData, state::VMCOptimizationState)
    # Get parameters
    n_site = data.modpara.nsite
    n_site2 = 2 * n_site
    n_sp_gauss_leg = max(1, data.modpara.nsp_gauss_leg)
    n_mp_trans = max(1, data.modpara.nmp_trans)
    n_qp_opt_trans = max(1, data.n_qp_opt_trans)

    # FSZ mode requires NSPGaussLeg = 1 (no spin projection)
    if n_sp_gauss_leg > 1
        @warn "FSZ mode requires NSPGaussLeg = 1, but got $n_sp_gauss_leg. Using 1."
        n_sp_gauss_leg = 1
    end

    # Calculate NQPFix and NQPFull
    n_qp_fix = n_sp_gauss_leg * n_mp_trans
    n_qp_full = n_qp_fix * n_qp_opt_trans

    # Use orbital index and sign matrices from parser (already built correctly)
    orbital_idx = data.orbital_idx_matrix
    orbital_sgn = data.orbital_sgn

    # Build slater array from orbital_terms values
    max_idx = 0
    for term in data.orbital_terms
        if term.idx > max_idx
            max_idx = term.idx
        end
    end
    slater = zeros(ComplexF64, max_idx + 1)
    for term in data.orbital_terms
        if term.idx >= 0 && term.idx <= max_idx
            slater[term.idx+1] = term.value
        end
    end

    # Build QPTrans matrices
    qp_trans, qp_trans_inv, qp_trans_sgn, qp_opt_trans, qp_opt_trans_sgn =
        build_qp_trans_matrices(data)

    # Get SlaterElm workspace from state
    # slater_elm is stored as 1D but accessed as 3D: [rsi + rsj * n_site2 + qp_idx * n_site2 * n_site2]
    # Size: n_site2 * n_site2 * n_qp_full
    slater_elm_flat = state.slater_matrix.slater_elm

    # Helper function to get/set slater_elm
    # C layout: slater_elm[qpidx * n_site2 * n_site2 + rsi * n_site2 + rsj]
    # Julia 1-based: slater_elm[(qp_idx - 1) * n_site2 * n_site2 + (rsi - 1) * n_site2 + rsj]
    function set_slater_elm!(elm, qp_idx, rsi, rsj, value)
        idx = (qp_idx - 1) * n_site2 * n_site2 + (rsi - 1) * n_site2 + rsj
        elm[idx] = value
    end

    # Main loop over quantum projections
    for qp_idx = 1:n_qp_full
        # Decompose qpidx: qpidx = optidx*NQPFix + NSPGaussLeg*mpidx + spidx
        qpidx_0 = qp_idx - 1  # 0-based index

        optidx = qpidx_0 ÷ n_qp_fix
        remainder = qpidx_0 % n_qp_fix
        mpidx = remainder ÷ n_sp_gauss_leg
        # spidx not used in FSZ (always 0)

        # Convert to 1-based for Julia arrays
        optidx_1 = optidx + 1
        mpidx_1 = mpidx + 1

        # Get translation mappings
        xqp_opt =
            optidx_1 <= length(qp_opt_trans) ? qp_opt_trans[optidx_1] :
            collect(0:(n_site-1))
        xqp_opt_sgn =
            optidx_1 <= length(qp_opt_trans_sgn) ? qp_opt_trans_sgn[optidx_1] :
            ones(Int, n_site)
        xqp = mpidx_1 <= length(qp_trans) ? qp_trans[mpidx_1] : collect(0:(n_site-1))
        xqp_sgn =
            mpidx_1 <= length(qp_trans_sgn) ? qp_trans_sgn[mpidx_1] : ones(Int, n_site)

        # FSZ: Loop over all site pairs (ri, rj) for both spins
        for ri = 0:(n_site-1)
            # Apply QPOptTrans: ri -> ori
            ori = ri < length(xqp_opt) ? xqp_opt[ri+1] : ri
            # Apply QPTrans: ori -> tri
            tri = ori < length(xqp) ? xqp[ori+1] : ori
            # Get sign
            sgni =
                (ri < length(xqp_opt_sgn) ? xqp_opt_sgn[ri+1] : 1) *
                (ori < length(xqp_sgn) ? xqp_sgn[ori+1] : 1)

            # FSZ: tri0 = tri (up), tri1 = tri + Nsite (down)
            tri0 = tri
            tri1 = tri + n_site

            # Row indices in SlaterElm (1-based)
            rsi0 = ri + 1         # up
            rsi1 = ri + n_site + 1 # down

            for rj = 0:(n_site-1)
                # Apply QPOptTrans: rj -> orj
                orj = rj < length(xqp_opt) ? xqp_opt[rj+1] : rj
                # Apply QPTrans: orj -> trj
                trj = orj < length(xqp) ? xqp[orj+1] : orj
                # Get sign
                sgnj =
                    (rj < length(xqp_opt_sgn) ? xqp_opt_sgn[rj+1] : 1) *
                    (orj < length(xqp_sgn) ? xqp_sgn[orj+1] : 1)

                # FSZ: trj0 = trj (up), trj1 = trj + Nsite (down)
                trj0 = trj
                trj1 = trj + n_site

                # Column indices in SlaterElm (1-based)
                rsj0 = rj + 1         # up
                rsj1 = rj + n_site + 1 # down

                # Get orbital indices and calculate Slater matrix elements
                # F_{IJ} - F_{JI} for each spin combination

                # up-up: sltE[rsi0][rsj0]
                idx_i0j0 = orbital_idx[tri0+1, trj0+1]
                idx_j0i0 = orbital_idx[trj0+1, tri0+1]
                sgn_i0j0 = orbital_sgn[tri0+1, trj0+1]
                sgn_j0i0 = orbital_sgn[trj0+1, tri0+1]
                slt_i0j0 =
                    (idx_i0j0 >= 0 && idx_i0j0 < length(slater)) ?
                    slater[idx_i0j0+1] * sgn_i0j0 * sgni * sgnj : 0.0
                slt_j0i0 =
                    (idx_j0i0 >= 0 && idx_j0i0 < length(slater)) ?
                    slater[idx_j0i0+1] * sgn_j0i0 * sgni * sgnj : 0.0
                set_slater_elm!(slater_elm_flat, qp_idx, rsi0, rsj0, slt_i0j0 - slt_j0i0)

                # up-down: sltE[rsi0][rsj1]
                idx_i0j1 = orbital_idx[tri0+1, trj1+1]
                idx_j1i0 = orbital_idx[trj1+1, tri0+1]
                sgn_i0j1 = orbital_sgn[tri0+1, trj1+1]
                sgn_j1i0 = orbital_sgn[trj1+1, tri0+1]
                slt_i0j1 =
                    (idx_i0j1 >= 0 && idx_i0j1 < length(slater)) ?
                    slater[idx_i0j1+1] * sgn_i0j1 * sgni * sgnj : 0.0
                slt_j1i0 =
                    (idx_j1i0 >= 0 && idx_j1i0 < length(slater)) ?
                    slater[idx_j1i0+1] * sgn_j1i0 * sgni * sgnj : 0.0
                set_slater_elm!(slater_elm_flat, qp_idx, rsi0, rsj1, slt_i0j1 - slt_j1i0)

                # down-up: sltE[rsi1][rsj0]
                idx_i1j0 = orbital_idx[tri1+1, trj0+1]
                idx_j0i1 = orbital_idx[trj0+1, tri1+1]
                sgn_i1j0 = orbital_sgn[tri1+1, trj0+1]
                sgn_j0i1 = orbital_sgn[trj0+1, tri1+1]
                slt_i1j0 =
                    (idx_i1j0 >= 0 && idx_i1j0 < length(slater)) ?
                    slater[idx_i1j0+1] * sgn_i1j0 * sgni * sgnj : 0.0
                slt_j0i1 =
                    (idx_j0i1 >= 0 && idx_j0i1 < length(slater)) ?
                    slater[idx_j0i1+1] * sgn_j0i1 * sgni * sgnj : 0.0
                set_slater_elm!(slater_elm_flat, qp_idx, rsi1, rsj0, slt_i1j0 - slt_j0i1)

                # down-down: sltE[rsi1][rsj1]
                idx_i1j1 = orbital_idx[tri1+1, trj1+1]
                idx_j1i1 = orbital_idx[trj1+1, tri1+1]
                sgn_i1j1 = orbital_sgn[tri1+1, trj1+1]
                sgn_j1i1 = orbital_sgn[trj1+1, tri1+1]
                slt_i1j1 =
                    (idx_i1j1 >= 0 && idx_i1j1 < length(slater)) ?
                    slater[idx_i1j1+1] * sgn_i1j1 * sgni * sgnj : 0.0
                slt_j1i1 =
                    (idx_j1i1 >= 0 && idx_j1i1 < length(slater)) ?
                    slater[idx_j1i1+1] * sgn_j1i1 * sgni * sgnj : 0.0
                set_slater_elm!(slater_elm_flat, qp_idx, rsi1, rsj1, slt_i1j1 - slt_j1i1)
            end
        end
    end

    return nothing
end
