"""
Orbital and QPTrans Utilities

Helper functions for building OrbitalSgn matrices and QPTrans mappings.
"""

"""
    build_orbital_sgn_matrix!(data::ExpertModeData)

Build OrbitalIdx and OrbitalSgn matrices from orbital_terms.
Equivalent to C's GetInfoOrbitalAntiParallel/GetInfoOrbitalGeneral/GetInfoOrbitalParallel.

# Note
- For iFlgOrbitalGeneral == 0: OrbitalIdx/OrbitalSgn[Nsite][Nsite]
- For iFlgOrbitalGeneral == 1: OrbitalIdx/OrbitalSgn[2*Nsite][2*Nsite]
"""
function build_orbital_sgn_matrix!(data::ExpertModeData)
    n_site = data.modpara.nsite

    # Determine APFlag: APFlag = 1 if NMPTrans < 0 (anti-periodic boundary condition)
    ap_flag = (data.modpara.nmp_trans < 0) ? 1 : 0

    if data.i_flg_orbital_general == 0
        # sz-conserved case: OrbitalIdx/OrbitalSgn[Nsite][Nsite]
        orbital_idx_matrix = zeros(Int, n_site, n_site)
        orbital_sgn = zeros(Int, n_site, n_site)

        # Build from orbital_terms
        for term in data.orbital_terms
            site1 = term.site1
            site2 = term.site2
            idx = term.idx
            sign = term.sign

            if 0 <= site1 < n_site && 0 <= site2 < n_site
                # C uses 0-based indexing, Julia uses 1-based
                orbital_idx_matrix[site1+1, site2+1] = idx
                orbital_sgn[site1+1, site2+1] = sign
            end
        end

        # If APFlag == 0, set all signs to 1 (C implementation behavior)
        if ap_flag == 0
            for i = 1:n_site
                for j = 1:n_site
                    orbital_sgn[i, j] = 1
                end
            end
        end

        data.orbital_idx_matrix = orbital_idx_matrix
        data.orbital_sgn = orbital_sgn
    else
        # fsz case: OrbitalIdx/OrbitalSgn[2*Nsite][2*Nsite]
        orbital_idx_matrix = zeros(Int, 2 * n_site, 2 * n_site)
        orbital_sgn = zeros(Int, 2 * n_site, 2 * n_site)

        # Determine the number of anti-parallel orbitals
        # Anti-parallel orbitals are those with idx < n_orbital_anti_parallel
        # Need to figure out where anti-parallel ends and parallel begins
        # This depends on whether OrbitalGeneral or Orbital+OrbitalParallel is used

        if data.i_flg_orbital_general == 1 &&
           data.i_flg_orbital_parallel == 0 &&
           data.i_flg_orbital_anti_parallel == 0
            # Pure OrbitalGeneral case: site indices already include spin information
            for term in data.orbital_terms
                site1 = term.site1
                site2 = term.site2
                idx = term.idx
                sign = term.sign

                if 0 <= site1 < 2 * n_site && 0 <= site2 < 2 * n_site
                    orbital_idx_matrix[site1+1, site2+1] = idx
                    orbital_sgn[site1+1, site2+1] = sign
                    if site1 != site2
                        orbital_idx_matrix[site2+1, site1+1] = idx
                        orbital_sgn[site2+1, site1+1] = -sign
                    end
                end
            end
        else
            # Orbital + OrbitalParallel case: need to map to correct matrix positions
            # First, find the boundary between anti-parallel and parallel orbitals
            # Anti-parallel orbitals have sequential indices starting from 0
            # Parallel orbitals have interleaved indices starting from n_anti_parallel

            n_anti_parallel = 0
            for term in data.orbital_terms
                # Anti-parallel orbitals: site1, site2 in [0, Nsite-1], idx < some threshold
                # Find the largest idx among terms that look like anti-parallel
                # Anti-parallel terms have site1 < Nsite and site2 < Nsite
                if term.site1 < n_site && term.site2 < n_site
                    # Check if this could be anti-parallel (sequential idx)
                    # vs parallel (interleaved idx after n_anti_parallel)
                    # We'll determine n_anti_parallel by finding the max+1 of sequential indices
                    if term.idx >= n_anti_parallel && (
                        n_anti_parallel == 0 ||
                        term.idx == n_anti_parallel ||
                        term.idx < n_anti_parallel + 10
                    )
                        # Looks like anti-parallel (sequential)
                        # Actually, let's just count anti-parallel as idx where idx is not interleaved
                        pass = true
                    end
                end
            end

            if data.i_flg_orbital_anti_parallel == 1 && data.i_flg_orbital_parallel == 1
                # We have both Orbital and OrbitalParallel. The parallel block
                # begins exactly at NArrayAP, recorded at parse time, matching
                # C's iNOrbitalAntiParallel offset (readdef.c GetInfoOrbitalParallel).
                n_anti_parallel = data.n_orbital_anti_parallel
            else
                n_anti_parallel = 0
            end

            # Build the matrix
            for term in data.orbital_terms
                site1 = term.site1
                site2 = term.site2
                idx = term.idx
                sign = term.sign

                if idx < n_anti_parallel
                    # Anti-parallel orbital: F_{i↑, j↓}
                    # Goes to matrix[i][j + Nsite] and matrix[j + Nsite][i]
                    all_i = site1  # up-spin row
                    all_j = site2 + n_site  # down-spin column

                    orbital_idx_matrix[all_i+1, all_j+1] = idx
                    orbital_sgn[all_i+1, all_j+1] = sign
                    # F_{JI} = -F_{IJ}
                    orbital_idx_matrix[all_j+1, all_i+1] = idx
                    orbital_sgn[all_j+1, all_i+1] = -sign
                else
                    # Parallel orbital
                    rel_idx = idx - n_anti_parallel  # Relative index in parallel section
                    is_down_down = (rel_idx % 2 == 1)  # Odd indices are down-down

                    if is_down_down
                        # down-down: F_{i↓, j↓} at matrix[i + Nsite][j + Nsite]
                        all_i = site1 + n_site
                        all_j = site2 + n_site
                    else
                        # up-up: F_{i↑, j↑} at matrix[i][j]
                        all_i = site1
                        all_j = site2
                    end

                    orbital_idx_matrix[all_i+1, all_j+1] = idx
                    orbital_sgn[all_i+1, all_j+1] = sign
                    # F_{JI} = -F_{IJ}
                    if all_i != all_j
                        orbital_idx_matrix[all_j+1, all_i+1] = idx
                        orbital_sgn[all_j+1, all_i+1] = -sign
                    end
                end
            end
        end

        # If APFlag == 0, set signs according to C implementation
        if ap_flag == 0
            for i = 1:(2*n_site)
                for j = (i+1):(2*n_site)
                    orbital_sgn[i, j] = 1
                    orbital_sgn[j, i] = -1
                end
            end
        end

        data.orbital_idx_matrix = orbital_idx_matrix
        data.orbital_sgn = orbital_sgn
    end
end

"""
    build_qp_trans_mappings!(data::ExpertModeData, file_path::String)

Build QPTrans, QPTransInv, QPTransSgn, QPOptTrans, and QPOptTransSgn mappings from qptransidx.def.
Equivalent to C's GetInfoTransSym.

# Format
- First NQPTrans lines: idx value (ParaQPTrans) - already parsed
- Then Nsite * NQPTrans lines: i j itmp itmpsgn (QPTrans indices)
  - i: original site index (0-based)
  - j: translated site index (0-based)
  - itmp: translated site index (itmp = j)
  - itmpsgn: sign (+1 or -1)
"""
function build_qp_trans_mappings!(data::ExpertModeData, file_path::String)
    n_site = data.modpara.nsite
    n_qp_trans = data.n_qp_trans

    if n_qp_trans <= 0 || n_site <= 0
        return
    end

    # Initialize arrays
    qp_trans = Vector{Vector{Int}}()
    qp_trans_inv = Vector{Vector{Int}}()
    qp_trans_sgn = Vector{Vector{Int}}()

    # Read file content
    content = read_def_file(file_path)
    lines = split(content, '\n')

    IGNORE_LINES_IN_DEF = 5

    # Skip header and ParaQPTrans lines
    line_idx = IGNORE_LINES_IN_DEF + 1
    for i = 1:n_qp_trans
        # Skip ParaQPTrans line
        while line_idx <= length(lines)
            line = clean_line(lines[line_idx])
            if !isempty(line)
                line_idx += 1
                break
            end
            line_idx += 1
        end
    end

    # Read QPTrans mappings: mpidx j itmp itmpsgn
    # Format: first column (i) is mpidx (translation operator index),
    #         second column (j) is original site index,
    #         third column (itmp) is translated site index
    # Initialize arrays for each mpidx
    for mpidx = 1:n_qp_trans
        trans = fill(-1, n_site)  # Initialize with -1 (invalid)
        trans_inv = fill(-1, n_site)
        trans_sgn = ones(Int, n_site)  # Default sign is +1

        push!(qp_trans, trans)
        push!(qp_trans_inv, trans_inv)
        push!(qp_trans_sgn, trans_sgn)
    end

    # Read all mapping lines
    while line_idx <= length(lines)
        line = clean_line(lines[line_idx])
        if isempty(line)
            line_idx += 1
            continue
        end

        tokens = split_def_line(line)
        # Support both 3-column (mpidx, j, itmp) and 4-column (mpidx, j, itmp, itmpsgn) formats
        if length(tokens) >= 3
            mpidx_0 = safe_parse_int(tokens[1], -1)  # Translation operator index (0-based)
            j = safe_parse_int(tokens[2], -1)  # Original site index (0-based)
            itmp = safe_parse_int(tokens[3], -1)  # Translated site index (0-based)
            itmpsgn = length(tokens) >= 4 ? safe_parse_int(tokens[4], 1) : 1  # Sign (default +1)

            # Validate indices
            if mpidx_0 >= 0 &&
               j >= 0 &&
               itmp >= 0 &&
               0 <= mpidx_0 < n_qp_trans &&
               0 <= j < n_site &&
               0 <= itmp < n_site
                mpidx_1 = mpidx_0 + 1  # Convert to 1-based for Julia
                # QPTrans[mpidx][j] = itmp (j is original, itmp is translated)
                qp_trans[mpidx_1][j+1] = itmp  # Convert to 1-based for Julia
                qp_trans_inv[mpidx_1][itmp+1] = j  # Inverse mapping
                qp_trans_sgn[mpidx_1][j+1] = itmpsgn
            end
        end
        line_idx += 1
    end

    # QPOptTrans: keep an already parsed OptTrans mapping. Otherwise install
    # the C FlagOptTrans<=0 identity mapping.
    qp_opt_trans = data.qp_opt_trans
    qp_opt_trans_sgn = data.qp_opt_trans_sgn
    if isempty(qp_opt_trans) || isempty(qp_opt_trans_sgn)
        n_qp_opt_trans = max(1, data.n_qp_opt_trans)
        qp_opt_trans = Vector{Vector{Int}}()
        qp_opt_trans_sgn = Vector{Vector{Int}}()

        for optidx = 1:n_qp_opt_trans
            opt_trans = collect(0:(n_site-1))  # Identity mapping (0-based)
            opt_trans_sgn = ones(Int, n_site)  # All signs are +1

            push!(qp_opt_trans, opt_trans)
            push!(qp_opt_trans_sgn, opt_trans_sgn)
        end
    end

    # Apply APFlag: if APFlag == 0, set all QPTransSgn to 1 (C implementation behavior)
    ap_flag = (data.modpara.nmp_trans < 0) ? 1 : 0
    if ap_flag == 0
        for mpidx = 1:n_qp_trans
            for j = 1:n_site
                qp_trans_sgn[mpidx][j] = 1
            end
        end
    end

    # Store in data
    data.qp_trans = qp_trans
    data.qp_trans_inv = qp_trans_inv
    data.qp_trans_sgn = qp_trans_sgn
    data.qp_opt_trans = qp_opt_trans
    data.qp_opt_trans_sgn = qp_opt_trans_sgn
end
