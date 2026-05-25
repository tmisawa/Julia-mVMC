"""
VMC Sampling Functions

Monte Carlo sampling for electron configurations.
"""

# Update type enum
@enum UpdateType HOPPING EXCHANGE LOCALSPINFLIP NONE

# NOTE: the legacy TimerOutputs-based profiler (MVMC_TIMER / TIMER_ENABLED[] /
# VMC_TIMER / @timeit and enable_timer!/reset_timer!/show_timer/get_timer) was
# removed in favour of the C-compatible CTimer (see c_timer.jl). The old path
# left `if TIMER_ENABLED[]` Ref reads in the hot loop even when disabled; CTimer
# is threaded as a Val-typed argument so the disabled path inlines to no-ops.
# MVMC_TIMER is now a deprecated alias for MVMC_C_TIMER (handled in
# run_para_opt_from_namelist).

# RNG helpers to match C's gen_rand32/genrand_real2 behavior
const RNG_REAL2_INV = 1.0 / 4294967296.0
@inline rng_rand32(rng::AbstractRNG) = rand(rng, UInt32)
@inline rng_real2(rng::AbstractRNG) = Float64(rng_rand32(rng)) * RNG_REAL2_INV
@inline function rng_mod(rng::AbstractRNG, n::Int)
    n <= 0 && return 0
    return Int(rng_rand32(rng) % UInt32(n))
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
    get_loc_spn_array(data::ExpertModeData) -> Vector{Int}

Get LocSpn array from ExpertModeData.
Returns array where LocSpn[ri] = 1 if site ri has local spin, 0 otherwise.
"""
function get_loc_spn_array(data::ExpertModeData)::Vector{Int}
    n_site = data.modpara.nsite
    loc_spn = zeros(Int, n_site)
    for term in data.locspin_terms
        if term.spin_value == 1 && 0 <= term.site < n_site
            loc_spn[term.site+1] = 1  # Convert to 1-based indexing
        end
    end
    return loc_spn
end

function debug_output_dir()::String
    return get(ENV, "MVMC_DEBUG_OUTPUT_DIR", "")
end

function debug_step_tag()::String
    return get(ENV, "MVMC_DEBUG_STEP", "001")
end

function dump_loc_spn_if_enabled(loc_spn::Vector{Int}, n_site::Int)
    if get(ENV, "MVMC_DEBUG_LOCSPN", "0") == "0"
        return
    end
    base_dir = debug_output_dir()
    if isempty(base_dir)
        base_dir = "."
    end
    step_tag = debug_step_tag()
    out_path = joinpath(base_dir, "locspn_step_" * step_tag * ".dat")
    if isfile(out_path)
        return
    end
    open(out_path, "w") do f
        @printf(f, "# n_site=%d\n", n_site)
        @printf(f, "LocSpn ")
        for i in 1:n_site
            @printf(f, "%d ", loc_spn[i])
        end
        @printf(f, "\n")
    end
end

function dump_elec_initial_if_enabled(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_spn::Vector{Int},
    n_elec::Int,
    n_site::Int,
)
    if get(ENV, "MVMC_DEBUG_ELEC_PRE", "0") == "0"
        return
    end
    base_dir = debug_output_dir()
    if isempty(base_dir)
        base_dir = "."
    end
    step_tag = debug_step_tag()
    out_path = joinpath(base_dir, "elec_initial_sample_step_" * step_tag * ".dat")
    # Always overwrite to allow dumping at multiple steps
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    open(out_path, "w") do f
        @printf(f, "# n_elec=%d n_site=%d n_size=%d n_site2=%d\n", n_elec, n_site, n_size, n_site2)
        @printf(f, "EleIdx ")
        for i in 1:n_size
            @printf(f, "%d ", ele_idx[i])
        end
        @printf(f, "\nEleCfg ")
        for i in 1:n_site2
            @printf(f, "%d ", ele_cfg[i])
        end
        @printf(f, "\nEleNum ")
        for i in 1:n_site2
            @printf(f, "%d ", ele_num[i])
        end
        if !isempty(ele_spn)
            @printf(f, "\nEleSpn ")
            for i in 1:n_size
                @printf(f, "%d ", ele_spn[i])
            end
        end
        @printf(f, "\n")
    end
end

"""
    init_loc_spn!(state::VMCOptimizationState, data::ExpertModeData)

Initialize the cached loc_spn array in the workspace.
This should be called once before sampling begins.
"""
function init_loc_spn!(state::VMCOptimizationState, data::ExpertModeData)
    loc_spn = state.workspace.loc_spn
    fill!(loc_spn, 0)
    for term in data.locspin_terms
        if term.spin_value == 1 && 0 <= term.site < length(loc_spn)
            loc_spn[term.site+1] = 1  # Convert to 1-based indexing
        end
    end
end

"""
    make_proj_cnt!(proj_cnt::Vector{Int}, ele_num::Vector{Int}, data::ExpertModeData)

Calculate projection counts from electron number configuration.
Equivalent to C's `MakeProjCnt()`.
"""
function make_proj_cnt!(proj_cnt::Vector{Int}, ele_num::Vector{Int}, data::ExpertModeData)
    n_site = data.modpara.nsite
    n_proj = length(proj_cnt)
    n_gutzwiller_idx = data.n_gutzwiller_idx
    n_jastrow = length(data.jastrow_terms)
    gutzwiller_idx = data.gutzwiller_idx

    # Initialize
    fill!(proj_cnt, 0)

    # Get up-spin and down-spin electron numbers (0-based indexing in C, 1-based in Julia)
    n0 = @view ele_num[1:n_site]  # up-spin
    n1 = @view ele_num[(n_site+1):(2*n_site)]  # down-spin

    # Gutzwiller factor
    # C implementation: for(ri=0;ri<nSite;ri++) { projCnt[ GutzwillerIdx[ri] ] += n0[ri]*n1[ri]; }
    if n_gutzwiller_idx > 0 && !isempty(gutzwiller_idx)
        for ri = 0:(n_site-1)
            if ri + 1 <= length(gutzwiller_idx)
                idx = gutzwiller_idx[ri+1]  # 1-based array access, idx is 0-based value
                if idx + 1 <= n_proj
                    proj_cnt[idx+1] += n0[ri+1] * n1[ri+1]
                end
            end
        end
    elseif n_gutzwiller_idx > 0
        # Fallback: assume all sites map to index 0 (NGutzwillerIdx=1 case)
        for ri = 0:(n_site-1)
            proj_cnt[1] += n0[ri+1] * n1[ri+1]
        end
    end

    # Jastrow factor exp(sum {v_ij * (ni-1) * (nj-1)})
    # C implementation: projCnt[ NGutzwillerIdx + JastrowIdx[ri][rj] ] += xi * xj
    # Use jastrow_idx matrix for proper index mapping (like C implementation)
    offset = n_gutzwiller_idx
    jastrow_idx_matrix = data.jastrow_idx

    if data.n_jastrow_idx > 0 && !isempty(jastrow_idx_matrix)
        # Use proper JastrowIdx array (C-style: iterate over all site pairs)
        for ri = 0:(n_site-1)
            xi = n0[ri+1] + n1[ri+1] - 1
            if xi == 0
                continue
            end
            for rj = (ri+1):(n_site-1)
                xj = n0[rj+1] + n1[rj+1] - 1
                # jastrow_idx_matrix is 1-based in Julia, stores 0-based idx values
                if ri + 1 <= size(jastrow_idx_matrix, 1) &&
                   rj + 1 <= size(jastrow_idx_matrix, 2)
                    idx = jastrow_idx_matrix[ri+1, rj+1]  # 0-based idx value
                    proj_idx = offset + idx + 1  # 1-based Julia array index
                    if proj_idx <= n_proj
                        proj_cnt[proj_idx] += xi * xj
                    end
                end
            end
        end
    else
        # Fallback: single Jastrow parameter (backward compatibility)
        for term in data.jastrow_terms
            ri = term.site1
            rj = term.site2
            if 0 <= ri < n_site && 0 <= rj < n_site
                xi = n0[ri+1] + n1[ri+1] - 1
                if xi != 0
                    xj = n0[rj+1] + n1[rj+1] - 1
                    # Use offset + 1 for the single Jastrow index
                    jastrow_idx = offset + 1
                    if jastrow_idx <= n_proj
                        proj_cnt[jastrow_idx] += xi * xj
                    end
                end
            end
        end
    end

    # TODO: Add Doublon-Holon correlation factors (2-site and 4-site)
    # This requires more complex logic from C code
end

"""
    update_proj_cnt!(ri::Int, rj::Int, s::Int, proj_cnt_new::Vector{Int},
                    proj_cnt_old::Vector{Int}, ele_num::Vector{Int}, data::ExpertModeData)

Update projection counts when electron with spin s hops from ri to rj.
Equivalent to C's `UpdateProjCnt()`.
"""
function update_proj_cnt!(
    ri::Int,
    rj::Int,
    s::Int,
    proj_cnt_new::Vector{Int},
    proj_cnt_old::Vector{Int},
    ele_num::Vector{Int},
    data::ExpertModeData,
)
    n_site = data.modpara.nsite
    n_proj = length(proj_cnt_new)
    n_gutzwiller_idx = data.n_gutzwiller_idx
    gutzwiller_idx = data.gutzwiller_idx

    # Copy old counts
    if proj_cnt_new !== proj_cnt_old
        copy!(proj_cnt_new, proj_cnt_old)
    end

    if ri == rj
        return
    end

    # Get up-spin and down-spin electron numbers
    n0 = @view ele_num[1:n_site]  # up-spin
    n1 = @view ele_num[(n_site+1):(2*n_site)]  # down-spin

    # Gutzwiller factor
    # C implementation:
    #   idx = GutzwillerIdx[ri]; projCntNew[idx] -= n0[ri]+n1[ri];
    #   idx = GutzwillerIdx[rj]; projCntNew[idx] += n0[rj]*n1[rj];
    # Note: The subtraction uses ADDITION (n0+n1), the addition uses MULTIPLICATION (n0*n1)
    if n_gutzwiller_idx > 0
        if !isempty(gutzwiller_idx) &&
           ri + 1 <= length(gutzwiller_idx) &&
           rj + 1 <= length(gutzwiller_idx)
            idx_ri = gutzwiller_idx[ri+1]  # 0-based idx value
            idx_rj = gutzwiller_idx[rj+1]
            if idx_ri + 1 <= n_proj
                proj_cnt_new[idx_ri+1] -= n0[ri+1] + n1[ri+1]
            end
            if idx_rj + 1 <= n_proj
                proj_cnt_new[idx_rj+1] += n0[rj+1] * n1[rj+1]
            end
        else
            # Fallback: assume all sites map to index 0
            proj_cnt_new[1] -= n0[ri+1] + n1[ri+1]
            proj_cnt_new[1] += n0[rj+1] * n1[rj+1]
        end
    end

    # Jastrow offset
    offset = n_gutzwiller_idx
    jastrow_idx_matrix = data.jastrow_idx

    # Jastrow factor
    # C implementation uses JastrowIdx[ri][rj] to get the Jastrow parameter index.
    if data.n_jastrow_idx > 0 && !isempty(jastrow_idx_matrix)
        # Helper function to get Jastrow index (symmetric: ri < rj)
        function get_jastrow_idx(ra::Int, rb::Int)
            if ra < rb
                return jastrow_idx_matrix[ra+1, rb+1]
            else
                return jastrow_idx_matrix[rb+1, ra+1]
            end
        end

        # Update [ri][rj] term
        # C: projCntNew[idx] += n0[ri]+n1[ri]-n0[rj]-n1[rj]+1;
        idx = get_jastrow_idx(ri, rj)
        proj_idx = offset + idx + 1
        if proj_idx <= n_proj
            proj_cnt_new[proj_idx] += (n0[ri+1] + n1[ri+1]) - (n0[rj+1] + n1[rj+1]) + 1
        end

        # Update [ri][rk] terms (rk != ri, rj)
        # C: projCntNew[idx] -= n0[rk]+n1[rk]-1;
        for rk = 0:(n_site-1)
            if rk == rj || rk == ri
                continue
            end
            idx = get_jastrow_idx(ri, rk)
            proj_idx = offset + idx + 1
            if proj_idx <= n_proj
                proj_cnt_new[proj_idx] -= (n0[rk+1] + n1[rk+1] - 1)
            end
        end

        # Update [rj][rk] terms (rk != ri, rj)
        # C: projCntNew[idx] += n0[rk]+n1[rk]-1;
        for rk = 0:(n_site-1)
            if rk == ri || rk == rj
                continue
            end
            idx = get_jastrow_idx(rj, rk)
            proj_idx = offset + idx + 1
            if proj_idx <= n_proj
                proj_cnt_new[proj_idx] += (n0[rk+1] + n1[rk+1] - 1)
            end
        end
    elseif length(data.jastrow_terms) > 0
        # Fallback: single Jastrow parameter (backward compatibility)
        jastrow_idx = offset + 1
        if jastrow_idx <= n_proj
            proj_cnt_new[jastrow_idx] += (n0[ri+1] + n1[ri+1]) - (n0[rj+1] + n1[rj+1]) + 1
            for rk = 0:(n_site-1)
                if rk == rj || rk == ri
                    continue
                end
                proj_cnt_new[jastrow_idx] -= (n0[rk+1] + n1[rk+1] - 1)
            end
            for rk = 0:(n_site-1)
                if rk == ri || rk == rj
                    continue
                end
                proj_cnt_new[jastrow_idx] += (n0[rk+1] + n1[rk+1] - 1)
            end
        end
    end

    # TODO: Add Doublon-Holon correlation factors update
end

"""
    has_rbm_terms(data::ExpertModeData) -> Bool

Return true when any RBM section is present.
"""
function has_rbm_terms(data::ExpertModeData)::Bool
    return !isempty(data.charge_rbm_phys_layer_terms) ||
           !isempty(data.spin_rbm_phys_layer_terms) ||
           !isempty(data.general_rbm_phys_layer_terms) ||
           !isempty(data.charge_rbm_hidden_layer_terms) ||
           !isempty(data.spin_rbm_hidden_layer_terms) ||
           !isempty(data.general_rbm_hidden_layer_terms) ||
           !isempty(data.charge_rbm_phys_hidden_terms) ||
           !isempty(data.spin_rbm_phys_hidden_terms) ||
           !isempty(data.general_rbm_phys_hidden_terms)
end

@inline function _rbm_nidx(terms)::Int
    return isempty(terms) ? 0 : (maximum(t.idx for t in terms) + 1)
end

@inline function _rbm_nneuron(
    modpara_n::Int,
    hidden_terms,
    phys_hidden_terms,
)::Int
    nh = isempty(hidden_terms) ? 0 : (maximum(t.site for t in hidden_terms) + 1)
    nph = isempty(phys_hidden_terms) ? 0 : (maximum(t.site2 for t in phys_hidden_terms) + 1)
    return max(max(modpara_n, nh), nph)
end

"""
    make_rbm_cnt(ele_num::Vector{Int}, data::ExpertModeData) -> Vector{ComplexF64}

Build RBM counters from electron occupation.
Equivalent to C's `MakeRBMCnt()`.
"""
function make_rbm_cnt(ele_num::Vector{Int}, data::ExpertModeData)::Vector{ComplexF64}
    if !has_rbm_terms(data)
        return ComplexF64[]
    end

    n_site = data.modpara.nsite

    n_charge_phys = _rbm_nidx(data.charge_rbm_phys_layer_terms)
    n_spin_phys = _rbm_nidx(data.spin_rbm_phys_layer_terms)
    n_general_phys = _rbm_nidx(data.general_rbm_phys_layer_terms)
    n_phys = n_charge_phys + n_spin_phys + n_general_phys

    n_charge_neuron = _rbm_nneuron(
        data.modpara.nneuron_charge,
        data.charge_rbm_hidden_layer_terms,
        data.charge_rbm_phys_hidden_terms,
    )
    n_spin_neuron = _rbm_nneuron(
        data.modpara.nneuron_spin,
        data.spin_rbm_hidden_layer_terms,
        data.spin_rbm_phys_hidden_terms,
    )
    n_general_neuron = _rbm_nneuron(
        data.modpara.nneuron_general,
        data.general_rbm_hidden_layer_terms,
        data.general_rbm_phys_hidden_terms,
    )

    rbm_cnt = zeros(ComplexF64, n_phys + n_charge_neuron + n_spin_neuron + n_general_neuron)
    hidden_offset = n_phys

    # Potential on physical layer.
    for term in data.charge_rbm_phys_layer_terms
        ri = term.site
        if 0 <= ri < n_site
            idx = term.idx + 1
            if 1 <= idx <= n_charge_phys
                rbm_cnt[idx] += ele_num[ri+1] + ele_num[ri+n_site+1] - 1
            end
        end
    end
    spin_phys_offset = n_charge_phys
    for term in data.spin_rbm_phys_layer_terms
        ri = term.site
        if 0 <= ri < n_site
            idx = spin_phys_offset + term.idx + 1
            if 1 <= idx <= n_charge_phys + n_spin_phys
                rbm_cnt[idx] += ele_num[ri+1] - ele_num[ri+n_site+1]
            end
        end
    end
    general_phys_offset = n_charge_phys + n_spin_phys
    for term in data.general_rbm_phys_layer_terms
        ri = term.site
        si = term.spin
        if 0 <= ri < n_site && (si == 0 || si == 1)
            rsi = ri + si * n_site
            idx = general_phys_offset + term.idx + 1
            if 1 <= idx <= n_phys
                rbm_cnt[idx] += 2 * ele_num[rsi+1] - 1
            end
        end
    end

    # Potential on hidden layer.
    for term in data.charge_rbm_hidden_layer_terms
        hi = term.site + 1
        idx = hidden_offset + hi
        if 1 <= hi <= n_charge_neuron && 1 <= idx <= length(rbm_cnt)
            rbm_cnt[idx] += term.value
        end
    end
    spin_hidden_offset = hidden_offset + n_charge_neuron
    for term in data.spin_rbm_hidden_layer_terms
        hi = term.site + 1
        idx = spin_hidden_offset + hi
        if 1 <= hi <= n_spin_neuron && 1 <= idx <= length(rbm_cnt)
            rbm_cnt[idx] += term.value
        end
    end
    general_hidden_offset = hidden_offset + n_charge_neuron + n_spin_neuron
    for term in data.general_rbm_hidden_layer_terms
        hi = term.site + 1
        idx = general_hidden_offset + hi
        if 1 <= hi <= n_general_neuron && 1 <= idx <= length(rbm_cnt)
            rbm_cnt[idx] += term.value
        end
    end

    # Coupling between physical and hidden layers.
    for term in data.charge_rbm_phys_hidden_terms
        ri = term.site1
        hi = term.site2 + 1
        idx = hidden_offset + hi
        if 0 <= ri < n_site && 1 <= hi <= n_charge_neuron && 1 <= idx <= length(rbm_cnt)
            xi = ele_num[ri+1] + ele_num[ri+n_site+1] - 1
            rbm_cnt[idx] += term.value * xi
        end
    end
    for term in data.spin_rbm_phys_hidden_terms
        ri = term.site1
        hi = term.site2 + 1
        idx = spin_hidden_offset + hi
        if 0 <= ri < n_site && 1 <= hi <= n_spin_neuron && 1 <= idx <= length(rbm_cnt)
            xi = ele_num[ri+1] - ele_num[ri+n_site+1]
            rbm_cnt[idx] += term.value * xi
        end
    end
    for term in data.general_rbm_phys_hidden_terms
        ri = term.site1
        si = term.spin
        hi = term.site2 + 1
        idx = general_hidden_offset + hi
        if 0 <= ri < n_site &&
           (si == 0 || si == 1) &&
           1 <= hi <= n_general_neuron &&
           1 <= idx <= length(rbm_cnt)
            rsi = ri + si * n_site
            xi = 2 * ele_num[rsi+1] - 1
            rbm_cnt[idx] += term.value * xi
        end
    end

    return rbm_cnt
end

"""
    update_rbm_cnt_hopping!(rbm_cnt_new, rbm_cnt_old, ri, rj, s, data)

Incrementally update RBM counters for one hopping move `ri,s -> rj,s`.
Equivalent to C's `UpdateRBMCnt()`.
"""
function update_rbm_cnt_hopping!(
    rbm_cnt_new::Vector{ComplexF64},
    rbm_cnt_old::Vector{ComplexF64},
    ri::Int,
    rj::Int,
    s::Int,
    data::ExpertModeData,
)
    n = min(length(rbm_cnt_new), length(rbm_cnt_old))
    @inbounds for i = 1:n
        rbm_cnt_new[i] = rbm_cnt_old[i]
    end
    ri == rj && return

    n_site = data.modpara.nsite
    rsi = ri + s * n_site
    rsj = rj + s * n_site

    n_charge_phys = _rbm_nidx(data.charge_rbm_phys_layer_terms)
    n_spin_phys = _rbm_nidx(data.spin_rbm_phys_layer_terms)
    n_general_phys = _rbm_nidx(data.general_rbm_phys_layer_terms)
    n_phys = n_charge_phys + n_spin_phys + n_general_phys

    n_charge_neuron = _rbm_nneuron(
        data.modpara.nneuron_charge,
        data.charge_rbm_hidden_layer_terms,
        data.charge_rbm_phys_hidden_terms,
    )
    n_spin_neuron = _rbm_nneuron(
        data.modpara.nneuron_spin,
        data.spin_rbm_hidden_layer_terms,
        data.spin_rbm_phys_hidden_terms,
    )
    n_general_neuron = _rbm_nneuron(
        data.modpara.nneuron_general,
        data.general_rbm_hidden_layer_terms,
        data.general_rbm_phys_hidden_terms,
    )

    # Physical layer
    for term in data.charge_rbm_phys_layer_terms
        idx = term.idx + 1
        idx > length(rbm_cnt_new) && continue
        if term.site == ri
            rbm_cnt_new[idx] -= 1
        elseif term.site == rj
            rbm_cnt_new[idx] += 1
        end
    end
    spin_phys_offset = n_charge_phys
    spin_delta = 1 - 2 * s
    for term in data.spin_rbm_phys_layer_terms
        idx = spin_phys_offset + term.idx + 1
        idx > length(rbm_cnt_new) && continue
        if term.site == ri
            rbm_cnt_new[idx] -= spin_delta
        elseif term.site == rj
            rbm_cnt_new[idx] += spin_delta
        end
    end
    general_phys_offset = n_charge_phys + n_spin_phys
    for term in data.general_rbm_phys_layer_terms
        idx = general_phys_offset + term.idx + 1
        idx > length(rbm_cnt_new) && continue
        r = term.site + term.spin * n_site
        if r == rsi
            rbm_cnt_new[idx] -= 2
        elseif r == rsj
            rbm_cnt_new[idx] += 2
        end
    end

    # Phys-hidden coupling
    charge_cnt_offset = n_phys
    for term in data.charge_rbm_phys_hidden_terms
        if term.site1 == ri || term.site1 == rj
            cnt_idx = charge_cnt_offset + term.site2 + 1
            cnt_idx > length(rbm_cnt_new) && continue
            if term.site1 == ri
                rbm_cnt_new[cnt_idx] -= term.value
            else
                rbm_cnt_new[cnt_idx] += term.value
            end
        end
    end
    spin_cnt_offset = n_phys + n_charge_neuron
    for term in data.spin_rbm_phys_hidden_terms
        if term.site1 == ri || term.site1 == rj
            cnt_idx = spin_cnt_offset + term.site2 + 1
            cnt_idx > length(rbm_cnt_new) && continue
            if term.site1 == ri
                rbm_cnt_new[cnt_idx] -= spin_delta * term.value
            else
                rbm_cnt_new[cnt_idx] += spin_delta * term.value
            end
        end
    end
    general_cnt_offset = n_phys + n_charge_neuron + n_spin_neuron
    for term in data.general_rbm_phys_hidden_terms
        cnt_idx = general_cnt_offset + term.site2 + 1
        cnt_idx > length(rbm_cnt_new) && continue
        r = term.site1 + term.spin * n_site
        if r == rsi
            rbm_cnt_new[cnt_idx] -= 2 * term.value
        elseif r == rsj
            rbm_cnt_new[cnt_idx] += 2 * term.value
        end
    end

    return
end

"""
    log_rbm_ratio(rbm_cnt_new, rbm_cnt_old, data) -> ComplexF64

Compute log RBM ratio from RBM counters.
Equivalent to C's `LogRBMRatio()`.
"""
function log_rbm_ratio(
    rbm_cnt_new::Vector{ComplexF64},
    rbm_cnt_old::Vector{ComplexF64},
    data::ExpertModeData,
)::ComplexF64
    if !has_rbm_terms(data)
        return 0.0 + 0.0im
    end

    n_charge_phys = _rbm_nidx(data.charge_rbm_phys_layer_terms)
    n_spin_phys = _rbm_nidx(data.spin_rbm_phys_layer_terms)
    n_general_phys = _rbm_nidx(data.general_rbm_phys_layer_terms)
    n_phys = n_charge_phys + n_spin_phys + n_general_phys

    phys_params = zeros(ComplexF64, n_phys)
    for term in data.charge_rbm_phys_layer_terms
        idx = term.idx + 1
        if 1 <= idx <= n_charge_phys
            phys_params[idx] = term.value
        end
    end
    for term in data.spin_rbm_phys_layer_terms
        idx = n_charge_phys + term.idx + 1
        if n_charge_phys + 1 <= idx <= n_charge_phys + n_spin_phys
            phys_params[idx] = term.value
        end
    end
    for term in data.general_rbm_phys_layer_terms
        idx = n_charge_phys + n_spin_phys + term.idx + 1
        if n_charge_phys + n_spin_phys + 1 <= idx <= n_phys
            phys_params[idx] = term.value
        end
    end

    z = 0.0 + 0.0im
    n = min(length(rbm_cnt_new), length(rbm_cnt_old))
    n_phys_eff = min(n_phys, n)

    # Physical layer part (C: z += RBM[idx] * (rbmCntNew[idx] - rbmCntOld[idx]))
    for i = 1:n_phys_eff
        z += phys_params[i] * (rbm_cnt_new[i] - rbm_cnt_old[i])
    end

    # Hidden layer part.
    # Match C's LogRBMRatio() block-product algorithm, including principal-log branch.
    n_hidden = max(0, n - n_phys_eff)
    if n_hidden == 0
        return z
    end

    block_size = max(1, data.modpara.nblock_size_rbm_ratio)
    n_blk = (n_hidden - 1) ÷ block_size + 1

    for iblk = 0:(n_blk-1)
        hist = iblk * block_size + 1
        hiend = min(hist + block_size - 1, n_hidden)
        zz = 1.0 + 0.0im

        for hi = hist:hiend
            idx = n_phys_eff + hi
            rbm_new = rbm_cnt_new[idx]
            rbm_old = rbm_cnt_old[idx]

            if real(rbm_new) <= 0.0
                rbm_new = -rbm_new
            end
            if real(rbm_old) <= 0.0
                rbm_old = -rbm_old
            end

            z += (rbm_new - rbm_old)
            zz *= (1.0 + exp(-2.0 * rbm_new)) / (1.0 + exp(-2.0 * rbm_old))
        end

        z += log(zz)
    end

    return z
end

@inline function _rbm_log_cosh_stable(z::ComplexF64)::ComplexF64
    # Match C LogRBMRatio branch handling: flip sign when Re(z) <= 0.
    zp = real(z) <= 0.0 ? -z : z
    return zp + log1p(exp(-2.0 * zp)) - log(2.0)
end

"""
    log_rbm_val(ele_num::Vector{Int}, data::ExpertModeData) -> ComplexF64

Compute RBM log-weight for current electron occupation.
Equivalent to C's `LogWeightRBM(MakeRBMCnt(ele_num))`.
"""
function log_rbm_val(ele_num::Vector{Int}, data::ExpertModeData)::ComplexF64
    if !has_rbm_terms(data)
        return 0.0 + 0.0im
    end

    n_site = data.modpara.nsite
    n_charge = max(0, data.modpara.nneuron_charge)
    n_spin = max(0, data.modpara.nneuron_spin)
    n_general = max(0, data.modpara.nneuron_general)

    # Hidden-layer arguments are accumulated per neuron type.
    hidden_charge = zeros(ComplexF64, n_charge)
    hidden_spin = zeros(ComplexF64, n_spin)
    hidden_general = zeros(ComplexF64, n_general)

    z = 0.0 + 0.0im

    # Potential on physical layer.
    for term in data.charge_rbm_phys_layer_terms
        ri = term.site
        if 0 <= ri < n_site
            z += term.value * (ele_num[ri+1] + ele_num[ri+n_site+1] - 1)
        end
    end
    for term in data.spin_rbm_phys_layer_terms
        ri = term.site
        if 0 <= ri < n_site
            z += term.value * (ele_num[ri+1] - ele_num[ri+n_site+1])
        end
    end
    for term in data.general_rbm_phys_layer_terms
        ri = term.site
        si = term.spin
        if 0 <= ri < n_site && (si == 0 || si == 1)
            rsi = ri + si * n_site
            z += term.value * (2 * ele_num[rsi+1] - 1)
        end
    end

    # Potential on hidden layer.
    for term in data.charge_rbm_hidden_layer_terms
        hi = term.site + 1
        if 1 <= hi <= n_charge
            hidden_charge[hi] += term.value
        end
    end
    for term in data.spin_rbm_hidden_layer_terms
        hi = term.site + 1
        if 1 <= hi <= n_spin
            hidden_spin[hi] += term.value
        end
    end
    for term in data.general_rbm_hidden_layer_terms
        hi = term.site + 1
        if 1 <= hi <= n_general
            hidden_general[hi] += term.value
        end
    end

    # Coupling between physical and hidden layers.
    for term in data.charge_rbm_phys_hidden_terms
        ri = term.site1
        hi = term.site2 + 1
        if 0 <= ri < n_site && 1 <= hi <= n_charge
            xi = ele_num[ri+1] + ele_num[ri+n_site+1] - 1
            hidden_charge[hi] += term.value * xi
        end
    end
    for term in data.spin_rbm_phys_hidden_terms
        ri = term.site1
        hi = term.site2 + 1
        if 0 <= ri < n_site && 1 <= hi <= n_spin
            xi = ele_num[ri+1] - ele_num[ri+n_site+1]
            hidden_spin[hi] += term.value * xi
        end
    end
    for term in data.general_rbm_phys_hidden_terms
        ri = term.site1
        si = term.spin
        hi = term.site2 + 1
        if 0 <= ri < n_site && (si == 0 || si == 1) && 1 <= hi <= n_general
            rsi = ri + si * n_site
            xi = 2 * ele_num[rsi+1] - 1
            hidden_general[hi] += term.value * xi
        end
    end

    # Sum log(cosh(.)) over all hidden neurons.
    for v in hidden_charge
        z += _rbm_log_cosh_stable(v)
    end
    for v in hidden_spin
        z += _rbm_log_cosh_stable(v)
    end
    for v in hidden_general
        z += _rbm_log_cosh_stable(v)
    end

    return z
end

"""
    log_proj_ratio(proj_cnt_new::Vector{Int}, proj_cnt_old::Vector{Int}, data::ExpertModeData) -> Float64

Calculate log of projection ratio.
Equivalent to C's `LogProjRatio()`.
"""
function log_proj_ratio(
    proj_cnt_new::Vector{Int},
    proj_cnt_old::Vector{Int},
    data::ExpertModeData,
)::Float64
    z = 0.0
    n_proj = min(length(proj_cnt_new), length(proj_cnt_old))

    # Collect all projection parameters
    proj_params = ComplexF64[]
    for term in data.gutzwiller_terms
        push!(proj_params, term.value)
    end
    for term in data.jastrow_terms
        push!(proj_params, term.value)
    end
    # TODO: Add Doublon-Holon terms

    for idx = 1:n_proj
        if idx <= length(proj_params)
            z += real(proj_params[idx]) * (proj_cnt_new[idx] - proj_cnt_old[idx])
        end
    end

    return z
end

"""
    log_proj_val(proj_cnt::Vector{Int}, data::ExpertModeData)::Float64

Calculate the log of the projection value: sum(Proj[idx] * projCnt[idx]).
Equivalent to C's `LogProjVal()`.
"""
function log_proj_val(proj_cnt::Vector{Int}, data::ExpertModeData)::Float64
    z = 0.0
    n_proj = length(proj_cnt)

    # Collect all projection parameters
    proj_params = ComplexF64[]
    for term in data.gutzwiller_terms
        push!(proj_params, term.value)
    end
    for term in data.jastrow_terms
        push!(proj_params, term.value)
    end

    for idx = 1:n_proj
        if idx <= length(proj_params)
            z += real(proj_params[idx]) * proj_cnt[idx]
        end
    end

    return z
end

"""
    update_ele_config!(mi::Int, ri::Int, rj::Int, s::Int,
                      ele_idx::Vector{Int}, ele_cfg::Vector{Int}, ele_num::Vector{Int},
                      n_site::Int, n_elec::Int)

Update electron configuration when electron mi with spin s hops from ri to rj.
Equivalent to C's `updateEleConfig()`.
"""
function update_ele_config!(
    mi::Int,
    ri::Int,
    rj::Int,
    s::Int,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    n_site::Int,
    n_elec::Int,
)
    msa = mi + s * n_elec
    rsa_old = ri + s * n_site
    rsa_new = rj + s * n_site

    ele_idx[msa+1] = rj  # Convert to 1-based indexing
    ele_cfg[rsa_old+1] = -1
    ele_cfg[rsa_new+1] = mi
    ele_num[rsa_old+1] = 0
    ele_num[rsa_new+1] = 1
end

"""
    revert_ele_config!(mi::Int, ri::Int, rj::Int, s::Int,
                      ele_idx::Vector{Int}, ele_cfg::Vector{Int}, ele_num::Vector{Int},
                      n_site::Int, n_elec::Int)

Revert electron configuration (opposite of update_ele_config!).
Equivalent to C's `revertEleConfig()`.
"""
function revert_ele_config!(
    mi::Int,
    ri::Int,
    rj::Int,
    s::Int,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    n_site::Int,
    n_elec::Int,
)
    msa = mi + s * n_elec
    rsa_old = ri + s * n_site
    rsa_new = rj + s * n_site

    ele_idx[msa+1] = ri
    ele_cfg[rsa_old+1] = mi
    ele_cfg[rsa_new+1] = -1
    ele_num[rsa_old+1] = 1
    ele_num[rsa_new+1] = 0
end

"""
    make_candidate_hopping(ele_idx::Vector{Int}, ele_cfg::Vector{Int},
                          n_site::Int, n_elec::Int, loc_spn::Vector{Int},
                          rng::AbstractRNG) -> (mi, ri, rj, s, reject_flag)

Generate candidate hopping move.
Equivalent to C's `makeCandidate_hopping()`.

Returns a tuple (mi, ri, rj, s, reject_flag) where reject_flag is 1 if rejected, 0 otherwise.
"""
function make_candidate_hopping(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    n_site::Int,
    n_elec::Int,
    loc_spn::Vector{Int},
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    reject_flag = 0  # FALSE
    icnt_max = n_site * n_site

    # Select random electron (do-while loop to match C implementation)
    # C: do { mi = gen_rand32() % Ne; r = genrand_real2(); s = (r<0.5) ? 0 : 1; ri = eleIdx[mi+s*Ne]; } while (LocSpn[ri] == 1);
    # Note: Added icnt check to prevent infinite loop when all electrons are on local spin sites
    mi_val = 0
    s_val = 0
    ri_val = 0
    icnt_electron = 0
    while true
        # Check for infinite loop (all electrons on local spin sites)
        if icnt_electron > icnt_max
            reject_flag = 1  # TRUE - reject this candidate
            return (mi_val, ri_val, 0, s_val, reject_flag)
        end
        icnt_electron += 1

        # Select random electron
        mi_val = rng_mod(rng, n_elec)
        # Select random spin
        r = rng_real2(rng)
        s_val = r < 0.5 ? 0 : 1
        ri_val = ele_idx[mi_val+s_val*n_elec+1]  # 1-based indexing

        # Check if site has local spin (exit loop if not)
        if loc_spn[ri_val+1] != 1
            break
        end
    end

    # Select random destination site
    icnt = 0
    rj_val = rng_mod(rng, n_site)
    while ele_cfg[rj_val+s_val*n_site+1] != -1 || loc_spn[rj_val+1] == 1
        if icnt > icnt_max
            reject_flag = 1  # TRUE
            break
        end
        rj_val = rng_mod(rng, n_site)
        icnt += 1
    end

    return (mi_val, ri_val, rj_val, s_val, reject_flag)
end

"""
    make_candidate_exchange(ele_idx::Vector{Int}, ele_cfg::Vector{Int},
                           n_site::Int, n_elec::Int, ele_num::Vector{Int},
                           rng::AbstractRNG) -> (mi, ri, mj, rj, s, t, reject_flag)

Generate candidate exchange move (two electrons exchange positions).
Equivalent to C's `makeCandidate_exchange()` in vmcmake.c.

Returns a tuple (mi, ri, mj, rj, s, t, reject_flag) where reject_flag is 1 if rejected, 0 otherwise.

# Note
The mi-th electron with spin s exchanges with the electron on site rj with spin 1-s.
This requires finding sites where only one spin orientation is occupied.
"""
function make_candidate_exchange(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    n_site::Int,
    n_elec::Int,
    ele_num::Vector{Int},
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    # C implementation: makeCandidate_exchange in vmcmake.c
    # The mi-th electron with spin s exchanges with the electron on site rj with spin 1-s

    # First, check if there exists a site with only one electron (not half-filled)
    # C: for(ri=0;ri<Nsite;ri++) { if((eleNum[ri]+eleNum[ri+Nsite]) == 1) { flag = 0; break; } }
    flag = 1  # TRUE (reject by default)
    for ri_check = 0:(n_site-1)
        # ele_num uses 1-based indexing: ele_num[ri+1] for spin 0, ele_num[ri+n_site+1] for spin 1
        if (ele_num[ri_check+1] + ele_num[ri_check+n_site+1]) == 1
            flag = 0  # FALSE (found a site that can exchange)
            break
        end
    end

    if flag == 1
        return (0, 0, 0, 0, 0, 0, 1)  # reject_flag = 1
    end

    # Select first electron (mi, s) such that the opposite spin site is empty
    # C: do { mi = gen_rand32()%Ne; s = ...; ri = eleIdx[mi+s*Ne]; } while (eleCfg[ri+(1-s)*Nsite] != -1);
    mi_val = 0
    s_val = 0
    ri_val = 0
    while true
        mi_val = rng_mod(rng, n_elec)
        r = rng_real2(rng)
        s_val = r < 0.5 ? 0 : 1
        ri_val = ele_idx[mi_val+s_val*n_elec+1]  # 1-based indexing
        # Check if opposite spin site is empty (can receive the other electron)
        # C: eleCfg[ri+(1-s)*Nsite] != -1
        opposite_spin = 1 - s_val
        if ele_cfg[ri_val+opposite_spin*n_site+1] == -1
            break
        end
    end

    # Select second electron (mj, t=1-s) such that the opposite spin site (s) is empty
    # C: do { mj = gen_rand32()%Ne; t = 1-s; rj = eleIdx[mj+t*Ne]; } while (eleCfg[rj+(1-t)*Nsite] != -1);
    t_val = 1 - s_val
    mj_val = 0
    rj_val = 0
    while true
        mj_val = rng_mod(rng, n_elec)
        rj_val = ele_idx[mj_val+t_val*n_elec+1]  # 1-based indexing
        # Check if opposite spin site (which is s_val) is empty
        # C: eleCfg[rj+(1-t)*Nsite] != -1, where (1-t) = s
        if ele_cfg[rj_val+s_val*n_site+1] == -1
            break
        end
    end

    return (mi_val, ri_val, mj_val, rj_val, s_val, t_val, 0)  # reject_flag = 0
end

"""
    get_update_type(n_ex_update_path::Int, i_flg_orbital_general::Int, rng::AbstractRNG; two_sz::Int = 0) -> UpdateType

Get update type based on path parameter.
Equivalent to C's `getUpdateType()`.

# Arguments
- `n_ex_update_path`: Update path parameter (NExUpdatePath)
- `i_flg_orbital_general`: Orbital general flag (0 = sz conserved, non-zero = fsz)
- `rng`: Random number generator
- `two_sz`: TwoSz parameter (-1 = non-conserved, else conserved)
"""
function get_update_type(
    n_ex_update_path::Int,
    i_flg_orbital_general::Int,
    rng::AbstractRNG = Random.GLOBAL_RNG;
    two_sz::Int = 0,
)::UpdateType
    if n_ex_update_path == 0
        return HOPPING
    elseif n_ex_update_path == 1
        r = rng_real2(rng)
        return r < 0.5 ? EXCHANGE : HOPPING
    elseif n_ex_update_path == 2
        if i_flg_orbital_general == 0
            return EXCHANGE
        else
            # FSZ mode: When TwoSz == -1, use EXCHANGE or LOCALSPINFLIP
            # When TwoSz != -1, use EXCHANGE only
            # This matches C implementation: vmcmake.c line 611-620
            if two_sz == -1
                r = rng_real2(rng)
                return r < 0.5 ? EXCHANGE : LOCALSPINFLIP
            else
                return EXCHANGE
            end
        end
    elseif n_ex_update_path == 3
        # For KondoGC
        r1 = rng_real2(rng)
        if r1 < 0.5
            return HOPPING
        else
            r2 = rng_real2(rng)
            return r2 < 0.5 ? EXCHANGE : LOCALSPINFLIP
        end
    end
    return NONE
end

# ============================================================================
# Initial Sample Generation
# ============================================================================

"""
    make_initial_sample!(ele_idx::Vector{Int}, ele_cfg::Vector{Int}, ele_num::Vector{Int},
                        ele_proj_cnt::Vector{Int}, data::ExpertModeData, state::VMCOptimizationState, rng::AbstractRNG) -> Int

Generate initial electron configuration sample.
Equivalent to C's `makeInitialSample()`.
Returns 0 on success, non-zero on error.
"""
function make_initial_sample!(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)::Int
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    loc_spn = get_loc_spn_array(data)
    dump_loc_spn_if_enabled(loc_spn, n_site)

    flag = 1
    loop = 0
    max_loops = 100

    while flag > 0 && loop < max_loops
        # Initialize
        fill!(ele_idx, -1)
        fill!(ele_cfg, -1)

        # Local spin sites
        for ri = 0:(n_site-1)
            if loc_spn[ri+1] == 1
                mi_val = -1
                si_val = -1
                while true
                    mi_val = rng_mod(rng, n_elec)
                    r = rng_real2(rng)
                    si_val = r < 0.5 ? 0 : 1
                    if ele_idx[mi_val+si_val*n_elec+1] == -1
                        break
                    end
                end
                ele_cfg[ri+si_val*n_site+1] = mi_val
                ele_idx[mi_val+si_val*n_elec+1] = ri
            end
        end

        # Itinerant electrons
        for si = 0:1
            for mi = 0:(n_elec-1)
                if ele_idx[mi+si*n_elec+1] == -1
                    ri_val = -1
                    while true
                        ri_val = rng_mod(rng, n_site)
                        if ele_cfg[ri_val+si*n_site+1] == -1 && loc_spn[ri_val+1] != 1
                            break
                        end
                    end
                    ele_cfg[ri_val+si*n_site+1] = mi
                    ele_idx[mi+si*n_elec+1] = ri_val
                end
            end
        end

        # Calculate EleNum
        for rsi = 0:(n_site2-1)
            ele_num[rsi+1] = ele_cfg[rsi+1] < 0 ? 0 : 1
        end

        # Calculate projection counts
        make_proj_cnt!(ele_proj_cnt, ele_num, data)

        # Calculate Pfaffian and check validity (similar to C implementation)
        # C: flag = CalculateMAll_fcmp(eleIdx, qpStart, qpEnd);
        # Get QP range from state
        n_qp_full = length(state.slater_matrix.pf_m)
        qp_start = 1
        qp_end = n_qp_full + 1

        # Call calculate_m_all_fcmp! to check Pfaffian validity
        flag = calculate_m_all_fcmp!(ele_idx, qp_start, qp_end, data, state)

        loop += 1
    end

    if loop >= max_loops
        @error "make_initial_sample!: Too many loops"
        return 1
    end

    if get(ENV, "MVMC_DEBUG_INIT_SAMPLE", "0") != "0"
        @info "make_initial_sample!: loops=$loop flag=$flag"
    end

    dump_elec_initial_if_enabled(
        ele_idx,
        ele_cfg,
        ele_num,
        Int[],
        n_elec,
        n_site,
    )

    return 0
end

# ============================================================================
# Pfaffian Calculation Functions (Stub Implementation)
# ============================================================================

"""
    calculate_m_all_fcmp!(ele_idx::Vector{Int}, qp_start::Int, qp_end::Int,
                         data::ExpertModeData, state::VMCOptimizationState) -> Int

Calculate Pfaffian and inverse matrix for all QP indices.
Equivalent to C's `CalculateMAll_fcmp()`.
Returns 0 on success, non-zero on error.

This function wraps PfaPack.jl's calculate_m_all_fcmp! to work with
MVMCOptimizers.jl's data structures.

# Note
- Converts 1D arrays in SlaterMatrixData to 3D array views for PfaPack.jl
- Handles the extra element in inv_m (n_size * n_size + 1 per QP) by using only
  the first n_size * n_size elements per QP
"""
function calculate_m_all_fcmp!(
    ele_idx::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)::Int
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec  # Total number of electrons (2*Ne)
    n_site2 = 2 * n_site
    n_qp_full = length(state.slater_matrix.pf_m)

    # Validate qp_start and qp_end (1-based indexing)
    if qp_start < 1 || qp_end > n_qp_full + 1 || qp_start >= qp_end
        @error "Invalid qp_start or qp_end: qp_start=$qp_start, qp_end=$qp_end, n_qp_full=$n_qp_full"
        return 1
    end

    # slater_elm: 1D array with C row-major layout
    # [qp_idx * n_site2 * n_site2 + rsi * n_site2 + rsj]
    # We pass the full 1D array and let PfaPack handle the offset

    # Calculate the offset for the qp_start range
    qp_offset = (qp_start - 1) * n_site2 * n_site2
    qp_num = qp_end - qp_start

    # Create a view of slater_elm for the relevant QP range
    slater_elm_subset = view(
        state.slater_matrix.slater_elm,
        (qp_offset+1):(qp_offset+qp_num*n_site2*n_site2),
    )

    # Use pre-allocated workspace arrays instead of allocating new ones
    # Create views for the required qp_num (workspace is allocated for n_qp_full)
    ws = state.workspace
    inv_m_temp = view(ws.inv_m_temp,:,:,(1:qp_num))
    pf_m_temp = view(ws.pf_m_temp, 1:qp_num)

    # Zero out the workspace views before use
    fill!(inv_m_temp, 0.0 + 0.0im)
    fill!(pf_m_temp, 0.0 + 0.0im)

    # Call PfaPack.jl's function
    # Note: PfaPack uses 1-based indexing for qp_start and qp_end
    # n_elec in PfaPack = n_size (2*Ne) in MVMCOptimizers
    # The subset views are passed directly (1-based relative to the subset)
    info = calculate_m_all_fcmp_pfapack!(
        ele_idx,
        slater_elm_subset,
        inv_m_temp,
        pf_m_temp,
        1,  # qp_start relative to the subset (always 1 for the subset)
        qp_num + 1,  # qp_end relative to the subset (exclusive)
        n_site,
        n_size,  # n_elec in PfaPack = total electrons (2*Ne)
        state.workspace.pfapack_workspace,  # Pre-allocated workspace
    )

    if info != 0
        return info
    end

    # Copy results back to state
    # pf_m: direct copy
    for qp = 1:qp_num
        state.slater_matrix.pf_m[qp_start+qp-1] = pf_m_temp[qp]
    end

    # inv_m: need to copy from 3D temp to 1D state array
    # C layout: InvM + qpidx*Nsize*Nsize, then InvM[msi*nsize + msj] for access
    # This is row-major storage: inv_m[qpidx][msi][msj] = inv_m[qpidx * n_size * n_size + msi * n_size + msj]
    # Note: inv_m_temp contains the result from calculate_m_all_child_fcmp!
    # which stores inv_m[msj, msi] (column-major for LTL decomposition)
    # After M_ZSCAL, we need to copy it in row-major layout for C compatibility
    # Since inv_m_temp[msj, msi, qp] contains the value, we copy it to row-major layout
    # Use copyto! for efficient vectorized copy instead of triple nested loops
    nsq = n_size * n_size
    @inbounds for qp = 1:qp_num
        qp_idx = qp_start + qp - 1  # 1-based
        dst_start = (qp_idx - 1) * nsq + 1
        dst_end = dst_start + nsq - 1
        if dst_end <= length(state.slater_matrix.inv_m)
            dst = view(state.slater_matrix.inv_m, dst_start:dst_end)
            src = view(inv_m_temp,:,:,qp)
            copyto!(dst, vec(src))
        end
    end

    return info
end

# Note: calculate_m_all_fsz! is defined later in this file (around line 2490)
# using direct julia_zsktf2!, utu2pfa, cimpl_utu2inv! calls.

"""
    calculate_new_pf_m2!(ma::Int, s::Int, pf_m_new::Vector{ComplexF64},
                        ele_idx::Vector{Int}, qp_start::Int, qp_end::Int,
                        data::ExpertModeData, state::VMCOptimizationState)

Calculate new Pfaffian after electron ma with spin s hops.
Equivalent to C's `CalculateNewPfM2()`.

This function calculates the new Pfaffian using the formula:
    pfMNew[qpidx] = -ratio * PfM[qpidx]
where ratio = sum_j invM[msa][msj] * sltE[rsa][rsj]

# Reference
- C implementation: mVMC/src/mVMC/pfupdate.c:78-116 (CalculateNewPfM2)
"""
function calculate_new_pf_m2!(
    ma::Int,
    s::Int,
    pf_m_new::Vector{ComplexF64},
    ele_idx::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    msa = ma + s * n_elec  # 0-based electron index
    rsa = ele_idx[msa+1] + s * n_site  # 0-based site index (ele_idx is 1-based)

    # Get array references
    slater_elm = state.slater_matrix.slater_elm
    inv_m_arr = state.slater_matrix.inv_m
    pf_m = state.slater_matrix.pf_m

    # Process each QP index
    # C implementation: for(qpidx=0;qpidx<qpNum;qpidx++) where qpNum = qpEnd - qpStart
    # C uses local index qpidx for InvM, but global index (qpidx+qpStart) for SlaterElm
    # Julia: qpidx is global (1-based), so we need to convert to local for InvM
    qp_num = qp_end - qp_start
    for local_qpidx = 0:(qp_num-1)
        qpidx = qp_start + local_qpidx  # Global index (1-based)
        if qpidx < 1 || qpidx > length(pf_m)
            continue
        end

        # C: sltE = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2
        # In C, qpidx is local (0-based), and qpStart is 0-based offset
        # So (qpidx+qpStart) is the global 0-based index
        # In Julia, qpidx is global (1-based), so we use (qpidx-1) for 0-based global index
        slt_offset = (qpidx - 1) * n_site2 * n_site2

        # InvM offset (global index) - state.slater_matrix.inv_m stores all QP indices
        # Note: In C, InvM uses local index when qpStart=0 (typical case),
        # but Julia always stores with global indices for consistency
        inv_offset = (qpidx - 1) * n_size * n_size

        # Calculate ratio = sum_j invM[msa][msj] * sltE[rsa][rsj]
        ratio = 0.0 + 0.0im
        for msj = 0:(n_size-1)
            # Calculate rsj from msj
            if msj < n_elec
                rsj = ele_idx[msj+1]  # up-spin, ele_idx is 1-based
            else
                rsj = ele_idx[msj+1] + n_site  # down-spin
            end

            # C: invM_a[msj] = InvM[inv_offset + msa*Nsize + msj]
            # C: sltE_a[rsj] = SlaterElm[slt_offset + rsa*Nsite2 + rsj]
            inv_m_a_msj = inv_m_arr[inv_offset+msa*n_size+msj+1]
            slt_e_a_rsj = slater_elm[slt_offset+rsa*n_site2+rsj+1]

            ratio += inv_m_a_msj * slt_e_a_rsj
        end

        # Update pfMNew: pfMNew[qpidx] = -ratio * PfM[qpidx]
        pf_m_new[qpidx] = -ratio * pf_m[qpidx]
    end
end

"""
    calculate_new_pf_m2_fsz!(ma::Int, s::Int, pf_m_new::Vector{ComplexF64},
                             ele_idx::Vector{Int}, ele_spn::Vector{Int},
                             qp_start::Int, qp_end::Int,
                             data::ExpertModeData, state::VMCOptimizationState)

Calculate new Pfaffian for FSZ mode after hopping electron ma to site ele_idx[ma] with spin s.

In FSZ mode, electron spins are tracked individually via `ele_spn` array.
The key difference from regular mode:
- `msa = ma` (no spin offset, all electrons in single array)
- `rsa = ele_idx[ma] + s * n_site` (use new spin s for hopped electron)
- For other electrons: `rsj = ele_idx[msj] + ele_spn[msj] * n_site` (use ele_spn)

# Arguments
- `ma::Int`: Electron index (0-based, as in C)
- `s::Int`: New spin of the electron (0 = up, 1 = down)
- `pf_m_new::Vector{ComplexF64}`: Output array for new Pfaffian values
- `ele_idx::Vector{Int}`: Electron site indices (1-based array, 0-based values)
- `ele_spn::Vector{Int}`: Electron spin indices (1-based array, values 0 or 1)
- `qp_start::Int`: Start QP index (1-based)
- `qp_end::Int`: End QP index (1-based, exclusive)

# Reference
- C implementation: mVMC/src/mVMC/pfupdate_fsz.c:35-71 (CalculateNewPfM_fsz)
"""
function calculate_new_pf_m2_fsz!(
    ma::Int,
    s::Int,
    pf_m_new::Vector{ComplexF64},
    ele_idx::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site

    # FSZ: msa = ma (no spin offset)
    msa = ma
    # FSZ: rsa = ele_idx[ma] + s * n_site (use new spin s)
    rsa = ele_idx[msa+1] + s * n_site

    # Get array references
    slater_elm = state.slater_matrix.slater_elm
    inv_m_arr = state.slater_matrix.inv_m
    pf_m = state.slater_matrix.pf_m

    # Process each QP index
    qp_num = qp_end - qp_start
    for local_qpidx = 0:(qp_num-1)
        qpidx = qp_start + local_qpidx
        if qpidx < 1 || qpidx > length(pf_m)
            continue
        end

        # SlaterElm offset (global index)
        slt_offset = (qpidx - 1) * n_site2 * n_site2

        # InvM offset (global index) - state.slater_matrix.inv_m stores all QP indices
        inv_offset = (qpidx - 1) * n_size * n_size

        # Calculate ratio = sum_j invM[msa][msj] * sltE[rsa][rsj]
        ratio = 0.0 + 0.0im
        for msj = 0:(n_size-1)
            # FSZ: rsj = ele_idx[msj] + ele_spn[msj] * n_site
            rsj = ele_idx[msj+1] + ele_spn[msj+1] * n_site

            inv_m_a_msj = inv_m_arr[inv_offset+msa*n_size+msj+1]
            slt_e_a_rsj = slater_elm[slt_offset+rsa*n_site2+rsj+1]

            ratio += inv_m_a_msj * slt_e_a_rsj
        end

        # Update pfMNew: pfMNew[qpidx] = -ratio * PfM[qpidx]
        pf_m_new[qpidx] = -ratio * pf_m[qpidx]
    end
end

"""
    calculate_log_ip_fcmp(pf_m::Vector{ComplexF64}, qp_start::Int, qp_end::Int,
                         data::ExpertModeData) -> ComplexF64

Calculate logarithm of inner product <phi|L|x>.
Equivalent to C's `CalculateLogIP_fcmp()`.

This function calculates the inner product:
    ip = sum(QPFullWeight[qpidx] * pfM[qpidx])
and returns log(ip).

# Arguments
- `pf_m::Vector{ComplexF64}`: Pfaffian values [qp_idx] (1-based indexing)
- `qp_start::Int`: Start QP index (1-based)
- `qp_end::Int`: End QP index (1-based, exclusive)
- `data::ExpertModeData`: Expert Mode data containing qp_weights

# Returns
- `ComplexF64`: Logarithm of inner product log(ip)

# Reference
- C implementation: mVMC/src/mVMC/qp.c:90-107 (CalculateLogIP_fcmp)

# Note
- MPI communication is skipped (single-process implementation)
- qp_weights must be initialized via init_qp_weight!(data) before calling this function
"""
function calculate_log_ip_fcmp(
    pf_m::Vector{ComplexF64},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
)::ComplexF64
    # Check if qp_weights is initialized
    if data.qp_weights === nothing
        @error "Quantum projection weights not initialized. Call init_qp_weight!(data) first."
        return log(1e-100)  # Return log(0) equivalent to avoid errors
    end

    qp_num = qp_end - qp_start
    qp_full_weight::Vector{ComplexF64} = data.qp_weights.qp_full_weight

    # Validate indices
    if qp_start < 1 || qp_end > length(qp_full_weight) + 1 || qp_start >= qp_end
        @error "Invalid qp_start or qp_end: qp_start=$qp_start, qp_end=$qp_end, n_qp_full=$(length(qp_full_weight))"
        return log(1e-100)
    end

    # qp_end is exclusive (like C's qpEnd), so we need at least qp_end - 1 elements
    if length(pf_m) < qp_end - 1
        @error "pf_m array too short: length=$(length(pf_m)), required=$(qp_end - 1)"
        return log(1e-100)
    end

    # Calculate inner product: ip = sum(QPFullWeight[qpidx] * pfM[qpidx])
    # Note: C code uses 0-based indexing, Julia uses 1-based
    # C: QPFullWeight[qpidx+qpStart] * pfM[qpidx] for qpidx in 0:qpNum-1
    # Julia: qp_full_weight[qp_start + qpidx - 1] * pf_m[qp_start + qpidx] for qpidx in 1:qp_num
    ip = 0.0 + 0.0im
    for qpidx = 1:qp_num
        # C index: qpidx + qpStart - 1 (0-based)
        # Julia index: qp_start + qpidx - 1 (1-based)
        qp_idx = qp_start + qpidx - 1

        if qp_idx <= length(qp_full_weight) && qp_idx <= length(pf_m)
            ip += qp_full_weight[qp_idx] * pf_m[qp_idx]
        end
    end

    # Return log(ip) - equivalent to C's clog(ip)
    # Match C implementation: return clog(ip) directly
    # Note: C uses clog() which handles complex logarithm, including log(0) cases
    return log(ip)
end

"""
    calculate_ip_fcmp(pf_m::Vector{ComplexF64}, qp_start::Int, qp_end::Int,
                     data::ExpertModeData)::ComplexF64

Calculate inner product <phi|L|x> directly (without taking logarithm).
Equivalent to C's `CalculateIP_fcmp()`.

C implementation: mVMC/src/mVMC/qp.c:110-127 (CalculateIP_fcmp)

Note
- MPI communication is skipped (single-process implementation)
- qp_weights must be initialized via init_qp_weight!(data) before calling this function
"""
function calculate_ip_fcmp(
    pf_m::Vector{ComplexF64},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
)::ComplexF64
    # Check if qp_weights is initialized
    if data.qp_weights === nothing
        @error "Quantum projection weights not initialized. Call init_qp_weight!(data) first."
        return 0.0 + 0.0im
    end

    qp_num = qp_end - qp_start
    qp_full_weight::Vector{ComplexF64} = data.qp_weights.qp_full_weight

    # Validate indices
    if qp_start < 1 || qp_end > length(qp_full_weight) + 1 || qp_start >= qp_end
        @error "Invalid qp_start or qp_end: qp_start=$qp_start, qp_end=$qp_end, n_qp_full=$(length(qp_full_weight))"
        return 0.0 + 0.0im
    end

    # qp_end is exclusive (like C's qpEnd), so we need at least qp_end - 1 elements
    if length(pf_m) < qp_end - 1
        @error "pf_m array too short: length=$(length(pf_m)), required=$(qp_end - 1)"
        return 0.0 + 0.0im
    end

    # Calculate inner product: ip = sum(QPFullWeight[qpidx] * pfM[qpidx])
    # Note: C code uses 0-based indexing, Julia uses 1-based
    # C: QPFullWeight[qpidx+qpStart] * pfM[qpidx] for qpidx in 0:qpNum-1
    # Julia: qp_full_weight[qp_start + qpidx - 1] * pf_m[qp_start + qpidx] for qpidx in 1:qp_num
    ip = 0.0 + 0.0im
    for qpidx = 1:qp_num
        # C index: qpidx + qpStart - 1 (0-based)
        # Julia index: qp_start + qpidx - 1 (1-based)
        qp_idx = qp_start + qpidx - 1

        if qp_idx <= length(qp_full_weight) && qp_idx <= length(pf_m)
            ip += qp_full_weight[qp_idx] * pf_m[qp_idx]
        end
    end

    # Return ip directly (without taking logarithm)
    return ip
end

"""
    update_m_all!(ma::Int, s::Int, ele_idx::Vector{Int}, qp_start::Int, qp_end::Int,
                 data::ExpertModeData, state::VMCOptimizationState)

Update inverse matrix and Pfaffian after electron ma with spin s hops.
Equivalent to C's `UpdateMAll()`.

This function uses the Sherman-Morrison formula to efficiently update
the inverse matrix and Pfaffian after a single electron hop.

# Reference
- C implementation: mVMC/src/mVMC/pfupdate.c:119-216 (UpdateMAll, updateMAll_child)
"""
function update_m_all!(
    ma::Int,
    s::Int,
    ele_idx::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    msa = ma + s * n_elec  # 0-based electron index
    rsa = ele_idx[msa+1] + s * n_site  # 0-based site index (ele_idx is 1-based)

    # Get array references
    slater_elm = state.slater_matrix.slater_elm
    inv_m_arr = state.slater_matrix.inv_m
    pf_m = state.slater_matrix.pf_m

    # Work arrays for each QP index
    vec1 = zeros(ComplexF64, n_size)
    vec2 = zeros(ComplexF64, n_size)
    slt_vec = zeros(ComplexF64, n_size)
    BLAS = LinearAlgebra.BLAS

    # Process each QP index
    for qpidx = qp_start:(qp_end-1)
        if qpidx < 1 || qpidx > length(pf_m)
            continue
        end

        # C: sltE = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2
        slt_offset = (qpidx - 1) * n_site2 * n_site2
        # C: invM = InvM + qpidx*Nsize*Nsize
        inv_offset = (qpidx - 1) * n_size * n_size

        # NOTE:
        # - state.slater_matrix.inv_m is stored in C row-major layout:
        #     A[i,j] is at inv_m_arr[inv_offset + i*n_size + j + 1] (0-based i,j)
        # - Reshaping the contiguous block as a Julia matrix yields Aᵀ (column-major view).
        #   We keep all BLAS ops in that Aᵀ view to avoid transposes/copies.
        @views inv_m_t =
            reshape(inv_m_arr[(inv_offset+1):(inv_offset+n_size*n_size)], n_size, n_size) # == Aᵀ

        # Build slt_vec[j] = sltE[a][r_s(j)] for electron-index ordering j = msj
        @inbounds for msj = 0:(n_size-1)
            rsj = msj < n_elec ? ele_idx[msj+1] : (ele_idx[msj+1] + n_site)
            slt_vec[msj+1] = slater_elm[slt_offset+rsa*n_site2+rsj+1]
        end

        # vec1 = -Aᵀ * slt_vec  (matches original loop: vec1[i] += -A[j,i]*slt_vec[j])
        mul!(vec1, inv_m_t, slt_vec, -1.0 + 0.0im, 0.0 + 0.0im)

        # Update Pfaffian
        # PfM[qpidx] *= -vec1[msa]
        tmp = vec1[msa+1]
        pf_m[qpidx] *= -tmp
        inv_vec1_a = -1.0 / tmp

        # vec2[i] = InvM[a,i] * (-1/vec1[a])
        # Since inv_m_t == Aᵀ, we have InvM[a,i] == inv_m_t[i,a].
        @views copyto!(vec2, inv_m_t[:, msa+1])
        BLAS.scal!(inv_vec1_a, vec2)

        # Update InvM using Sherman-Morrison formula
        # invM[i][j] += vec1[i] * vec2[j] - vec1[j] * vec2[i]
        # Work in Aᵀ view (inv_m_t). Original update is on A; for B=Aᵀ we apply the transposed update:
        #   B += (vec1*vec2ᵀ - vec2*vec1ᵀ)ᵀ = -vec1*vec2ᵀ + vec2*vec1ᵀ
        BLAS.geru!(-1.0 + 0.0im, vec1, vec2, inv_m_t)
        BLAS.geru!(1.0 + 0.0im, vec2, vec1, inv_m_t)

        # Additional corrections from original code:
        #   A[:,a] -= vec2  <=>  B[a,:] -= vec2
        #   A[a,:] += vec2  <=>  B[:,a] += vec2
        @inbounds for j = 1:n_size
            inv_m_t[msa+1, j] -= vec2[j]
        end
        @views BLAS.axpy!(1.0 + 0.0im, vec2, inv_m_t[:, msa+1])
    end
end

"""
    calculate_new_pf_m_two2!(ma::Int, s::Int, mb::Int, t::Int, pf_m_new::Vector{ComplexF64},
                            ele_idx::Vector{Int}, qp_start::Int, qp_end::Int,
                            data::ExpertModeData, state::VMCOptimizationState)

Calculate new Pfaffian after two electrons exchange positions.
The ma-th electron with spin s hops, then the mb-th electron with spin t hops.
Equivalent to C's `CalculateNewPfMTwo2_fcmp()`.

# Reference
- C implementation: mVMC/src/mVMC/pfupdate_two_fcmp.c:73-184 (CalculateNewPfMTwo2_fcmp, calculateNewPfMTwo_child_fcmp)
"""
function calculate_new_pf_m_two2!(
    ma::Int,
    s::Int,
    mb::Int,
    t::Int,
    pf_m_new::Vector{ComplexF64},
    ele_idx::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    msa = ma + s * n_elec  # 0-based electron index
    msb = mb + t * n_elec  # 0-based electron index

    # If same electron, use single hop calculation
    if msa == msb
        calculate_new_pf_m2!(mb, t, pf_m_new, ele_idx, qp_start, qp_end, data, state)
        return
    end

    # Work arrays
    vec_a = zeros(ComplexF64, n_size)
    vec_b = zeros(ComplexF64, n_size)

    # Get array references
    slater_elm = state.slater_matrix.slater_elm
    inv_m_arr = state.slater_matrix.inv_m
    pf_m = state.slater_matrix.pf_m

    # Process each QP index
    # C: for(qpidx=0;qpidx<qpNum;qpidx++) where qpNum = qpEnd - qpStart
    # C uses local index qpidx for InvM, but global index (qpidx+qpStart) for SlaterElm
    # Julia: qp_start is 1-based, so we need to convert to local for InvM
    qp_num = qp_end - qp_start
    for local_qpidx = 0:(qp_num-1)
        qpidx = qp_start + local_qpidx  # Global index (1-based)
        if qpidx < 1 || qpidx > length(pf_m)
            continue
        end

        rsa = ele_idx[msa+1] + s * n_site  # 0-based site index
        rsb = ele_idx[msb+1] + t * n_site  # 0-based site index

        # C: sltE = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2
        # In C, qpidx is local and (qpidx+qpStart) is global 0-based
        # In Julia, qpidx is global 1-based, so (qpidx-1) is global 0-based
        slt_offset = (qpidx - 1) * n_site2 * n_site2

        # InvM offset (global index) - state.slater_matrix.inv_m stores all QP indices
        inv_offset = (qpidx - 1) * n_size * n_size

        # Calculate vec_a[i] = sltE[rsa][rsi] and vec_b[i] = sltE[rsb][rsi]
        # C: vec_a[msi] = sltE_a[rsi] = sltE[rsa*Nsite2 + rsi]
        for msi = 0:(n_size-1)
            if msi < n_elec
                rsi = ele_idx[msi+1]  # up-spin (0-based)
            else
                rsi = ele_idx[msi+1] + n_site  # down-spin (0-based)
            end
            # C: sltE_a[rsi] = SlaterElm[slt_offset + rsa*Nsite2 + rsi]
            vec_a[msi+1] = slater_elm[slt_offset+rsa*n_site2+rsi+1]
            vec_b[msi+1] = slater_elm[slt_offset+rsb*n_site2+rsi+1]
        end
        vec_ba = vec_b[msa+1]

        # Get invM_a, invM_b, invM_ab
        # C: invM_a = invM + msa*Nsize, invM_b = invM + msb*Nsize
        # C: invM_a[msi] = InvM[inv_offset + msa*Nsize + msi]

        # Calculate p_a, p_b, q_a, q_b
        p_a = p_b = q_a = q_b = 0.0 + 0.0im
        for msi = 0:(n_size-1)
            # C: invM_a[msi] = InvM[inv_offset + msa*Nsize + msi]
            inv_m_a_msi = inv_m_arr[inv_offset+msa*n_size+msi+1]
            inv_m_b_msi = inv_m_arr[inv_offset+msb*n_size+msi+1]

            p_a += inv_m_a_msi * vec_a[msi+1]
            p_b += inv_m_b_msi * vec_a[msi+1]
            q_a += inv_m_a_msi * vec_b[msi+1]
            q_b += inv_m_b_msi * vec_b[msi+1]
        end

        # invM_ab = invM_a[msb]
        inv_m_ab = inv_m_arr[inv_offset+msa*n_size+msb+1]

        # Calculate bMa = sum_i vec_b[i] * (sum_j invM[i][j] * vec_a[j])
        bMa = 0.0 + 0.0im
        for msi = 0:(n_size-1)
            tmp = 0.0 + 0.0im
            for msj = 0:(n_size-1)
                # C: invM_i[msj] = InvM[inv_offset + msi*Nsize + msj]
                tmp += inv_m_arr[inv_offset+msi*n_size+msj+1] * vec_a[msj+1]
            end
            bMa += vec_b[msi+1] * tmp
        end

        # Calculate ratio = PfMNew / PfMOld
        ratio = inv_m_ab * vec_ba + inv_m_ab * bMa + p_a * q_b - p_b * q_a

        # Update pfMNew
        pf_m_new[qpidx] = ratio * pf_m[qpidx]
    end
end

"""
    calculate_new_pf_m_two_fsz!(ma::Int, s::Int, mb::Int, t::Int, pf_m_new::Vector{ComplexF64},
                                ele_idx::Vector{Int}, ele_spn::Vector{Int},
                                qp_start::Int, qp_end::Int,
                                data::ExpertModeData, state::VMCOptimizationState)

Calculate new Pfaffian for FSZ mode after two electrons hop.
The ma-th electron hops to ele_idx[ma] with spin s.
The mb-th electron hops to ele_idx[mb] with spin t.

In FSZ mode:
- msa = ma, msb = mb (no spin offset)
- rsa = ele_idx[ma] + s * n_site, rsb = ele_idx[mb] + t * n_site
- For other electrons: rsi = ele_idx[msi] + ele_spn[msi] * n_site

# Arguments
- `ma`, `mb`: Electron indices (0-based)
- `s`, `t`: New spins of the electrons (0 = up, 1 = down)
- `pf_m_new`: Output array for new Pfaffian values
- `ele_idx`: Electron site indices (1-based array, 0-based values)
- `ele_spn`: Electron spin indices (1-based array, values 0 or 1)
- `qp_start`, `qp_end`: QP index range (1-based)

# Reference
- C implementation: mVMC/src/mVMC/pfupdate_two_fsz.c:50-189 (CalculateNewPfMTwo_fsz, calculateNewPfMTwo_child_fsz)
"""
function calculate_new_pf_m_two_fsz!(
    ma::Int,
    s::Int,
    mb::Int,
    t::Int,
    pf_m_new::Vector{ComplexF64},
    ele_idx::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site

    # FSZ: msa = ma, msb = mb (no spin offset)
    msa = ma
    msb = mb

    # If same electron, use single-electron update
    if msa == msb
        calculate_new_pf_m2_fsz!(
            mb,
            t,
            pf_m_new,
            ele_idx,
            ele_spn,
            qp_start,
            qp_end,
            data,
            state,
        )
        return
    end

    # FSZ: rsa, rsb use new spins s, t
    rsa = ele_idx[msa+1] + s * n_site
    rsb = ele_idx[msb+1] + t * n_site

    # Get array references
    slater_elm = state.slater_matrix.slater_elm
    inv_m_arr = state.slater_matrix.inv_m
    pf_m = state.slater_matrix.pf_m

    # Allocate work vectors
    vec_a = zeros(ComplexF64, n_size)
    vec_b = zeros(ComplexF64, n_size)

    # Process each QP index
    qp_num = qp_end - qp_start
    for local_qpidx = 0:(qp_num-1)
        qpidx = qp_start + local_qpidx
        if qpidx < 1 || qpidx > length(pf_m)
            continue
        end

        # SlaterElm offset (global index)
        slt_offset = (qpidx - 1) * n_site2 * n_site2

        # InvM offset (global index) - state.slater_matrix.inv_m stores all QP indices
        inv_offset = (qpidx - 1) * n_size * n_size

        # Fill vec_a and vec_b from SlaterElm
        # FSZ: rsi = ele_idx[msi] + ele_spn[msi] * n_site
        for msi = 0:(n_size-1)
            rsi = ele_idx[msi+1] + ele_spn[msi+1] * n_site
            vec_a[msi+1] = slater_elm[slt_offset+rsa*n_site2+rsi+1]
            vec_b[msi+1] = slater_elm[slt_offset+rsb*n_site2+rsi+1]
        end
        vec_ba = vec_b[msa+1]

        # Calculate p_a, p_b, q_a, q_b
        p_a = p_b = q_a = q_b = 0.0 + 0.0im
        for msi = 0:(n_size-1)
            inv_m_a_msi = inv_m_arr[inv_offset+msa*n_size+msi+1]
            inv_m_b_msi = inv_m_arr[inv_offset+msb*n_size+msi+1]

            p_a += inv_m_a_msi * vec_a[msi+1]
            p_b += inv_m_b_msi * vec_a[msi+1]
            q_a += inv_m_a_msi * vec_b[msi+1]
            q_b += inv_m_b_msi * vec_b[msi+1]
        end

        # invM_ab = invM_a[msb]
        inv_m_ab = inv_m_arr[inv_offset+msa*n_size+msb+1]

        # Calculate bMa = sum_i vec_b[i] * (sum_j invM[i][j] * vec_a[j])
        bMa = 0.0 + 0.0im
        for msi = 0:(n_size-1)
            tmp = 0.0 + 0.0im
            for msj = 0:(n_size-1)
                tmp += inv_m_arr[inv_offset+msi*n_size+msj+1] * vec_a[msj+1]
            end
            bMa += vec_b[msi+1] * tmp
        end

        # Calculate ratio = PfMNew / PfMOld
        ratio = inv_m_ab * vec_ba + inv_m_ab * bMa + p_a * q_b - p_b * q_a

        # Update pfMNew
        pf_m_new[qpidx] = ratio * pf_m[qpidx]
    end
end

"""
    update_m_all_two!(ma::Int, s::Int, mb::Int, t::Int, ra_old::Int, rb_old::Int,
                     ele_idx::Vector{Int}, qp_start::Int, qp_end::Int,
                     data::ExpertModeData, state::VMCOptimizationState)

Update inverse matrix and Pfaffian after two electrons exchange positions.
The ma-th electron with spin s hops from ra_old, and the mb-th electron with spin t hops from rb_old.
Equivalent to C's `UpdateMAllTwo_fcmp()`.

# Reference
- C implementation: mVMC/src/mVMC/pfupdate_two_fcmp.c:189-344 (UpdateMAllTwo_fcmp, updateMAllTwo_child_fcmp)
"""
function update_m_all_two!(
    ma::Int,
    s::Int,
    mb::Int,
    t::Int,
    ra_old::Int,
    rb_old::Int,
    ele_idx::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    msa = ma + s * n_elec  # 0-based electron index
    msb = mb + t * n_elec  # 0-based electron index
    rsa = ele_idx[msa+1] + s * n_site  # 0-based site index (new position)
    rsb = ele_idx[msb+1] + t * n_site  # 0-based site index (new position)
    rsa_old = ra_old + s * n_site  # 0-based site index (old position)
    rsb_old = rb_old + t * n_site  # 0-based site index (old position)

    # Get array references
    slater_elm = state.slater_matrix.slater_elm
    inv_m_arr = state.slater_matrix.inv_m
    pf_m = state.slater_matrix.pf_m

    # Work arrays
    vec_p = zeros(ComplexF64, n_size)
    vec_q = zeros(ComplexF64, n_size)
    vec_s = zeros(ComplexF64, n_size)
    vec_t = zeros(ComplexF64, n_size)

    # Process each QP index
    for qpidx = qp_start:(qp_end-1)
        if qpidx < 1 || qpidx > length(pf_m)
            continue
        end

        # C: sltE = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2
        slt_offset = (qpidx - 1) * n_site2 * n_site2
        # C: invM = InvM + qpidx*Nsize*Nsize
        inv_offset = (qpidx - 1) * n_size * n_size

        # Get m_old_ab = sltE[rsa_old][rsb_old]
        # C: sltE[rsa_old*Nsite2 + rsb_old]
        m_old_ab = slater_elm[slt_offset+rsa_old*n_site2+rsb_old+1]

        # Get invM_ab = invM[msa][msb]
        inv_m_ab = inv_m_arr[inv_offset+msa*n_size+msb+1]

        # Initialize vecP, vecQ
        fill!(vec_p, 0.0 + 0.0im)
        fill!(vec_q, 0.0 + 0.0im)

        # vecS[i] = sltE[rsa][rsi], vecT[i] = sltE[rsb][rsi]
        for msi = 0:(n_size-1)
            if msi < n_elec
                rsi = ele_idx[msi+1]  # up-spin
            else
                rsi = ele_idx[msi+1] + n_site  # down-spin
            end
            # C: sltE_a[rsi] = sltE[rsa*Nsite2 + rsi]
            vec_s[msi+1] = slater_elm[slt_offset+rsa*n_site2+rsi+1]
            vec_t[msi+1] = slater_elm[slt_offset+rsb*n_site2+rsi+1]
        end
        # Set vecS[b] = mOld_ab (old (a,b) element)
        vec_s[msb+1] = m_old_ab

        # Calculate vecP[i] = sum_j invM[i][j] * vecS[j]
        # Calculate vecQ[i] = sum_j invM[i][j] * vecT[j]
        for msi = 0:(n_size-1)
            for msj = 0:(n_size-1)
                # C: invM_i[msj] = invM[msi*Nsize + msj]
                inv_m_ij = inv_m_arr[inv_offset+msi*n_size+msj+1]
                vec_p[msi+1] += inv_m_ij * vec_s[msj+1]
                vec_q[msi+1] += inv_m_ij * vec_t[msj+1]
            end
        end

        # Update Pfaffian
        bMa = 0.0 + 0.0im
        for msi = 0:(n_size-1)
            bMa += vec_t[msi+1] * vec_p[msi+1]
        end
        ratio =
            inv_m_ab * vec_t[msa+1] + inv_m_ab * bMa + vec_p[msa+1] * vec_q[msb+1] -
            vec_p[msb+1] * vec_q[msa+1]
        pf_m[qpidx] *= ratio

        # Set coefficients
        a = -vec_p[msa+1]
        b = vec_p[msb+1]
        c = vec_q[msa+1]
        d = -vec_q[msb+1]
        e = -bMa - vec_t[msa+1]
        # f = invM[a][b]
        f = inv_m_arr[inv_offset+msa*n_size+msb+1]

        det = a * d - b * c - e * f
        inv_det = 1.0 / det

        # Calculate vecS[i] = invM[a][i] / det, vecT[i] = invM[b][i] / det
        for msi = 0:(n_size-1)
            vec_s[msi+1] = inv_det * inv_m_arr[inv_offset+msa*n_size+msi+1]
            vec_t[msi+1] = inv_det * inv_m_arr[inv_offset+msb*n_size+msi+1]
        end

        # Update InvM
        for msi = 0:(n_size-1)
            msi_next = msi + 1
            p_i = vec_p[msi_next]
            q_i = vec_q[msi_next]
            s_i = vec_s[msi_next]
            t_i = vec_t[msi_next]

            for msj = 0:(n_size-1)
                msj_next = msj + 1
                p_j = vec_p[msj_next]
                q_j = vec_q[msj_next]
                s_j = vec_s[msj_next]
                t_j = vec_t[msj_next]

                # C: invM[msi*Nsize + msj]
                inv_m_arr[inv_offset+msi*n_size+msj+1] +=
                    a * (q_i * t_j - q_j * t_i) +
                    b * (q_i * s_j - q_j * s_i) +
                    c * (p_i * t_j - p_j * t_i) +
                    d * (p_i * s_j - p_j * s_i) +
                    e * det * (s_i * t_j - s_j * t_i) +
                    f * inv_det * (p_i * q_j - q_i * p_j)
            end
            inv_m_arr[inv_offset+msi*n_size+msa+1] += -c * t_i - d * s_i - f * inv_det * q_i
            inv_m_arr[inv_offset+msi*n_size+msb+1] += -a * t_i - b * s_i + f * inv_det * p_i
        end

        for msj = 0:(n_size-1)
            msj_next = msj + 1
            p_j = vec_p[msj_next]
            q_j = vec_q[msj_next]
            s_j = vec_s[msj_next]
            t_j = vec_t[msj_next]

            inv_m_arr[inv_offset+msa*n_size+msj+1] += c * t_j + d * s_j + f * inv_det * q_j
            inv_m_arr[inv_offset+msb*n_size+msj+1] += a * t_j + b * s_j - f * inv_det * p_j
        end

        inv_m_arr[inv_offset+msa*n_size+msb+1] += f * inv_det
        inv_m_arr[inv_offset+msb*n_size+msa+1] -= f * inv_det
    end
end

"""
    update_m_all_two_real!(ma::Int, s::Int, mb::Int, t::Int, ra_old::Int, rb_old::Int,
                           ele_idx::Vector{Int}, qp_start::Int, qp_end::Int,
                           data::ExpertModeData, state::VMCOptimizationState)

Update inverse matrix and Pfaffian after two electrons exchange positions (real version).
The ma-th electron with spin s hops from ra_old, and the mb-th electron with spin t hops from rb_old.
Equivalent to C's `UpdateMAllTwo_real()`.

This is an O(N²) Sherman-Morrison update, much faster than full O(N³) recalculation.

# Reference
- C implementation: mVMC/src/mVMC/pfupdate_two_real.c:189-344 (UpdateMAllTwo_real, updateMAllTwo_child_real)
"""
function update_m_all_two_real!(
    ma::Int,
    s::Int,
    mb::Int,
    t::Int,
    ra_old::Int,
    rb_old::Int,
    ele_idx::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    msa = ma + s * n_elec  # 0-based electron index
    msb = mb + t * n_elec  # 0-based electron index
    rsa = ele_idx[msa+1] + s * n_site  # 0-based site index (new position)
    rsb = ele_idx[msb+1] + t * n_site  # 0-based site index (new position)
    rsa_old = ra_old + s * n_site  # 0-based site index (old position)
    # Note: C code has a bug here using raOld for both, but we follow the algorithm
    rsb_old = ra_old + t * n_site  # 0-based site index (matching C: rsbOld = raOld + t*Nsite)

    # Get array references
    slater_elm_real = state.slater_matrix.slater_elm_real
    inv_m_real = state.slater_matrix.inv_m_real
    pf_m_real = state.slater_matrix.pf_m_real

    # Work arrays
    vec_p = zeros(Float64, n_size)
    vec_q = zeros(Float64, n_size)
    vec_s = zeros(Float64, n_size)
    vec_t = zeros(Float64, n_size)

    # Process each QP index
    for qpidx = qp_start:(qp_end-1)
        if qpidx < 1 || qpidx > length(pf_m_real)
            continue
        end

        # C: sltE = SlaterElm_real + (qpidx+qpStart)*Nsite2*Nsite2
        slt_offset = (qpidx - 1) * n_site2 * n_site2
        # C: invM = InvM_real + qpidx*Nsize*Nsize
        inv_offset = (qpidx - 1) * n_size * n_size

        # Get m_old_ab = sltE[rsa_old][rsb_old]
        m_old_ab = slater_elm_real[slt_offset+rsa_old*n_site2+rsb_old+1]

        # Get invM_ab = invM[msa][msb]
        inv_m_ab = inv_m_real[inv_offset+msa*n_size+msb+1]

        # Initialize vecP, vecQ
        fill!(vec_p, 0.0)
        fill!(vec_q, 0.0)

        # vecS[i] = sltE[rsa][rsi], vecT[i] = sltE[rsb][rsi]
        @inbounds for msi = 0:(n_size-1)
            if msi < n_elec
                rsi = ele_idx[msi+1]  # up-spin
            else
                rsi = ele_idx[msi+1] + n_site  # down-spin
            end
            vec_s[msi+1] = slater_elm_real[slt_offset+rsa*n_site2+rsi+1]
            vec_t[msi+1] = slater_elm_real[slt_offset+rsb*n_site2+rsi+1]
        end
        # Set vecS[b] = mOld_ab (old (a,b) element)
        vec_s[msb+1] = m_old_ab

        # Calculate vecP[i] = sum_j invM[i][j] * vecS[j]
        # Calculate vecQ[i] = sum_j invM[i][j] * vecT[j]
        @turbo for msi = 0:(n_size-1)
            msi_next = msi + 1
            for msj = 0:(n_size-1)
                msj_next = msj + 1
                inv_m_ij = inv_m_real[inv_offset+msi*n_size+msj+1]
                vec_p[msi_next] += inv_m_ij * vec_s[msj_next]
                vec_q[msi_next] += inv_m_ij * vec_t[msj_next]
            end
        end

        # Update Pfaffian
        bMa = 0.0
        @inbounds for msi = 0:(n_size-1)
            bMa += vec_t[msi+1] * vec_p[msi+1]
        end
        ratio =
            inv_m_ab * vec_t[msa+1] + inv_m_ab * bMa + vec_p[msa+1] * vec_q[msb+1] -
            vec_p[msb+1] * vec_q[msa+1]
        pf_m_real[qpidx] *= ratio

        # Set coefficients
        a = -vec_p[msa+1]
        b = vec_p[msb+1]
        c = vec_q[msa+1]
        d = -vec_q[msb+1]
        e = -bMa - vec_t[msa+1]
        f = inv_m_real[inv_offset+msa*n_size+msb+1]

        det = a * d - b * c - e * f
        inv_det = 1.0 / det

        # Calculate vecS[i] = invM[a][i] / det, vecT[i] = invM[b][i] / det
        @inbounds for msi = 0:(n_size-1)
            vec_s[msi+1] = inv_det * inv_m_real[inv_offset+msa*n_size+msi+1]
            vec_t[msi+1] = inv_det * inv_m_real[inv_offset+msb*n_size+msi+1]
        end

        # Update InvM
        @turbo for msi = 0:(n_size-1)
            msi_next = msi + 1
            p_i = vec_p[msi_next]
            q_i = vec_q[msi_next]
            s_i = vec_s[msi_next]
            t_i = vec_t[msi_next]

            for msj = 0:(n_size-1)
                msj_next = msj + 1
                p_j = vec_p[msj_next]
                q_j = vec_q[msj_next]
                s_j = vec_s[msj_next]
                t_j = vec_t[msj_next]

                inv_m_real[inv_offset+msi*n_size+msj+1] +=
                    a * (q_i * t_j - q_j * t_i) +
                    b * (q_i * s_j - q_j * s_i) +
                    c * (p_i * t_j - p_j * t_i) +
                    d * (p_i * s_j - p_j * s_i) +
                    e * det * (s_i * t_j - s_j * t_i) +
                    f * inv_det * (p_i * q_j - q_i * p_j)
            end
            inv_m_real[inv_offset+msi*n_size+msa+1] +=
                -c * t_i - d * s_i - f * inv_det * q_i
            inv_m_real[inv_offset+msi*n_size+msb+1] +=
                -a * t_i - b * s_i + f * inv_det * p_i
        end

        @inbounds for msj = 0:(n_size-1)
            p_j = vec_p[msj+1]
            q_j = vec_q[msj+1]
            s_j = vec_s[msj+1]
            t_j = vec_t[msj+1]

            inv_m_real[inv_offset+msa*n_size+msj+1] += c * t_j + d * s_j + f * inv_det * q_j
            inv_m_real[inv_offset+msb*n_size+msj+1] += a * t_j + b * s_j - f * inv_det * p_j
        end

        inv_m_real[inv_offset+msa*n_size+msb+1] += f * inv_det
        inv_m_real[inv_offset+msb*n_size+msa+1] -= f * inv_det
    end
end

# ============================================================================
# Burn Sample Functions
# ============================================================================

"""
    copy_from_burn_sample!(tmp_ele_idx::Vector{Int}, tmp_ele_cfg::Vector{Int},
                          tmp_ele_num::Vector{Int}, tmp_ele_proj_cnt::Vector{Int},
                          state::VMCOptimizationState)

Copy electron configuration from burn sample to temporary arrays.
Equivalent to C's `copyFromBurnSample()`.
"""
function copy_from_burn_sample!(
    tmp_ele_idx::Vector{Int},
    tmp_ele_cfg::Vector{Int},
    tmp_ele_num::Vector{Int},
    tmp_ele_proj_cnt::Vector{Int},
    state::VMCOptimizationState,
)
    # C implementation: BurnEleIdx contains all data in one array
    # BurnEleIdx[0:Nsize-1] = eleIdx
    # BurnEleCfg = BurnEleIdx + Nsize (points to BurnEleIdx[Nsize])
    # BurnEleNum = BurnEleCfg + 2*Nsite (points to BurnEleIdx[Nsize+2*Nsite])
    # BurnEleProjCnt = BurnEleNum + 2*Nsite (points to BurnEleIdx[Nsize+2*Nsite+2*Nsite])
    # In Julia, we use separate arrays, so copy from burn_ele_idx to tmp arrays
    n_size = length(tmp_ele_idx)
    n_site2 = length(tmp_ele_cfg)
    n_proj = length(tmp_ele_proj_cnt)

    # Copy from burn_ele_idx (which contains all data in C, but we use separate arrays)
    # For Julia, burn_ele_idx is a combined storage array
    burn_combined = state.electron_config.burn_ele_idx
    if length(burn_combined) >= n_size + n_site2 + n_site2 + n_proj
        # Copy eleIdx
        for i = 1:n_size
            tmp_ele_idx[i] = burn_combined[i]
        end
        # Copy eleCfg
        for i = 1:n_site2
            tmp_ele_cfg[i] = burn_combined[n_size+i]
        end
        # Copy eleNum
        for i = 1:n_site2
            tmp_ele_num[i] = burn_combined[n_size+n_site2+i]
        end
        # Copy eleProjCnt
        for i = 1:n_proj
            tmp_ele_proj_cnt[i] = burn_combined[n_size+n_site2+n_site2+i]
        end
    else
        # Fallback: use separate arrays if combined storage is not available
        copy!(
            tmp_ele_idx,
            state.electron_config.burn_ele_idx[1:min(
                n_size,
                length(state.electron_config.burn_ele_idx),
            )],
        )
        if !isempty(state.electron_config.burn_ele_cfg)
            copy!(tmp_ele_cfg, state.electron_config.burn_ele_cfg)
        end
        if !isempty(state.electron_config.burn_ele_num)
            copy!(tmp_ele_num, state.electron_config.burn_ele_num)
        end
        if !isempty(state.electron_config.burn_ele_proj_cnt)
            copy!(tmp_ele_proj_cnt, state.electron_config.burn_ele_proj_cnt)
        end
    end
end

"""
    copy_to_burn_sample!(tmp_ele_idx::Vector{Int}, tmp_ele_cfg::Vector{Int},
                        tmp_ele_num::Vector{Int}, tmp_ele_proj_cnt::Vector{Int},
                        state::VMCOptimizationState)

Copy electron configuration from temporary arrays to burn sample.
Equivalent to C's `copyToBurnSample()`.
"""
function copy_to_burn_sample!(
    tmp_ele_idx::Vector{Int},
    tmp_ele_cfg::Vector{Int},
    tmp_ele_num::Vector{Int},
    tmp_ele_proj_cnt::Vector{Int},
    state::VMCOptimizationState,
)
    # C implementation: BurnEleIdx contains all data in one array
    # Copy to burn_ele_idx (combined storage array)
    n_size = length(tmp_ele_idx)
    n_site2 = length(tmp_ele_cfg)
    n_proj = length(tmp_ele_proj_cnt)

    burn_combined = state.electron_config.burn_ele_idx
    if length(burn_combined) >= n_size + n_site2 + n_site2 + n_proj
        # Copy eleIdx
        for i = 1:n_size
            burn_combined[i] = tmp_ele_idx[i]
        end
        # Copy eleCfg
        for i = 1:n_site2
            burn_combined[n_size+i] = tmp_ele_cfg[i]
        end
        # Copy eleNum
        for i = 1:n_site2
            burn_combined[n_size+n_site2+i] = tmp_ele_num[i]
        end
        # Copy eleProjCnt
        for i = 1:n_proj
            burn_combined[n_size+n_site2+n_site2+i] = tmp_ele_proj_cnt[i]
        end
    else
        # Fallback: use separate arrays if combined storage is not available
        if length(state.electron_config.burn_ele_idx) >= n_size
            copy!(state.electron_config.burn_ele_idx[1:n_size], tmp_ele_idx)
        end
        if !isempty(state.electron_config.burn_ele_cfg) &&
           length(state.electron_config.burn_ele_cfg) >= n_site2
            copy!(state.electron_config.burn_ele_cfg, tmp_ele_cfg)
        end
        if !isempty(state.electron_config.burn_ele_num) &&
           length(state.electron_config.burn_ele_num) >= n_site2
            copy!(state.electron_config.burn_ele_num, tmp_ele_num)
        end
        if !isempty(state.electron_config.burn_ele_proj_cnt) &&
           length(state.electron_config.burn_ele_proj_cnt) >= n_proj
            copy!(state.electron_config.burn_ele_proj_cnt, tmp_ele_proj_cnt)
        end
    end
end

# ============================================================================
# Main Sampling Functions
# ============================================================================

"""
    vmc_make_sample!(data::ExpertModeData, state::VMCOptimizationState, rng::AbstractRNG)

VMC sampling (complex version, sz-conserved).
Equivalent to C's `VMCMakeSample()`.

This function performs Monte Carlo sampling using the Metropolis-Hastings algorithm
with electron hopping and exchange updates. It uses PfaPack.jl for Pfaffian
calculations.

# Implementation details
- Uses PfaPack.jl's `calculate_m_all_fcmp!` for Pfaffian and inverse matrix calculations
- Implements Metropolis-Hastings algorithm with electron hopping and exchange updates
- Periodically recalculates Pfaffian from scratch to maintain numerical stability
- Supports both hopping and exchange update types based on `n_ex_update_path` parameter

# Note
- Local spin flip updates are not yet implemented (TODO)
- Burn-in sample storage is not yet implemented (TODO)
"""
function vmc_make_sample!(
    data::ExpertModeData,
    state::VMCOptimizationState,
    rng::AbstractRNG = Random.GLOBAL_RNG,
    c_timer::CTimer = CTIMER_DISABLED,
)

    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_proj = length(data.gutzwiller_terms) + length(data.jastrow_terms)
    n_vmc_sample = data.modpara.nvmc_sample
    n_vmc_warmup = data.modpara.nvmc_warmup
    n_vmc_interval = data.modpara.nvmc_interval
    # Get NExUpdatePath from modpara (0: hopping only, 1: hopping+exchange, 2: exchange only, 3: KondoGC)
    n_ex_update_path = data.modpara.nex_update_path

    # Debug: Verify state.electron_config.ele_idx size
    n_size = 2 * n_elec
    expected_ele_idx_size = n_vmc_sample * n_size
    actual_ele_idx_size = length(state.electron_config.ele_idx)
    @debug "vmc_make_sample!: n_vmc_sample=$n_vmc_sample, n_size=$n_size, expected_ele_idx_size=$expected_ele_idx_size, actual_ele_idx_size=$actual_ele_idx_size"
    if actual_ele_idx_size != expected_ele_idx_size
        @error "vmc_make_sample!: ele_idx size mismatch! expected=$expected_ele_idx_size, actual=$actual_ele_idx_size. This will cause BoundsError when saving samples."
    end

    # Get temporary arrays
    tmp_ele_idx = state.electron_config.tmp_ele_idx
    tmp_ele_cfg = state.electron_config.tmp_ele_cfg
    tmp_ele_num = state.electron_config.tmp_ele_num
    tmp_ele_proj_cnt = state.electron_config.tmp_ele_proj_cnt

    # Use cached loc_spn from workspace (initialize if needed)
    ws = state.workspace
    if all(x -> x == 0, ws.loc_spn) && !isempty(data.locspin_terms)
        init_loc_spn!(state, data)
    end
    loc_spn = ws.loc_spn

    # Initialize sample
    # burn_flag is passed as a parameter (managed by vmc_para_opt!)
    # Ensure counter array is large enough
    if length(state.electron_config.counter) < 11
        old_len = length(state.electron_config.counter)
        resize!(state.electron_config.counter, 11)
        # CRITICAL: Initialize new elements to zero (resize! leaves them uninitialized!)
        for i = (old_len+1):11
            state.electron_config.counter[i] = 0
        end
    end
    burn_flag = state.electron_config.counter[11] != 0  # Use counter[11] as burn_flag storage
    if !burn_flag
        result = make_initial_sample!(
            tmp_ele_idx,
            tmp_ele_cfg,
            tmp_ele_num,
            tmp_ele_proj_cnt,
            data,
            state,
            rng,
        )
        if result != 0
            @error "Failed to generate initial sample"
            return
        end

        @debug "After make_initial_sample!: tmp_ele_idx[1:10] = $(tmp_ele_idx[1:min(10, length(tmp_ele_idx))])"
        if all(x -> x == 0 || x == -1, tmp_ele_idx)
            @error "tmp_ele_idx is invalid after make_initial_sample! (all zeros or -1)"
            @error "  tmp_ele_idx[1:min(20, length(tmp_ele_idx))] = $(tmp_ele_idx[1:min(20, length(tmp_ele_idx))])"
            return
        end
    else
        # Copy from burn sample (from previous step)
        copy_from_burn_sample!(
            tmp_ele_idx,
            tmp_ele_cfg,
            tmp_ele_num,
            tmp_ele_proj_cnt,
            state,
        )
    end

    # Calculate Pfaffian (CalculateMAll_fcmp)
    n_qp_full = length(state.slater_matrix.pf_m)
    qp_start = 1
    qp_end = n_qp_full + 1

    # Retry loop for initial Pfaffian calculation (similar to C implementation)
    max_retries = 100
    retry_count = 0
    result = 1
    while result != 0 && retry_count < max_retries
        # Debug: Check tmp_ele_idx before calculate_m_all_fcmp!
        if retry_count == 0 && all(x -> x == 0 || x == -1, tmp_ele_idx)
            @error "tmp_ele_idx is invalid before calculate_m_all_fcmp! (all zeros or -1). First 10 values: $(tmp_ele_idx[1:min(10, length(tmp_ele_idx))])"
        end

        result = calculate_m_all_fcmp!(tmp_ele_idx, qp_start, qp_end, data, state)
        if result != 0
            # Regenerate sample if Pfaffian calculation fails
            retry_count += 1
            if retry_count <= 3 || retry_count % 10 == 0
                @warn "calculate_m_all_fcmp! failed with result=$result, retry_count=$retry_count. Regenerating sample..."
            end
            result_init = make_initial_sample!(
                tmp_ele_idx,
                tmp_ele_cfg,
                tmp_ele_num,
                tmp_ele_proj_cnt,
                data,
                state,
                rng,
            )
            if result_init != 0
                @error "Failed to regenerate initial sample after Pfaffian calculation failure"
                return
            end
        end
    end

    if result != 0
        @error "Failed to calculate initial Pfaffian after $max_retries retries, result=$result. tmp_ele_idx first 10 values: $(tmp_ele_idx[1:min(10, length(tmp_ele_idx))])"
        return
    end

    # Calculate initial log inner product
    log_ip_old = calculate_log_ip_fcmp(state.slater_matrix.pf_m, qp_start, qp_end, data)
    use_rbm = has_rbm_terms(data)
    rbm_cnt_old = use_rbm ? make_rbm_cnt(tmp_ele_num, data) : ComplexF64[]
    rbm_cnt_new = use_rbm ? similar(rbm_cnt_old) : ComplexF64[]

    if !isfinite(real(log_ip_old)) || !isfinite(imag(log_ip_old))
        @warn "Initial log inner product is not finite, regenerating sample"
        result = make_initial_sample!(
            tmp_ele_idx,
            tmp_ele_cfg,
            tmp_ele_num,
            tmp_ele_proj_cnt,
            data,
            state,
            rng,
        )
        if result != 0
            @error "Failed to regenerate initial sample"
            return
        end
        result = calculate_m_all_fcmp!(tmp_ele_idx, qp_start, qp_end, data, state)
        if result != 0
            @error "Failed to recalculate Pfaffian"
            return
        end
        log_ip_old = calculate_log_ip_fcmp(state.slater_matrix.pf_m, qp_start, qp_end, data)
        rbm_cnt_old = use_rbm ? make_rbm_cnt(tmp_ele_num, data) : ComplexF64[]
        if use_rbm && length(rbm_cnt_new) != length(rbm_cnt_old)
            rbm_cnt_new = similar(rbm_cnt_old)
        end
        burn_flag = false
    end

    # Main sampling loop
    n_out_step = burn_flag ? n_vmc_sample + 1 : n_vmc_warmup + n_vmc_sample
    n_in_step = n_vmc_interval * n_site

    # Reset counters
    fill!(state.electron_config.counter, 0)

    n_accept = 0
    # Use pre-allocated workspace arrays instead of allocating new ones
    proj_cnt_new = ws.proj_cnt_new
    pf_m_new = ws.pf_m_new
    fill!(proj_cnt_new, 0)  # Reset workspace arrays
    fill!(pf_m_new, 0.0 + 0.0im)

    for out_step = 0:(n_out_step-1)
        for in_step = 0:(n_in_step-1)
            update_type = get_update_type(n_ex_update_path, data.i_flg_orbital_general, rng)

            if update_type == HOPPING
                state.electron_config.counter[1] += 1

                mi, ri, rj, s, reject_flag = make_candidate_hopping(
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    n_site,
                    n_elec,
                    loc_spn,
                    rng,
                )

                if reject_flag != 0
                    continue
                end

                # Update electron configuration
                update_ele_config!(
                    mi,
                    ri,
                    rj,
                    s,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    n_site,
                    n_elec,
                )

                # Update projection counts
                update_proj_cnt!(
                    ri,
                    rj,
                    s,
                    proj_cnt_new,
                    tmp_ele_proj_cnt,
                    tmp_ele_num,
                    data,
                )

                # Calculate new Pfaffian (CalculateNewPfM2)
                # Use pre-allocated workspace array (reset before use)
                fill!(pf_m_new, 0.0 + 0.0im)
                calculate_new_pf_m2!(
                    mi,
                    s,
                    pf_m_new,
                    tmp_ele_idx,
                    qp_start,
                    qp_end,
                    data,
                    state,
                )

                # Calculate log inner product (CalculateLogIP_fcmp)
                log_ip_new = calculate_log_ip_fcmp(pf_m_new, qp_start, qp_end, data)
                if use_rbm
                    update_rbm_cnt_hopping!(rbm_cnt_new, rbm_cnt_old, ri, rj, s, data)
                end

                # Metropolis acceptance/rejection
                x = log_proj_ratio(proj_cnt_new, tmp_ele_proj_cnt, data)
                if use_rbm
                    x += real(log_rbm_ratio(rbm_cnt_new, rbm_cnt_old, data))
                end
                w = exp(2.0 * real(x + log_ip_new - log_ip_old))
                if !isfinite(w)
                    w = -1.0  # Should be rejected
                end

                r_metro = rng_real2(rng)
                if w > r_metro
                    # Accept
                    update_m_all!(mi, s, tmp_ele_idx, qp_start, qp_end, data, state)
                    copy!(tmp_ele_proj_cnt, proj_cnt_new)
                    state.slater_matrix.pf_m .= pf_m_new
                    log_ip_old = log_ip_new
                    if use_rbm
                        copy!(rbm_cnt_old, rbm_cnt_new)
                    end
                    n_accept += 1
                    state.electron_config.counter[2] += 1
                else
                    # Reject
                    revert_ele_config!(
                        mi,
                        ri,
                        rj,
                        s,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        n_site,
                        n_elec,
                    )
                end

            elseif update_type == EXCHANGE
                # Exchange update: two electrons exchange positions
                mi, ri, mj, rj, s, t, reject_flag = make_candidate_exchange(
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    n_site,
                    n_elec,
                    tmp_ele_num,
                    rng,
                )

                if reject_flag != 0
                    continue
                end

                # Debug: Validate exchange candidate (s != t should always be true for valid exchange)
                if s == t
                    @error "EXCHANGE BUG: s == t (both spins are $(s)). This should have been rejected!"
                    @error "  mi=$(mi), ri=$(ri), mj=$(mj), rj=$(rj), s=$(s), t=$(t)"
                    continue
                end

                # Store old positions
                ri_old = ri
                rj_old = rj

                # Debug: Check ele_num before exchange for half-filling
                n0_ri_before = tmp_ele_num[ri+1]
                n1_ri_before = tmp_ele_num[ri+n_site+1]
                n0_rj_before = tmp_ele_num[rj+1]
                n1_rj_before = tmp_ele_num[rj+n_site+1]

                # Update electron configuration (first electron)
                update_ele_config!(
                    mi,
                    ri,
                    rj,
                    s,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    n_site,
                    n_elec,
                )
                update_proj_cnt!(
                    ri,
                    rj,
                    s,
                    proj_cnt_new,
                    tmp_ele_proj_cnt,
                    tmp_ele_num,
                    data,
                )

                # Update electron configuration (second electron)
                update_ele_config!(
                    mj,
                    rj,
                    ri,
                    t,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    n_site,
                    n_elec,
                )
                update_proj_cnt!(rj, ri, t, proj_cnt_new, proj_cnt_new, tmp_ele_num, data)

                # Debug: Check ele_num after exchange for half-filling
                n0_ri_after = tmp_ele_num[ri+1]
                n1_ri_after = tmp_ele_num[ri+n_site+1]
                n0_rj_after = tmp_ele_num[rj+1]
                n1_rj_after = tmp_ele_num[rj+n_site+1]

                # Verify half-filling is preserved at affected sites
                if (n0_ri_after + n1_ri_after != 1) || (n0_rj_after + n1_rj_after != 1)
                    @error "EXCHANGE BUG: Half-filling violated after exchange!"
                    @error "  mi=$(mi), ri=$(ri), mj=$(mj), rj=$(rj), s=$(s), t=$(t)"
                    @error "  Before: ri(n0=$n0_ri_before,n1=$n1_ri_before), rj(n0=$n0_rj_before,n1=$n1_rj_before)"
                    @error "  After:  ri(n0=$n0_ri_after,n1=$n1_ri_after), rj(n0=$n0_rj_after,n1=$n1_rj_after)"
                end

                # Calculate new Pfaffian (CalculateNewPfMTwo2_fcmp)
                pf_m_new = zeros(ComplexF64, length(state.slater_matrix.pf_m))
                calculate_new_pf_m_two2!(
                    mi,
                    s,
                    mj,
                    t,
                    pf_m_new,
                    tmp_ele_idx,
                    qp_start,
                    qp_end,
                    data,
                    state,
                )

                # Calculate log inner product (CalculateLogIP_fcmp)
                log_ip_new = calculate_log_ip_fcmp(pf_m_new, qp_start, qp_end, data)
                if use_rbm
                    update_rbm_cnt_hopping!(rbm_cnt_new, rbm_cnt_old, ri, rj, s, data)
                    update_rbm_cnt_hopping!(rbm_cnt_new, rbm_cnt_new, rj, ri, t, data)
                end

                # Metropolis acceptance/rejection
                x = log_proj_ratio(proj_cnt_new, tmp_ele_proj_cnt, data)
                if use_rbm
                    x += real(log_rbm_ratio(rbm_cnt_new, rbm_cnt_old, data))
                end
                w = exp(2.0 * real(x + log_ip_new - log_ip_old))
                if !isfinite(w)
                    w = -1.0  # Should be rejected
                end

                r_metro = rng_real2(rng)
                if w > r_metro
                    # Accept
                    update_m_all_two!(
                        mi,
                        s,
                        mj,
                        t,
                        ri_old,
                        rj_old,
                        tmp_ele_idx,
                        qp_start,
                        qp_end,
                        data,
                        state,
                    )
                    tmp_ele_proj_cnt .= proj_cnt_new
                    log_ip_old = log_ip_new
                    if use_rbm
                        copy!(rbm_cnt_old, rbm_cnt_new)
                    end
                    n_accept += 1
                else
                    # Reject: revert electron configuration
                    revert_ele_config!(
                        mj,
                        rj,
                        ri,
                        t,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        n_site,
                        n_elec,
                    )
                    revert_ele_config!(
                        mi,
                        ri,
                        rj,
                        s,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        n_site,
                        n_elec,
                    )
                end
            end

            # Recalculate Pfaffian periodically
            if n_accept > n_site
                # Recalculate PfM and InvM from scratch
                result = calculate_m_all_fcmp!(tmp_ele_idx, qp_start, qp_end, data, state)
                if result == 0
                    log_ip_old = calculate_log_ip_fcmp(
                        state.slater_matrix.pf_m,
                        qp_start,
                        qp_end,
                        data,
                    )
                    if use_rbm
                        rbm_cnt_old = make_rbm_cnt(tmp_ele_num, data)
                        if length(rbm_cnt_new) != length(rbm_cnt_old)
                            rbm_cnt_new = similar(rbm_cnt_old)
                        end
                    end
                end
                n_accept = 0
            end
        end

        # Save sample
        # C implementation: saveEleConfig is called for each sample after warmup
        # C: if(outStep >= nOutStep-NVMCSample) where nOutStep = NVMCWarmUp+NVMCSample
        # So: outStep >= NVMCWarmUp, sample = outStep - NVMCWarmUp
        save_start = burn_flag ? 1 : n_vmc_warmup
        if out_step >= save_start
            sample = out_step - save_start
            if sample >= 0 && sample < n_vmc_sample
                # Validate tmp_ele_idx before saving
                if all(x -> x == 0, tmp_ele_idx)
                    @error "tmp_ele_idx is all zeros at sample=$sample, out_step=$out_step, n_vmc_warmup=$n_vmc_warmup. Skipping save for this sample."
                    @error "  This indicates the electron configuration was corrupted during sampling. tmp_ele_idx[1:min(10, length(tmp_ele_idx))] = $(tmp_ele_idx[1:min(10, length(tmp_ele_idx))])"
                    continue  # Skip saving this invalid sample
                end

                if all(x -> x < 0, tmp_ele_idx)
                    @error "tmp_ele_idx is all negative at sample=$sample, out_step=$out_step. Skipping save for this sample."
                    @error "  tmp_ele_idx[1:min(10, length(tmp_ele_idx))] = $(tmp_ele_idx[1:min(10, length(tmp_ele_idx))])"
                    continue  # Skip saving this invalid sample
                end

                # Save electron configuration
                n_size = 2 * n_elec
                n_site2 = 2 * n_site
                offset_idx = sample * n_size
                offset_cfg = sample * n_site2
                offset_num = sample * n_site2
                offset_proj = sample * n_proj

                # Debug: Check ele_idx size before saving
                ele_idx_size = length(state.electron_config.ele_idx)
                expected_size = n_vmc_sample * n_size
                if ele_idx_size != expected_size
                    @error "ele_idx size mismatch at sample=$sample: expected=$expected_size, actual=$ele_idx_size. offset_idx=$offset_idx, max_index=$(offset_idx + n_size)"
                end

                # Debug: Check tmp_ele_idx before saving
                if sample == 0
                    @debug "Before saving sample=$sample: tmp_ele_idx[1:10] = $(tmp_ele_idx[1:min(10, length(tmp_ele_idx))]), offset_idx=$offset_idx, ele_idx_size=$ele_idx_size"
                end

                # Check bounds before accessing
                if offset_idx + n_size > ele_idx_size
                    @error "BoundsError: offset_idx + n_size = $(offset_idx + n_size) > ele_idx_size = $ele_idx_size at sample=$sample, out_step=$out_step"
                    continue
                end

                for i = 1:n_size
                    state.electron_config.ele_idx[offset_idx+i] = tmp_ele_idx[i]
                end

                # Debug: Check saved ele_idx after saving
                if sample == 0
                    saved_ele_idx =
                        state.electron_config.ele_idx[(offset_idx+1):(offset_idx+min(
                            10,
                            n_size,
                        ))]
                    @debug "After saving sample=$sample: saved_ele_idx[1:10] = $saved_ele_idx"
                end

                # Debug: Check tmp_ele_num consistency before saving
                # Note: Half-filling check (n0 + n1 == 1 for each site) only applies to pure spin models
                # (e.g., Heisenberg) where nlocspin > 0 and there are no conduction electrons.
                # For itinerant electron models (Hubbard: nlocspin == 0) or Kondo (ncond > 0),
                # sites can be empty or doubly occupied, so disable this check.
                if data.modpara.nlocspin > 0 && data.modpara.ncond <= 0
                    n_violations = 0
                    violation_sites = Int[]
                    for ri = 0:(n_site-1)
                        n0_i = tmp_ele_num[ri+1]
                        n1_i = tmp_ele_num[ri+n_site+1]
                        if n0_i + n1_i != 1
                            n_violations += 1
                            push!(violation_sites, ri)
                        end
                    end
                    if n_violations > 0
                        @error "Sample $sample: tmp_ele_num violates half-filling at $n_violations sites: $violation_sites"
                        @error "  tmp_ele_num n0[0:$(n_site-1)] = $(tmp_ele_num[1:n_site])"
                        @error "  tmp_ele_num n1[0:$(n_site-1)] = $(tmp_ele_num[n_site+1:2*n_site])"
                        @error "  tmp_ele_idx = $(tmp_ele_idx[1:n_size])"
                        # Dump detailed info for debugging
                        for vi in violation_sites[1:min(5, length(violation_sites))]
                            n0_vi = tmp_ele_num[vi+1]
                            n1_vi = tmp_ele_num[vi+n_site+1]
                            cfg0_vi = tmp_ele_cfg[vi+1]
                            cfg1_vi = tmp_ele_cfg[vi+n_site+1]
                            @error "  Site $vi: n0=$n0_vi, n1=$n1_vi, cfg0=$cfg0_vi, cfg1=$cfg1_vi (n0+n1=$(n0_vi+n1_vi))"
                        end
                    end
                end

                # Debug: Verify consistency between tmp_ele_num and tmp_ele_cfg
                cfg_num_mismatch = 0
                for rsi = 0:(n_site2-1)
                    expected_num = tmp_ele_cfg[rsi+1] < 0 ? 0 : 1
                    if tmp_ele_num[rsi+1] != expected_num
                        cfg_num_mismatch += 1
                    end
                end
                if cfg_num_mismatch > 0
                    @error "Sample $sample: tmp_ele_num and tmp_ele_cfg inconsistent at $cfg_num_mismatch positions!"
                    @error "  tmp_ele_cfg[1:8] = $(tmp_ele_cfg[1:min(8, n_site2)])"
                    @error "  tmp_ele_num[1:8] = $(tmp_ele_num[1:min(8, n_site2)])"
                end

                if sample == 0
                    @debug "Sample 0: Before saving, tmp_ele_num[1:$n_site2] = $(tmp_ele_num[1:n_site2])"
                end

                for i = 1:n_site2
                    state.electron_config.ele_cfg[offset_cfg+i] = tmp_ele_cfg[i]
                    state.electron_config.ele_num[offset_num+i] = tmp_ele_num[i]
                end

                # Debug: Verify saved ele_num for sample 0
                if sample == 0
                    saved_ele_num =
                        state.electron_config.ele_num[(offset_num+1):(offset_num+n_site2)]
                    @debug "Sample 0: After saving, saved_ele_num = $saved_ele_num"
                    if saved_ele_num != tmp_ele_num[1:n_site2]
                        @error "Sample 0: Mismatch between tmp_ele_num and saved ele_num!"
                    end
                end
                for i = 1:n_proj
                    state.electron_config.ele_proj_cnt[offset_proj+i] = tmp_ele_proj_cnt[i]
                end
            end
        end
    end

    # Debug: Verify that samples were saved
    n_saved_samples = 0
    if length(state.electron_config.ele_idx) > 0
        n_size = 2 * n_elec
        for s = 0:(n_vmc_sample-1)
            offset = s * n_size
            sample_ele_idx = state.electron_config.ele_idx[(offset+1):(offset+n_size)]
            if !all(x -> x == 0, sample_ele_idx)
                n_saved_samples += 1
            end
        end
        if n_saved_samples == 0
            @error "After vmc_make_sample!, no samples were saved. n_vmc_warmup=$n_vmc_warmup, n_vmc_sample=$n_vmc_sample, n_out_step=$n_out_step"
        end
    end

    # Copy to burn sample (for next step)
    copy_to_burn_sample!(tmp_ele_idx, tmp_ele_cfg, tmp_ele_num, tmp_ele_proj_cnt, state)
    # Set burn_flag for next step
    if length(state.electron_config.counter) >= 11
        state.electron_config.counter[11] = 1
    else
        # Extend counter array if needed
        resize!(
            state.electron_config.counter,
            max(11, length(state.electron_config.counter)),
        )
        state.electron_config.counter[11] = 1
    end
end

# ============================================================================
# FSZ (Free Spin Z) Mode Implementation
# ============================================================================

"""
    calculate_m_all_fsz!(ele_idx::Vector{Int}, ele_spn::Vector{Int},
                         qp_start::Int, qp_end::Int,
                         data::ExpertModeData, state::VMCOptimizationState)::Int

Calculate Pfaffian and inverse matrix for FSZ mode.
Equivalent to C's `CalculateMAll_fsz()`.

In FSZ mode, each electron's spin is tracked in ele_spn array instead of
being determined by electron index (msi < Ne → spin 0, else spin 1).

# Arguments
- `ele_idx`: Electron site indices [mi] (1-based site index)
- `ele_spn`: Electron spin values [mi] (0 or 1)
- `qp_start`, `qp_end`: QP index range (1-based, exclusive end)
- `data`: Expert mode data
- `state`: VMC optimization state

# Returns
- `info::Int`: 0 on success, >0 on error
"""
function calculate_m_all_fsz!(
    ele_idx::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)::Int
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec  # Nsize = 2*Ne
    n_site2 = 2 * n_site
    n_qp_full = length(state.slater_matrix.pf_m)

    # Validate indices
    if qp_start < 1 || qp_end > n_qp_full + 1 || qp_start >= qp_end
        @error "Invalid qp range: qp_start=$qp_start, qp_end=$qp_end, n_qp_full=$n_qp_full"
        return 1
    end

    qp_num = qp_end - qp_start

    # Use workspace arrays
    ws = state.workspace
    inv_m_temp = view(ws.inv_m_temp,:,:,(1:qp_num))
    pf_m_temp = view(ws.pf_m_temp, 1:qp_num)

    # Get thread-local workspace from pfapack_workspace
    thread_ws = get_thread_workspace(ws.pfapack_workspace)

    # Process each QP index
    for qp_local = 1:qp_num
        qp_idx = qp_start + qp_local - 1

        # Get SlaterElm for this QP (1D array with C row-major layout)
        slater_offset = (qp_idx - 1) * n_site2 * n_site2

        # Construct invM from SlaterElm using ele_spn
        # C: invM[msi*Nsize + msj] = -sltE[rsi*Nsite2 + rsj]
        # For LAPACK column-major: invM(msj, msi) = -sltE[rsi][rsj]
        # Julia column-major: inv_m[msj, msi] = -sltE[rsi][rsj]
        # (Same storage as non-FSZ calculate_m_all_child_fcmp!)
        inv_m = view(inv_m_temp,:,:,qp_local)
        @inbounds for msi = 1:n_size
            # ele_idx[msi] is 0-based (0 to n_site-1)
            rsi = ele_idx[msi] + ele_spn[msi] * n_site  # 0-based
            for msj = 1:n_size
                rsj = ele_idx[msj] + ele_spn[msj] * n_site  # 0-based
                # Index into slater_elm (1-based): rsi * n_site2 + rsj + 1
                slater_idx = slater_offset + rsi * n_site2 + rsj + 1
                inv_m[msj, msi] = -state.slater_matrix.slater_elm[slater_idx]
            end
        end

        # LTL decomposition using thread workspace
        iwork = thread_ws.iwork
        work_vec = thread_ws.v_t_complex
        buf_m = thread_ws.buf_m_complex

        # Call the LTL factorization (zsktf2) - only needs matrix and pivot array
        julia_zsktf2!(inv_m, iwork)

        # Calculate Pfaffian using utu2pfa(n, a, lda, ipiv)
        pfaff = utu2pfa(n_size, inv_m, n_size, iwork)
        if !isfinite(real(pfaff)) || !isfinite(imag(pfaff))
            return qp_idx
        end
        pf_m_temp[qp_local] = pfaff

        # Calculate inverse matrix
        # cimpl_utu2inv!(n, a, lda, ipiv, work, buf, ldb)
        cimpl_utu2inv!(n_size, inv_m, n_size, iwork, work_vec, buf_m, n_size)

        # Negate (InvM = -InvM)
        @inbounds for i = 1:n_size
            for j = 1:n_size
                inv_m[i, j] = -inv_m[i, j]
            end
        end
    end

    # Copy results back to state
    @inbounds for qp_local = 1:qp_num
        qp_global = qp_start + qp_local - 1
        state.slater_matrix.pf_m[qp_global] = pf_m_temp[qp_local]
    end

    # Copy inv_m to state using same method as non-FSZ calculate_m_all_fcmp!
    # Use vec() to copy in column-major order (consistent with non-FSZ)
    nsq = n_size * n_size
    @inbounds for qp_local = 1:qp_num
        qp_global = qp_start + qp_local - 1
        dst_start = (qp_global - 1) * nsq + 1
        dst_end = dst_start + nsq - 1
        if dst_end <= length(state.slater_matrix.inv_m)
            dst = view(state.slater_matrix.inv_m, dst_start:dst_end)
            src = view(inv_m_temp,:,:,qp_local)
            copyto!(dst, vec(src))
        end
    end

    return 0
end

"""
    calculate_new_pf_m_fsz!(ma::Int, s::Int, pf_m_new::Vector{ComplexF64},
                            ele_idx::Vector{Int}, ele_spn::Vector{Int},
                            qp_start::Int, qp_end::Int,
                            data::ExpertModeData, state::VMCOptimizationState)

Calculate new Pfaffian after electron ma moves to site with spin s.
Equivalent to C's `CalculateNewPfM2_fsz()`.

# Arguments
- `ma`: Electron index (1-based)
- `s`: New spin (0 or 1)
- `pf_m_new`: Output array for new Pfaffian values
- `ele_idx`, `ele_spn`: Current electron configuration
- `qp_start`, `qp_end`: QP index range
"""
function calculate_new_pf_m_fsz!(
    ma::Int,
    s::Int,
    pf_m_new::Vector{ComplexF64},
    ele_idx::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    qp_num = qp_end - qp_start

    # ma is 0-based electron index
    # msa is 1-based Julia array index for ele_idx/ele_spn access
    msa = ma + 1
    # rsa = eleIdx[ma] + s*Nsite (new position, 0-based site-spin index)
    # ele_idx stores 0-based site indices
    rsa = ele_idx[msa] + s * n_site

    @inbounds for qp_local = 1:qp_num
        qp_idx = qp_start + qp_local - 1

        # Get SlaterElm row for new position
        # rsa is 0-based (0 to n_site2-1), so use rsa directly for offset calculation
        slater_offset = (qp_idx - 1) * n_site2 * n_site2 + rsa * n_site2

        # Get InvM row for electron ma (0-based)
        # C uses: invM + ma*Nsize
        inv_m_offset = (qp_idx - 1) * n_size * n_size + ma * n_size

        # Calculate ratio = sum_j invM[msa][msj] * sltE[rsa][rsj]
        ratio = zero(ComplexF64)
        for msj = 1:n_size
            # rsj is 0-based site-spin index
            rsj = ele_idx[msj] + ele_spn[msj] * n_site
            inv_m_val = state.slater_matrix.inv_m[inv_m_offset+msj]
            # slater_elm is 1-based, so add 1 to 0-based rsj
            slater_val = state.slater_matrix.slater_elm[slater_offset+rsj+1]
            ratio += inv_m_val * slater_val
        end

        # pfMNew[qp] = -ratio * PfM[qp]
        pf_m_new[qp_local] = -ratio * state.slater_matrix.pf_m[qp_idx]
    end
end

"""
    update_m_all_fsz!(ma::Int, s::Int, ele_idx::Vector{Int}, ele_spn::Vector{Int},
                      qp_start::Int, qp_end::Int,
                      data::ExpertModeData, state::VMCOptimizationState)

Update InvM and PfM after electron ma moves with spin s.
Equivalent to C's `UpdateMAll_fsz()`.
"""
function update_m_all_fsz!(
    ma::Int,
    s::Int,
    ele_idx::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    qp_num = qp_end - qp_start

    # ma is 0-based electron index
    # msa is 1-based Julia array index for ele_idx/ele_spn access
    msa = ma + 1
    # rsa = eleIdx[ma] + s*Nsite (0-based site-spin index)
    rsa = ele_idx[msa] + s * n_site

    # Workspace vectors
    vec1 = zeros(ComplexF64, n_size)
    vec2 = zeros(ComplexF64, n_size)

    @inbounds for qp_local = 1:qp_num
        qp_idx = qp_start + qp_local - 1

        # SlaterElm offset for row rsa (rsa is 0-based)
        slater_offset = (qp_idx - 1) * n_site2 * n_site2 + rsa * n_site2

        # InvM offset
        inv_m_base = (qp_idx - 1) * n_size * n_size

        # Calculate vec1[i] = sum_j invM[i][j] * sltE[a][j]
        # Note: invM[i][j] = -invM[j][i] (skew-symmetric)
        fill!(vec1, zero(ComplexF64))
        for msj = 1:n_size
            # rsj is 0-based site-spin index
            rsj = ele_idx[msj] + ele_spn[msj] * n_site
            # slater_elm is 1-based array, so add 1 to 0-based rsj
            slater_aj = state.slater_matrix.slater_elm[slater_offset+rsj+1]

            # invM_j = invM + msj*Nsize (column j)
            for msi = 1:n_size
                # invM[msj][msi] at offset inv_m_base + (msj-1)*n_size + msi
                inv_m_ji = state.slater_matrix.inv_m[inv_m_base+(msj-1)*n_size+msi]
                vec1[msi] += -inv_m_ji * slater_aj
            end
        end

        # Update Pfaffian
        tmp = vec1[msa]
        state.slater_matrix.pf_m[qp_idx] *= -tmp
        inv_vec1_a = -one(ComplexF64) / tmp

        # Calculate vec2[i] = -InvM[a][i]/vec1[a]
        for msi = 1:n_size
            inv_m_ai = state.slater_matrix.inv_m[inv_m_base+(msa-1)*n_size+msi]
            vec2[msi] = inv_m_ai * inv_vec1_a
        end

        # Update InvM
        for msi = 1:n_size
            vec1_i = vec1[msi]
            vec2_i = vec2[msi]
            for msj = 1:n_size
                idx = inv_m_base + (msi-1)*n_size + msj
                state.slater_matrix.inv_m[idx] += vec1_i * vec2[msj] - vec1[msj] * vec2_i
            end
            # invM_i[msa] -= vec2_i
            state.slater_matrix.inv_m[inv_m_base+(msi-1)*n_size+msa] -= vec2_i
        end

        # invM_a[msj] += vec2[msj]
        for msj = 1:n_size
            state.slater_matrix.inv_m[inv_m_base+(msa-1)*n_size+msj] += vec2[msj]
        end
    end
end

"""
    make_initial_sample_fsz!(ele_idx::Vector{Int}, ele_cfg::Vector{Int},
                             ele_num::Vector{Int}, ele_proj_cnt::Vector{Int},
                             ele_spn::Vector{Int},
                             qp_start::Int, qp_end::Int,
                             data::ExpertModeData, state::VMCOptimizationState,
                             rng::AbstractRNG)::Int

Initialize electron configuration for FSZ mode.
Equivalent to C's `makeInitialSample_fsz()`.

# Returns
- 0 on success
"""
function make_initial_sample_fsz!(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
    rng::AbstractRNG,
)::Int
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    two_sz = data.modpara.two_sz
    loc_spn = get_loc_spn_array(data)
    dump_loc_spn_if_enabled(loc_spn, n_site)

    max_loops = 100
    loop = 0

    while true
        # Initialize arrays
        fill!(ele_idx, -1)
        fill!(ele_spn, -1)
        fill!(ele_cfg, -1)

        # Determine initial Sz configuration
        # If TwoSz == -1, start with Sz = 0
        tmp_two_sz = (two_sz == -1) ? 0 : div(two_sz, 2)

        # Assign spins to electrons
        # Electrons 0 to Ne+tmp_two_sz-1 get spin 0 (up)
        # Electrons Ne+tmp_two_sz to Nsize-1 get spin 1 (down)
        # mi is 0-based electron index, stored in 1-based Julia array
        for mi = 0:(n_size-1)
            if mi < n_elec + tmp_two_sz
                ele_spn[mi+1] = 0  # up spin
            else
                ele_spn[mi+1] = 1  # down spin
            end
        end

        # Place local spins first (ri is 0-based site index)
        for ri = 0:(n_site-1)
            if loc_spn[ri+1] == 1
                # Find an unplaced electron
                mi = 0
                while true
                    mi = rng_mod(rng, n_size)  # 0-based electron index
                    if ele_idx[mi+1] == -1
                        break
                    end
                end
                si = ele_spn[mi+1]  # 0 or 1
                # ele_cfg uses 1-based Julia indexing, stores 0-based electron index
                ele_cfg[ri+1+si*n_site] = mi
                # ele_idx uses 1-based Julia indexing, stores 0-based site index
                ele_idx[mi+1] = ri
            end
        end

        # Place itinerant electrons
        for mi = 0:(n_size-1)
            if ele_idx[mi+1] == -1
                si = ele_spn[mi+1]
                # Find an empty itinerant site
                while true
                    ri = rng_mod(rng, n_site)  # 0-based site index
                    if ele_cfg[ri+1+si*n_site] == -1 && loc_spn[ri+1] == 0
                        ele_cfg[ri+1+si*n_site] = mi
                        ele_idx[mi+1] = ri
                        break
                    end
                end
            end
        end

        # Set EleNum (1-based Julia indexing)
        for rsi = 1:n_site2
            ele_num[rsi] = (ele_cfg[rsi] == -1) ? 0 : 1
        end

        # Calculate projection counts
        make_proj_cnt!(ele_proj_cnt, ele_num, data)

        # Calculate Pfaffian
        flag = calculate_m_all_fsz!(ele_idx, ele_spn, qp_start, qp_end, data, state)

        if flag == 0
            break
        end

        loop += 1
        if loop > max_loops
            @error "makeInitialSample_fsz: Too many loops"
            return 1
        end
    end

    dump_elec_initial_if_enabled(
        ele_idx,
        ele_cfg,
        ele_num,
        ele_spn,
        n_elec,
        n_site,
    )

    return 0
end

"""
    make_initial_sample_fsz_real!(ele_idx::Vector{Int}, ele_cfg::Vector{Int},
                                  ele_num::Vector{Int}, ele_proj_cnt::Vector{Int},
                                  ele_spn::Vector{Int},
                                  qp_start::Int, qp_end::Int,
                                  data::ExpertModeData, state::VMCOptimizationState,
                                  rng::AbstractRNG)::Int

Initialize electron configuration for FSZ real mode.
Equivalent to C's `makeInitialSample_fsz_real()`.
"""
function make_initial_sample_fsz_real!(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
    rng::AbstractRNG,
)::Int
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    two_sz = data.modpara.two_sz
    loc_spn = get_loc_spn_array(data)
    dump_loc_spn_if_enabled(loc_spn, n_site)

    max_loops = 100
    loop = 0

    while true
        fill!(ele_idx, -1)
        fill!(ele_spn, -1)
        fill!(ele_cfg, -1)

        tmp_two_sz = (two_sz == -1) ? 0 : div(two_sz, 2)
        for mi = 0:(n_size-1)
            if mi < n_elec + tmp_two_sz
                ele_spn[mi+1] = 0
            else
                ele_spn[mi+1] = 1
            end
        end

        for ri = 0:(n_site-1)
            if loc_spn[ri+1] == 1
                mi = 0
                while true
                    mi = rng_mod(rng, n_size)
                    if ele_idx[mi+1] == -1
                        break
                    end
                end
                si = ele_spn[mi+1]
                ele_cfg[ri+1+si*n_site] = mi
                ele_idx[mi+1] = ri
            end
        end

        for mi = 0:(n_size-1)
            if ele_idx[mi+1] == -1
                si = ele_spn[mi+1]
                while true
                    ri = rng_mod(rng, n_site)
                    if ele_cfg[ri+1+si*n_site] == -1 && loc_spn[ri+1] == 0
                        ele_cfg[ri+1+si*n_site] = mi
                        ele_idx[mi+1] = ri
                        break
                    end
                end
            end
        end

        for rsi = 1:n_site2
            ele_num[rsi] = (ele_cfg[rsi] == -1) ? 0 : 1
        end

        make_proj_cnt!(ele_proj_cnt, ele_num, data)

        flag = calculate_m_all_fsz_real!(ele_idx, ele_spn, qp_start, qp_end, data, state)
        if flag == 0
            break
        end

        loop += 1
        if loop > max_loops
            @error "makeInitialSample_fsz_real: Too many loops"
            return 1
        end
    end

    dump_elec_initial_if_enabled(
        ele_idx,
        ele_cfg,
        ele_num,
        ele_spn,
        n_elec,
        n_site,
    )

    return 0
end

"""
    update_ele_config_fsz!(mi::Int, org_r::Int, dst_r::Int, org_spn::Int, dst_spn::Int,
                           ele_idx::Vector{Int}, ele_cfg::Vector{Int},
                           ele_num::Vector{Int}, ele_spn::Vector{Int}, n_site::Int)

Update electron configuration for FSZ mode.
Electron mi (0-based) moves from (org_r, org_spn) to (dst_r, dst_spn).
org_r and dst_r are 0-based site indices.
"""
function update_ele_config_fsz!(
    mi::Int,
    org_r::Int,
    dst_r::Int,
    org_spn::Int,
    dst_spn::Int,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_spn::Vector{Int},
    n_site::Int,
)
    # ele_idx, ele_spn use 1-based Julia indexing, mi is 0-based
    ele_idx[mi+1] = dst_r
    ele_spn[mi+1] = dst_spn

    # ele_cfg, ele_num use 1-based Julia indexing, org_r/dst_r are 0-based
    ele_cfg[org_r+1+org_spn*n_site] = -1
    ele_cfg[dst_r+1+dst_spn*n_site] = mi  # Store 0-based electron index

    ele_num[org_r+1+org_spn*n_site] = 0
    ele_num[dst_r+1+dst_spn*n_site] = 1
end

"""
    revert_ele_config_fsz!(mi::Int, org_r::Int, dst_r::Int, org_spn::Int, dst_spn::Int,
                           ele_idx::Vector{Int}, ele_cfg::Vector{Int},
                           ele_num::Vector{Int}, ele_spn::Vector{Int}, n_site::Int)

Revert electron configuration for FSZ mode.
mi is 0-based electron index, org_r/dst_r are 0-based site indices.
"""
function revert_ele_config_fsz!(
    mi::Int,
    org_r::Int,
    dst_r::Int,
    org_spn::Int,
    dst_spn::Int,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_spn::Vector{Int},
    n_site::Int,
)
    # ele_idx, ele_spn use 1-based Julia indexing
    ele_idx[mi+1] = org_r
    ele_spn[mi+1] = org_spn

    # ele_cfg uses 1-based Julia indexing
    ele_cfg[org_r+1+org_spn*n_site] = mi  # Store 0-based electron index
    ele_cfg[dst_r+1+dst_spn*n_site] = -1

    # ele_num uses 1-based Julia indexing
    ele_num[org_r+1+org_spn*n_site] = 1
    ele_num[dst_r+1+dst_spn*n_site] = 0
end

"""
    vmc_make_sample_fsz!(data::ExpertModeData, state::VMCOptimizationState, rng::AbstractRNG)

VMC sampling (complex version, fsz).
Equivalent to C's `VMCMakeSample_fsz()`.

Implements Monte Carlo sampling with free spin Z (non-conserved total Sz).
"""
function vmc_make_sample_fsz!(
    data::ExpertModeData,
    state::VMCOptimizationState,
    rng::AbstractRNG = Random.GLOBAL_RNG,
    c_timer::CTimer = CTIMER_DISABLED,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    n_proj = length(state.electron_config.tmp_ele_proj_cnt)
    n_vmc_warmup = data.modpara.nvmc_warmup
    n_vmc_interval = data.modpara.nvmc_interval
    n_vmc_sample = data.modpara.nvmc_sample
    n_ex_update_path = data.modpara.nex_update_path
    two_sz = data.modpara.two_sz
    i_flg_orbital_general = data.i_flg_orbital_general
    loc_spn = get_loc_spn_array(data)

    # QP range
    n_qp_full = length(state.slater_matrix.pf_m)
    qp_start = 1
    qp_end = n_qp_full + 1

    # Get burn flag
    burn_flag =
        length(state.electron_config.counter) >= 11 ? state.electron_config.counter[11] : 0

    # Temporary arrays
    tmp_ele_idx = state.electron_config.tmp_ele_idx
    tmp_ele_cfg = state.electron_config.tmp_ele_cfg
    tmp_ele_num = state.electron_config.tmp_ele_num
    tmp_ele_proj_cnt = state.electron_config.tmp_ele_proj_cnt
    tmp_ele_spn = state.electron_config.tmp_ele_spn

    # Workspace for new Pfaffian
    pf_m_new = zeros(ComplexF64, n_qp_full)
    proj_cnt_new = zeros(Int, n_proj)

    # Initialize or restore configuration
    if burn_flag == 0
        make_initial_sample_fsz!(
            tmp_ele_idx,
            tmp_ele_cfg,
            tmp_ele_num,
            tmp_ele_proj_cnt,
            tmp_ele_spn,
            qp_start,
            qp_end,
            data,
            state,
            rng,
        )
    else
        copy_from_burn_sample_fsz!(
            tmp_ele_idx,
            tmp_ele_cfg,
            tmp_ele_num,
            tmp_ele_proj_cnt,
            tmp_ele_spn,
            state,
        )
    end

    # Calculate initial Pfaffian
    calculate_m_all_fsz!(tmp_ele_idx, tmp_ele_spn, qp_start, qp_end, data, state)
    log_ip_old = calculate_log_ip_fcmp(state.slater_matrix.pf_m, qp_start, qp_end, data)

    # Validate initial configuration
    if !isfinite(real(log_ip_old)) || !isfinite(imag(log_ip_old))
        @warn "VMCMakeSample_fsz: remakeSample logIpOld not finite"
        make_initial_sample_fsz!(
            tmp_ele_idx,
            tmp_ele_cfg,
            tmp_ele_num,
            tmp_ele_proj_cnt,
            tmp_ele_spn,
            qp_start,
            qp_end,
            data,
            state,
            rng,
        )
        calculate_m_all_fsz!(tmp_ele_idx, tmp_ele_spn, qp_start, qp_end, data, state)
        log_ip_old = calculate_log_ip_fcmp(state.slater_matrix.pf_m, qp_start, qp_end, data)
        burn_flag = 0
    end

    # Sampling loops
    n_out_step = (burn_flag == 0) ? n_vmc_warmup + n_vmc_sample : n_vmc_sample + 1
    n_in_step = n_vmc_interval * n_site
    n_accept = 0

    # Reset counters
    if length(state.electron_config.counter) < 6
        resize!(state.electron_config.counter, 6)
    end
    fill!(state.electron_config.counter, 0)

    for out_step = 0:(n_out_step-1)
        for in_step = 0:(n_in_step-1)
            update_type = get_update_type(
                n_ex_update_path,
                i_flg_orbital_general,
                rng;
                two_sz = two_sz,
            )

            if update_type == HOPPING
                # Hopping update
                state.electron_config.counter[1] += 1

                # For FSZ with TwoSz == -1, allow spin flips (C: rand<0.5 => hopping)
                if two_sz == -1 && rng_real2(rng) < 0.5
                    # Standard hopping (spin conserved or FSZ hopping)
                    mi, ri, rj, s, t, reject_flag = make_candidate_hopping_fsz(
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        tmp_ele_spn,
                        loc_spn,
                        n_site,
                        n_size,
                        two_sz,
                        rng,
                    )
                else
                    # Local spin flip for conduction electrons
                    mi, ri, rj, s, t, reject_flag =
                        make_candidate_local_spin_flip_conduction(
                            tmp_ele_idx,
                            tmp_ele_cfg,
                            tmp_ele_num,
                            tmp_ele_spn,
                            loc_spn,
                            n_site,
                            n_size,
                            rng,
                        )
                end

                if reject_flag
                    continue
                end

                # Update configuration temporarily
                update_ele_config_fsz!(
                    mi,
                    ri,
                    rj,
                    s,
                    t,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    tmp_ele_spn,
                    n_site,
                )

                # Update projection counts (ri and rj are already 0-based)
                if s == t
                    update_proj_cnt!(
                        ri,
                        rj,
                        s,
                        proj_cnt_new,
                        tmp_ele_proj_cnt,
                        tmp_ele_num,
                        data,
                    )
                else
                    update_proj_cnt_fsz!(
                        ri,
                        rj,
                        s,
                        t,
                        proj_cnt_new,
                        tmp_ele_proj_cnt,
                        tmp_ele_num,
                        data,
                    )
                end

                # Calculate new Pfaffian
                calculate_new_pf_m_fsz!(
                    mi,
                    t,
                    pf_m_new,
                    tmp_ele_idx,
                    tmp_ele_spn,
                    qp_start,
                    qp_end,
                    data,
                    state,
                )

                # Calculate new log IP
                log_ip_new = calculate_log_ip_fcmp(pf_m_new, 1, n_qp_full + 1, data)

                # Metropolis criterion
                x = log_proj_ratio(proj_cnt_new, tmp_ele_proj_cnt, data)
                w = exp(2.0 * (x + real(log_ip_new - log_ip_old)))
                if !isfinite(w)
                    w = -1.0
                end

                r = rng_real2(rng)
                accept = w > r
                if accept
                    # Accept
                    update_m_all_fsz!(
                        mi,
                        t,
                        tmp_ele_idx,
                        tmp_ele_spn,
                        qp_start,
                        qp_end,
                        data,
                        state,
                    )
                    copyto!(tmp_ele_proj_cnt, proj_cnt_new)
                    log_ip_old = log_ip_new
                    n_accept += 1
                    state.electron_config.counter[2] += 1
                else
                    # Reject
                    revert_ele_config_fsz!(
                        mi,
                        ri,
                        rj,
                        s,
                        t,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        tmp_ele_spn,
                        n_site,
                    )
                end

            elseif update_type == EXCHANGE
                # Exchange update
                state.electron_config.counter[3] += 1

                mi, ri, rj, s, reject_flag = make_candidate_exchange_fsz(
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    tmp_ele_spn,
                    loc_spn,
                    n_site,
                    n_size,
                    rng,
                )

                if reject_flag
                    continue
                end

                t = 1 - s
                # rj is 0-based, ele_cfg uses 1-based Julia indexing
                mj = tmp_ele_cfg[rj+1+t*n_site]

                # First hop: mi (ri, s) -> (rj, s)
                update_ele_config_fsz!(
                    mi,
                    ri,
                    rj,
                    s,
                    s,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    tmp_ele_spn,
                    n_site,
                )
                # ri and rj are already 0-based
                update_proj_cnt!(
                    ri,
                    rj,
                    s,
                    proj_cnt_new,
                    tmp_ele_proj_cnt,
                    tmp_ele_num,
                    data,
                )

                # Second hop: mj (rj, t) -> (ri, t)
                update_ele_config_fsz!(
                    mj,
                    rj,
                    ri,
                    t,
                    t,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    tmp_ele_spn,
                    n_site,
                )
                update_proj_cnt!(rj, ri, t, proj_cnt_new, proj_cnt_new, tmp_ele_num, data)

                # Calculate new Pfaffian (two-electron update)
                calculate_new_pf_m_two_fsz!(
                    mi,
                    s,
                    mj,
                    t,
                    pf_m_new,
                    tmp_ele_idx,
                    tmp_ele_spn,
                    qp_start,
                    qp_end,
                    data,
                    state,
                )

                log_ip_new = calculate_log_ip_fcmp(pf_m_new, 1, n_qp_full + 1, data)

                x = log_proj_ratio(proj_cnt_new, tmp_ele_proj_cnt, data)
                w = exp(2.0 * (x + real(log_ip_new - log_ip_old)))
                if !isfinite(w)
                    w = -1.0
                end

                r = rng_real2(rng)
                accept = w > r
                if accept
                    # Accept
                    update_m_all_two_fsz!(
                        mi,
                        s,
                        mj,
                        t,
                        ri,
                        rj,
                        tmp_ele_idx,
                        tmp_ele_spn,
                        qp_start,
                        qp_end,
                        data,
                        state,
                    )
                    copyto!(tmp_ele_proj_cnt, proj_cnt_new)
                    log_ip_old = log_ip_new
                    n_accept += 1
                    state.electron_config.counter[4] += 1
                else
                    # Reject
                    revert_ele_config_fsz!(
                        mj,
                        rj,
                        ri,
                        t,
                        t,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        tmp_ele_spn,
                        n_site,
                    )
                    revert_ele_config_fsz!(
                        mi,
                        ri,
                        rj,
                        s,
                        s,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        tmp_ele_spn,
                        n_site,
                    )
                end

            elseif update_type == LOCALSPINFLIP
                # Local spin flip for localized spins
                state.electron_config.counter[5] += 1

                mi, ri, rj, s, t, reject_flag = make_candidate_local_spin_flip_localspin(
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    tmp_ele_spn,
                    loc_spn,
                    n_site,
                    n_size,
                    rng,
                )

                if reject_flag
                    continue
                end

                update_ele_config_fsz!(
                    mi,
                    ri,
                    rj,
                    s,
                    t,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    tmp_ele_spn,
                    n_site,
                )
                # ri and rj are already 0-based
                update_proj_cnt_fsz!(
                    ri,
                    rj,
                    s,
                    t,
                    proj_cnt_new,
                    tmp_ele_proj_cnt,
                    tmp_ele_num,
                    data,
                )

                calculate_new_pf_m_fsz!(
                    mi,
                    t,
                    pf_m_new,
                    tmp_ele_idx,
                    tmp_ele_spn,
                    qp_start,
                    qp_end,
                    data,
                    state,
                )

                log_ip_new = calculate_log_ip_fcmp(pf_m_new, 1, n_qp_full + 1, data)

                x = log_proj_ratio(proj_cnt_new, tmp_ele_proj_cnt, data)
                w = exp(2.0 * (x + real(log_ip_new - log_ip_old)))
                if !isfinite(w)
                    w = -1.0
                end

                r = rng_real2(rng)
                accept = w > r
                if accept
                    update_m_all_fsz!(
                        mi,
                        t,
                        tmp_ele_idx,
                        tmp_ele_spn,
                        qp_start,
                        qp_end,
                        data,
                        state,
                    )
                    copyto!(tmp_ele_proj_cnt, proj_cnt_new)
                    log_ip_old = log_ip_new
                    n_accept += 1
                    state.electron_config.counter[6] += 1
                else
                    revert_ele_config_fsz!(
                        mi,
                        ri,
                        rj,
                        s,
                        t,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        tmp_ele_spn,
                        n_site,
                    )
                end
            end

            # Recalculate if too many accepts
            if n_accept > n_site
                calculate_m_all_fsz!(
                    tmp_ele_idx,
                    tmp_ele_spn,
                    qp_start,
                    qp_end,
                    data,
                    state,
                )
                log_ip_old =
                    calculate_log_ip_fcmp(state.slater_matrix.pf_m, qp_start, qp_end, data)
                n_accept = 0
            end
        end

        # Save sample
        if out_step >= n_out_step - n_vmc_sample
            sample = out_step - (n_out_step - n_vmc_sample)
            save_ele_config_fsz!(
                sample,
                log_ip_old,
                tmp_ele_idx,
                tmp_ele_cfg,
                tmp_ele_num,
                tmp_ele_proj_cnt,
                tmp_ele_spn,
                data,
                state,
            )
        end
    end

    # Copy to burn sample
    copy_to_burn_sample_fsz!(
        tmp_ele_idx,
        tmp_ele_cfg,
        tmp_ele_num,
        tmp_ele_proj_cnt,
        tmp_ele_spn,
        state,
    )

    # Set burn flag
    if length(state.electron_config.counter) >= 11
        state.electron_config.counter[11] = 1
    else
        resize!(state.electron_config.counter, 11)
        state.electron_config.counter[11] = 1
    end
end

# FSZ helper functions

function make_candidate_hopping_fsz(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_spn::Vector{Int},
    loc_spn::Vector{Int},
    n_site::Int,
    n_size::Int,
    two_sz::Int,
    rng::AbstractRNG,
)
    icnt_max = n_site * n_site

    # Select an itinerant electron
    # mi is 0-based electron index (but stored in 1-based Julia array)
    # ri is 0-based site index
    mi = 0
    s = 0
    ri = 0
    while true
        mi = rng_mod(rng, n_size)  # 0-based electron index
        s = ele_spn[mi+1]             # ele_spn uses 1-based Julia indexing
        ri = ele_idx[mi+1]            # ele_idx uses 1-based Julia indexing, stores 0-based site
        if loc_spn[ri+1] == 0         # loc_spn uses 1-based Julia indexing
            break
        end
    end

    # Find empty destination
    # rj is 0-based site index
    rj = 0
    t = 0
    icnt = 0
    flag = false
    while true
        rj = rng_mod(rng, n_site)  # 0-based site index
        if two_sz == -1
            t = rng_real2(rng) < 0.5 ? 0 : 1
        else
            t = s  # Conserve Sz
        end
        if icnt > icnt_max
            flag = true
            break
        end
        icnt += 1
        # ele_cfg uses 1-based Julia indexing
        if ele_cfg[rj+1+t*n_site] == -1 && loc_spn[rj+1] == 0
            break
        end
    end

    return mi, ri, rj, s, t, flag
end

function make_candidate_local_spin_flip_conduction(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_spn::Vector{Int},
    loc_spn::Vector{Int},
    n_site::Int,
    n_size::Int,
    rng::AbstractRNG,
)
    icnt_max = n_site * n_site
    icnt = 0
    flag = false

    # mi is 0-based electron index, ri/rj are 0-based site indices
    mi = 0
    s = 0
    t = 0
    ri = 0
    rj = 0

    while true
        mi = rng_mod(rng, n_size)  # 0-based electron index
        s = ele_spn[mi+1]             # 1-based Julia array
        t = 1 - s
        ri = ele_idx[mi+1]            # 0-based site index
        rj = ri  # Local spin flip
        if icnt > icnt_max
            flag = true
            break
        end
        icnt += 1
        # 1-based Julia array indexing
        if loc_spn[ri+1] == 0 && ele_cfg[ri+1+t*n_site] == -1
            break
        end
    end

    return mi, ri, rj, s, t, flag
end

function make_candidate_local_spin_flip_localspin(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_spn::Vector{Int},
    loc_spn::Vector{Int},
    n_site::Int,
    n_size::Int,
    rng::AbstractRNG,
)
    # Select a local spin
    # mi is 0-based electron index, ri is 0-based site index
    mi = 0
    s = 0
    t = 0
    ri = 0
    flag = false

    while true
        mi = rng_mod(rng, n_size)  # 0-based electron index
        s = ele_spn[mi+1]             # 1-based Julia array
        t = 1 - s
        ri = ele_idx[mi+1]            # 0-based site index
        if loc_spn[ri+1] == 1         # 1-based Julia array
            break
        end
    end
    rj = ri

    return mi, ri, rj, s, t, flag
end

function make_candidate_exchange_fsz(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_spn::Vector{Int},
    loc_spn::Vector{Int},
    n_site::Int,
    n_size::Int,
    rng::AbstractRNG,
)
    # Check if exchange is possible
    # ele_num uses 1-based Julia indexing
    flag = true
    spn_0 = 0
    spn_1 = 0
    for ri = 1:n_site
        if ele_num[ri] + ele_num[ri+n_site] == 1
            if spn_0 == 0
                spn_0 = 2 * ele_num[ri] - 1
            else
                spn_1 = 2 * ele_num[ri] - 1
            end
            if spn_0 * spn_1 < 0
                flag = false
                break
            end
        end
    end

    if flag
        return 0, 0, 0, 0, true
    end

    # Select electron for exchange
    # mi is 0-based electron index, ri is 0-based site index
    mi = 0
    s = 0
    ri = 0
    while true
        mi = rng_mod(rng, n_size)  # 0-based electron index
        s = ele_spn[mi+1]             # 1-based Julia array
        ri = ele_idx[mi+1]            # 0-based site index
        # ele_cfg uses 1-based Julia indexing
        if ele_cfg[ri+1+(1-s)*n_site] == -1
            break
        end
    end

    t = 1 - s
    mj = 0
    rj = 0
    while true
        mj = rng_mod(rng, n_size)  # 0-based electron index
        rj = ele_idx[mj+1]            # 0-based site index
        # ele_cfg uses 1-based Julia indexing
        if ele_cfg[rj+1+(1-t)*n_site] == -1 && ele_spn[mj+1] == t
            break
        end
    end

    return mi, ri, rj, s, false
end

function update_proj_cnt_fsz!(
    ri::Int,
    rj::Int,
    s::Int,
    t::Int,
    proj_cnt_new::Vector{Int},
    proj_cnt_old::Vector{Int},
    ele_num::Vector{Int},
    data::ExpertModeData,
)
    n_site = data.modpara.nsite
    n_proj = length(proj_cnt_new)
    n_gutzwiller_idx = data.n_gutzwiller_idx
    gutzwiller_idx = data.gutzwiller_idx

    # Copy old counts
    if proj_cnt_new !== proj_cnt_old
        copy!(proj_cnt_new, proj_cnt_old)
    end

    if ri == rj
        return
    end

    # Get up-spin and down-spin electron numbers
    n0 = @view ele_num[1:n_site]  # up-spin
    n1 = @view ele_num[(n_site+1):(2*n_site)]  # down-spin

    # Gutzwiller factor (C: UpdateProjCnt_fsz)
    if n_gutzwiller_idx > 0
        if !isempty(gutzwiller_idx) &&
           ri + 1 <= length(gutzwiller_idx) &&
           rj + 1 <= length(gutzwiller_idx)
            idx_ri = gutzwiller_idx[ri+1]  # 0-based idx value
            idx_rj = gutzwiller_idx[rj+1]
            if idx_ri + 1 <= n_proj
                proj_cnt_new[idx_ri+1] -= n0[ri+1] + n1[ri+1]
            end
            if idx_rj + 1 <= n_proj
                proj_cnt_new[idx_rj+1] += n0[rj+1] * n1[rj+1]
            end
        else
            proj_cnt_new[1] -= n0[ri+1] + n1[ri+1]
            proj_cnt_new[1] += n0[rj+1] * n1[rj+1]
        end
    end

    # Jastrow factor
    offset = n_gutzwiller_idx
    jastrow_idx_matrix = data.jastrow_idx
    if data.n_jastrow_idx > 0 && !isempty(jastrow_idx_matrix)
        function get_jastrow_idx(ra::Int, rb::Int)
            if ra < rb
                return jastrow_idx_matrix[ra+1, rb+1]
            else
                return jastrow_idx_matrix[rb+1, ra+1]
            end
        end

        idx = get_jastrow_idx(ri, rj)
        proj_idx = offset + idx + 1
        if proj_idx <= n_proj
            proj_cnt_new[proj_idx] += (n0[ri+1] + n1[ri+1]) - (n0[rj+1] + n1[rj+1]) + 1
        end

        for rk = 0:(n_site-1)
            if rk == rj || rk == ri
                continue
            end
            idx = get_jastrow_idx(ri, rk)
            proj_idx = offset + idx + 1
            if proj_idx <= n_proj
                proj_cnt_new[proj_idx] -= (n0[rk+1] + n1[rk+1] - 1)
            end
        end

        for rk = 0:(n_site-1)
            if rk == ri || rk == rj
                continue
            end
            idx = get_jastrow_idx(rj, rk)
            proj_idx = offset + idx + 1
            if proj_idx <= n_proj
                proj_cnt_new[proj_idx] += (n0[rk+1] + n1[rk+1] - 1)
            end
        end
    elseif length(data.jastrow_terms) > 0
        jastrow_idx = offset + 1
        if jastrow_idx <= n_proj
            proj_cnt_new[jastrow_idx] += (n0[ri+1] + n1[ri+1]) - (n0[rj+1] + n1[rj+1]) + 1
            for rk = 0:(n_site-1)
                if rk == rj || rk == ri
                    continue
                end
                proj_cnt_new[jastrow_idx] -= (n0[rk+1] + n1[rk+1] - 1)
            end
            for rk = 0:(n_site-1)
                if rk == ri || rk == rj
                    continue
                end
                proj_cnt_new[jastrow_idx] += (n0[rk+1] + n1[rk+1] - 1)
            end
        end
    end

    # Doublon-Holon terms are not yet implemented (no precomputed DH index tables).
end

function copy_from_burn_sample_fsz!(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    ele_spn::Vector{Int},
    state::VMCOptimizationState,
)
    n_size = length(ele_idx)
    n_site2 = length(ele_cfg)
    n_proj = length(ele_proj_cnt)

    burn = state.electron_config.burn_ele_idx
    offset = 0
    copyto!(ele_idx, 1, burn, offset + 1, n_size)
    offset += n_size
    copyto!(ele_cfg, 1, burn, offset + 1, n_site2)
    offset += n_site2
    copyto!(ele_num, 1, burn, offset + 1, n_site2)
    offset += n_site2
    copyto!(ele_proj_cnt, 1, burn, offset + 1, n_proj)
    offset += n_proj
    copyto!(ele_spn, 1, burn, offset + 1, n_size)
end

function copy_to_burn_sample_fsz!(
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    ele_spn::Vector{Int},
    state::VMCOptimizationState,
)
    n_size = length(ele_idx)
    n_site2 = length(ele_cfg)
    n_proj = length(ele_proj_cnt)

    burn = state.electron_config.burn_ele_idx
    if length(burn) < n_size + n_site2 + n_site2 + n_proj + n_size
        resize!(burn, n_size + n_site2 + n_site2 + n_proj + n_size)
    end

    offset = 0
    copyto!(burn, offset + 1, ele_idx, 1, n_size)
    offset += n_size
    copyto!(burn, offset + 1, ele_cfg, 1, n_site2)
    offset += n_site2
    copyto!(burn, offset + 1, ele_num, 1, n_site2)
    offset += n_site2
    copyto!(burn, offset + 1, ele_proj_cnt, 1, n_proj)
    offset += n_proj
    copyto!(burn, offset + 1, ele_spn, 1, n_size)
end

function save_ele_config_fsz!(
    sample::Int,
    log_ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    ele_spn::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_size = length(ele_idx)
    n_site2 = length(ele_cfg)
    n_proj = length(ele_proj_cnt)

    # Save to state arrays
    # ele_idx already stores 0-based site indices
    offset = sample * n_size
    copyto!(state.electron_config.ele_idx, offset + 1, ele_idx, 1, n_size)

    offset = sample * n_site2
    copyto!(state.electron_config.ele_cfg, offset + 1, ele_cfg, 1, n_site2)
    copyto!(state.electron_config.ele_num, offset + 1, ele_num, 1, n_site2)

    offset = sample * n_proj
    copyto!(state.electron_config.ele_proj_cnt, offset + 1, ele_proj_cnt, 1, n_proj)

    offset = sample * n_size
    copyto!(state.electron_config.ele_spn, offset + 1, ele_spn, 1, n_size)

    # Note: log_sq_pf_full_slater is not used in the Julia implementation
    # The reweighting is handled differently
end

function update_m_all_two_fsz!(
    mi::Int,
    s::Int,
    mj::Int,
    t::Int,
    ri::Int,
    rj::Int,
    ele_idx::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    # Full recalculation for two-electron update
    calculate_m_all_fsz!(ele_idx, ele_spn, qp_start, qp_end, data, state)
end

"""
    calculate_m_all_real!(ele_idx::Vector{Int}, qp_start::Int, qp_end::Int,
                          data::ExpertModeData, state::VMCOptimizationState)::Int

Calculate M matrices using real arithmetic.
Equivalent to C's `CalculateMAll_real()`.

Uses real implementation directly (not complex version).
This matches the C implementation which uses `SlaterElm_real` and real LAPACK functions.
"""
function calculate_m_all_real!(
    ele_idx::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)::Int
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec  # Total number of electrons (2*Ne)
    n_site2 = 2 * n_site
    n_qp_full = length(state.slater_matrix.pf_m_real)

    # Validate qp_start and qp_end (1-based indexing)
    if qp_start < 1 || qp_end > n_qp_full + 1 || qp_start >= qp_end
        @error "Invalid qp_start or qp_end: qp_start=$qp_start, qp_end=$qp_end, n_qp_full=$n_qp_full"
        return 1
    end

    # slater_elm_real: 1D array with C row-major layout
    # [qp_idx * n_site2 * n_site2 + rsi * n_site2 + rsj]
    # We pass the full 1D array and let PfaPack handle the offset

    # Calculate the offset for the qp_start range
    qp_offset = (qp_start - 1) * n_site2 * n_site2
    qp_num = qp_end - qp_start

    # Create a view of slater_elm_real for the relevant QP range
    slater_elm_real_subset = view(
        state.slater_matrix.slater_elm_real,
        (qp_offset+1):(qp_offset+qp_num*n_site2*n_site2),
    )

    # Use pre-allocated workspace arrays instead of allocating new ones
    # Create views for the required qp_num (workspace is allocated for n_qp_full)
    ws = state.workspace
    inv_m_real_temp = view(ws.inv_m_real_temp,:,:,(1:qp_num))
    pf_m_real_temp = view(ws.pf_m_real_temp, 1:qp_num)

    # Zero out the workspace views before use
    fill!(inv_m_real_temp, 0.0)
    fill!(pf_m_real_temp, 0.0)

    # Call PfaPack's calculate_m_all_real!
    # Note: We use the imported function from PfaPack
    # PfaPack expects n_elec to be total number of electrons (2*Ne), not Ne
    info = calculate_m_all_real_pfapack!(
        ele_idx,
        slater_elm_real_subset,
        inv_m_real_temp,
        pf_m_real_temp,
        1,  # qp_start (1-based, relative to subset)
        qp_num + 1,  # qp_end (1-based, exclusive, relative to subset)
        n_site,
        n_size,  # Pass n_size (2*Ne) instead of n_elec (Ne)
        ws.pfapack_workspace,  # Pre-allocated workspace
    )

    if info != 0
        return info
    end

    # Copy results back to state.slater_matrix
    # pf_m_real: 1D array indexed by qp (1-based)
    @inbounds for qp_local = 1:qp_num
        qp_global = qp_start + qp_local - 1
        if qp_global <= length(state.slater_matrix.pf_m_real)
            state.slater_matrix.pf_m_real[qp_global] = pf_m_real_temp[qp_local]
        end
    end

    # inv_m_real: 1D array with C row-major layout
    # The column-major inv_m_real_temp can be copied directly to the 1D array
    # because the memory layout is equivalent when reshaped
    # Use copyto! for efficient vectorized copy instead of triple nested loops
    nsq = n_size * n_size
    @inbounds for qp_local = 1:qp_num
        qp_global = qp_start + qp_local - 1
        dst_start = (qp_global - 1) * nsq + 1
        dst_end = dst_start + nsq - 1
        if dst_end <= length(state.slater_matrix.inv_m_real)
            dst = view(state.slater_matrix.inv_m_real, dst_start:dst_end)
            src = view(inv_m_real_temp,:,:,qp_local)
            copyto!(dst, vec(src))
        end
    end

    return info
end

"""
    calculate_m_all_fsz_real!(ele_idx::Vector{Int}, ele_spn::Vector{Int},
                              qp_start::Int, qp_end::Int,
                              data::ExpertModeData, state::VMCOptimizationState)::Int

Calculate M matrices using real arithmetic for FSZ (free Sz).
Equivalent to C's `CalculateMAll_fsz_real()`.
"""
function calculate_m_all_fsz_real!(
    ele_idx::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)::Int
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_qp_full = length(state.slater_matrix.pf_m_real)

    if qp_start < 1 || qp_end > n_qp_full + 1 || qp_start >= qp_end
        @error "Invalid qp_start or qp_end: qp_start=$qp_start, qp_end=$qp_end, n_qp_full=$n_qp_full"
        return 1
    end

    info = calculate_m_all_fsz!(ele_idx, ele_spn, qp_start, qp_end, data, state)
    if info != 0
        return info
    end

    qp_num = qp_end - qp_start
    @inbounds for qp_local = 1:qp_num
        qp_global = qp_start + qp_local - 1
        if qp_global <= length(state.slater_matrix.pf_m_real)
            state.slater_matrix.pf_m_real[qp_global] = real(state.slater_matrix.pf_m[qp_global])
        end
    end

    nsq = n_size * n_size
    @inbounds for qp_local = 1:qp_num
        qp_global = qp_start + qp_local - 1
        dst_start = (qp_global - 1) * nsq + 1
        dst_end = dst_start + nsq - 1
        if dst_end <= length(state.slater_matrix.inv_m_real)
            dst = view(state.slater_matrix.inv_m_real, dst_start:dst_end)
            src = view(state.slater_matrix.inv_m, dst_start:dst_end)
            for i = 1:nsq
                dst[i] = real(src[i])
            end
        end
    end

    return 0
end

"""
    calculate_log_ip_real(pf_m_real::Vector{Float64}, qp_start::Int, qp_end::Int,
                          data::ExpertModeData)::Float64

Calculate logarithm of inner product using real arithmetic.
Equivalent to C's `CalculateLogIP_real()`.
"""
function calculate_log_ip_real(
    pf_m_real::Vector{Float64},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
)::Float64
    if data.qp_weights === nothing
        @error "Quantum projection weights not initialized. Call init_qp_weight!(data) first."
        return log(1e-100)
    end

    qp_num = qp_end - qp_start
    qp_full_weight::Vector{ComplexF64} = data.qp_weights.qp_full_weight

    if qp_start < 1 || qp_end > length(qp_full_weight) + 1 || qp_start >= qp_end
        @error "Invalid qp_start or qp_end: qp_start=$qp_start, qp_end=$qp_end"
        return log(1e-100)
    end

    # Calculate inner product: ip = sum(QPFullWeight[qpidx] * pfM_real[qpidx])
    ip = 0.0
    for qpidx = 1:qp_num
        qp_idx = qp_start + qpidx - 1
        if qp_idx <= length(qp_full_weight) && qp_idx <= length(pf_m_real)
            ip += real(qp_full_weight[qp_idx]) * pf_m_real[qp_idx]
        end
    end

    return log(abs(ip) + 1e-100)
end

"""
    calculate_ip_real(pf_m_real::Vector{Float64}, qp_start::Int, qp_end::Int,
                     data::ExpertModeData)::Float64

Calculate inner product using real arithmetic.
Equivalent to C's `CalculateIP_real()`.
"""
function calculate_ip_real(
    pf_m_real::Vector{Float64},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
)::Float64
    if data.qp_weights === nothing
        @error "Quantum projection weights not initialized. Call init_qp_weight!(data) first."
        return 0.0
    end

    qp_num = qp_end - qp_start
    qp_full_weight::Vector{ComplexF64} = data.qp_weights.qp_full_weight

    if qp_start < 1 || qp_end > length(qp_full_weight) + 1 || qp_start >= qp_end
        @error "Invalid qp_start or qp_end: qp_start=$qp_start, qp_end=$qp_end"
        return 0.0
    end

    # Calculate inner product: ip = sum(QPFullWeight[qpidx] * pfM_real[qpidx])
    ip = 0.0
    @inbounds for qpidx = 1:qp_num
        qp_idx = qp_start + qpidx - 1
        if qp_idx <= length(qp_full_weight) && qp_idx <= length(pf_m_real)
            ip += real(qp_full_weight[qp_idx]) * pf_m_real[qp_idx]
        end
    end

    return ip
end

"""
    calculate_new_pf_m2_real!(ma::Int, s::Int, pf_m_new_real::Vector{Float64},
                              ele_idx::Vector{Int}, qp_start::Int, qp_end::Int,
                              data::ExpertModeData, state::VMCOptimizationState)

Calculate new Pfaffian using real arithmetic.
Equivalent to C's `CalculateNewPfM2_real()`.
"""
function calculate_new_pf_m2_real!(
    ma::Int,
    s::Int,
    pf_m_new_real::Vector{Float64},
    ele_idx::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    msa = ma + s * n_elec  # 0-based electron index
    rsa = ele_idx[msa+1] + s * n_site  # 0-based site index

    # Use slater_elm_real and inv_m_real
    slater_elm_real = state.slater_matrix.slater_elm_real
    inv_m_real = state.slater_matrix.inv_m_real
    pf_m_real = state.slater_matrix.pf_m_real

    for qpidx = qp_start:(qp_end-1)
        if qpidx < 1 || qpidx > length(pf_m_real)
            continue
        end

        # Calculate ratio = sum_j invM_real[msa][msj] * sltE_real[rsa][rsj]
        # C: invM_a = InvM_real + qpidx*Nsize*Nsize + msa*Nsize
        #    invM_a[msj] = InvM_real[(qpidx+qpStart)*Nsize*Nsize + msa*Nsize + msj]
        # Julia: inv_m_real is stored as [qp * n_size * n_size + msi * n_size + msj] (1-based)
        #        So inv_m_real[(qpidx-1)*n_size*n_size + msa*n_size + msj + 1]
        #        But msa is 0-based, so we need: inv_m_real[(qpidx-1)*n_size*n_size + (msa+1-1)*n_size + msj + 1]
        #        = inv_m_real[(qpidx-1)*n_size*n_size + msa*n_size + msj + 1]
        ratio = 0.0
        inv_base = (qpidx - 1) * n_size * n_size + msa * n_size
        slt_base = (qpidx - 1) * n_site2 * n_site2 + rsa * n_site2
        @inbounds for msj = 0:(n_size-1)
            rsj = if msj < n_elec
                ele_idx[msj+1]
            else
                ele_idx[msj+1] + n_site
            end

            # Linear index for inv_m_real: inv_base + msj + 1 (msj is 0-based, +1 for 1-based indexing)
            inv_idx = inv_base + msj + 1
            # Linear index for slater_elm_real: slt_base + rsj + 1 (rsj is 0-based, +1 for 1-based indexing)
            slt_idx = slt_base + rsj + 1

            if inv_idx <= length(inv_m_real) && slt_idx <= length(slater_elm_real)
                ratio += inv_m_real[inv_idx] * slater_elm_real[slt_idx]
            end
        end

        # pfMNew_real[qpidx] = -ratio * PfM_real[qpidx]
        pf_m_new_real[qpidx] = -ratio * pf_m_real[qpidx]
    end
end

"""
    calculate_new_pf_m2_fsz_real!(ma::Int, s::Int, pf_m_new_real::Vector{Float64},
                                  ele_idx::Vector{Int}, ele_spn::Vector{Int},
                                  qp_start::Int, qp_end::Int,
                                  data::ExpertModeData, state::VMCOptimizationState)

Calculate new Pfaffian using real arithmetic for FSZ.
Equivalent to C's `CalculateNewPfM2_fsz_real()`.
"""
function calculate_new_pf_m2_fsz_real!(
    ma::Int,
    s::Int,
    pf_m_new_real::Vector{Float64},
    ele_idx::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    msa = ma
    rsa = ele_idx[msa+1] + s * n_site

    slater_elm_real = state.slater_matrix.slater_elm_real
    inv_m_real = state.slater_matrix.inv_m_real
    pf_m_real = state.slater_matrix.pf_m_real

    for qpidx = qp_start:(qp_end-1)
        if qpidx < 1 || qpidx > length(pf_m_real)
            continue
        end
        inv_base = (qpidx - 1) * n_size * n_size + msa * n_size
        slt_base = (qpidx - 1) * n_site2 * n_site2 + rsa * n_site2

        ratio = 0.0
        @inbounds for msj = 0:(n_size-1)
            rsj = ele_idx[msj+1] + ele_spn[msj+1] * n_site
            inv_idx = inv_base + msj + 1
            slt_idx = slt_base + rsj + 1
            ratio += inv_m_real[inv_idx] * slater_elm_real[slt_idx]
        end

        pf_m_new_real[qpidx] = -ratio * pf_m_real[qpidx]
    end
end

"""
    calculate_new_pf_m_two2_real!(ma::Int, s::Int, mb::Int, t::Int, pf_m_new_real::Vector{Float64},
                                  ele_idx::Vector{Int}, qp_start::Int, qp_end::Int,
                                  data::ExpertModeData, state::VMCOptimizationState)

Calculate new Pfaffian for two-electron hopping using real arithmetic.
Equivalent to C's `CalculateNewPfMTwo2_real()`.
"""
function calculate_new_pf_m_two2_real!(
    ma::Int,
    s::Int,
    mb::Int,
    t::Int,
    pf_m_new_real::Vector{Float64},
    ele_idx::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    msa = ma + s * n_elec  # 0-based electron index
    msb = mb + t * n_elec  # 0-based electron index

    # If same electron, use single hop calculation
    if msa == msb
        calculate_new_pf_m2_real!(
            mb,
            t,
            pf_m_new_real,
            ele_idx,
            qp_start,
            qp_end,
            data,
            state,
        )
        return
    end

    # Use slater_elm_real and inv_m_real
    slater_elm_real = state.slater_matrix.slater_elm_real
    inv_m_real = state.slater_matrix.inv_m_real
    pf_m_real = state.slater_matrix.pf_m_real

    # Work arrays
    vec_a = zeros(Float64, n_size)
    vec_b = zeros(Float64, n_size)

    qp_num = qp_end - qp_start
    for local_qpidx = 0:(qp_num-1)
        qpidx = qp_start + local_qpidx  # Global index (1-based)
        if qpidx < 1 || qpidx > length(pf_m_real)
            continue
        end

        rsa = ele_idx[msa+1] + s * n_site  # 0-based site index
        rsb = ele_idx[msb+1] + t * n_site  # 0-based site index

        slt_offset = (qpidx - 1) * n_site2 * n_site2
        # inv_m_real uses global index: [qp * n_size * n_size + msi * n_size + msj]
        inv_offset = (qpidx - 1) * n_size * n_size

        # Calculate vec_a[i] = sltE_real[rsa][rsi] and vec_b[i] = sltE_real[rsb][rsi]
        @turbo for msi = 0:(n_size-1)
            #=
            if msi < n_elec
                rsi = ele_idx[msi + 1]  # up-spin (0-based)
            else
                rsi = ele_idx[msi + 1] + n_site  # down-spin (0-based)
            end
            =#
            rsi = ifelse(msi < n_elec, ele_idx[msi+1], ele_idx[msi+1] + n_site)

            vec_a[msi+1] = slater_elm_real[slt_offset+rsa*n_site2+rsi+1]
            vec_b[msi+1] = slater_elm_real[slt_offset+rsb*n_site2+rsi+1]
        end
        vec_ba = vec_b[msa+1]

        # Calculate p_a, p_b, q_a, q_b
        p_a = p_b = q_a = q_b = 0.0
        @inbounds for msi = 0:(n_size-1)
            inv_m_a_msi = inv_m_real[inv_offset+msa*n_size+msi+1]
            inv_m_b_msi = inv_m_real[inv_offset+msb*n_size+msi+1]

            p_a += inv_m_a_msi * vec_a[msi+1]
            p_b += inv_m_b_msi * vec_a[msi+1]
            q_a += inv_m_a_msi * vec_b[msi+1]
            q_b += inv_m_b_msi * vec_b[msi+1]
        end

        # invM_ab = invM_a[msb]
        inv_m_ab = inv_m_real[inv_offset+msa*n_size+msb+1]

        # Calculate bMa = sum_i vec_b[i] * (sum_j invM[i][j] * vec_a[j])
        bMa = 0.0
        @turbo for msi = 0:(n_size-1)
            tmp = 0.0
            for msj = 0:(n_size-1)
                tmp += inv_m_real[inv_offset+msi*n_size+msj+1] * vec_a[msj+1]
            end
            bMa += vec_b[msi+1] * tmp
        end

        # Calculate ratio = PfMNew / PfMOld
        ratio = inv_m_ab * vec_ba + inv_m_ab * bMa + p_a * q_b - p_b * q_a

        # Update pfMNew
        pf_m_new_real[qpidx] = ratio * pf_m_real[qpidx]
    end
end

"""
    calculate_new_pf_m_two2_fsz_real!(ma::Int, s::Int, mb::Int, t::Int,
                                      pf_m_new_real::Vector{Float64},
                                      ele_idx::Vector{Int}, ele_spn::Vector{Int},
                                      qp_start::Int, qp_end::Int,
                                      data::ExpertModeData, state::VMCOptimizationState)

Calculate new Pfaffian for two-electron hopping (FSZ, real).
Equivalent to C's `CalculateNewPfMTwo2_fsz_real()`.
"""
function calculate_new_pf_m_two2_fsz_real!(
    ma::Int,
    s::Int,
    mb::Int,
    t::Int,
    pf_m_new_real::Vector{Float64},
    ele_idx::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    msa = ma
    msb = mb

    if msa == msb
        calculate_new_pf_m2_fsz_real!(
            mb,
            t,
            pf_m_new_real,
            ele_idx,
            ele_spn,
            qp_start,
            qp_end,
            data,
            state,
        )
        return
    end

    slater_elm_real = state.slater_matrix.slater_elm_real
    inv_m_real = state.slater_matrix.inv_m_real
    pf_m_real = state.slater_matrix.pf_m_real

    vec_a = Vector{Float64}(undef, n_size)
    vec_b = Vector{Float64}(undef, n_size)

    qp_num = qp_end - qp_start
    for local_qpidx = 0:(qp_num-1)
        qpidx = qp_start + local_qpidx
        if qpidx < 1 || qpidx > length(pf_m_real)
            continue
        end
        rsa = ele_idx[msa+1] + s * n_site
        rsb = ele_idx[msb+1] + t * n_site

        slt_offset = (qpidx - 1) * n_site2 * n_site2
        inv_offset = (qpidx - 1) * n_size * n_size

        @inbounds for msi = 0:(n_size-1)
            rsi = ele_idx[msi+1] + ele_spn[msi+1] * n_site
            vec_a[msi+1] = slater_elm_real[slt_offset+rsa*n_site2+rsi+1]
            vec_b[msi+1] = slater_elm_real[slt_offset+rsb*n_site2+rsi+1]
        end
        vec_ba = vec_b[msa+1]

        p_a = p_b = q_a = q_b = bMa = 0.0
        @inbounds for msi = 0:(n_size-1)
            inv_m_ai = inv_m_real[inv_offset+msa*n_size+msi+1]
            inv_m_bi = inv_m_real[inv_offset+msb*n_size+msi+1]
            vec_ai = vec_a[msi+1]
            vec_bi = vec_b[msi+1]

            p_a += inv_m_ai * vec_ai
            p_b += inv_m_bi * vec_ai
            q_a += inv_m_ai * vec_bi
            q_b += inv_m_bi * vec_bi
        end

        @inbounds for msi = 0:(n_size-1)
            tmp = 0.0
            base_i = inv_offset + msi * n_size
            for msj = 0:(n_size-1)
                tmp += inv_m_real[base_i+msj+1] * vec_a[msj+1]
            end
            bMa += vec_b[msi+1] * tmp
        end

        inv_m_ab = inv_m_real[inv_offset+msa*n_size+msb+1]
        ratio = inv_m_ab * vec_ba + inv_m_ab * bMa + p_a * q_b - p_b * q_a
        pf_m_new_real[qpidx] = ratio * pf_m_real[qpidx]
    end
end

"""
    update_m_all_real!(ma::Int, s::Int, ele_idx::Vector{Int}, qp_start::Int, qp_end::Int,
                       data::ExpertModeData, state::VMCOptimizationState)

Update inverse matrix and Pfaffian using direct real arithmetic (Sherman-Morrison O(N²)).
Equivalent to C's `UpdateMAll_real()`.

This is a direct real implementation that avoids the overhead of complex arithmetic.
The algorithm uses Sherman-Morrison formula for efficient O(N²) inverse matrix update.

# Reference
- C implementation: mVMC/src/mVMC/pfupdate_real.c:147-220 (updateMAll_child_real)
"""
function update_m_all_real!(
    ma::Int,
    s::Int,
    ele_idx::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    msa = ma + s * n_elec  # 0-based electron index
    rsa = ele_idx[msa+1] + s * n_site  # 0-based site index (ele_idx is 1-based)

    # Get array references
    slater_elm_real = state.slater_matrix.slater_elm_real
    inv_m_real = state.slater_matrix.inv_m_real
    pf_m_real = state.slater_matrix.pf_m_real

    # Work arrays
    vec1 = zeros(Float64, n_size)
    vec2 = zeros(Float64, n_size)

    # Process each QP index
    for qpidx = qp_start:(qp_end-1)
        if qpidx < 1 || qpidx > length(pf_m_real)
            continue
        end

        # C: sltE = SlaterElm_real + (qpidx+qpStart)*Nsite2*Nsite2
        slt_offset = (qpidx - 1) * n_site2 * n_site2
        # C: invM = InvM_real + qpidx*Nsize*Nsize
        inv_offset = (qpidx - 1) * n_size * n_size

        # Initialize vec1 to zero
        fill!(vec1, 0.0)

        # Calculate vec1[i] = sum_j (-invM[j][i] * sltE[a][rsj])
        # Note: invM[i][j] = -invM[j][i] (skew-symmetric)
        # C: vec1[msi] += -invM_j[msi] * sltE_aj
        @inbounds for msj = 0:(n_size-1)
            # rsj = eleIdx[msj] + (msj/Ne)*Nsite
            rsj = msj < n_elec ? ele_idx[msj+1] : (ele_idx[msj+1] + n_site)
            slt_e_aj = slater_elm_real[slt_offset+rsa*n_site2+rsj+1]

            for msi = 0:(n_size-1)
                # invM_j[msi] = invM[msi][msj] stored at inv_offset + msj*n_size + msi
                # But C stores row-major: invM[j][i] at invM + j*Nsize + i
                inv_m_ji = inv_m_real[inv_offset+msj*n_size+msi+1]
                vec1[msi+1] += -inv_m_ji * slt_e_aj
            end
        end

        # Update Pfaffian
        # C: PfM_real[qpidx] *= -vec1[msa]
        tmp = vec1[msa+1]
        pf_m_real[qpidx] *= -tmp
        inv_vec1_a = -1.0 / tmp

        # Calculate vec2[i] = invM[a][i] * (-1/vec1[a])
        # C: vec2[msi] = invM_a[msi] * invVec1_a
        @inbounds for msi = 0:(n_size-1)
            inv_m_ai = inv_m_real[inv_offset+msa*n_size+msi+1]
            vec2[msi+1] = inv_m_ai * inv_vec1_a
        end

        # Update InvM using Sherman-Morrison formula
        # C: invM_i[msj] += vec1_i * vec2[msj] - vec1[msj] * vec2_i
        @inbounds for msi = 0:(n_size-1)
            vec1_i = vec1[msi+1]
            vec2_i = vec2[msi+1]

            for msj = 0:(n_size-1)
                inv_m_real[inv_offset+msi*n_size+msj+1] +=
                    vec1_i * vec2[msj+1] - vec1[msj+1] * vec2_i
            end

            # C: invM_i[msa] -= vec2_i
            inv_m_real[inv_offset+msi*n_size+msa+1] -= vec2_i
        end

        # C: invM_a[msj] += vec2[msj]
        @inbounds for msj = 0:(n_size-1)
            inv_m_real[inv_offset+msa*n_size+msj+1] += vec2[msj+1]
        end
    end
end

"""
    update_m_all_fsz_real!(ma::Int, s::Int, ele_idx::Vector{Int}, ele_spn::Vector{Int},
                           qp_start::Int, qp_end::Int,
                           data::ExpertModeData, state::VMCOptimizationState)

Update inverse matrix and Pfaffian using real arithmetic for FSZ.
Equivalent to C's `UpdateMAll_fsz_real()`.
"""
function update_m_all_fsz_real!(
    ma::Int,
    s::Int,
    ele_idx::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    msa = ma
    rsa = ele_idx[msa+1] + s * n_site

    slater_elm_real = state.slater_matrix.slater_elm_real
    inv_m_real = state.slater_matrix.inv_m_real
    pf_m_real = state.slater_matrix.pf_m_real

    vec1 = Vector{Float64}(undef, n_size)
    vec2 = Vector{Float64}(undef, n_size)

    for qpidx = qp_start:(qp_end-1)
        if qpidx < 1 || qpidx > length(pf_m_real)
            continue
        end
        inv_offset = (qpidx - 1) * n_size * n_size
        slt_offset = (qpidx - 1) * n_site2 * n_site2 + rsa * n_site2

        fill!(vec1, 0.0)
        @inbounds for msj = 0:(n_size-1)
            rsj = ele_idx[msj+1] + ele_spn[msj+1] * n_site
            slt_e_aj = slater_elm_real[slt_offset+rsj+1]
            inv_base = inv_offset + msj * n_size
            for msi = 0:(n_size-1)
                vec1[msi+1] += -inv_m_real[inv_base+msi+1] * slt_e_aj
            end
        end

        tmp = vec1[msa+1]
        pf_m_real[qpidx] *= -tmp
        inv_vec1_a = -1.0 / tmp

        @inbounds for msi = 0:(n_size-1)
            inv_m_ai = inv_m_real[inv_offset+msa*n_size+msi+1]
            vec2[msi+1] = inv_m_ai * inv_vec1_a
        end

        @inbounds for msi = 0:(n_size-1)
            vec1_i = vec1[msi+1]
            vec2_i = vec2[msi+1]
            base_i = inv_offset + msi * n_size
            for msj = 0:(n_size-1)
                inv_m_real[base_i+msj+1] +=
                    vec1_i * vec2[msj+1] - vec1[msj+1] * vec2_i
            end
            inv_m_real[base_i+msa+1] -= vec2_i
        end

        inv_base_a = inv_offset + msa * n_size
        @inbounds for msj = 0:(n_size-1)
            inv_m_real[inv_base_a+msj+1] += vec2[msj+1]
        end
    end
end

"""
    update_m_all_two_fsz_real!(ma::Int, s::Int, mb::Int, t::Int,
                               ra_old::Int, rb_old::Int,
                               ele_idx::Vector{Int}, ele_spn::Vector{Int},
                               qp_start::Int, qp_end::Int,
                               data::ExpertModeData, state::VMCOptimizationState)

Update inverse matrix and Pfaffian for two-electron hopping (FSZ, real).
Equivalent to C's `UpdateMAllTwo_fsz_real()`.
"""
function update_m_all_two_fsz_real!(
    ma::Int,
    s::Int,
    mb::Int,
    t::Int,
    ra_old::Int,
    rb_old::Int,
    ele_idx::Vector{Int},
    ele_spn::Vector{Int},
    qp_start::Int,
    qp_end::Int,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    msa = ma
    msb = mb
    rsa = ele_idx[msa+1] + s * n_site
    rsb = ele_idx[msb+1] + t * n_site
    rsa_old = ra_old + s * n_site
    rsb_old = ra_old + t * n_site  # match C behavior

    slater_elm_real = state.slater_matrix.slater_elm_real
    inv_m_real = state.slater_matrix.inv_m_real
    pf_m_real = state.slater_matrix.pf_m_real

    vecP = Vector{Float64}(undef, n_size)
    vecQ = Vector{Float64}(undef, n_size)
    vecS = Vector{Float64}(undef, n_size)
    vecT = Vector{Float64}(undef, n_size)

    for qpidx = qp_start:(qp_end-1)
        if qpidx < 1 || qpidx > length(pf_m_real)
            continue
        end
        slt_offset = (qpidx - 1) * n_site2 * n_site2
        inv_offset = (qpidx - 1) * n_size * n_size
        sltE_a_base = slt_offset + rsa * n_site2
        sltE_b_base = slt_offset + rsb * n_site2
        m_old_ab = slater_elm_real[slt_offset+rsa_old*n_site2+rsb_old+1]

        fill!(vecP, 0.0)
        fill!(vecQ, 0.0)

        @inbounds for msi = 0:(n_size-1)
            rsi = ele_idx[msi+1] + ele_spn[msi+1] * n_site
            vecS[msi+1] = slater_elm_real[sltE_a_base+rsi+1]
            vecT[msi+1] = slater_elm_real[sltE_b_base+rsi+1]
        end
        vecS[msb+1] = m_old_ab

        @inbounds for msi = 0:(n_size-1)
            base_i = inv_offset + msi * n_size
            for msj = 0:(n_size-1)
                vecP[msi+1] += inv_m_real[base_i+msj+1] * vecS[msj+1]
                vecQ[msi+1] += inv_m_real[base_i+msj+1] * vecT[msj+1]
            end
        end

        bMa = 0.0
        @inbounds for msi = 0:(n_size-1)
            bMa += vecT[msi+1] * vecP[msi+1]
        end
        inv_m_ab = inv_m_real[inv_offset+msa*n_size+msb+1]
        ratio =
            inv_m_ab * vecT[msa+1] + inv_m_ab * bMa +
            vecP[msa+1] * vecQ[msb+1] - vecP[msb+1] * vecQ[msa+1]
        pf_m_real[qpidx] *= ratio

        a = -vecP[msa+1]
        b = vecP[msb+1]
        c = vecQ[msa+1]
        d = -vecQ[msb+1]
        e = -bMa - vecT[msa+1]
        f = inv_m_real[inv_offset+msa*n_size+msb+1]
        det = a * d - b * c - e * f
        inv_det = 1.0 / det

        @inbounds for msi = 0:(n_size-1)
            vecS[msi+1] = inv_det * inv_m_real[inv_offset+msa*n_size+msi+1]
            vecT[msi+1] = inv_det * inv_m_real[inv_offset+msb*n_size+msi+1]
        end

        @inbounds for msi = 0:(n_size-1)
            base_i = inv_offset + msi * n_size
            p_i = vecP[msi+1]
            q_i = vecQ[msi+1]
            s_i = vecS[msi+1]
            t_i = vecT[msi+1]
            for msj = 0:(n_size-1)
                p_j = vecP[msj+1]
                q_j = vecQ[msj+1]
                s_j = vecS[msj+1]
                t_j = vecT[msj+1]
                inv_m_real[base_i+msj+1] +=
                    a * (q_i * t_j - q_j * t_i) +
                    b * (q_i * s_j - q_j * s_i) +
                    c * (p_i * t_j - p_j * t_i) +
                    d * (p_i * s_j - p_j * s_i) +
                    e * det * (s_i * t_j - s_j * t_i) +
                    f * inv_det * (p_i * q_j - q_i * p_j)
            end
            inv_m_real[base_i+msa+1] += -c * t_i - d * s_i - f * inv_det * q_i
            inv_m_real[base_i+msb+1] += -a * t_i - b * s_i + f * inv_det * p_i
        end

        inv_base_a = inv_offset + msa * n_size
        inv_base_b = inv_offset + msb * n_size
        @inbounds for msj = 0:(n_size-1)
            p_j = vecP[msj+1]
            q_j = vecQ[msj+1]
            s_j = vecS[msj+1]
            t_j = vecT[msj+1]
            inv_m_real[inv_base_a+msj+1] += c * t_j + d * s_j + f * inv_det * q_j
            inv_m_real[inv_base_b+msj+1] += a * t_j + b * s_j - f * inv_det * p_j
        end
        inv_m_real[inv_base_a+msb+1] += f * inv_det
        inv_m_real[inv_base_b+msa+1] -= f * inv_det
    end
end

"""
    vmc_make_sample_real!(data::ExpertModeData, state::VMCOptimizationState, rng::AbstractRNG)

VMC sampling (real version, sz-conserved).
Equivalent to C's `VMCMakeSample_real()`.

Uses real-valued arrays (slater_elm_real, inv_m_real, pf_m_real) for computation.
This is more efficient than the complex version for real-valued problems.
"""
function vmc_make_sample_real!(
    data::ExpertModeData,
    state::VMCOptimizationState,
    rng::AbstractRNG = Random.GLOBAL_RNG,
    c_timer::CTimer = CTIMER_DISABLED,
)
    # Get parameters
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    n_proj = length(data.gutzwiller_terms) + length(data.jastrow_terms)
    n_qp_full = length(state.slater_matrix.pf_m)
    n_vmc_warmup = data.modpara.nvmc_warmup
    n_vmc_sample = data.modpara.nvmc_sample

    # Debug: Verify state.electron_config.ele_idx size
    expected_ele_idx_size = n_vmc_sample * n_size
    actual_ele_idx_size = length(state.electron_config.ele_idx)
    @debug "vmc_make_sample_real!: n_vmc_sample=$n_vmc_sample, n_size=$n_size, expected_ele_idx_size=$expected_ele_idx_size, actual_ele_idx_size=$actual_ele_idx_size"
    if actual_ele_idx_size != expected_ele_idx_size
        @error "vmc_make_sample_real!: ele_idx size mismatch! expected=$expected_ele_idx_size, actual=$actual_ele_idx_size. This will cause BoundsError when saving samples."
    end
    n_vmc_interval = data.modpara.nvmc_interval
    n_ex_update_path = data.modpara.nex_update_path

    # Get temporary arrays from state
    tmp_ele_idx = state.electron_config.tmp_ele_idx
    tmp_ele_cfg = state.electron_config.tmp_ele_cfg
    tmp_ele_num = state.electron_config.tmp_ele_num
    tmp_ele_proj_cnt = state.electron_config.tmp_ele_proj_cnt

    # Use cached loc_spn from workspace (initialize if needed)
    ws = state.workspace
    if all(x -> x == 0, ws.loc_spn) && !isempty(data.locspin_terms)
        init_loc_spn!(state, data)
    end
    loc_spn = ws.loc_spn

    # Initialize random sample
    # burn_flag is passed as a parameter (managed by vmc_para_opt!)
    # Ensure counter array is large enough
    if length(state.electron_config.counter) < 11
        old_len = length(state.electron_config.counter)
        resize!(state.electron_config.counter, 11)
        # CRITICAL: Initialize new elements to zero (resize! leaves them uninitialized!)
        for i = (old_len+1):11
            state.electron_config.counter[i] = 0
        end
    end
    burn_flag = state.electron_config.counter[11] != 0  # Use counter[11] as burn_flag storage
    # [30] makeInitialSample: initial/burn sample + initial CalculateMAll + logIP.
    # (early `return`s below are fatal-error exits; leaving [30] open is fine then.)
    ctimer_start!(c_timer, 30)
    if !burn_flag
        info = make_initial_sample!(
            tmp_ele_idx,
            tmp_ele_cfg,
            tmp_ele_num,
            tmp_ele_proj_cnt,
            data,
            state,
            rng,
        )
        if info != 0
            @error "make_initial_sample! failed with info=$info"
            return
        end

        # Debug: Check tmp_ele_idx after make_initial_sample!
        @debug "vmc_make_sample_real!: After make_initial_sample!, tmp_ele_idx[1:10] = $(tmp_ele_idx[1:min(10, length(tmp_ele_idx))])"
        if all(x -> x == 0 || x == -1, tmp_ele_idx)
            @error "vmc_make_sample_real!: tmp_ele_idx is invalid after make_initial_sample! (all zeros or -1)"
            return
        end
    else
        # Copy from burn sample (from previous step)
        copy_from_burn_sample!(
            tmp_ele_idx,
            tmp_ele_cfg,
            tmp_ele_num,
            tmp_ele_proj_cnt,
            state,
        )
    end

    # Copy slater_elm to slater_elm_real
    # Optimized: use broadcast assignment for better performance
    n_copy = min(
        length(state.slater_matrix.slater_elm),
        length(state.slater_matrix.slater_elm_real),
    )
    @inbounds @simd for i = 1:n_copy
        state.slater_matrix.slater_elm_real[i] = real(state.slater_matrix.slater_elm[i])
    end

    # Initialize QP range (single-process: use all)
    qp_start = 1
    qp_end = n_qp_full + 1

    # Calculate M matrices using real arithmetic
    # Retry loop for initial Pfaffian calculation (similar to complex version)
    max_retries = 100
    retry_count = 0
    info = 1
    while info != 0 && retry_count < max_retries
        info = calculate_m_all_real!(tmp_ele_idx, qp_start, qp_end, data, state)
        if info != 0
            # Regenerate sample if Pfaffian calculation fails
            retry_count += 1
            if retry_count <= 3 || retry_count % 10 == 0
                @warn "calculate_m_all_real! failed with info=$info, retry_count=$retry_count. Regenerating sample..."
            end
            result_init = make_initial_sample!(
                tmp_ele_idx,
                tmp_ele_cfg,
                tmp_ele_num,
                tmp_ele_proj_cnt,
                data,
                state,
                rng,
            )
            if result_init != 0
                @error "Failed to regenerate initial sample after Pfaffian calculation failure"
                return
            end
        end
    end

    if info != 0
        @error "calculate_m_all_real! failed after $max_retries retries with info=$info"
        return
    end

    # Calculate initial log(ip) using real version
    log_ip_old =
        calculate_log_ip_real(state.slater_matrix.pf_m_real, qp_start, qp_end, data)

    if !isfinite(log_ip_old)
        @warn "Initial logIpOld is not finite, remaking sample"
        info = make_initial_sample!(
            tmp_ele_idx,
            tmp_ele_cfg,
            tmp_ele_num,
            tmp_ele_proj_cnt,
            data,
            state,
            rng,
        )
        if info != 0
            return
        end
        info = calculate_m_all_real!(tmp_ele_idx, qp_start, qp_end, data, state)
        if info != 0
            return
        end
        log_ip_old =
            calculate_log_ip_real(state.slater_matrix.pf_m_real, qp_start, qp_end, data)
    end
    ctimer_stop!(c_timer, 30)

    # Sampling loop
    n_out_step = burn_flag ? n_vmc_sample + 1 : n_vmc_warmup + n_vmc_sample
    n_in_step = n_vmc_interval * n_site

    # Reset counters
    fill!(state.electron_config.counter, 0)

    n_accept = 0
    # Use pre-allocated workspace arrays instead of allocating new ones
    proj_cnt_new = ws.proj_cnt_new
    pf_m_new_real = ws.pf_m_new_real
    fill!(proj_cnt_new, 0)  # Reset workspace arrays
    fill!(pf_m_new_real, 0.0)
    n_saved_samples = 0

    for out_step = 0:(n_out_step-1)
        for in_step = 0:(n_in_step-1)
            update_type = get_update_type(n_ex_update_path, data.i_flg_orbital_general, rng)

            if update_type == HOPPING
                state.electron_config.counter[1] += 1

                # [31] make candidate (closed before the reject `continue`)
                ctimer_start!(c_timer, 31)
                mi, ri, rj, s, reject_flag = make_candidate_hopping(
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    n_site,
                    n_elec,
                    loc_spn,
                    rng,
                )
                ctimer_stop!(c_timer, 31)

                if reject_flag != 0
                    continue
                end

                # [32] hopping update (children [60]-[63]); starts after the
                # reject `continue` so the timer is never left open.
                ctimer_start!(c_timer, 32)

                # [60] UpdateProjCnt: electron-config + projection update
                ctimer_start!(c_timer, 60)
                update_ele_config!(
                    mi,
                    ri,
                    rj,
                    s,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    n_site,
                    n_elec,
                )
                update_proj_cnt!(
                    ri,
                    rj,
                    s,
                    proj_cnt_new,
                    tmp_ele_proj_cnt,
                    tmp_ele_num,
                    data,
                )
                ctimer_stop!(c_timer, 60)

                # [61] CalculateNewPfM2
                ctimer_start!(c_timer, 61)
                calculate_new_pf_m2_real!(
                    mi,
                    s,
                    pf_m_new_real,
                    tmp_ele_idx,
                    qp_start,
                    qp_end,
                    data,
                    state,
                )
                ctimer_stop!(c_timer, 61)

                # [62] CalculateLogIP
                ctimer_start!(c_timer, 62)
                log_ip_new = calculate_log_ip_real(pf_m_new_real, qp_start, qp_end, data)
                ctimer_stop!(c_timer, 62)

                # Metropolis acceptance/rejection
                x = log_proj_ratio(proj_cnt_new, tmp_ele_proj_cnt, data)
                w = exp(2.0 * (x + log_ip_new - log_ip_old))
                if !isfinite(w)
                    w = -1.0  # Should be rejected
                end

                if w > rng_real2(rng)
                    # Accept
                    # [63] UpdateMAll
                    ctimer_start!(c_timer, 63)
                    update_m_all_real!(
                        mi,
                        s,
                        tmp_ele_idx,
                        qp_start,
                        qp_end,
                        data,
                        state,
                    )
                    ctimer_stop!(c_timer, 63)
                    copy!(tmp_ele_proj_cnt, proj_cnt_new)
                    state.slater_matrix.pf_m_real .= pf_m_new_real
                    log_ip_old = log_ip_new
                    n_accept += 1
                    state.electron_config.counter[2] += 1
                else
                    # Reject
                    revert_ele_config!(
                        mi,
                        ri,
                        rj,
                        s,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        n_site,
                        n_elec,
                    )
                end
                ctimer_stop!(c_timer, 32)

            elseif update_type == EXCHANGE
                # Exchange update: two electrons exchange positions
                # Uses O(N²) Sherman-Morrison update instead of O(N³) full recalculation

                # [31] make candidate (closed before the reject `continue`)
                ctimer_start!(c_timer, 31)
                mi, ri, mj, rj, s, t, reject_flag = make_candidate_exchange(
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    n_site,
                    n_elec,
                    tmp_ele_num,
                    rng,
                )
                ctimer_stop!(c_timer, 31)

                if reject_flag != 0
                    continue
                end

                # [33] exchange update (children [65]-[68]); starts after the
                # reject `continue` so the timer is never left open.
                ctimer_start!(c_timer, 33)

                # Store old positions for Sherman-Morrison update
                ri_old = ri
                rj_old = rj

                # [65] UpdateProjCnt: both electrons' config + projection update
                ctimer_start!(c_timer, 65)
                # Update electron configuration (first electron)
                update_ele_config!(
                    mi,
                    ri,
                    rj,
                    s,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    n_site,
                    n_elec,
                )
                update_proj_cnt!(
                    ri,
                    rj,
                    s,
                    proj_cnt_new,
                    tmp_ele_proj_cnt,
                    tmp_ele_num,
                    data,
                )

                # Update electron configuration (second electron)
                update_ele_config!(
                    mj,
                    rj,
                    ri,
                    t,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    n_site,
                    n_elec,
                )
                update_proj_cnt!(rj, ri, t, proj_cnt_new, proj_cnt_new, tmp_ele_num, data)
                ctimer_stop!(c_timer, 65)

                # [66] CalculateNewPfMTwo2: O(N²) algorithm (no inv_m update yet)
                ctimer_start!(c_timer, 66)
                calculate_new_pf_m_two2_real!(
                    mi,
                    s,
                    mj,
                    t,
                    pf_m_new_real,
                    tmp_ele_idx,
                    qp_start,
                    qp_end,
                    data,
                    state,
                )
                ctimer_stop!(c_timer, 66)

                # [67] CalculateLogIP
                ctimer_start!(c_timer, 67)
                log_ip_new = calculate_log_ip_real(pf_m_new_real, qp_start, qp_end, data)
                ctimer_stop!(c_timer, 67)

                # Metropolis acceptance/rejection
                x = log_proj_ratio(proj_cnt_new, tmp_ele_proj_cnt, data)
                w = exp(2.0 * (x + log_ip_new - log_ip_old))
                if !isfinite(w)
                    w = -1.0  # Should be rejected
                end

                if w > rng_real2(rng)
                    # Accept: update inv_m using O(N²) Sherman-Morrison
                    # [68] UpdateMAllTwo
                    ctimer_start!(c_timer, 68)
                    update_m_all_two_real!(
                        mi,
                        s,
                        mj,
                        t,
                        ri_old,
                        rj_old,
                        tmp_ele_idx,
                        qp_start,
                        qp_end,
                        data,
                        state,
                    )
                    ctimer_stop!(c_timer, 68)
                    tmp_ele_proj_cnt .= proj_cnt_new
                    log_ip_old = log_ip_new
                    n_accept += 1
                else
                    # Reject: just revert electron configuration (no matrix update needed!)
                    revert_ele_config!(
                        mj,
                        rj,
                        ri,
                        t,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        n_site,
                        n_elec,
                    )
                    revert_ele_config!(
                        mi,
                        ri,
                        rj,
                        s,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        n_site,
                        n_elec,
                    )
                    # No need to recalculate M matrix - inv_m was never modified!
                end
                ctimer_stop!(c_timer, 33)
            end

            # [34] recal PfM and InvM: periodic full recompute + logIP update
            if n_accept > n_site
                ctimer_start!(c_timer, 34)
                result = calculate_m_all_real!(tmp_ele_idx, qp_start, qp_end, data, state)
                if result == 0
                    log_ip_old = calculate_log_ip_real(
                        state.slater_matrix.pf_m_real,
                        qp_start,
                        qp_end,
                        data,
                    )
                end
                n_accept = 0
                ctimer_stop!(c_timer, 34)
            end
        end

        # Save sample
        # C: if (outStep >= nOutStep - NVMCSample) { sample = outStep - (nOutStep - NVMCSample); }
        if out_step >= n_out_step - n_vmc_sample
            sample = out_step - (n_out_step - n_vmc_sample)
            if sample >= 0 && sample < n_vmc_sample
                # Validate tmp_ele_idx before saving
                if all(x -> x == 0, tmp_ele_idx)
                    @error "tmp_ele_idx is all zeros at sample=$sample"
                    continue
                end

                # [35] save electron config (after the validation `continue`)
                ctimer_start!(c_timer, 35)
                # Save electron configuration
                offset_idx = sample * n_size
                offset_cfg = sample * n_site2
                offset_num = sample * n_site2
                offset_proj = sample * n_proj

                # Debug: Check ele_idx size before saving
                ele_idx_size = length(state.electron_config.ele_idx)
                expected_size = n_vmc_sample * n_size
                if ele_idx_size != expected_size
                    @error "vmc_make_sample_real!: ele_idx size mismatch at sample=$sample: expected=$expected_size, actual=$ele_idx_size. offset_idx=$offset_idx, max_index=$(offset_idx + n_size)"
                end

                for i = 1:n_size
                    if offset_idx + i <= length(state.electron_config.ele_idx)
                        state.electron_config.ele_idx[offset_idx+i] = tmp_ele_idx[i]
                    else
                        @error "vmc_make_sample_real!: BoundsError: offset_idx + i = $(offset_idx + i) > ele_idx_size = $ele_idx_size at sample=$sample, i=$i"
                    end
                end
                for i = 1:n_site2
                    if offset_cfg + i <= length(state.electron_config.ele_cfg)
                        state.electron_config.ele_cfg[offset_cfg+i] = tmp_ele_cfg[i]
                    end
                    if offset_num + i <= length(state.electron_config.ele_num)
                        state.electron_config.ele_num[offset_num+i] = tmp_ele_num[i]
                    end
                end
                for i = 1:n_proj
                    if offset_proj + i <= length(state.electron_config.ele_proj_cnt)
                        state.electron_config.ele_proj_cnt[offset_proj+i] =
                            tmp_ele_proj_cnt[i]
                    end
                end

                n_saved_samples += 1
                ctimer_stop!(c_timer, 35)
            end
        end
    end

    if n_saved_samples == 0
        @error "vmc_make_sample_real!: no samples saved"
    end

    # Copy to burn sample (for next step)
    copy_to_burn_sample!(tmp_ele_idx, tmp_ele_cfg, tmp_ele_num, tmp_ele_proj_cnt, state)
    # Set burn_flag for next step
    if length(state.electron_config.counter) >= 11
        state.electron_config.counter[11] = 1
    else
        # Extend counter array if needed
        resize!(
            state.electron_config.counter,
            max(11, length(state.electron_config.counter)),
        )
        state.electron_config.counter[11] = 1
    end
end

"""
    vmc_make_sample_fsz_real!(data::ExpertModeData, state::VMCOptimizationState, rng::AbstractRNG)

VMC sampling (real version, fsz).
Equivalent to C's `VMCMakeSample_fsz_real()`.

# Note
This is a stub implementation. Full implementation requires:
- fsz-specific electron configuration handling
- Real-valued Pfaffian calculation
"""
function vmc_make_sample_fsz_real!(
    data::ExpertModeData,
    state::VMCOptimizationState,
    rng::AbstractRNG = Random.GLOBAL_RNG,
    c_timer::CTimer = CTIMER_DISABLED,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    n_proj = length(state.electron_config.tmp_ele_proj_cnt)
    n_vmc_warmup = data.modpara.nvmc_warmup
    n_vmc_interval = data.modpara.nvmc_interval
    n_vmc_sample = data.modpara.nvmc_sample
    n_ex_update_path = data.modpara.nex_update_path
    two_sz = data.modpara.two_sz
    i_flg_orbital_general = data.i_flg_orbital_general
    loc_spn = get_loc_spn_array(data)

    n_qp_full = length(state.slater_matrix.pf_m_real)
    qp_start = 1
    qp_end = n_qp_full + 1

    if length(state.electron_config.counter) < 11
        resize!(state.electron_config.counter, 11)
        fill!(state.electron_config.counter, 0)
    end
    burn_flag = state.electron_config.counter[11] != 0

    tmp_ele_idx = state.electron_config.tmp_ele_idx
    tmp_ele_cfg = state.electron_config.tmp_ele_cfg
    tmp_ele_num = state.electron_config.tmp_ele_num
    tmp_ele_proj_cnt = state.electron_config.tmp_ele_proj_cnt
    tmp_ele_spn = state.electron_config.tmp_ele_spn

    ws = state.workspace
    proj_cnt_new = ws.proj_cnt_new
    pf_m_new_real = ws.pf_m_new_real
    fill!(proj_cnt_new, 0)
    fill!(pf_m_new_real, 0.0)

    # Ensure slater_elm_real is in sync
    n_copy = min(
        length(state.slater_matrix.slater_elm),
        length(state.slater_matrix.slater_elm_real),
    )
    @inbounds @simd for i = 1:n_copy
        state.slater_matrix.slater_elm_real[i] = real(state.slater_matrix.slater_elm[i])
    end

    if !burn_flag
        info = make_initial_sample_fsz_real!(
            tmp_ele_idx,
            tmp_ele_cfg,
            tmp_ele_num,
            tmp_ele_proj_cnt,
            tmp_ele_spn,
            qp_start,
            qp_end,
            data,
            state,
            rng,
        )
        if info != 0
            @error "make_initial_sample_fsz_real! failed with info=$info"
            return
        end
    else
        copy_from_burn_sample_fsz!(
            tmp_ele_idx,
            tmp_ele_cfg,
            tmp_ele_num,
            tmp_ele_proj_cnt,
            tmp_ele_spn,
            state,
        )
    end

    calculate_m_all_fsz_real!(tmp_ele_idx, tmp_ele_spn, qp_start, qp_end, data, state)
    log_ip_old = calculate_log_ip_real(state.slater_matrix.pf_m_real, qp_start, qp_end, data)

    if !isfinite(log_ip_old)
        @warn "VMCMakeSample_fsz_real: remakeSample logIpOld not finite"
        info = make_initial_sample_fsz_real!(
            tmp_ele_idx,
            tmp_ele_cfg,
            tmp_ele_num,
            tmp_ele_proj_cnt,
            tmp_ele_spn,
            qp_start,
            qp_end,
            data,
            state,
            rng,
        )
        if info != 0
            return
        end
        calculate_m_all_fsz_real!(tmp_ele_idx, tmp_ele_spn, qp_start, qp_end, data, state)
        log_ip_old = calculate_log_ip_real(
            state.slater_matrix.pf_m_real,
            qp_start,
            qp_end,
            data,
        )
        burn_flag = false
    end

    n_out_step = burn_flag ? n_vmc_sample + 1 : n_vmc_warmup + n_vmc_sample
    n_in_step = n_vmc_interval * n_site
    n_accept = 0

    if length(state.electron_config.counter) < 6
        resize!(state.electron_config.counter, 6)
    end
    fill!(state.electron_config.counter, 0)

    for out_step = 0:(n_out_step-1)
        for in_step = 0:(n_in_step-1)
            update_type = get_update_type(
                n_ex_update_path,
                i_flg_orbital_general,
                rng;
                two_sz = two_sz,
            )

            if update_type == HOPPING
                flag_hop = false
                if two_sz == -1
                    if rng_real2(rng) < 0.5
                        flag_hop = true
                        state.electron_config.counter[1] += 1
                        mi, ri, rj, s, t, reject_flag = make_candidate_hopping_fsz(
                            tmp_ele_idx,
                            tmp_ele_cfg,
                            tmp_ele_num,
                            tmp_ele_spn,
                            loc_spn,
                            n_site,
                            n_size,
                            two_sz,
                            rng,
                        )
                    else
                        state.electron_config.counter[5] += 1
                        mi, ri, rj, s, t, reject_flag =
                            make_candidate_local_spin_flip_conduction(
                                tmp_ele_idx,
                                tmp_ele_cfg,
                                tmp_ele_num,
                                tmp_ele_spn,
                                loc_spn,
                                n_site,
                                n_size,
                                rng,
                            )
                    end
                else
                    flag_hop = true
                    state.electron_config.counter[1] += 1
                    mi, ri, rj, s, t, reject_flag = make_candidate_hopping_fsz(
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        tmp_ele_spn,
                        loc_spn,
                        n_site,
                        n_size,
                        two_sz,
                        rng,
                    )
                end

                if reject_flag
                    continue
                end

                update_ele_config_fsz!(
                    mi,
                    ri,
                    rj,
                    s,
                    t,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    tmp_ele_spn,
                    n_site,
                )

                if s == t
                    update_proj_cnt!(
                        ri,
                        rj,
                        s,
                        proj_cnt_new,
                        tmp_ele_proj_cnt,
                        tmp_ele_num,
                        data,
                    )
                else
                    update_proj_cnt_fsz!(
                        ri,
                        rj,
                        s,
                        t,
                        proj_cnt_new,
                        tmp_ele_proj_cnt,
                        tmp_ele_num,
                        data,
                    )
                end

                calculate_new_pf_m2_fsz_real!(
                    mi,
                    t,
                    pf_m_new_real,
                    tmp_ele_idx,
                    tmp_ele_spn,
                    qp_start,
                    qp_end,
                    data,
                    state,
                )

                log_ip_new =
                    calculate_log_ip_real(pf_m_new_real, qp_start, qp_end, data)

                x = log_proj_ratio(proj_cnt_new, tmp_ele_proj_cnt, data)
                w = exp(2.0 * (x + (log_ip_new - log_ip_old)))
                if !isfinite(w)
                    w = -1.0
                end

                if w > rng_real2(rng)
                    update_m_all_fsz_real!(
                        mi,
                        t,
                        tmp_ele_idx,
                        tmp_ele_spn,
                        qp_start,
                        qp_end,
                        data,
                        state,
                    )
                    copyto!(tmp_ele_proj_cnt, proj_cnt_new)
                    log_ip_old = log_ip_new
                    n_accept += 1
                    if flag_hop
                        state.electron_config.counter[2] += 1
                    else
                        state.electron_config.counter[6] += 1
                    end
                else
                    revert_ele_config_fsz!(
                        mi,
                        ri,
                        rj,
                        s,
                        t,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        tmp_ele_spn,
                        n_site,
                    )
                end

            elseif update_type == EXCHANGE
                state.electron_config.counter[3] += 1
                mi, ri, rj, s, reject_flag = make_candidate_exchange_fsz(
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    tmp_ele_spn,
                    loc_spn,
                    n_site,
                    n_size,
                    rng,
                )
                if reject_flag
                    continue
                end

                t = 1 - s
                mj = tmp_ele_cfg[rj+1+t*n_site]

                update_ele_config_fsz!(
                    mi,
                    ri,
                    rj,
                    s,
                    s,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    tmp_ele_spn,
                    n_site,
                )
                update_proj_cnt!(ri, rj, s, proj_cnt_new, tmp_ele_proj_cnt, tmp_ele_num, data)

                update_ele_config_fsz!(
                    mj,
                    rj,
                    ri,
                    t,
                    t,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    tmp_ele_spn,
                    n_site,
                )
                update_proj_cnt!(rj, ri, t, proj_cnt_new, proj_cnt_new, tmp_ele_num, data)

                calculate_new_pf_m_two2_fsz_real!(
                    mi,
                    s,
                    mj,
                    t,
                    pf_m_new_real,
                    tmp_ele_idx,
                    tmp_ele_spn,
                    qp_start,
                    qp_end,
                    data,
                    state,
                )

                log_ip_new =
                    calculate_log_ip_real(pf_m_new_real, qp_start, qp_end, data)

                x = log_proj_ratio(proj_cnt_new, tmp_ele_proj_cnt, data)
                w = exp(2.0 * (x + (log_ip_new - log_ip_old)))
                if !isfinite(w)
                    w = -1.0
                end

                if w > rng_real2(rng)
                    update_m_all_two_fsz_real!(
                        mi,
                        s,
                        mj,
                        t,
                        ri,
                        rj,
                        tmp_ele_idx,
                        tmp_ele_spn,
                        qp_start,
                        qp_end,
                        data,
                        state,
                    )
                    copyto!(tmp_ele_proj_cnt, proj_cnt_new)
                    log_ip_old = log_ip_new
                    n_accept += 1
                    state.electron_config.counter[4] += 1
                else
                    revert_ele_config_fsz!(
                        mj,
                        rj,
                        ri,
                        t,
                        t,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        tmp_ele_spn,
                        n_site,
                    )
                    revert_ele_config_fsz!(
                        mi,
                        ri,
                        rj,
                        s,
                        s,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        tmp_ele_spn,
                        n_site,
                    )
                end

            elseif update_type == LOCALSPINFLIP
                state.electron_config.counter[5] += 1
                mi, ri, rj, s, t, reject_flag = make_candidate_local_spin_flip_localspin(
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    tmp_ele_spn,
                    loc_spn,
                    n_site,
                    n_size,
                    rng,
                )

                if reject_flag
                    continue
                end

                update_ele_config_fsz!(
                    mi,
                    ri,
                    rj,
                    s,
                    t,
                    tmp_ele_idx,
                    tmp_ele_cfg,
                    tmp_ele_num,
                    tmp_ele_spn,
                    n_site,
                )
                update_proj_cnt_fsz!(
                    ri,
                    rj,
                    s,
                    t,
                    proj_cnt_new,
                    tmp_ele_proj_cnt,
                    tmp_ele_num,
                    data,
                )

                calculate_new_pf_m2_fsz_real!(
                    mi,
                    t,
                    pf_m_new_real,
                    tmp_ele_idx,
                    tmp_ele_spn,
                    qp_start,
                    qp_end,
                    data,
                    state,
                )

                log_ip_new =
                    calculate_log_ip_real(pf_m_new_real, qp_start, qp_end, data)

                x = log_proj_ratio(proj_cnt_new, tmp_ele_proj_cnt, data)
                w = exp(2.0 * (x + (log_ip_new - log_ip_old)))
                if !isfinite(w)
                    w = -1.0
                end

                if w > rng_real2(rng)
                    update_m_all_fsz_real!(
                        mi,
                        t,
                        tmp_ele_idx,
                        tmp_ele_spn,
                        qp_start,
                        qp_end,
                        data,
                        state,
                    )
                    copyto!(tmp_ele_proj_cnt, proj_cnt_new)
                    log_ip_old = log_ip_new
                    n_accept += 1
                    state.electron_config.counter[6] += 1
                else
                    revert_ele_config_fsz!(
                        mi,
                        ri,
                        rj,
                        s,
                        t,
                        tmp_ele_idx,
                        tmp_ele_cfg,
                        tmp_ele_num,
                        tmp_ele_spn,
                        n_site,
                    )
                end
            end

            if n_accept > n_site
                calculate_m_all_fsz_real!(
                    tmp_ele_idx,
                    tmp_ele_spn,
                    qp_start,
                    qp_end,
                    data,
                    state,
                )
                log_ip_old = calculate_log_ip_real(
                    state.slater_matrix.pf_m_real,
                    qp_start,
                    qp_end,
                    data,
                )
                n_accept = 0
            end
        end

        if out_step >= n_out_step - n_vmc_sample
            sample = out_step - (n_out_step - n_vmc_sample)
            save_ele_config_fsz!(
                sample,
                log_ip_old,
                tmp_ele_idx,
                tmp_ele_cfg,
                tmp_ele_num,
                tmp_ele_proj_cnt,
                tmp_ele_spn,
                data,
                state,
            )
        end
    end

    copy_to_burn_sample_fsz!(
        tmp_ele_idx,
        tmp_ele_cfg,
        tmp_ele_num,
        tmp_ele_proj_cnt,
        tmp_ele_spn,
        state,
    )

    state.electron_config.counter[11] = 1
end

"""
    vmc_bf_make_sample!(data::ExpertModeData, state::VMCOptimizationState, rng::AbstractRNG)

VMC sampling with Back Flow (complex version).
Equivalent to C's `VMC_BF_MakeSample()`.

BackFlow is not supported in Julia-mVMC v0.1: calling this function raises
an error. Inputs that drive `n_proj_bf > 0` should not be passed to v0.1.
"""
function vmc_bf_make_sample!(
    ::ExpertModeData,
    ::VMCOptimizationState,
    ::AbstractRNG = Random.GLOBAL_RNG,
)
    error("BackFlow is not supported in Julia-mVMC v0.1. Remove BackFlow keywords from namelist.def, " *
          "or fall back to the C reference at https://github.com/issp-center-dev/mVMC.")
end

"""
    vmc_bf_make_sample_real!(data::ExpertModeData, state::VMCOptimizationState, rng::AbstractRNG)

VMC sampling with Back Flow (real version).
Equivalent to C's `VMC_BF_MakeSample_real()`.

BackFlow is not supported in Julia-mVMC v0.1: calling this function raises
an error. Inputs that drive `n_proj_bf > 0` should not be passed to v0.1.
"""
function vmc_bf_make_sample_real!(
    ::ExpertModeData,
    ::VMCOptimizationState,
    ::AbstractRNG = Random.GLOBAL_RNG,
)
    error("BackFlow is not supported in Julia-mVMC v0.1. Remove BackFlow keywords from namelist.def, " *
          "or fall back to the C reference at https://github.com/issp-center-dev/mVMC.")
end
