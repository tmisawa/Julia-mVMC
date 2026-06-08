# VMC Main Calculation Functions
# Calculate energy expectation values and SR optimization quantities.

# ============================================================================
# Physical Quantity Initialization
# ============================================================================

# Import green function calculation
include("green_func_calc.jl")

"""
    clear_phys_quantity!(state::VMCOptimizationState)

Clear all physical quantities for a new calculation.
Equivalent to C's `clearPhysQuantity()`.
"""

function clear_phys_quantity!(state::VMCOptimizationState)
    # Clear energy data
    state.energy.wc = 0.0 + 0.0im
    state.energy.etot = 0.0 + 0.0im
    state.energy.etot2 = 0.0 + 0.0im
    state.energy.sztot = 0.0 + 0.0im
    state.energy.sztot2 = 0.0 + 0.0im

    # Clear SR optimization data
    fill!(state.sr_opt.sr_opt_oo, 0.0 + 0.0im)
    fill!(state.sr_opt.sr_opt_ho, 0.0 + 0.0im)
    fill!(state.sr_opt.sr_opt_o, 0.0 + 0.0im)

    if !isempty(state.sr_opt.sr_opt_oo_real)
        fill!(state.sr_opt.sr_opt_oo_real, 0.0)
        fill!(state.sr_opt.sr_opt_ho_real, 0.0)
        fill!(state.sr_opt.sr_opt_o_real, 0.0)
    end
end

# ============================================================================
# Projection Ratio
# ============================================================================

"""
    proj_ratio(proj_cnt_new::Vector{Int}, proj_cnt_old::Vector{Int}, data::ExpertModeData) -> Float64

Calculate projection ratio exp(z) where z is the log projection ratio.
Equivalent to C's `ProjRatio()`.
"""
function proj_ratio(
    proj_cnt_new::Vector{Int},
    proj_cnt_old::Vector{Int},
    data::ExpertModeData,
)::Float64
    return exp(log_proj_ratio(proj_cnt_new, proj_cnt_old, data))
end

"""
    proj_rbm_ratio(... ) -> ComplexF64

Combined projection/RBM ratio for local Green-function updates.
Equivalent to C's ProjRatio * RBMRatio.
"""
function proj_rbm_ratio(
    proj_cnt_new::Vector{Int},
    proj_cnt_old::Vector{Int},
    ele_num_new::Vector{Int},
    ele_num_old::Vector{Int},
    data::ExpertModeData,
)::ComplexF64
    z = ComplexF64(proj_ratio(proj_cnt_new, proj_cnt_old, data), 0.0)
    if has_rbm_terms(data)
        rbm_cnt_new = make_rbm_cnt(ele_num_new, data)
        rbm_cnt_old = make_rbm_cnt(ele_num_old, data)
        z *= exp(log_rbm_ratio(rbm_cnt_new, rbm_cnt_old, data))
    end
    return z
end

# ============================================================================
# Green Functions
# ============================================================================

"""
    green_func1(ri::Int, rj::Int, s::Int, ip::ComplexF64,
               ele_idx::Vector{Int}, ele_cfg::Vector{Int}, ele_num::Vector{Int},
               ele_proj_cnt::Vector{Int}, data::ExpertModeData,
               state::VMCOptimizationState) -> ComplexF64

Calculate 1-body Green function <ψ|c†_{ri,s} c_{rj,s}|x> / <ψ|x>.
Equivalent to C's `GreenFunc1()`.

# Arguments
- `ri`, `rj`: Site indices (0-based, as in C)
- `s`: Spin (0 = up, 1 = down)
- `ip`: Inner product <ψ|x>
- Other arguments: Configuration data
"""
function green_func1(
    ri::Int,
    rj::Int,
    s::Int,
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState;
    all_complex::Bool = true,
)::ComplexF64
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_qp_full = get_n_qp_full(data)
    n_proj = length(ele_proj_cnt)

    # Simple cases
    if ri == rj
        return ComplexF64(ele_num[ri+1+s*n_site])
    end
    if ele_num[ri+1+s*n_site] == 1 || ele_num[rj+1+s*n_site] == 0
        return 0.0 + 0.0im
    end

    # Get electron index at site rj with spin s
    mj = ele_cfg[rj+1+s*n_site]  # 0-based in C
    mj < 0 && return 0.0 + 0.0im
    msj = mj + s * n_elec  # 0-based index in ele_idx
    rsi = ri + s * n_site
    rsj = rj + s * n_site

    # Allocate temporary arrays
    my_ele_idx = copy(ele_idx)
    my_ele_num = copy(ele_num)
    proj_cnt_new = zeros(Int, n_proj)

    # Hopping: move electron from rj to ri
    my_ele_idx[msj+1] = ri  # Update position (1-based Julia indexing)
    my_ele_num[rsj+1] = 0
    my_ele_num[rsi+1] = 1

    # Update projection count
    update_proj_cnt!(rj, ri, s, proj_cnt_new, ele_proj_cnt, my_ele_num, data)

    # Calculate projection/RBM ratio
    z = proj_rbm_ratio(proj_cnt_new, ele_proj_cnt, my_ele_num, ele_num, data)

    # Calculate new Pfaffian and inner product (use real version if !all_complex)
    if !all_complex
        pf_m_new_real = zeros(Float64, n_qp_full)
        calculate_new_pf_m2_real!(
            mj,
            s,
            pf_m_new_real,
            my_ele_idx,
            1,
            n_qp_full + 1,
            data,
            state,
        )
        ip_new_real = calculate_ip_real(pf_m_new_real, 1, n_qp_full + 1, data)
        ip_new = ComplexF64(ip_new_real, 0.0)
    else
        pf_m_new = zeros(ComplexF64, n_qp_full)
        calculate_new_pf_m2!(mj, s, pf_m_new, my_ele_idx, 1, n_qp_full + 1, data, state)
        ip_new = calculate_ip_fcmp(pf_m_new, 1, n_qp_full + 1, data)
    end
    z *= ip_new

    result = conj(z / ip)


    return result
end

"""
    green_func1_fsz(ri::Int, rj::Int, s::Int, ip::ComplexF64,
                    ele_idx::Vector{Int}, ele_cfg::Vector{Int}, ele_num::Vector{Int},
                    ele_proj_cnt::Vector{Int}, ele_spn::Vector{Int},
                    data::ExpertModeData, state::VMCOptimizationState;
                    all_complex::Bool = true) -> ComplexF64

Calculate 1-body Green function for FSZ mode: <ψ|c†_{ri,s} c_{rj,s}|x> / <ψ|x>.

In FSZ mode:
- Electron spins are tracked individually via `ele_spn` array
- `mj = ele_cfg[rj + s*n_site]` (electron index at site rj with spin s)
- After hopping: `ele_idx[mj] = ri`, `ele_spn[mj] = s`

# Arguments
- `ri`, `rj`: Site indices (0-based)
- `s`: Spin (0 = up, 1 = down)
- `ip`: Inner product <ψ|x>
- `ele_spn`: Electron spin array (1-based, values 0 or 1)
- Other arguments: Configuration data

# Reference
- C implementation: mVMC/src/mVMC/locgrn_fsz.c:50-84 (GreenFunc1_fsz)
"""
function green_func1_fsz(
    ri::Int,
    rj::Int,
    s::Int,
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    ele_spn::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState;
    all_complex::Bool = true,
)::ComplexF64
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_qp_full = get_n_qp_full(data)
    n_proj = length(ele_proj_cnt)

    # Simple cases
    if ri == rj
        return ComplexF64(ele_num[ri+1+s*n_site])
    end
    if ele_num[ri+1+s*n_site] == 1 || ele_num[rj+1+s*n_site] == 0
        return 0.0 + 0.0im
    end

    # FSZ: mj = ele_cfg[rj + s*n_site] (electron index at site rj with spin s)
    mj = ele_cfg[rj+1+s*n_site]  # 0-based electron index
    # FSZ: msj = mj (no spin offset in FSZ)
    msj = mj
    rsi = ri + s * n_site
    rsj = rj + s * n_site

    # Allocate temporary arrays
    my_ele_idx = copy(ele_idx)
    my_ele_spn = copy(ele_spn)  # FSZ: copy ele_spn too
    my_ele_num = copy(ele_num)
    proj_cnt_new = zeros(Int, n_proj)

    # Hopping: move electron from rj to ri with spin s
    # FSZ: update both ele_idx and ele_spn
    my_ele_idx[msj+1] = ri  # New site (0-based)
    my_ele_spn[msj+1] = s   # FSZ: update spin
    my_ele_num[rsj+1] = 0
    my_ele_num[rsi+1] = 1

    # Update projection count
    update_proj_cnt!(rj, ri, s, proj_cnt_new, ele_proj_cnt, my_ele_num, data)

    # Calculate projection/RBM ratio
    z = proj_rbm_ratio(proj_cnt_new, ele_proj_cnt, my_ele_num, ele_num, data)

    # Calculate new Pfaffian using FSZ version
    pf_m_new = zeros(ComplexF64, n_qp_full)
    calculate_new_pf_m2_fsz!(
        mj,
        s,
        pf_m_new,
        my_ele_idx,
        my_ele_spn,
        1,
        n_qp_full + 1,
        data,
        state,
    )
    ip_new = calculate_ip_fcmp(pf_m_new, 1, n_qp_full + 1, data)
    z *= ip_new

    result = conj(z / ip)

    return result
end

"""
    green_func1_fsz2(ri::Int, rj::Int, s::Int, t::Int, ip::ComplexF64,
                     ele_idx::Vector{Int}, ele_cfg::Vector{Int}, ele_num::Vector{Int},
                     ele_proj_cnt::Vector{Int}, ele_spn::Vector{Int},
                     data::ExpertModeData, state::VMCOptimizationState;
                     all_complex::Bool = true) -> ComplexF64

Calculate 1-body Green function for FSZ mode with different spins: <ψ|c†_{ri,s} c_{rj,t}|x> / <ψ|x>.

This function handles the case where the creation and annihilation operators have different spins (s != t).

# Arguments
- `ri`: Site index for creation operator (0-based)
- `rj`: Site index for annihilation operator (0-based)
- `s`: Spin of creation operator (0 = up, 1 = down)
- `t`: Spin of annihilation operator (0 = up, 1 = down)
- `ip`: Inner product <ψ|x>
- Other arguments: Configuration data

# Reference
- C implementation: mVMC/src/mVMC/locgrn_fsz.c:86-120 (GreenFunc1_fsz2)
"""
function green_func1_fsz2(
    ri::Int,
    rj::Int,
    s::Int,
    t::Int,
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    ele_spn::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState;
    all_complex::Bool = true,
)::ComplexF64
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_qp_full = get_n_qp_full(data)
    n_proj = length(ele_proj_cnt)

    # For s == t, fall back to green_func1_fsz
    if s == t
        return green_func1_fsz(
            ri,
            rj,
            s,
            ip,
            ele_idx,
            ele_cfg,
            ele_num,
            ele_proj_cnt,
            ele_spn,
            data,
            state;
            all_complex = all_complex,
        )
    end

    rsi = ri + s * n_site  # Target site with spin s
    rtj = rj + t * n_site  # Source site with spin t

    # Check occupation: need electron at (rj, t) and vacancy at (ri, s)
    if ele_num[rsi+1] == 1 || ele_num[rtj+1] == 0
        return 0.0 + 0.0im
    end

    # Get electron index at (rj, t)
    mj = ele_cfg[rtj+1]

    # Allocate temporary arrays
    my_ele_idx = copy(ele_idx)
    my_ele_spn = copy(ele_spn)
    my_ele_num = copy(ele_num)
    proj_cnt_new = zeros(Int, n_proj)

    # Hopping: move electron from (rj, t) to (ri, s) - spin flip!
    my_ele_idx[mj+1] = ri
    my_ele_spn[mj+1] = s
    my_ele_num[rtj+1] = 0
    my_ele_num[rsi+1] = 1

    # Update projection count for spin-flip hopping
    update_proj_cnt_fsz!(rj, ri, t, s, proj_cnt_new, ele_proj_cnt, my_ele_num, data)

    # Calculate projection/RBM ratio
    z = proj_rbm_ratio(proj_cnt_new, ele_proj_cnt, my_ele_num, ele_num, data)

    # Calculate new Pfaffian using FSZ version (note: ele_spn[mj] = s after hopping)
    pf_m_new = zeros(ComplexF64, n_qp_full)
    calculate_new_pf_m2_fsz!(
        mj,
        s,
        pf_m_new,
        my_ele_idx,
        my_ele_spn,
        1,
        n_qp_full + 1,
        data,
        state,
    )
    ip_new = calculate_ip_fcmp(pf_m_new, 1, n_qp_full + 1, data)
    z *= ip_new

    result = conj(z / ip)

    return result
end

"""
    green_func2_fsz(ri::Int, rj::Int, rk::Int, rl::Int, s::Int, t::Int, ip::ComplexF64,
                    ele_idx::Vector{Int}, ele_cfg::Vector{Int}, ele_num::Vector{Int},
                    ele_proj_cnt::Vector{Int}, ele_spn::Vector{Int},
                    data::ExpertModeData, state::VMCOptimizationState;
                    all_complex::Bool = true) -> ComplexF64

Calculate 2-body Green function for FSZ mode: <ψ|c†_{ri,s} c_{rj,s} c†_{rk,t} c_{rl,t}|x> / <ψ|x>.

This function handles many special cases (when indices coincide) and the general case
where two electrons simultaneously hop.

# Arguments
- `ri`, `rj`, `rk`, `rl`: Site indices (0-based)
- `s`, `t`: Spins (0 = up, 1 = down)
- `ip`: Inner product <ψ|x>
- `ele_spn`: Electron spin array (1-based, values 0 or 1)

# Reference
- C implementation: mVMC/src/mVMC/locgrn_fsz.c:128-216 (GreenFunc2_fsz)
"""
function green_func2_fsz(
    ri::Int,
    rj::Int,
    rk::Int,
    rl::Int,
    s::Int,
    t::Int,
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    ele_spn::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState;
    all_complex::Bool = true,
)::ComplexF64
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_qp_full = get_n_qp_full(data)
    n_proj = length(ele_proj_cnt)

    rsi = ri + s * n_site
    rsj = rj + s * n_site
    rtk = rk + t * n_site
    rtl = rl + t * n_site

    # Handle special cases when s == t
    if s == t
        if rk == rl  # CisAjsNks
            if ele_num[rtk+1] == 0
                return 0.0 + 0.0im
            else
                return green_func1_fsz(
                    ri,
                    rj,
                    s,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    ele_spn,
                    data,
                    state;
                    all_complex = all_complex,
                )
            end
        elseif rj == rl
            return 0.0 + 0.0im  # CisAjsCksAjs (j!=k)
        elseif ri == rl  # AjsCksNis
            if ele_num[rsi+1] == 0
                return 0.0 + 0.0im
            elseif rj == rk
                return ComplexF64(1.0 - ele_num[rsj+1])
            else
                return -green_func1_fsz(
                    rk,
                    rj,
                    s,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    ele_spn,
                    data,
                    state;
                    all_complex = all_complex,
                )
            end
        elseif rj == rk  # CisAls(1-Njs)
            if ele_num[rsj+1] == 1
                return 0.0 + 0.0im
            elseif ri == rl
                return ComplexF64(ele_num[rsi+1])
            else
                return green_func1_fsz(
                    ri,
                    rl,
                    s,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    ele_spn,
                    data,
                    state;
                    all_complex = all_complex,
                )
            end
        elseif ri == rk
            return 0.0 + 0.0im  # CisAjsCisAls (i!=j)
        elseif ri == rj  # NisCksAls (i!=k,l)
            if ele_num[rsi+1] == 0
                return 0.0 + 0.0im
            else
                return green_func1_fsz(
                    rk,
                    rl,
                    s,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    ele_spn,
                    data,
                    state;
                    all_complex = all_complex,
                )
            end
        end
    else  # s != t
        if rk == rl  # CisAjsNkt
            if ele_num[rtk+1] == 0
                return 0.0 + 0.0im
            elseif ri == rj
                return ComplexF64(ele_num[rsi+1])
            else
                return green_func1_fsz(
                    ri,
                    rj,
                    s,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    ele_spn,
                    data,
                    state;
                    all_complex = all_complex,
                )
            end
        elseif ri == rj  # NisCktAlt
            if ele_num[rsi+1] == 0
                return 0.0 + 0.0im
            else
                return green_func1_fsz(
                    rk,
                    rl,
                    t,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    ele_spn,
                    data,
                    state;
                    all_complex = all_complex,
                )
            end
        end
    end

    # General case: check occupation conditions
    if ele_num[rsi+1] == 1 ||
       ele_num[rsj+1] == 0 ||
       ele_num[rtk+1] == 1 ||
       ele_num[rtl+1] == 0
        return 0.0 + 0.0im
    end

    # Get electron indices
    mj = ele_cfg[rj+1+s*n_site]  # Electron at (rj, s)
    ml = ele_cfg[rl+1+t*n_site]  # Electron at (rl, t)
    # FSZ: msj = mj, mtl = ml (no spin offset)
    msj = mj
    mtl = ml

    # Allocate temporary arrays
    my_ele_idx = copy(ele_idx)
    my_ele_spn = copy(ele_spn)
    my_ele_num = copy(ele_num)
    proj_cnt_new = zeros(Int, n_proj)

    # First hopping: electron ml from (rl, t) to (rk, t)
    my_ele_idx[mtl+1] = rk
    my_ele_spn[mtl+1] = t
    my_ele_num[rtl+1] = 0
    my_ele_num[rtk+1] = 1
    update_proj_cnt!(rl, rk, t, proj_cnt_new, ele_proj_cnt, my_ele_num, data)

    # Second hopping: electron mj from (rj, s) to (ri, s)
    my_ele_idx[msj+1] = ri
    my_ele_spn[msj+1] = s
    my_ele_num[rsj+1] = 0
    my_ele_num[rsi+1] = 1
    proj_cnt_new2 = zeros(Int, n_proj)
    update_proj_cnt!(rj, ri, s, proj_cnt_new2, proj_cnt_new, my_ele_num, data)

    # Calculate projection ratio
    z = proj_rbm_ratio(proj_cnt_new2, ele_proj_cnt, my_ele_num, ele_num, data)

    # Calculate new Pfaffian for two-electron hopping
    pf_m_new = zeros(ComplexF64, n_qp_full)
    calculate_new_pf_m_two_fsz!(
        ml,
        t,
        mj,
        s,
        pf_m_new,
        my_ele_idx,
        my_ele_spn,
        1,
        n_qp_full + 1,
        data,
        state,
    )
    ip_new = calculate_ip_fcmp(pf_m_new, 1, n_qp_full + 1, data)
    z *= ip_new

    result = conj(z / ip)

    return result
end

"""
    green_func2_fsz2(ri::Int, rj::Int, rk::Int, rl::Int, s::Int, t::Int, u::Int, v::Int, ip::ComplexF64,
                     ele_idx::Vector{Int}, ele_cfg::Vector{Int}, ele_num::Vector{Int},
                     ele_proj_cnt::Vector{Int}, ele_spn::Vector{Int},
                     data::ExpertModeData, state::VMCOptimizationState;
                     all_complex::Bool = true) -> ComplexF64

Calculate 2-body Green function for FSZ mode with all 4 different spins:
<ψ|c†_{ri,s} c_{rj,t} c†_{rk,u} c_{rl,v}|x> / <ψ|x>.

This function handles the most general case where all 4 operators can have different spins.
Used for Sz-non-conserving terms in InterAll.

# Arguments
- `ri`, `rj`, `rk`, `rl`: Site indices (0-based)
- `s`, `t`, `u`, `v`: Spins for each operator (0 = up, 1 = down)
- `ip`: Inner product <ψ|x>

# Reference
- C implementation: mVMC/src/mVMC/locgrn_fsz.c:220-342 (GreenFunc2_fsz2)
"""
function green_func2_fsz2(
    ri::Int,
    rj::Int,
    rk::Int,
    rl::Int,
    s::Int,
    t::Int,
    u::Int,
    v::Int,
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    ele_spn::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState;
    all_complex::Bool = true,
)::ComplexF64
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_qp_full = get_n_qp_full(data)
    n_proj = length(ele_proj_cnt)

    # Combined site-spin indices (FSZ)
    XI = ri + s * n_site
    XJ = rj + t * n_site
    XK = rk + u * n_site
    XL = rl + v * n_site

    # Handle all special cases (following C implementation exactly)

    if XI == XJ
        if XJ == XK
            if XK == XL
                # Case #1: I=J=K=L
                return ComplexF64(ele_num[XI+1])
            else
                # Case #2: I=J=K != L
                if ele_num[XI+1] == 1
                    return 0.0 + 0.0im
                else
                    return green_func1_fsz2(
                        rk,
                        rl,
                        u,
                        v,
                        ip,
                        ele_idx,
                        ele_cfg,
                        ele_num,
                        ele_proj_cnt,
                        ele_spn,
                        data,
                        state;
                        all_complex = all_complex,
                    )
                end
            end
        elseif XJ == XL
            # Case #3: I=J=L != K
            return 0.0 + 0.0im
        elseif XK == XL
            # Case #4: I=J != K=L
            return ComplexF64(ele_num[XI+1] * ele_num[XK+1])
        else
            # Case #5: I=J, K!=L
            if ele_num[XI+1] == 0
                return 0.0 + 0.0im
            else
                return green_func1_fsz2(
                    rk,
                    rl,
                    u,
                    v,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    ele_spn,
                    data,
                    state;
                    all_complex = all_complex,
                )
            end
        end
    elseif XI == XK
        if XJ == XL
            # Case #6: I=K != J=L
            return 0.0 + 0.0im
        elseif XK == XL
            # Case #7: I=K=L != J
            return 0.0 + 0.0im
        else
            # Case #8: I=K != L != J
            return 0.0 + 0.0im
        end
    elseif XI == XL
        if XJ == XK
            # Case #9: I=L != J=K
            return ComplexF64(ele_num[XI+1] * (1 - ele_num[XJ+1]))
        else
            # Case #10: I=L != J != K
            if ele_num[XI+1] == 0
                return 0.0 + 0.0im
            else
                return -green_func1_fsz2(
                    rk,
                    rj,
                    u,
                    t,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    ele_spn,
                    data,
                    state;
                    all_complex = all_complex,
                )
            end
        end
    elseif XJ == XK
        if XK == XL
            # Case #11: I != J=K=L
            if ele_num[XJ+1] == 0
                return 0.0 + 0.0im
            else
                return green_func1_fsz2(
                    ri,
                    rj,
                    s,
                    t,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    ele_spn,
                    data,
                    state;
                    all_complex = all_complex,
                )
            end
        else
            # Case #12: I != J=K != L
            if ele_num[XJ+1] == 1
                return 0.0 + 0.0im
            else
                return green_func1_fsz2(
                    ri,
                    rl,
                    s,
                    v,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    ele_spn,
                    data,
                    state;
                    all_complex = all_complex,
                )
            end
        end
    elseif XJ == XL
        # Case #13: I != J=L != K
        return 0.0 + 0.0im
    elseif XK == XL
        # Case #14: I != J != K=L
        if ele_num[XK+1] == 0
            return 0.0 + 0.0im
        else
            return green_func1_fsz2(
                ri,
                rj,
                s,
                t,
                ip,
                ele_idx,
                ele_cfg,
                ele_num,
                ele_proj_cnt,
                ele_spn,
                data,
                state;
                all_complex = all_complex,
            )
        end
    end

    # General case: no pairs exist
    # Check occupation: need electron at XJ and XL, vacancy at XI and XK
    if ele_num[XI+1] == 1 || ele_num[XJ+1] == 0 || ele_num[XK+1] == 1 || ele_num[XL+1] == 0
        return 0.0 + 0.0im
    end

    # Get electron indices
    mj = ele_cfg[XJ+1]  # Electron at (rj, t)
    ml = ele_cfg[XL+1]  # Electron at (rl, v)

    # Allocate temporary arrays
    my_ele_idx = copy(ele_idx)
    my_ele_spn = copy(ele_spn)
    my_ele_num = copy(ele_num)
    proj_cnt_new = zeros(Int, n_proj)
    proj_cnt_new2 = zeros(Int, n_proj)

    # First hopping: electron ml from (rl, v) to (rk, u) - may include spin flip
    my_ele_idx[ml+1] = rk
    my_ele_spn[ml+1] = u
    my_ele_num[XL+1] = 0
    my_ele_num[XK+1] = 1
    update_proj_cnt_fsz!(rl, rk, v, u, proj_cnt_new, ele_proj_cnt, my_ele_num, data)

    # Second hopping: electron mj from (rj, t) to (ri, s) - may include spin flip
    my_ele_idx[mj+1] = ri
    my_ele_spn[mj+1] = s
    my_ele_num[XJ+1] = 0
    my_ele_num[XI+1] = 1
    update_proj_cnt_fsz!(rj, ri, t, s, proj_cnt_new2, proj_cnt_new, my_ele_num, data)

    # Calculate projection ratio
    z = proj_rbm_ratio(proj_cnt_new2, ele_proj_cnt, my_ele_num, ele_num, data)

    # Calculate new Pfaffian for two-electron hopping with spin flips
    pf_m_new = zeros(ComplexF64, n_qp_full)
    calculate_new_pf_m_two_fsz!(
        ml,
        u,
        mj,
        s,
        pf_m_new,
        my_ele_idx,
        my_ele_spn,
        1,
        n_qp_full + 1,
        data,
        state,
    )
    ip_new = calculate_ip_fcmp(pf_m_new, 1, n_qp_full + 1, data)
    z *= ip_new

    result = conj(z / ip)

    return result
end

"""
    green_func2(ri::Int, rj::Int, rk::Int, rl::Int, s::Int, t::Int, ip::ComplexF64,
               ele_idx::Vector{Int}, ele_cfg::Vector{Int}, ele_num::Vector{Int},
               ele_proj_cnt::Vector{Int}, data::ExpertModeData,
               state::VMCOptimizationState) -> ComplexF64

Calculate 2-body Green function <ψ|c†_{ri,s} c_{rj,s} c†_{rk,t} c_{rl,t}|x> / <ψ|x>.
Equivalent to C's `GreenFunc2()`.

# Arguments
- `ri`, `rj`, `rk`, `rl`: Site indices (0-based)
- `s`, `t`: Spins (0 = up, 1 = down)
- `ip`: Inner product <ψ|x>
"""
function green_func2(
    ri::Int,
    rj::Int,
    rk::Int,
    rl::Int,
    s::Int,
    t::Int,
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState;
    all_complex::Bool = true,
)::ComplexF64
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_qp_full = get_n_qp_full(data)
    n_proj = length(ele_proj_cnt)

    rsi = ri + s * n_site
    rsj = rj + s * n_site
    rtk = rk + t * n_site
    rtl = rl + t * n_site

    # Handle special cases (same logic as C code)
    if s == t
        if rk == rl  # CisAjsNks
            if ele_num[rtk+1] == 0
                return 0.0 + 0.0im
            else
                return green_func1(
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
            end
        elseif rj == rl
            return 0.0 + 0.0im  # CisAjsCksAjs (j != k)
        elseif ri == rl  # AjsCksNis
            if ele_num[rsi+1] == 0
                return 0.0 + 0.0im
            elseif rj == rk
                return ComplexF64(1 - ele_num[rsj+1])
            else
                return -green_func1(
                    rk,
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
            end
        elseif rj == rk  # CisAls(1-Njs)
            if ele_num[rsj+1] == 1
                return 0.0 + 0.0im
            elseif ri == rl
                return ComplexF64(ele_num[rsi+1])
            else
                return green_func1(
                    ri,
                    rl,
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
            end
        elseif ri == rk
            return 0.0 + 0.0im  # CisAjsCisAls (i != j)
        elseif ri == rj  # NisCksAls (i != k, l)
            if ele_num[rsi+1] == 0
                return 0.0 + 0.0im
            else
                return green_func1(
                    rk,
                    rl,
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
            end
        end
    else  # s != t
        if rk == rl  # CisAjsNkt
            if ele_num[rtk+1] == 0
                return 0.0 + 0.0im
            elseif ri == rj
                return ComplexF64(ele_num[rsi+1])
            else
                return green_func1(
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
            end
        elseif ri == rj  # NisCktAlt
            if ele_num[rsi+1] == 0
                return 0.0 + 0.0im
            else
                return green_func1(
                    rk,
                    rl,
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
            end
        end
    end

    # General case: two-body term
    if ele_num[rsi+1] == 1 ||
       ele_num[rsj+1] == 0 ||
       ele_num[rtk+1] == 1 ||
       ele_num[rtl+1] == 0
        return 0.0 + 0.0im
    end

    mj = ele_cfg[rj+1+s*n_site]
    ml = ele_cfg[rl+1+t*n_site]
    (mj < 0 || ml < 0) && return 0.0 + 0.0im
    msj = mj + s * n_elec
    mtl = ml + t * n_elec

    # Allocate temporary arrays
    my_ele_idx = copy(ele_idx)
    my_ele_num = copy(ele_num)
    proj_cnt_new = zeros(Int, n_proj)

    # First hopping: move electron from rl to rk (spin t)
    my_ele_idx[mtl+1] = rk
    my_ele_num[rtl+1] = 0
    my_ele_num[rtk+1] = 1
    update_proj_cnt!(rl, rk, t, proj_cnt_new, ele_proj_cnt, my_ele_num, data)

    # Second hopping: move electron from rj to ri (spin s)
    my_ele_idx[msj+1] = ri
    my_ele_num[rsj+1] = 0
    my_ele_num[rsi+1] = 1
    update_proj_cnt!(rj, ri, s, proj_cnt_new, proj_cnt_new, my_ele_num, data)

    # Calculate projection ratio
    z = proj_rbm_ratio(proj_cnt_new, ele_proj_cnt, my_ele_num, ele_num, data)

    # Calculate new Pfaffian and inner product (use real version if !all_complex)
    if !all_complex
        pf_m_new_real = zeros(Float64, n_qp_full)
        calculate_new_pf_m_two2_real!(
            ml,
            t,
            mj,
            s,
            pf_m_new_real,
            my_ele_idx,
            1,
            n_qp_full + 1,
            data,
            state,
        )
        ip_new_real = calculate_ip_real(pf_m_new_real, 1, n_qp_full + 1, data)
        ip_new = ComplexF64(ip_new_real, 0.0)
    else
        pf_m_new = zeros(ComplexF64, n_qp_full)
        calculate_new_pf_m_two2!(
            ml,
            t,
            mj,
            s,
            pf_m_new,
            my_ele_idx,
            1,
            n_qp_full + 1,
            data,
            state,
        )
        ip_new = calculate_ip_fcmp(pf_m_new, 1, n_qp_full + 1, data)
    end
    z *= ip_new

    result = conj(z / ip)


    return result
end


# ============================================================================
# Hamiltonian Calculation
# ============================================================================

"""
    calculate_hamiltonian(ip::ComplexF64, ele_idx::Vector{Int}, ele_cfg::Vector{Int},
                         ele_num::Vector{Int}, ele_proj_cnt::Vector{Int},
                         data::ExpertModeData, state::VMCOptimizationState) -> ComplexF64

Calculate local energy (Hamiltonian expectation value).
Equivalent to C's `CalculateHamiltonian()`.
"""
function calculate_hamiltonian(
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState;
    all_complex::Bool = true,
    c_timer::CTimer = CTIMER_DISABLED,
)::ComplexF64
    n_site = data.modpara.nsite

    # Get up-spin and down-spin electron numbers (0-based indexing in C)
    n0 = @view ele_num[1:n_site]  # up-spin
    n1 = @view ele_num[(n_site+1):(2*n_site)]  # down-spin

    e = 0.0 + 0.0im

    # [70] CalHamiltonian0: diagonal terms (CoulombIntra/CoulombInter/Hund)
    ctimer_start!(c_timer, 70)

    # Debug: Track energy contributions
    e_coulomb_intra = 0.0 + 0.0im
    e_coulomb_inter = 0.0 + 0.0im
    e_hund = 0.0 + 0.0im
    e_transfer = 0.0 + 0.0im
    e_pair_hop = 0.0 + 0.0im
    e_exchange = 0.0 + 0.0im

    # CoulombIntra: sum_i U_i * n0[i] * n1[i]
    for term in data.coulomb_intra_terms
        ri = term.site
        if 0 <= ri < n_site
            contrib = term.value * n0[ri+1] * n1[ri+1]
            e_coulomb_intra += contrib
            e += contrib
        end
    end

    # CoulombInter: sum_{i,j} V_{ij} * (n0[i] + n1[i]) * (n0[j] + n1[j])
    for term in data.coulomb_inter_terms
        ri = term.site1
        rj = term.site2
        if 0 <= ri < n_site && 0 <= rj < n_site
            contrib = term.value * (n0[ri+1] + n1[ri+1]) * (n0[rj+1] + n1[rj+1])
            e_coulomb_inter += contrib
            e += contrib
        end
    end

    # HundCoupling: -sum_{i,j} J_{ij} * (n0[i]*n0[j] + n1[i]*n1[j])
    for term in data.hund_terms
        ri = term.site1
        rj = term.site2
        if 0 <= ri < n_site && 0 <= rj < n_site
            contrib = -term.value * (n0[ri+1] * n0[rj+1] + n1[ri+1] * n1[rj+1])
            e_hund += contrib
            e += contrib
        end
    end
    ctimer_stop!(c_timer, 70)

    # [71] CalHamiltonian1: one-body Transfer terms
    ctimer_start!(c_timer, 71)
    # Transfer: -sum_{i,j,s} t_{ij} * <c†_{i,s} c_{j,s}>
    for term in data.transfer_terms
        ri = term.site1
        rj = term.site2
        if 0 <= ri < n_site && 0 <= rj < n_site
            # Determine spin based on term
            if term.spin == :up || term.spin == :both
                contrib =
                    -term.value * green_func1(
                        ri,
                        rj,
                        0,
                        ip,
                        ele_idx,
                        ele_cfg,
                        ele_num,
                        ele_proj_cnt,
                        data,
                        state;
                        all_complex = all_complex,
                    )
                e_transfer += contrib
                e += contrib
            end
            if term.spin == :down || term.spin == :both
                contrib =
                    -term.value * green_func1(
                        ri,
                        rj,
                        1,
                        ip,
                        ele_idx,
                        ele_cfg,
                        ele_num,
                        ele_proj_cnt,
                        data,
                        state;
                        all_complex = all_complex,
                    )
                e_transfer += contrib
                e += contrib
            end
        end
    end
    ctimer_stop!(c_timer, 71)

    # [72] CalHamiltonian2: two-body terms (PairHopping / Exchange / InterAll)
    ctimer_start!(c_timer, 72)
    # PairHopping: sum_{i,j} P_{ij} * <c†_{i,↑} c_{j,↑} c†_{i,↓} c_{j,↓}>
    for term in data.pair_hop_terms
        ri = term.site1
        rj = term.site2
        if 0 <= ri < n_site && 0 <= rj < n_site
            contrib =
                term.value * green_func2(
                    ri,
                    rj,
                    ri,
                    rj,
                    0,
                    1,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    data,
                    state;
                    all_complex = all_complex,
                )
            e_pair_hop += contrib
            e += contrib
        end
    end

    # ExchangeCoupling: sum_{i,j} J_{ij} * (<c†_{i,↑} c_{j,↑} c†_{j,↓} c_{i,↓}> + <c†_{i,↓} c_{j,↓} c†_{j,↑} c_{i,↑}>)
    n_exchange_zero = 0
    for term in data.exchange_terms
        ri = term.site1
        rj = term.site2
        if 0 <= ri < n_site && 0 <= rj < n_site
            tmp = green_func2(
                ri,
                rj,
                rj,
                ri,
                0,
                1,
                ip,
                ele_idx,
                ele_cfg,
                ele_num,
                ele_proj_cnt,
                data,
                state;
                all_complex = all_complex,
            )
            tmp += green_func2(
                ri,
                rj,
                rj,
                ri,
                1,
                0,
                ip,
                ele_idx,
                ele_cfg,
                ele_num,
                ele_proj_cnt,
                data,
                state;
                all_complex = all_complex,
            )
            if abs(tmp) < 1e-15
                n_exchange_zero += 1
            end
            contrib = term.value * tmp
            e_exchange += contrib
            e += contrib
        end
    end



    # InterAll: sum_{i,j,k,l,s,t} V_{ijkl} * <c†_{i,s} c_{j,s} c†_{k,t} c_{l,t}>
    # Note: InterAllTerm.sites = [ri, rj, rk, rl] (4 sites)
    # In C code: InterAll[idx][0]=ri, [2]=rj, [3]=s, [4]=rk, [6]=rl, [7]=t
    # For now, we assume spins are determined by context or use default (0, 1)
    # TODO: Parse spin information from InterAll file if available
    for term in data.inter_all_terms
        if length(term.sites) >= 4
            ri = term.sites[1]
            rj = term.sites[2]
            rk = term.sites[3]
            rl = term.sites[4]
            # Default: assume first pair is up-spin, second pair is down-spin
            # This may need adjustment based on actual file format
            s = 0  # up-spin for first pair
            t = 1  # down-spin for second pair
            if 0 <= ri < n_site && 0 <= rj < n_site && 0 <= rk < n_site && 0 <= rl < n_site
                e +=
                    term.value * green_func2(
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
            end
        end
    end
    ctimer_stop!(c_timer, 72)

    return e
end

"""
    calculate_hamiltonian_fsz(ip::ComplexF64, ele_idx::Vector{Int}, ele_cfg::Vector{Int},
                              ele_num::Vector{Int}, ele_proj_cnt::Vector{Int}, ele_spn::Vector{Int},
                              data::ExpertModeData, state::VMCOptimizationState;
                              all_complex::Bool = true) -> ComplexF64

Calculate local energy (Hamiltonian expectation value) for FSZ mode.
Uses green_func1_fsz and green_func2_fsz for Green function calculations.

# Arguments
- `ele_spn`: Electron spin array (1-based, values 0 or 1)

# Reference
- C implementation: mVMC/src/mVMC/calham_fsz.c:49-197 (CalculateHamiltonian_fsz)
"""
function calculate_hamiltonian_fsz(
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_cfg::Vector{Int},
    ele_num::Vector{Int},
    ele_proj_cnt::Vector{Int},
    ele_spn::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState;
    all_complex::Bool = true,
    c_timer::CTimer = CTIMER_DISABLED,
)::ComplexF64
    n_site = data.modpara.nsite

    # Get up-spin and down-spin electron numbers
    n0 = @view ele_num[1:n_site]  # up-spin
    n1 = @view ele_num[(n_site+1):(2*n_site)]  # down-spin

    e = 0.0 + 0.0im

    # [70] CalHamiltonian0: diagonal terms (CoulombIntra/CoulombInter/Hund)
    ctimer_start!(c_timer, 70)
    # CoulombIntra: sum_i U_i * n0[i] * n1[i]
    for term in data.coulomb_intra_terms
        ri = term.site
        if 0 <= ri < n_site
            e += term.value * n0[ri+1] * n1[ri+1]
        end
    end

    # CoulombInter: sum_{i,j} V_{ij} * (n0[i] + n1[i]) * (n0[j] + n1[j])
    for term in data.coulomb_inter_terms
        ri = term.site1
        rj = term.site2
        if 0 <= ri < n_site && 0 <= rj < n_site
            e += term.value * (n0[ri+1] + n1[ri+1]) * (n0[rj+1] + n1[rj+1])
        end
    end

    # HundCoupling: -sum_{i,j} J_{ij} * (n0[i]*n0[j] + n1[i]*n1[j])
    for term in data.hund_terms
        ri = term.site1
        rj = term.site2
        if 0 <= ri < n_site && 0 <= rj < n_site
            e += -term.value * (n0[ri+1] * n0[rj+1] + n1[ri+1] * n1[rj+1])
        end
    end
    ctimer_stop!(c_timer, 70)

    # [71] CalHamiltonian1: one-body Transfer terms
    ctimer_start!(c_timer, 71)
    # Transfer: -sum_{i,j,s} t_{ij} * <c†_{i,s} c_{j,s}>
    for term in data.transfer_terms
        ri = term.site1
        rj = term.site2
        if 0 <= ri < n_site && 0 <= rj < n_site
            s = term.spin1
            t = term.spin2
            if s == t
                e +=
                    -term.value * green_func1_fsz(
                        ri,
                        rj,
                        s,
                        ip,
                        ele_idx,
                        ele_cfg,
                        ele_num,
                        ele_proj_cnt,
                        ele_spn,
                        data,
                        state;
                        all_complex = all_complex,
                    )
            else
                e +=
                    -term.value * green_func1_fsz2(
                        ri,
                        rj,
                        s,
                        t,
                        ip,
                        ele_idx,
                        ele_cfg,
                        ele_num,
                        ele_proj_cnt,
                        ele_spn,
                        data,
                        state;
                        all_complex = all_complex,
                    )
            end
        end
    end
    ctimer_stop!(c_timer, 71)

    # [72] CalHamiltonian2: two-body terms (PairHopping / Exchange / InterAll)
    ctimer_start!(c_timer, 72)
    # PairHopping: sum_{i,j} P_{ij} * <c†_{i,↑} c_{j,↑} c†_{i,↓} c_{j,↓}>
    for term in data.pair_hop_terms
        ri = term.site1
        rj = term.site2
        if 0 <= ri < n_site && 0 <= rj < n_site
            e +=
                term.value * green_func2_fsz(
                    ri,
                    rj,
                    ri,
                    rj,
                    0,
                    1,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    ele_spn,
                    data,
                    state;
                    all_complex = all_complex,
                )
        end
    end

    # ExchangeCoupling: sum_{i,j} J_{ij} * (<c†_{i,↑} c_{j,↑} c†_{j,↓} c_{i,↓}> + <c†_{i,↓} c_{j,↓} c†_{j,↑} c_{i,↑}>)
    for term in data.exchange_terms
        ri = term.site1
        rj = term.site2
        if 0 <= ri < n_site && 0 <= rj < n_site
            tmp = green_func2_fsz(
                ri,
                rj,
                rj,
                ri,
                0,
                1,
                ip,
                ele_idx,
                ele_cfg,
                ele_num,
                ele_proj_cnt,
                ele_spn,
                data,
                state;
                all_complex = all_complex,
            )
            tmp += green_func2_fsz(
                ri,
                rj,
                rj,
                ri,
                1,
                0,
                ip,
                ele_idx,
                ele_cfg,
                ele_num,
                ele_proj_cnt,
                ele_spn,
                data,
                state;
                all_complex = all_complex,
            )
            e += term.value * tmp
        end
    end

    # InterAll: sum_{i,j,k,l,s,t,u,v} V_{ijkl} * <c†_{i,s} c_{j,t} c†_{k,u} c_{l,v}>
    # C implementation stores 8 indices per term: site0, spin0, site1, spin1, site2, spin2, site3, spin3
    # The term represents: value * <c†_{site0,spin0} c_{site1,spin1} c†_{site2,spin2} c_{site3,spin3}>
    for term in data.inter_all_terms
        ri = term.site0
        s = term.spin0
        rj = term.site1
        t = term.spin1
        rk = term.site2
        u = term.spin2
        rl = term.site3
        v = term.spin3

        # Validate site indices
        if 0 <= ri < n_site && 0 <= rj < n_site && 0 <= rk < n_site && 0 <= rl < n_site
            # C code: if(s==t && u==v) uses GreenFunc2_fsz, else uses GreenFunc2_fsz2
            if s == t && u == v
                e +=
                    term.value * green_func2_fsz(
                        ri,
                        rj,
                        rk,
                        rl,
                        s,
                        u,
                        ip,
                        ele_idx,
                        ele_cfg,
                        ele_num,
                        ele_proj_cnt,
                        ele_spn,
                        data,
                        state;
                        all_complex = all_complex,
                    )
            else
                # Need to implement green_func2_fsz2 for non-conserved Sz terms
                # For now, use the general case with all 4 spins
                e +=
                    term.value * green_func2_fsz2(
                        ri,
                        rj,
                        rk,
                        rl,
                        s,
                        t,
                        u,
                        v,
                        ip,
                        ele_idx,
                        ele_cfg,
                        ele_num,
                        ele_proj_cnt,
                        ele_spn,
                        data,
                        state;
                        all_complex = all_complex,
                    )
            end
        end
    end
    ctimer_stop!(c_timer, 72)

    return e
end

# ============================================================================
# SR Optimization Quantities
# ============================================================================

"""
    slater_elm_diff_fcmp!(sr_opt_o::AbstractVector{ComplexF64}, ip::ComplexF64,
                          ele_idx::Vector{Int}, data::ExpertModeData,
                          state::VMCOptimizationState)

Calculate Slater element derivative: Tr[M^{-1} * ∂M/∂f_{ij}].
Equivalent to C's `SlaterElmDiff_fcmp()`.
"""
function slater_elm_diff_fcmp!(
    sr_opt_o::AbstractVector{ComplexF64},
    ip::ComplexF64,
    ele_idx::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_qp_full = get_n_qp_full(data)

    # Use n_orbital_idx (number of unique orbital parameters) instead of length(orbital_terms)
    # C implementation: NSlater = NOrbitalIdx
    n_slater = if data.modpara.n_orbital_idx > 0
        data.modpara.n_orbital_idx
    elseif !isempty(data.orbital_terms)
        maximum(t.idx for t in data.orbital_terms) + 1
    else
        0
    end

    # NMPTrans from modpara.def (absolute value, as C implementation does)
    n_mp_trans = abs(data.modpara.nmp_trans)
    n_sp_gauss_leg = data.modpara.nsp_gauss_leg

    if n_slater == 0 || n_qp_full == 0
        return
    end

    inv_ip = 1.0 / ip

    # Get inverse matrix view
    n_size_sq = n_size * n_size
    inv_m_flat::Vector{ComplexF64} = state.slater_matrix.inv_m

    # Initialize buffer for accumulation
    buffer = zeros(ComplexF64, n_qp_full * n_slater)

    # Check if QPTrans is available
    if isempty(data.qp_trans) || isempty(data.qp_trans_sgn)
        @warn "QPTrans mappings not found, skipping SlaterElmDiff calculation"
        return
    end

    # Build transOrbIdx and transOrbSgn arrays (pre-compute orbital indices)
    n_trans = n_mp_trans * data.n_qp_opt_trans
    trans_orb_idx = zeros(Int, n_trans * n_size * n_size)
    trans_orb_sgn = zeros(Int, n_trans * n_size * n_size)

    # Get orbital index mapping
    orbital_idx = data.orbital_terms

    # Pre-compute transOrbIdx and transOrbSgn
    for qpidx = 0:(n_trans-1)
        optidx = qpidx ÷ n_mp_trans
        mpidx = qpidx % n_mp_trans

        if mpidx + 1 > n_mp_trans || optidx + 1 > data.n_qp_opt_trans
            continue
        end

        # Get QP translation arrays
        xqp_opt =
            optidx + 1 <= length(data.qp_opt_trans) ? data.qp_opt_trans[optidx+1] : Int[]
        xqp_opt_sgn =
            optidx + 1 <= length(data.qp_opt_trans_sgn) ? data.qp_opt_trans_sgn[optidx+1] :
            Int[]
        xqp = mpidx + 1 <= length(data.qp_trans) ? data.qp_trans[mpidx+1] : Int[]
        xqp_sgn =
            mpidx + 1 <= length(data.qp_trans_sgn) ? data.qp_trans_sgn[mpidx+1] : Int[]

        if isempty(xqp) || isempty(xqp_sgn)
            continue
        end

        t_orb_idx_base = qpidx * n_size * n_size
        t_orb_sgn_base = qpidx * n_size * n_size

        for msi = 0:(n_size-1)
            ri = ele_idx[msi+1]  # Site index (0-based, stored in ele_idx)
            if ri < 0 || ri >= n_site
                continue
            end

            # Apply QPOptTrans and QPTrans
            # ele_idx contains 0-based site indices, so ri is 0-based
            # xqp_opt and xqp_opt_sgn are 1-based arrays (Julia indexing)
            ori = isempty(xqp_opt) || (ri + 1) > length(xqp_opt) ? ri : xqp_opt[ri+1]  # 0-based
            tri = (ori + 1) > length(xqp) ? ori : xqp[ori+1]  # 0-based
            sgni =
                (
                    isempty(xqp_opt_sgn) || (ri + 1) > length(xqp_opt_sgn) ? 1 :
                    xqp_opt_sgn[ri+1]
                ) * ((ori + 1) > length(xqp_sgn) ? 1 : xqp_sgn[ori+1])

            t_orb_idx_i = t_orb_idx_base + msi * n_size
            t_orb_sgn_i = t_orb_sgn_base + msi * n_size

            for msj = 0:(n_size-1)
                rj = ele_idx[msj+1]  # Site index (0-based)
                if rj < 0 || rj >= n_site
                    continue
                end

                orj = isempty(xqp_opt) || (rj + 1) > length(xqp_opt) ? rj : xqp_opt[rj+1]
                trj = (orj + 1) > length(xqp) ? orj : xqp[orj+1]
                sgnj =
                    (
                        isempty(xqp_opt_sgn) || (rj + 1) > length(xqp_opt_sgn) ? 1 :
                        xqp_opt_sgn[rj+1]
                    ) * ((orj + 1) > length(xqp_sgn) ? 1 : xqp_sgn[orj+1])

                # Find orbital index using orbital_idx_matrix (like C's OrbitalIdx[tri][trj])
                if data.orbital_idx_matrix !== nothing &&
                   tri + 1 <= size(data.orbital_idx_matrix, 1) &&
                   trj + 1 <= size(data.orbital_idx_matrix, 2)
                    orbidx = data.orbital_idx_matrix[tri+1, trj+1]
                else
                    # Fallback to search if matrix not available
                    orbidx = find_orbital_idx(tri, trj, orbital_idx, n_site)
                end
                trans_orb_idx[t_orb_idx_i+msj+1] = orbidx

                # Calculate sign including OrbitalSgn
                sgn = sgni * sgnj
                if data.orbital_sgn !== nothing &&
                   tri + 1 <= size(data.orbital_sgn, 1) &&
                   trj + 1 <= size(data.orbital_sgn, 2)
                    sgn *= data.orbital_sgn[tri+1, trj+1]
                end
                trans_orb_sgn[t_orb_sgn_i+msj+1] = sgn
            end
        end
    end

    # Get spin projection parameters
    sp_gl_cos_sin::Vector{ComplexF64},
    sp_gl_cos_cos::Vector{ComplexF64},
    sp_gl_sin_sin::Vector{ComplexF64} = get_spin_proj_params(data)

    # Loop over QP indices and accumulate buffer
    for qpidx = 0:(n_qp_full-1)
        mpidx = qpidx ÷ n_sp_gauss_leg
        spidx = qpidx % n_sp_gauss_leg

        if mpidx + 1 > n_mp_trans || spidx + 1 > length(sp_gl_cos_sin)
            continue
        end

        # Spin rotation factors
        cs = state.slater_matrix.pf_m[qpidx+1] * sp_gl_cos_sin[spidx+1]
        cc = state.slater_matrix.pf_m[qpidx+1] * sp_gl_cos_cos[spidx+1]
        ss = state.slater_matrix.pf_m[qpidx+1] * sp_gl_sin_sin[spidx+1]

        # Get inverse matrix for this QP index
        # inv_m_flat is stored in row-major (C) layout: InvM[msi*n_size + msj]
        # qpidx offset: qpidx * n_size * n_size
        inv_m_base = qpidx * n_size_sq

        # Get pre-computed orbital indices for this mpidx
        t_orb_idx_base = mpidx * n_size * n_size
        t_orb_sgn_base = mpidx * n_size * n_size

        # Accumulate contributions
        # inv_m_flat is stored in row-major (C) layout: InvM[msi*n_size + msj]
        # Access: inv_m_flat[inv_m_base + msi * n_size + msj + 1]
        for msi = 0:(n_elec-1)  # up-spin electrons
            t_orb_idx_i = t_orb_idx_base + msi * n_size
            t_orb_sgn_i = t_orb_sgn_base + msi * n_size

            for msj = 0:(n_elec-1)  # up-up
                orbidx = trans_orb_idx[t_orb_idx_i+msj+1]
                if orbidx >= 0 && orbidx < n_slater
                    sgn = trans_orb_sgn[t_orb_sgn_i+msj+1]
                    # C: InvM[msi][msj] = InvM[msi*n_size + msj]
                    inv_idx = inv_m_base + msi * n_size + msj + 1
                    buffer[qpidx*n_slater+orbidx+1] += (inv_m_flat[inv_idx] * cs) * sgn
                end
            end

            for msj = n_elec:(n_size-1)  # up-down
                orbidx = trans_orb_idx[t_orb_idx_i+msj+1]
                if orbidx >= 0 && orbidx < n_slater
                    sgn = trans_orb_sgn[t_orb_sgn_i+msj+1]
                    inv_idx = inv_m_base + msi * n_size + msj + 1
                    buffer[qpidx*n_slater+orbidx+1] -= (inv_m_flat[inv_idx] * cc) * sgn
                end
            end
        end

        for msi = n_elec:(n_size-1)  # down-spin electrons
            t_orb_idx_i = t_orb_idx_base + msi * n_size
            t_orb_sgn_i = t_orb_sgn_base + msi * n_size

            for msj = 0:(n_elec-1)  # down-up
                orbidx = trans_orb_idx[t_orb_idx_i+msj+1]
                if orbidx >= 0 && orbidx < n_slater
                    sgn = trans_orb_sgn[t_orb_sgn_i+msj+1]
                    inv_idx = inv_m_base + msi * n_size + msj + 1
                    buffer[qpidx*n_slater+orbidx+1] += (inv_m_flat[inv_idx] * ss) * sgn
                end
            end

            for msj = n_elec:(n_size-1)  # down-down
                orbidx = trans_orb_idx[t_orb_idx_i+msj+1]
                if orbidx >= 0 && orbidx < n_slater
                    sgn = trans_orb_sgn[t_orb_sgn_i+msj+1]
                    inv_idx = inv_m_base + msi * n_size + msj + 1
                    buffer[qpidx*n_slater+orbidx+1] -= (inv_m_flat[inv_idx] * cs) * sgn
                end
            end
        end
    end

    # Accumulate results weighted by QPFullWeight
    if data.qp_weights === nothing
        return
    end
    qp_full_weight::Vector{ComplexF64} = data.qp_weights.qp_full_weight

    # Note: This function receives a view of sr_opt_o that starts AFTER projection parameters
    # So we don't need to add n_proj offset here - the view already starts at the right position
    # sr_opt_o[1] here corresponds to the full array's sr_opt_o[2*n_proj + 3] (the first Slater parameter)

    for orbidx = 0:(n_slater-1)
        # Index in the view: 2*orbidx + 1 for real part, 2*orbidx + 2 for imaginary part
        real_idx = 2 * orbidx + 1
        imag_idx = 2 * orbidx + 2

        # Ensure indices are within bounds
        if real_idx > length(sr_opt_o) || imag_idx > length(sr_opt_o)
            @warn "sr_opt_o index out of bounds: real_idx=$real_idx, imag_idx=$imag_idx, length=$(length(sr_opt_o)), orbidx=$orbidx, n_slater=$n_slater" maxlog=1
            continue
        end

        sr_opt_o[real_idx] = 0.0 + 0.0im  # real part index
        sr_opt_o[imag_idx] = 0.0 + 0.0im  # imaginary part index

        @inbounds for qpidx = 0:(n_qp_full-1)
            if qpidx + 1 <= length(qp_full_weight)
                tmp = qp_full_weight[qpidx+1] * buffer[qpidx*n_slater+orbidx+1]
                sr_opt_o[real_idx] += tmp
                sr_opt_o[imag_idx] += tmp * im
            end
        end

        sr_opt_o[real_idx] *= inv_ip
        sr_opt_o[imag_idx] *= inv_ip
    end
end

"""
    slater_elm_diff_fsz!(sr_opt_o::AbstractVector{ComplexF64}, ip::ComplexF64,
                         ele_idx::Vector{Int}, ele_spn::Vector{Int},
                         data::ExpertModeData, state::VMCOptimizationState)

Calculate Slater element derivative for FSZ mode: Tr[M^{-1} * ∂M/∂f_{ij}].
Uses ele_spn array to determine electron spins.

Equivalent to C's `SlaterElmDiff_fsz()`.

# Reference
- C implementation: mVMC/src/mVMC/slater_fsz.c:119-245 (SlaterElmDiff_fsz)
"""
function slater_elm_diff_fsz!(
    sr_opt_o::AbstractVector{ComplexF64},
    ip::ComplexF64,
    ele_idx::Vector{Int},
    ele_spn::Vector{Int},
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_qp_full = get_n_qp_full(data)

    n_slater = if data.modpara.n_orbital_idx > 0
        data.modpara.n_orbital_idx
    elseif !isempty(data.orbital_terms)
        maximum(t.idx for t in data.orbital_terms) + 1
    else
        0
    end

    n_mp_trans = abs(data.modpara.nmp_trans)
    n_sp_gauss_leg = data.modpara.nsp_gauss_leg

    if n_slater == 0 || n_qp_full == 0
        return
    end

    inv_ip = 1.0 / ip

    n_size_sq = n_size * n_size
    inv_m_flat::Vector{ComplexF64} = state.slater_matrix.inv_m

    buffer = zeros(ComplexF64, n_qp_full * n_slater)

    if isempty(data.qp_trans) || isempty(data.qp_trans_sgn)
        @warn "QPTrans mappings not found, skipping SlaterElmDiff_fsz calculation"
        return
    end

    # Build transOrbIdx and transOrbSgn arrays using ele_spn
    n_trans = n_mp_trans * data.n_qp_opt_trans
    trans_orb_idx = zeros(Int, n_trans * n_size * n_size)
    trans_orb_sgn = zeros(Int, n_trans * n_size * n_size)

    orbital_idx = data.orbital_terms

    # Pre-compute transOrbIdx and transOrbSgn (FSZ version)
    for qpidx = 0:(n_trans-1)
        optidx = qpidx ÷ n_mp_trans
        mpidx = qpidx % n_mp_trans

        if mpidx + 1 > n_mp_trans || optidx + 1 > data.n_qp_opt_trans
            continue
        end

        xqp_opt =
            optidx + 1 <= length(data.qp_opt_trans) ? data.qp_opt_trans[optidx+1] : Int[]
        xqp_opt_sgn =
            optidx + 1 <= length(data.qp_opt_trans_sgn) ? data.qp_opt_trans_sgn[optidx+1] :
            Int[]
        xqp = mpidx + 1 <= length(data.qp_trans) ? data.qp_trans[mpidx+1] : Int[]
        xqp_sgn =
            mpidx + 1 <= length(data.qp_trans_sgn) ? data.qp_trans_sgn[mpidx+1] : Int[]

        if isempty(xqp) || isempty(xqp_sgn)
            continue
        end

        t_orb_idx_base = qpidx * n_size * n_size
        t_orb_sgn_base = qpidx * n_size * n_size

        for msi = 0:(n_size-1)
            ri = ele_idx[msi+1]
            si = ele_spn[msi+1]  # FSZ: get spin from ele_spn
            if ri < 0 || ri >= n_site
                continue
            end

            ori = isempty(xqp_opt) || (ri + 1) > length(xqp_opt) ? ri : xqp_opt[ri+1]
            tri = (ori + 1) > length(xqp) ? ori : xqp[ori+1]
            tri = tri + si * n_site  # FSZ: add spin offset
            sgni =
                (
                    isempty(xqp_opt_sgn) || (ri + 1) > length(xqp_opt_sgn) ? 1 :
                    xqp_opt_sgn[ri+1]
                ) * ((ori + 1) > length(xqp_sgn) ? 1 : xqp_sgn[ori+1])

            t_orb_idx_i = t_orb_idx_base + msi * n_size
            t_orb_sgn_i = t_orb_sgn_base + msi * n_size

            for msj = 0:(n_size-1)
                rj = ele_idx[msj+1]
                sj = ele_spn[msj+1]  # FSZ: get spin from ele_spn
                if rj < 0 || rj >= n_site
                    continue
                end

                orj = isempty(xqp_opt) || (rj + 1) > length(xqp_opt) ? rj : xqp_opt[rj+1]
                trj = (orj + 1) > length(xqp) ? orj : xqp[orj+1]
                trj = trj + sj * n_site  # FSZ: add spin offset
                sgnj =
                    (
                        isempty(xqp_opt_sgn) || (rj + 1) > length(xqp_opt_sgn) ? 1 :
                        xqp_opt_sgn[rj+1]
                    ) * ((orj + 1) > length(xqp_sgn) ? 1 : xqp_sgn[orj+1])

                # Find orbital index
                if data.orbital_idx_matrix !== nothing &&
                   tri + 1 <= size(data.orbital_idx_matrix, 1) &&
                   trj + 1 <= size(data.orbital_idx_matrix, 2)
                    orbidx = data.orbital_idx_matrix[tri+1, trj+1]
                else
                    orbidx = find_orbital_idx(tri, trj, orbital_idx, n_site)
                end
                trans_orb_idx[t_orb_idx_i+msj+1] = orbidx

                sgn = sgni * sgnj
                if data.orbital_sgn !== nothing &&
                   tri + 1 <= size(data.orbital_sgn, 1) &&
                   trj + 1 <= size(data.orbital_sgn, 2)
                    sgn *= data.orbital_sgn[tri+1, trj+1]
                end
                trans_orb_sgn[t_orb_sgn_i+msj+1] = sgn
            end
        end
    end

    # FSZ: no spin projection, use PfM directly
    for qpidx = 0:(n_qp_full-1)
        mpidx = qpidx ÷ n_sp_gauss_leg

        if mpidx + 1 > n_mp_trans
            continue
        end

        # FSZ: cs = PfM[qpidx] (no spin projection factor)
        cs = state.slater_matrix.pf_m[qpidx+1]

        inv_m_base = qpidx * n_size_sq
        t_orb_idx_base = mpidx * n_size * n_size
        t_orb_sgn_base = mpidx * n_size * n_size

        # FSZ: simple loop without spin block separation
        for msi = 0:(n_size-1)
            t_orb_idx_i = t_orb_idx_base + msi * n_size
            t_orb_sgn_i = t_orb_sgn_base + msi * n_size

            for msj = 0:(n_size-1)
                orbidx = trans_orb_idx[t_orb_idx_i+msj+1]
                if orbidx >= 0 && orbidx < n_slater
                    sgn = trans_orb_sgn[t_orb_sgn_i+msj+1]
                    inv_idx = inv_m_base + msi * n_size + msj + 1
                    # FSZ: Tr(X^{-1}*dX/df) = -1.0 * invM * cs * sgn
                    buffer[qpidx*n_slater+orbidx+1] += -1.0 * inv_m_flat[inv_idx] * cs * sgn
                end
            end
        end
    end

    # Accumulate results
    if data.qp_weights === nothing
        return
    end
    qp_full_weight::Vector{ComplexF64} = data.qp_weights.qp_full_weight

    for orbidx = 0:(n_slater-1)
        real_idx = 2 * orbidx + 1
        imag_idx = 2 * orbidx + 2

        if real_idx > length(sr_opt_o) || imag_idx > length(sr_opt_o)
            continue
        end

        sr_opt_o[real_idx] = 0.0 + 0.0im
        sr_opt_o[imag_idx] = 0.0 + 0.0im

        @inbounds for qpidx = 0:(n_qp_full-1)
            if qpidx + 1 <= length(qp_full_weight)
                tmp = qp_full_weight[qpidx+1] * buffer[qpidx*n_slater+orbidx+1]
                sr_opt_o[real_idx] += tmp
                sr_opt_o[imag_idx] += tmp * im
            end
        end

        sr_opt_o[real_idx] *= inv_ip
        sr_opt_o[imag_idx] *= inv_ip
    end
end

function opt_trans_diff!(
    sr_opt_o::AbstractVector{ComplexF64},
    ip::ComplexF64,
    data::ExpertModeData,
    state::VMCOptimizationState,
)
    n_opt_trans = MVMCExpertModeParsers.count_opt_trans_parameters(data)
    n_opt_trans == 0 && return
    data.qp_weights === nothing && return

    n_qp_fix = length(data.qp_weights.qp_fix_weight)
    n_qp_fix > 0 || return
    pf_values = state.slater_matrix.pf_m

    for optidx = 1:n_opt_trans
        acc = 0.0 + 0.0im
        base = (optidx - 1) * n_qp_fix
        for j = 1:n_qp_fix
            pf_idx = base + j
            if j <= length(data.qp_weights.qp_fix_weight) && pf_idx <= length(pf_values)
                acc += data.qp_weights.qp_fix_weight[j] * pf_values[pf_idx]
            end
        end

        val = acc / ip
        real_idx = 2 * (optidx - 1) + 1
        imag_idx = real_idx + 1
        if imag_idx <= length(sr_opt_o)
            sr_opt_o[real_idx] = val
            sr_opt_o[imag_idx] = im * val
        end
    end
end

"""
    find_orbital_idx(tri::Int, trj::Int, orbital_terms::Vector, n_site::Int) -> Int

Find the orbital index for the given translated site pair.
Returns -1 if not found.
"""
function find_orbital_idx(tri::Int, trj::Int, orbital_terms::Vector, n_site::Int)::Int
    for term in orbital_terms
        if term.site1 == tri && term.site2 == trj
            return term.idx  # Return the orbital parameter index from term.idx
        end
    end
    return -1
end

"""
    get_spin_proj_params(data::ExpertModeData) -> Tuple{Vector{ComplexF64}, Vector{ComplexF64}, Vector{ComplexF64}}

Get spin projection parameters (SPGLCosSin, SPGLCosCos, SPGLSinSin).
"""
function get_spin_proj_params(
    data::ExpertModeData,
)::Tuple{Vector{ComplexF64},Vector{ComplexF64},Vector{ComplexF64}}
    # Use the properly initialized values from qp_weights
    if data.qp_weights !== nothing
        qp_w = data.qp_weights
        if !isempty(qp_w.spgl_cos_sin) &&
           !isempty(qp_w.spgl_cos_cos) &&
           !isempty(qp_w.spgl_sin_sin)
            sp_gl_cos_sin::Vector{ComplexF64} = qp_w.spgl_cos_sin
            sp_gl_cos_cos::Vector{ComplexF64} = qp_w.spgl_cos_cos
            sp_gl_sin_sin::Vector{ComplexF64} = qp_w.spgl_sin_sin
            return sp_gl_cos_sin, sp_gl_cos_cos, sp_gl_sin_sin
        end
    end

    # Fallback: For simple case with NSPGaussLeg=1, use identity
    n_sp = data.modpara.nsp_gauss_leg
    if n_sp <= 0
        n_sp = 1
    end

    sp_gl_cos_sin = fill(ComplexF64(1.0), n_sp)
    sp_gl_cos_cos = fill(ComplexF64(1.0), n_sp)
    sp_gl_sin_sin = fill(ComplexF64(1.0), n_sp)

    return sp_gl_cos_sin, sp_gl_cos_cos, sp_gl_sin_sin
end

"""
    set_projection_diff!(sr_opt_o::AbstractVector{ComplexF64}, ele_proj_cnt::Vector{Int}, n_proj::Int)

Set projection factor derivatives in sr_opt_o.
Equivalent to C's projection calculation in VMCMainCal.
"""
function set_projection_diff!(
    sr_opt_o::AbstractVector{ComplexF64},
    ele_proj_cnt::Vector{Int},
    n_proj::Int,
)
    # srOptO[0] = 1.0 + 0.0*I
    # srOptO[1] = 0.0 + 0.0*I
    sr_opt_o[1] = 1.0 + 0.0im  # real part
    sr_opt_o[2] = 0.0 + 0.0im  # imaginary part

    for i = 0:(n_proj-1)
        if i < length(ele_proj_cnt)
            sr_opt_o[(i+1)*2+1] = ComplexF64(ele_proj_cnt[i+1])  # even: real
            sr_opt_o[(i+1)*2+2] = 0.0 + 0.0im  # odd: complex part
        end
    end
end

"""
    set_rbm_diff!(sr_opt_o, rbm_cnt, ele_num, data)

Set RBM derivatives in `sr_opt_o` (view starting at C offset `SROptO + 2*NProj + 2`).
Equivalent to C's `RBMDiff()`.
"""
function set_rbm_diff!(
    sr_opt_o::AbstractVector{ComplexF64},
    rbm_cnt::Vector{ComplexF64},
    ele_num::Vector{Int},
    data::ExpertModeData,
)
    n_site = data.modpara.nsite
    n_site2 = 2 * n_site

    n_charge_phys = _rbm_nidx(data.charge_rbm_phys_layer_terms)
    n_spin_phys = _rbm_nidx(data.spin_rbm_phys_layer_terms)
    n_general_phys = _rbm_nidx(data.general_rbm_phys_layer_terms)
    n_phys = n_charge_phys + n_spin_phys + n_general_phys

    n_charge_hidden = _rbm_nidx(data.charge_rbm_hidden_layer_terms)
    n_spin_hidden = _rbm_nidx(data.spin_rbm_hidden_layer_terms)
    n_general_hidden = _rbm_nidx(data.general_rbm_hidden_layer_terms)
    n_hidden = n_charge_hidden + n_spin_hidden + n_general_hidden

    n_charge_ph = _rbm_nidx(data.charge_rbm_phys_hidden_terms)
    n_spin_ph = _rbm_nidx(data.spin_rbm_phys_hidden_terms)
    n_general_ph = _rbm_nidx(data.general_rbm_phys_hidden_terms)
    n_ph = n_charge_ph + n_spin_ph + n_general_ph

    n_rbm = n_phys + n_hidden + n_ph
    if n_rbm == 0
        return
    end

    if length(sr_opt_o) < 2 * n_rbm
        @warn "set_rbm_diff!: sr_opt_o too small for RBM derivatives (need $(2*n_rbm), got $(length(sr_opt_o)))" maxlog=1
        return
    end

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

    rbm_cnt_len = length(rbm_cnt)
    n_eff_phys = min(n_phys, rbm_cnt_len)

    # Physical-layer derivatives: O[idx] = rbmCnt[idx]
    @inbounds for idx0 = 0:(n_eff_phys-1)
        ctmp = rbm_cnt[idx0+1]
        sr_opt_o[2*idx0+1] = ctmp
        sr_opt_o[2*idx0+2] = im * ctmp
    end

    n0 = @view ele_num[1:n_site]
    n1 = @view ele_num[(n_site+1):n_site2]

    # RBM parameter offsets (0-based in C layout)
    hidden_param_offset = n_phys
    ph_param_offset = n_phys + n_hidden

    # RBM counter offsets (1-based in rbm_cnt)
    charge_cnt_offset = n_phys
    spin_cnt_offset = n_phys + n_charge_neuron
    general_cnt_offset = n_phys + n_charge_neuron + n_spin_neuron

    # Hidden-layer parameter derivatives.
    for term in data.charge_rbm_hidden_layer_terms
        hi = term.site
        idx = term.idx
        if !(0 <= hi < n_charge_neuron && 0 <= idx < n_charge_hidden)
            continue
        end
        cnt_idx = charge_cnt_offset + hi + 1
        cnt_idx > rbm_cnt_len && continue
        ctmp = tanh(rbm_cnt[cnt_idx])
        param_idx0 = hidden_param_offset + idx
        sr_opt_o[2*param_idx0+1] += ctmp
        sr_opt_o[2*param_idx0+2] += im * ctmp
    end
    for term in data.spin_rbm_hidden_layer_terms
        hi = term.site
        idx = term.idx
        if !(0 <= hi < n_spin_neuron && 0 <= idx < n_spin_hidden)
            continue
        end
        cnt_idx = spin_cnt_offset + hi + 1
        cnt_idx > rbm_cnt_len && continue
        ctmp = tanh(rbm_cnt[cnt_idx])
        param_idx0 = hidden_param_offset + n_charge_hidden + idx
        sr_opt_o[2*param_idx0+1] += ctmp
        sr_opt_o[2*param_idx0+2] += im * ctmp
    end
    for term in data.general_rbm_hidden_layer_terms
        hi = term.site
        idx = term.idx
        if !(0 <= hi < n_general_neuron && 0 <= idx < n_general_hidden)
            continue
        end
        cnt_idx = general_cnt_offset + hi + 1
        cnt_idx > rbm_cnt_len && continue
        ctmp = tanh(rbm_cnt[cnt_idx])
        param_idx0 = hidden_param_offset + n_charge_hidden + n_spin_hidden + idx
        sr_opt_o[2*param_idx0+1] += ctmp
        sr_opt_o[2*param_idx0+2] += im * ctmp
    end

    # Phys-hidden coupling derivatives.
    for term in data.charge_rbm_phys_hidden_terms
        ri = term.site1
        hi = term.site2
        idx = term.idx
        if !(0 <= ri < n_site && 0 <= hi < n_charge_neuron && 0 <= idx < n_charge_ph)
            continue
        end
        cnt_idx = charge_cnt_offset + hi + 1
        cnt_idx > rbm_cnt_len && continue
        xi = n0[ri+1] + n1[ri+1] - 1
        ctmp = xi * tanh(rbm_cnt[cnt_idx])
        param_idx0 = ph_param_offset + idx
        sr_opt_o[2*param_idx0+1] += ctmp
        sr_opt_o[2*param_idx0+2] += im * ctmp
    end
    for term in data.spin_rbm_phys_hidden_terms
        ri = term.site1
        hi = term.site2
        idx = term.idx
        if !(0 <= ri < n_site && 0 <= hi < n_spin_neuron && 0 <= idx < n_spin_ph)
            continue
        end
        cnt_idx = spin_cnt_offset + hi + 1
        cnt_idx > rbm_cnt_len && continue
        xi = n0[ri+1] - n1[ri+1]
        ctmp = xi * tanh(rbm_cnt[cnt_idx])
        param_idx0 = ph_param_offset + n_charge_ph + idx
        sr_opt_o[2*param_idx0+1] += ctmp
        sr_opt_o[2*param_idx0+2] += im * ctmp
    end
    for term in data.general_rbm_phys_hidden_terms
        ri = term.site1
        si = term.spin
        hi = term.site2
        idx = term.idx
        if !(0 <= ri < n_site &&
             (si == 0 || si == 1) &&
             0 <= hi < n_general_neuron &&
             0 <= idx < n_general_ph)
            continue
        end
        cnt_idx = general_cnt_offset + hi + 1
        cnt_idx > rbm_cnt_len && continue
        rsi = ri + si * n_site
        xi = 2 * ele_num[rsi+1] - 1
        ctmp = xi * tanh(rbm_cnt[cnt_idx])
        param_idx0 = ph_param_offset + n_charge_ph + n_spin_ph + idx
        sr_opt_o[2*param_idx0+1] += ctmp
        sr_opt_o[2*param_idx0+2] += im * ctmp
    end

    return
end

"""
    calculate_oo!(sr_opt_oo::Vector{ComplexF64}, sr_opt_ho::Vector{ComplexF64},
                 sr_opt_o::Vector{ComplexF64}, w::Float64, e::ComplexF64, sr_opt_size::Int)

Calculate <O†O> and <HO> matrices/vectors.
Equivalent to C's `calculateOO()`.
"""
function calculate_oo!(
    sr_opt_oo::Vector{ComplexF64},
    sr_opt_ho::Vector{ComplexF64},
    sr_opt_o::Vector{ComplexF64},
    w::Float64,
    e::ComplexF64,
    sr_opt_size::Int,
    ;
    threaded::Bool = false,
)
    size_2 = 2 * sr_opt_size

    # Update <O> and <HO>
    if vmc_inner_threading_enabled(size_2, threaded)
        Base.Threads.@threads :static for j = 0:(size_2-1)
            @inbounds begin
                tmp = w * sr_opt_o[j+1]
                sr_opt_oo[0*size_2+j+1] += tmp  # First row: <O>
                sr_opt_ho[j+1] += e * tmp  # <HO>
            end
        end
    else
        @inbounds @simd for j = 0:(size_2-1)
            tmp = w * sr_opt_o[j+1]
            sr_opt_oo[0*size_2+j+1] += tmp  # First row: <O>
            sr_opt_ho[j+1] += e * tmp  # <HO>
        end
    end

    # Update <O†O>
    if vmc_inner_threading_enabled(size_2 - 2, threaded)
        Base.Threads.@threads :static for i = 2:(size_2-1)
            @inbounds for j = 0:(size_2-1)
                sr_opt_oo[i*size_2+j+1] += w * sr_opt_o[j+1] * conj(sr_opt_o[i+1])
            end
        end
    else
        @inbounds for i = 2:(size_2-1)
            for j = 0:(size_2-1)
                sr_opt_oo[i*size_2+j+1] += w * sr_opt_o[j+1] * conj(sr_opt_o[i+1])
            end
        end
    end
end

"""
    calculate_oo_real!(sr_opt_oo::Vector{Float64}, sr_opt_ho::Vector{Float64},
                      sr_opt_o::Vector{Float64}, w::Float64, e::Float64, sr_opt_size::Int)

Calculate <O†O> and <HO> matrices/vectors (real version).
Equivalent to C's `calculateOO_real()`.
"""
function calculate_oo_real!(
    sr_opt_oo::Vector{Float64},
    sr_opt_ho::Vector{Float64},
    sr_opt_o::Vector{Float64},
    w::Float64,
    e::Float64,
    sr_opt_size::Int,
    ;
    threaded::Bool = false,
)
    we = w * e

    # OO[i][j] += w * O[i] * O[j]
    # C implementation: M_DGER with lda=srOptSize (column-major storage)
    # BLAS column-major: A[i,j] = A[i + lda*j] for i in 0..m-1, j in 0..n-1
    # So OO[i,j] is at index i + srOptSize*j (0-based) or (i-1) + srOptSize*(j-1) + 1 (1-based)
    lda = sr_opt_size
    if vmc_inner_threading_enabled(sr_opt_size, threaded)
        Base.Threads.@threads :static for i = 1:sr_opt_size
            @inbounds for j = 1:sr_opt_size
                # Column-major indexing: idx = row + lda * column (0-based)
                idx = (i - 1) + lda * (j - 1) + 1  # 1-based Julia index
                sr_opt_oo[idx] += w * sr_opt_o[i] * sr_opt_o[j]
            end
        end
    else
        @inbounds for i = 1:sr_opt_size
            for j = 1:sr_opt_size
                # Column-major indexing: idx = row + lda * column (0-based)
                idx = (i - 1) + lda * (j - 1) + 1  # 1-based Julia index
                sr_opt_oo[idx] += w * sr_opt_o[i] * sr_opt_o[j]
            end
        end
    end

    # HO[i] += w * e * O[i]
    if vmc_inner_threading_enabled(sr_opt_size, threaded)
        Base.Threads.@threads :static for i = 1:sr_opt_size
            @inbounds sr_opt_ho[i] += we * sr_opt_o[i]
        end
    else
        @inbounds @simd for i = 1:sr_opt_size
            sr_opt_ho[i] += we * sr_opt_o[i]
        end
    end
end

"""
    calculate_oo_store!(sr_opt_oo::Vector{ComplexF64}, sr_opt_ho::Vector{ComplexF64},
                       sr_opt_o_store::Vector{ComplexF64}, w::Float64, e::ComplexF64,
                       sample::Int, sr_opt_size::Int, sqrtw::Float64)

Store sample O and accumulate HO.
Equivalent to C's sample storage in VMCMainCal.
"""
function calculate_oo_store!(
    sr_opt_oo::Vector{ComplexF64},
    sr_opt_ho::Vector{ComplexF64},
    sr_opt_o_store::Vector{ComplexF64},
    sr_opt_o::Vector{ComplexF64},
    w::Float64,
    e::ComplexF64,
    sample::Int,
    sr_opt_size::Int,
    ;
    threaded::Bool = false,
)
    we = w * e
    sqrtw = sqrt(w)
    size_2 = 2 * sr_opt_size

    if vmc_inner_threading_enabled(size_2, threaded)
        Base.Threads.@threads :static for i = 0:(size_2-1)
            @inbounds begin
                # Store sqrt(w) * O for later matrix multiplication
                sr_opt_o_store[i+sample*size_2+1] = sqrtw * sr_opt_o[i+1]
                # Accumulate HO
                sr_opt_ho[i+1] += we * sr_opt_o[i+1]
            end
        end
    else
        @inbounds @simd for i = 0:(size_2-1)
            # Store sqrt(w) * O for later matrix multiplication
            sr_opt_o_store[i+sample*size_2+1] = sqrtw * sr_opt_o[i+1]
            # Accumulate HO
            sr_opt_ho[i+1] += we * sr_opt_o[i+1]
        end
    end
end

"""
    finalize_oo_store!(sr_opt_oo::Vector{ComplexF64}, sr_opt_o_store::Vector{ComplexF64},
                       sr_opt_size::Int, sample_size::Int)

Finalize <O†O> calculation from stored samples using matrix multiplication.
Equivalent to C's `calculateOO_Store()`.
"""
function finalize_oo_store!(
    sr_opt_oo::Vector{ComplexF64},
    sr_opt_o_store::Vector{ComplexF64},
    sr_opt_size::Int,
    sample_size::Int,
    ;
    threaded::Bool = false,
)
    size_2 = 2 * sr_opt_size

    # sr_opt_oo = sr_opt_o_store' * sr_opt_o_store
    # Reshape for matrix multiplication
    O_store = reshape(sr_opt_o_store[1:(size_2*sample_size)], size_2, sample_size)

    # OO = O * O^H (Hermitian)
    if vmc_inner_threading_enabled(size_2, threaded)
        Base.Threads.@threads :static for i = 1:size_2
            @inbounds for j = 1:size_2
                sum_val = 0.0 + 0.0im
                for s = 1:sample_size
                    sum_val += O_store[i, s] * conj(O_store[j, s])
                end
                sr_opt_oo[(i-1)*size_2+j] = sum_val
            end
        end
    else
        @inbounds for i = 1:size_2
            for j = 1:size_2
                sum_val = 0.0 + 0.0im
                for s = 1:sample_size
                    sum_val += O_store[i, s] * conj(O_store[j, s])
                end
                sr_opt_oo[(i-1)*size_2+j] = sum_val
            end
        end
    end
end

function vmc_sample_chunks(n_samples::Integer, n_chunks::Integer)
    n = Int(n_samples)
    c = Int(n_chunks)
    c >= 1 || error("n_chunks must be >= 1, got $c")
    n >= 0 || error("n_samples must be >= 0, got $n")
    chunk_size = div(n, c)
    remainder = rem(n, c)
    chunks = Vector{UnitRange{Int}}(undef, c)
    first_sample = 0
    for chunk = 1:c
        len = chunk_size + (chunk <= remainder ? 1 : 0)
        last_sample = first_sample + len - 1
        chunks[chunk] = first_sample:last_sample
        first_sample += len
    end
    return chunks
end

function make_vmc_main_cal_worker_state(
    data::ExpertModeData,
    parent::VMCOptimizationState,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_proj = MVMCExpertModeParsers.projection_layout(data).n_proj
    n_para = parent.sr_opt.sr_opt_size - 1
    n_qp_full = length(parent.slater_matrix.pf_m)
    n_vmc_sample = data.modpara.nvmc_sample
    all_complex = get_all_complex_flag(data)
    use_fsz = !isempty(parent.electron_config.ele_spn)

    worker = VMCOptimizationState(
        n_site,
        n_elec,
        n_proj,
        n_para,
        n_qp_full,
        n_vmc_sample,
        all_complex,
        use_fsz,
    )

    # Saved sample arrays are read-only in VMCMainCal; avoid duplicating them.
    # Scratch arrays and counters remain worker-local.
    worker.electron_config.ele_idx = parent.electron_config.ele_idx
    worker.electron_config.ele_cfg = parent.electron_config.ele_cfg
    worker.electron_config.ele_num = parent.electron_config.ele_num
    worker.electron_config.ele_proj_cnt = parent.electron_config.ele_proj_cnt
    worker.electron_config.ele_spn = parent.electron_config.ele_spn
    worker.phys_quantities = parent.phys_quantities
    worker.slater_matrix.slater_elm = parent.slater_matrix.slater_elm
    worker.slater_matrix.slater_elm_real = parent.slater_matrix.slater_elm_real
    copyto!(worker.slater_matrix.inv_m, parent.slater_matrix.inv_m)
    copyto!(worker.slater_matrix.pf_m, parent.slater_matrix.pf_m)
    if !isempty(worker.slater_matrix.inv_m_real)
        copyto!(worker.slater_matrix.inv_m_real, parent.slater_matrix.inv_m_real)
        copyto!(worker.slater_matrix.pf_m_real, parent.slater_matrix.pf_m_real)
    end
    return worker
end

# ============================================================================
# Main Calculation Functions
# ============================================================================

"""
    vmc_main_cal!(data::ExpertModeData, state::VMCOptimizationState)

Main calculation: energy expectation values and SR optimization quantities (sz-conserved).
Equivalent to C's `VMCMainCal()`.
"""
function vmc_main_cal!(
    data::ExpertModeData,
    state::VMCOptimizationState,
    c_timer::CTimer = CTIMER_DISABLED;
    requested_threads::Integer = Base.Threads.nthreads(),
    use_store::Bool = true,
)
    n_site = data.modpara.nsite
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_site2 = 2 * n_site
    n_qp_full = get_n_qp_full(data)
    n_vmc_sample = data.modpara.nvmc_sample
    n_proj = MVMCExpertModeParsers.projection_layout(data).n_proj
    n_rbm = has_rbm_terms(data) ? MVMCExpertModeParsers.count_rbm_parameters(data) : 0
    n_opt_trans = MVMCExpertModeParsers.count_opt_trans_parameters(data)

    # IMPORTANT: n_slater should be NOrbitalIdx (number of unique orbital parameters),
    # NOT length(orbital_terms) (total number of site pairs)
    # C implementation: NSlater = NOrbitalIdx
    n_slater = if data.modpara.n_orbital_idx > 0
        data.modpara.n_orbital_idx
    elseif !isempty(data.orbital_terms)
        maximum(t.idx for t in data.orbital_terms) + 1
    else
        0
    end
    sr_opt_size = state.sr_opt.sr_opt_size

    # Get all_complex flag from data (not from empty arrays, as both are now always allocated)
    all_complex = get_all_complex_flag(data)

    # NVMCCalMode: 0 = optimization, 1 = measurement
    # We only implement optimization mode (0)
    nvmc_cal_mode = data.modpara.vmc_calc_mode

    # Clear physical quantities ([24] cal)
    ctimer_start!(c_timer, 24)
    clear_phys_quantity!(state)
    ctimer_stop!(c_timer, 24)

    # Debug: Check slater_elm at start of vmc_main_cal!
    slater_elm_max = maximum(abs.(state.slater_matrix.slater_elm))
    @debug "vmc_main_cal!: slater_elm max=$(slater_elm_max), size=$(length(state.slater_matrix.slater_elm))"
    if slater_elm_max < 1e-14
        @error "vmc_main_cal!: slater_elm is all zeros! This indicates update_slater_elm_fcmp! was not called or failed."
    end

    function process_sample_range!(
        sample_range,
        worker_state::VMCOptimizationState,
        local_acc::VMCThreadAccumulator,
        c_timer::CTimer,
        allow_inner_threads::Bool,
    )
        @debug "VMCMainCal worker sample range" sample_range length_ele_idx=length(worker_state.electron_config.ele_idx) length_ele_cfg=length(worker_state.electron_config.ele_cfg) length_ele_num=length(worker_state.electron_config.ele_num) length_ele_proj_cnt=length(worker_state.electron_config.ele_proj_cnt)

        # Sample loop
        for sample_value in sample_range
            sample = Int(sample_value)
            sample_n_site2 = n_site2
            # Get electron configuration for this sample
            ele_idx_start = sample * n_size + 1
            ele_idx_end = ele_idx_start + n_size - 1
            if ele_idx_end > length(worker_state.electron_config.ele_idx)
                @error "VMCMainCal: sample ele_idx slice is out of bounds" sample ele_idx_start ele_idx_end saved_ele_idx_len=length(worker_state.electron_config.ele_idx)
                continue
            end
            ele_idx = Vector{Int}(undef, n_size)
            copyto!(ele_idx, 1, worker_state.electron_config.ele_idx, ele_idx_start, n_size)

            # Validate ele_idx: check if it's all zeros (uninitialized)
            if all(x -> x == 0, ele_idx)
                @error "VMCMainCal: sample=$sample has invalid ele_idx (all zeros). This sample was not properly saved by vmc_make_sample!. Skipping this sample."
                @error "  ele_idx_start=$ele_idx_start, ele_idx_end=$ele_idx_end, n_size=$n_size"
                @error "  ele_idx[1:min(10, length(ele_idx))] = $(ele_idx[1:min(10, length(ele_idx))])"
                continue
            end

            # Validate ele_idx: check if it contains invalid values (all -1 or negative)
            if all(x -> x < 0, ele_idx)
                @error "VMCMainCal: sample=$sample has invalid ele_idx (all negative). This sample was not properly initialized. Skipping this sample."
                @error "  ele_idx[1:min(10, length(ele_idx))] = $(ele_idx[1:min(10, length(ele_idx))])"
                continue
            end

            ele_cfg_start = sample * sample_n_site2 + 1
            ele_cfg_end = ele_cfg_start + sample_n_site2 - 1
            if ele_cfg_end > length(worker_state.electron_config.ele_cfg)
                @error "VMCMainCal: sample ele_cfg slice is out of bounds" sample ele_cfg_start ele_cfg_end n_site2=sample_n_site2 saved_ele_cfg_len=length(worker_state.electron_config.ele_cfg)
                continue
            end
            ele_cfg = Vector{Int}(undef, sample_n_site2)
            copyto!(ele_cfg, 1, worker_state.electron_config.ele_cfg, ele_cfg_start, sample_n_site2)

            ele_num_start = sample * sample_n_site2 + 1
            ele_num_end = ele_num_start + sample_n_site2 - 1
            if ele_num_end > length(worker_state.electron_config.ele_num)
                @error "VMCMainCal: sample ele_num slice is out of bounds" sample ele_num_start ele_num_end n_site2=sample_n_site2 saved_ele_num_len=length(worker_state.electron_config.ele_num)
                continue
            end
            ele_num = Vector{Int}(undef, sample_n_site2)
            copyto!(ele_num, 1, worker_state.electron_config.ele_num, ele_num_start, sample_n_site2)

            # Debug: Validate ele_num for half-filling on localized spin sites only
            # Note: Half-filling check (n0[ri] + n1[ri] == 1) only applies to
            # sites marked as localized spins (LocSpn[ri] == 1).
            # For itinerant sites (LocSpn[ri] == 0), sites can be empty or doubly occupied.
            n0 = @view ele_num[1:n_site]
            n1 = @view ele_num[(n_site+1):(2*n_site)]
            loc_spn = get_loc_spn_array(data)
            half_filling_violations = 0
            for ri = 1:n_site
                if loc_spn[ri] == 1 && n0[ri] + n1[ri] != 1
                    # Only check localized spin sites
                    half_filling_violations += 1
                end
            end
            if half_filling_violations > 0
                @warn "VMCMainCal: sample=$sample violates half-filling at $half_filling_violations localized spin sites"
                @warn "  n0[1:min(8,$n_site)] = $(n0[1:min(8, n_site)])"
                @warn "  n1[1:min(8,$n_site)] = $(n1[1:min(8, n_site)])"
                @warn "  loc_spn[1:min(8,$n_site)] = $(loc_spn[1:min(8, n_site)])"
                @warn "  ele_idx = $(ele_idx[1:min(10, length(ele_idx))])"
            end

            ele_proj_cnt_start = sample * n_proj + 1
            ele_proj_cnt_end =
                min(ele_proj_cnt_start + n_proj - 1, length(worker_state.electron_config.ele_proj_cnt))
            if n_proj > 0 && ele_proj_cnt_end >= ele_proj_cnt_start
                ele_proj_cnt =
                    worker_state.electron_config.ele_proj_cnt[ele_proj_cnt_start:ele_proj_cnt_end]
            else
                ele_proj_cnt = Int[]
            end

            # Calculate M inverse and Pfaffian for this sample ([40] CalculateMAll;
            # wraps the call site and the real->complex copy, as in C's vmccal.c).
            # Use real version if !all_complex, matching C implementation
            ctimer_start!(c_timer, 40)
            info = 0
            if !all_complex
                if allow_inner_threads
                    info = calculate_m_all_real!(
                        ele_idx,
                        1,
                        n_qp_full + 1,
                        data,
                        worker_state;
                        threaded = true,
                    )
                else
                    info = calculate_m_all_real!(
                        ele_idx,
                        1,
                        n_qp_full + 1,
                        data,
                        worker_state;
                        threaded = false,
                    )
                end
                # Copy inv_m_real to inv_m (needed for SlaterElmDiff_fcmp)
                # C does: for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++) InvM[tmp_i]=InvM_real[tmp_i];
                # which copies both InvM_real and PfM_real (since PfM = InvM + NQPFull*Nsize*Nsize)
                n_size_sq = n_size * n_size
                copy_real_to_complex!(
                    worker_state.slater_matrix.inv_m,
                    worker_state.slater_matrix.inv_m_real,
                    n_qp_full * n_size_sq;
                    threaded = allow_inner_threads,
                )
                # Also copy pf_m_real to pf_m (needed for SlaterElmDiff_fcmp)
                # In C, PfM and InvM are contiguous, so the above copy covers both
                # In Julia, they are separate arrays, so we need explicit copy
                copy_real_to_complex!(
                    worker_state.slater_matrix.pf_m,
                    worker_state.slater_matrix.pf_m_real,
                    n_qp_full;
                    threaded = allow_inner_threads,
                )
            else
                if allow_inner_threads
                    info = calculate_m_all_fcmp!(
                        ele_idx,
                        1,
                        n_qp_full + 1,
                        data,
                        worker_state;
                        threaded = true,
                    )
                else
                    info = calculate_m_all_fcmp!(
                        ele_idx,
                        1,
                        n_qp_full + 1,
                        data,
                        worker_state;
                        threaded = false,
                    )
                end
            end
            ctimer_stop!(c_timer, 40)

            if info != 0
                @warn "VMCMainCal: sample=$sample info=$info (CalculateMAll)"
                continue
            end


            # Debug: Check slater_elm values
            if sample == 0 || sample == 1
                # Check up-spin part (rsi=0-15): slater_elm[rsi * n_site2 + rsj + 1] for rsi=0-15
                # Check down-spin part (rsi=16-31): slater_elm[rsi * n_site2 + rsj + 1] for rsi=16-31
                # For qp_idx=1 (first QP), offset = 0
                qp_offset = 0 * sample_n_site2 * sample_n_site2
                @debug "VMCMainCal: sample=$sample, slater_elm size=$(length(worker_state.slater_matrix.slater_elm)), n_qp_full=$(n_qp_full)"
                # Check all 4 blocks of slater_elm
                @debug "VMCMainCal: sample=$sample, slater_elm blocks:"
                @debug "  up-up [rsi=0, rsj=0-4]: $(worker_state.slater_matrix.slater_elm[qp_offset + 1:qp_offset + min(5, sample_n_site2)])"
                @debug "  up-down [rsi=0, rsj=16-20]: $(worker_state.slater_matrix.slater_elm[qp_offset + 16 + 1:qp_offset + min(21, sample_n_site2)])"
                @debug "  down-up [rsi=16, rsj=0-4]: $(worker_state.slater_matrix.slater_elm[qp_offset + 16 * sample_n_site2 + 1:qp_offset + 16 * sample_n_site2 + min(5, sample_n_site2)])"
                @debug "  down-down [rsi=16, rsj=16-20]: $(worker_state.slater_matrix.slater_elm[qp_offset + 16 * sample_n_site2 + 16 + 1:qp_offset + 16 * sample_n_site2 + min(21, sample_n_site2)])"
                @debug "  down-down [rsi=17, rsj=16-20]: $(worker_state.slater_matrix.slater_elm[qp_offset + 17 * sample_n_site2 + 16 + 1:qp_offset + 17 * sample_n_site2 + min(21, sample_n_site2)])"
            end

            # Calculate inner product (use real version if !all_complex, matching C implementation)
            if !all_complex
                ip_real =
                    calculate_ip_real(worker_state.slater_matrix.pf_m_real, 1, n_qp_full + 1, data)
                ip = ComplexF64(ip_real, 0.0)  # Convert to ComplexF64 for calculate_hamiltonian
            else
                ip = calculate_ip_fcmp(worker_state.slater_matrix.pf_m, 1, n_qp_full + 1, data)
            end

            # Debug: Check for zero or extreme ip values
            if abs(ip) < 1e-100
                @warn "VMCMainCal: sample=$sample has zero ip: ip=$ip, pf_m[1]=$(!all_complex ? worker_state.slater_matrix.pf_m_real[1] : worker_state.slater_matrix.pf_m[1])" maxlog=5
            elseif abs(ip) > 1e100
                @warn "VMCMainCal: sample=$sample has extreme ip: ip=$ip, pf_m[1]=$(!all_complex ? worker_state.slater_matrix.pf_m_real[1] : worker_state.slater_matrix.pf_m[1])" maxlog=5
            end

            # Calculate weight (currently w = 1.0 for unbiased sampling)
            w = 1.0

            if !isfinite(w)
                @warn "VMCMainCal: sample=$sample w=$w"
                continue
            end

            # Calculate energy
            # C implementation uses CalculateHamiltonian_real(creal(ip), ...) for real mode
            # For real mode, pass real part of ip to calculate_hamiltonian
            # Note: calculate_hamiltonian uses green_func1/green_func2 which need to handle real mode
            # [41] LocEnergyCal
            ctimer_start!(c_timer, 41)
            if !all_complex
                # Real mode: use real ip and pass all_complex flag
                e = calculate_hamiltonian(
                    ComplexF64(ip_real, 0.0),
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    data,
                    worker_state;
                    all_complex = false,
                    c_timer = c_timer,
                )
            else
                e = calculate_hamiltonian(
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                    data,
                    worker_state;
                    all_complex = true,
                    c_timer = c_timer,
                )
            end
            ctimer_stop!(c_timer, 41)

            if !isfinite(real(e) + imag(e))
                @warn "VMCMainCal: sample=$sample e=$e ip=$ip"
                continue
            end

            # Debug: Track zero energy cases
            if abs(real(e)) < 1e-15 && sample < 5
                @debug "VMCMainCal: Zero energy at sample=$sample, ip=$ip, pf_m[1]=$(worker_state.slater_matrix.pf_m[1]), ele_idx[1:8]=$(ele_idx[1:min(8, length(ele_idx))])"
            end

            # Debug: Check energy values
            if sample == 0 || sample == 1
                @debug "VMCMainCal: sample=$sample, e=$(real(e)), ip=$ip, wc_before=$(local_acc.energy.wc)"
            end

            # Accumulate energy statistics
            accumulate_energy!(local_acc.energy, w, e)

            # Calculate Green's functions (only in measurement mode)
            if nvmc_cal_mode == 1 && worker_state.phys_quantities !== nothing
                calculate_green_func!(
                    data,
                    worker_state,
                    local_acc.phys,
                    w,
                    ip,
                    ele_idx,
                    ele_cfg,
                    ele_num,
                    ele_proj_cnt,
                )
            end

            # Calculate SR optimization quantities (only in optimization mode)
            if nvmc_cal_mode == 0
                sr_opt_o = local_acc.sr_opt.sr_opt_o
                fill!(sr_opt_o, 0.0 + 0.0im)

                # Set projection factor derivatives
                set_projection_diff!(sr_opt_o, ele_proj_cnt, n_proj)

                # Set RBM derivatives
                if n_rbm > 0
                    rbm_offset = 2 * n_proj + 2
                    rbm_end = rbm_offset + 2 * n_rbm
                    if rbm_end <= length(sr_opt_o)
                        rbm_cnt = make_rbm_cnt(ele_num, data)
                        rbm_view = @view sr_opt_o[(rbm_offset+1):rbm_end]
                        set_rbm_diff!(rbm_view, rbm_cnt, ele_num, data)
                    else
                        @warn "RBM derivative range out of bounds: rbm_end=$rbm_end, length(sr_opt_o)=$(length(sr_opt_o))" maxlog=1
                    end
                end

                # Set Slater element derivatives ([42] ReturnSlaterElmDiff)
                if n_slater > 0
                    slater_offset = 2 * (n_proj + n_rbm) + 2  # After projection + RBM factors
                    slater_end = slater_offset + 2 * n_slater
                    slater_view = @view sr_opt_o[(slater_offset+1):slater_end]
                    ctimer_start!(c_timer, 42)
                    slater_elm_diff_fcmp!(slater_view, ip, ele_idx, data, worker_state)
                    ctimer_stop!(c_timer, 42)
                end
                if n_opt_trans > 0
                    opt_offset = 2 * (n_proj + n_rbm + n_slater) + 2
                    opt_end = opt_offset + 2 * n_opt_trans
                    if opt_end <= length(sr_opt_o)
                        opt_trans_diff!(@view(sr_opt_o[(opt_offset+1):opt_end]), ip, data, worker_state)
                    end
                end

                # Accumulate OO and HO ([43] calculate OO and HO)
                ctimer_start!(c_timer, 43)
                if all_complex
                    if use_store
                        calculate_oo_store!(
                            local_acc.sr_opt.sr_opt_oo,
                            local_acc.sr_opt.sr_opt_ho,
                            local_acc.sr_opt.sr_opt_o_store,
                            sr_opt_o,
                            w,
                            e,
                            sample,
                            sr_opt_size,
                            threaded = allow_inner_threads,
                        )
                    else
                        calculate_oo!(
                            local_acc.sr_opt.sr_opt_oo,
                            local_acc.sr_opt.sr_opt_ho,
                            sr_opt_o,
                            w,
                            e,
                            sr_opt_size,
                            threaded = allow_inner_threads,
                        )
                    end
                else
                    # Convert to real and calculate
                    sr_opt_o_real = local_acc.sr_opt.sr_opt_o_real
                    for i = 1:sr_opt_size
                        sr_opt_o_real[i] = real(sr_opt_o[2*(i-1)+1])
                    end

                    calculate_oo_real!(
                        local_acc.sr_opt.sr_opt_oo_real,
                        local_acc.sr_opt.sr_opt_ho_real,
                        sr_opt_o_real,
                        w,
                        real(e),
                        sr_opt_size,
                        threaded = allow_inner_threads,
                    )
                end
                ctimer_stop!(c_timer, 43)
            end
        end

        return nothing
    end

    config = VMCThreadConfig(n_vmc_sample; requested_threads = requested_threads)
    if vmc_threading_enabled(config)
        local_accs = make_thread_accumulators(state, config, c_timer)
        worker_states = [
            make_vmc_main_cal_worker_state(data, state)
            for _ = 1:config.effective_threads
        ]
        chunks = vmc_sample_chunks(n_vmc_sample, config.effective_threads)

        Base.Threads.@threads :static for worker_id = 1:config.effective_threads
            process_sample_range!(
                chunks[worker_id],
                worker_states[worker_id],
                local_accs[worker_id],
                local_accs[worker_id].timer,
                false,
            )
        end

        # Finalize OO from stored samples ([45] multiply store OO)
        if nvmc_cal_mode == 0 && use_store && all_complex
            ctimer_start!(c_timer, 45)
            for local_acc in local_accs
                finalize_oo_store!(
                    local_acc.sr_opt.sr_opt_oo,
                    local_acc.sr_opt.sr_opt_o_store,
                    sr_opt_size,
                    n_vmc_sample,
                    threaded = true,
                )
            end
            ctimer_stop!(c_timer, 45)
        end

        clear_sropt_store!(state.sr_opt)
        merge_thread_accumulators!(state, c_timer, local_accs)
    else
        local_acc = VMCThreadAccumulator(state, c_timer)
        process_sample_range!(
            0:(n_vmc_sample-1),
            state,
            local_acc,
            c_timer,
            true,
        )

        # Finalize OO from stored samples ([45] multiply store OO)
        if nvmc_cal_mode == 0 && use_store && all_complex
            ctimer_start!(c_timer, 45)
            finalize_oo_store!(
                local_acc.sr_opt.sr_opt_oo,
                local_acc.sr_opt.sr_opt_o_store,
                sr_opt_size,
                n_vmc_sample,
                threaded = true,
            )
            ctimer_stop!(c_timer, 45)
        end

        clear_sropt_store!(state.sr_opt)
        merge_thread_accumulator!(state, c_timer, local_acc)
    end

    # Debug: Check etot at the end of vmc_main_cal!
    if abs(state.energy.etot) < 1e-15 && real(state.energy.wc) > 0
        @warn "VMCMainCal: etot is zero at end of sampling. wc=$(state.energy.wc), n_vmc_sample=$n_vmc_sample"
        # Check the first sample's ele_num
        ele_num_start = 0 * n_site2 + 1
        ele_num_end = ele_num_start + n_site2 - 1
        sample_ele_num = state.electron_config.ele_num[ele_num_start:ele_num_end]
        n0 = sample_ele_num[1:n_site]
        n1 = sample_ele_num[(n_site+1):(2*n_site)]
        @warn "  Sample 0 ele_num: n0=$(n0[1:min(8,length(n0))]), n1=$(n1[1:min(8,length(n1))])"
        # Check exchange term conditions
        n_exchange_possible = 0
        for term in data.exchange_terms
            ri = term.site1
            rj = term.site2
            if 0 <= ri < n_site && 0 <= rj < n_site && ri != rj
                # Check condition for green_func2(ri, rj, rj, ri, 0, 1, ...)
                # Needs: n0[ri]==0, n0[rj]==1, n1[rj]==0, n1[ri]==1
                if n0[ri+1] == 0 && n0[rj+1] == 1 && n1[rj+1] == 0 && n1[ri+1] == 1
                    n_exchange_possible += 1
                end
            end
        end
        @warn "  Number of possible exchange terms (sample 0): $n_exchange_possible / $(length(data.exchange_terms))"
    end
end

"""
    calculate_sz_fsz(ele_num, n_site) -> Float64

Sample-local Sz = sum_i (n_up[i] - n_down[i]) for the fsz code path.
Mirrors C's `CalculateSz_fsz` in mVMC/src/mVMC/calham_fsz.c.
"""
function calculate_sz_fsz(ele_num::AbstractVector{<:Integer}, n_site::Integer)
    sz = 0
    @inbounds for ri = 1:n_site
        sz += ele_num[ri] - ele_num[ri+n_site]
    end
    return Float64(sz)
end

"""
    vmc_main_cal_fsz!(data::ExpertModeData, state::VMCOptimizationState)

Main calculation (fsz version) for Full Spin Z mode.
Uses FSZ-specific Green function and Pfaffian calculations where electrons
have individually tracked spins.

Equivalent to C's `VMCMainCal_fsz()`.

# Reference
- C implementation: mVMC/src/mVMC/vmccal_fsz.c:36-248 (VMCMainCal_fsz)
"""
function vmc_main_cal_fsz!(
    data::ExpertModeData,
    state::VMCOptimizationState,
    c_timer::CTimer = CTIMER_DISABLED;
    requested_threads::Integer = Base.Threads.nthreads(),
)
    n_site = data.modpara.nsite
    n_site2 = 2 * n_site
    n_elec = data.modpara.nelec
    n_size = 2 * n_elec
    n_vmc_sample = data.modpara.nvmc_sample
    n_qp_full = get_n_qp_full(data)
    nvmc_cal_mode = data.modpara.vmc_calc_mode
    n_proj = MVMCExpertModeParsers.projection_layout(data).n_proj
    n_slater = data.modpara.n_orbital_idx
    n_opt_trans = MVMCExpertModeParsers.count_opt_trans_parameters(data)

    # Determine if using all complex or real mode
    # FSZ typically uses complex mode
    all_complex = true
    orbital_complex = any(term -> term.is_complex, data.orbital_terms)
    if data.modpara.complex_flag != 0 || orbital_complex
        all_complex = true
    end

    # Clear all physical quantities (energy and SR optimization data)
    # This is equivalent to C's clearPhysQuantity() call in VMCMainCal_fsz ([24] cal)
    ctimer_start!(c_timer, 24)
    clear_phys_quantity!(state)
    ctimer_stop!(c_timer, 24)

    function process_sample_range_fsz!(
        sample_range,
        worker_state::VMCOptimizationState,
        local_acc::VMCThreadAccumulator,
        c_timer::CTimer,
        allow_inner_threads::Bool,
    )
        @debug "VMCMainCal_fsz worker sample range" sample_range length_ele_idx=length(worker_state.electron_config.ele_idx) length_ele_cfg=length(worker_state.electron_config.ele_cfg) length_ele_num=length(worker_state.electron_config.ele_num) length_ele_proj_cnt=length(worker_state.electron_config.ele_proj_cnt) length_ele_spn=length(worker_state.electron_config.ele_spn)

    # Process each sample
    for sample in sample_range
        # Extract electron configuration for this sample
        ele_idx_start = sample * n_size + 1
        ele_idx_end = ele_idx_start + n_size - 1
        ele_idx = worker_state.electron_config.ele_idx[ele_idx_start:ele_idx_end]

        # Validate ele_idx
        if all(x -> x == 0, ele_idx) || all(x -> x < 0, ele_idx)
            continue
        end

        ele_cfg_start = sample * n_site2 + 1
        ele_cfg_end = ele_cfg_start + n_site2 - 1
        ele_cfg = worker_state.electron_config.ele_cfg[ele_cfg_start:ele_cfg_end]

        ele_num_start = sample * n_site2 + 1
        ele_num_end = ele_num_start + n_site2 - 1
        ele_num = worker_state.electron_config.ele_num[ele_num_start:ele_num_end]

        ele_proj_cnt_start = sample * n_proj + 1
        ele_proj_cnt_end =
            min(ele_proj_cnt_start + n_proj - 1, length(worker_state.electron_config.ele_proj_cnt))
        if n_proj > 0 && ele_proj_cnt_end >= ele_proj_cnt_start
            ele_proj_cnt =
                worker_state.electron_config.ele_proj_cnt[ele_proj_cnt_start:ele_proj_cnt_end]
        else
            ele_proj_cnt = Int[]
        end

        # FSZ: Extract ele_spn for this sample
        ele_spn_start = sample * n_size + 1
        ele_spn_end = ele_spn_start + n_size - 1
        ele_spn = worker_state.electron_config.ele_spn[ele_spn_start:ele_spn_end]

        # Calculate M inverse and Pfaffian using FSZ version ([40] CalculateMAll)
        ctimer_start!(c_timer, 40)
        info = calculate_m_all_fsz!(ele_idx, ele_spn, 1, n_qp_full + 1, data, worker_state)
        ctimer_stop!(c_timer, 40)

        if info != 0
            @warn "VMCMainCal_fsz: sample=$sample info=$info (CalculateMAll_fsz)"
            continue
        end

        # Calculate inner product
        ip = calculate_ip_fcmp(worker_state.slater_matrix.pf_m, 1, n_qp_full + 1, data)

        if abs(ip) < 1e-100
            @warn "VMCMainCal_fsz: sample=$sample has zero ip" maxlog=5
        end

        # Weight (currently w = 1.0)
        w = 1.0

        if !isfinite(w)
            continue
        end

        # Calculate energy using FSZ version ([41] LocEnergyCal)
        ctimer_start!(c_timer, 41)
        e = calculate_hamiltonian_fsz(
            ip,
            ele_idx,
            ele_cfg,
            ele_num,
            ele_proj_cnt,
            ele_spn,
            data,
            worker_state;
            all_complex = all_complex,
            c_timer = c_timer,
        )
        ctimer_stop!(c_timer, 41)

        if !isfinite(real(e) + imag(e))
            @warn "VMCMainCal_fsz: sample=$sample e=$e ip=$ip"
            continue
        end

        # Accumulate energy statistics. For fsz mode the magnetisation is not
        # constrained, so we also accumulate <Sz> and <Sz^2> following the
        # C reference (mVMC/src/mVMC/calham_fsz.c:CalculateSz_fsz +
        # vmccal_fsz.c:144-148).
        sz_local = calculate_sz_fsz(ele_num, n_site)
        accumulate_energy!(local_acc.energy, w, e; sz = sz_local)

        # Calculate Green's functions (only in measurement mode)
        # C implementation: VMCMainCal_fsz calls CalculateGreenFunc_fsz when NVMCCalMode==1
        if nvmc_cal_mode == 1 && worker_state.phys_quantities !== nothing
            # TODO: Implement FSZ-specific green function calculation
            # For now, use standard version (may not be accurate for spin-flip terms)
            calculate_green_func!(
                data,
                worker_state,
                local_acc.phys,
                w,
                ip,
                ele_idx,
                ele_cfg,
                ele_num,
                ele_proj_cnt,
            )
        end

        # Calculate SR optimization quantities (only in optimization mode)
        if nvmc_cal_mode == 0
            sr_opt_o = local_acc.sr_opt.sr_opt_o
            fill!(sr_opt_o, 0.0 + 0.0im)

            # Set projection factor derivatives
            set_projection_diff!(sr_opt_o, ele_proj_cnt, n_proj)

            # Set Slater element derivatives (using FSZ version if available)
            # TODO: Implement slater_elm_diff_fsz for full FSZ support
            if n_slater > 0
                slater_offset = 2 * n_proj + 2
                slater_end = slater_offset + 2 * n_slater
                slater_view = @view sr_opt_o[(slater_offset+1):slater_end]
                # Use FSZ version for Slater element derivative ([42] ReturnSlaterElmDiff)
                ctimer_start!(c_timer, 42)
                slater_elm_diff_fsz!(slater_view, ip, ele_idx, ele_spn, data, worker_state)
                ctimer_stop!(c_timer, 42)
            end
            if n_opt_trans > 0
                opt_offset = 2 * (n_proj + n_slater) + 2
                opt_end = opt_offset + 2 * n_opt_trans
                if opt_end <= length(sr_opt_o)
                    opt_trans_diff!(@view(sr_opt_o[(opt_offset+1):opt_end]), ip, data, worker_state)
                end
            end

            # Accumulate OO and HO ([43] calculate OO and HO)
            ctimer_start!(c_timer, 43)
            calculate_oo!(
                local_acc.sr_opt.sr_opt_oo,
                local_acc.sr_opt.sr_opt_ho,
                sr_opt_o,
                w,
                e,
                length(sr_opt_o) ÷ 2,
                threaded = allow_inner_threads,
            )
            ctimer_stop!(c_timer, 43)
        end
    end

        return nothing
    end

    config = VMCThreadConfig(n_vmc_sample; requested_threads = requested_threads)
    if vmc_threading_enabled(config)
        local_accs = make_thread_accumulators(state, config, c_timer)
        worker_states = [
            make_vmc_main_cal_worker_state(data, state)
            for _ = 1:config.effective_threads
        ]
        chunks = vmc_sample_chunks(n_vmc_sample, config.effective_threads)

        Base.Threads.@threads :static for worker_id = 1:config.effective_threads
            process_sample_range_fsz!(
                chunks[worker_id],
                worker_states[worker_id],
                local_accs[worker_id],
                local_accs[worker_id].timer,
                false,
            )
        end

        clear_sropt_store!(state.sr_opt)
        merge_thread_accumulators!(state, c_timer, local_accs)
    else
        local_acc = VMCThreadAccumulator(state, c_timer)
        process_sample_range_fsz!(0:(n_vmc_sample-1), state, local_acc, c_timer, true)

        clear_sropt_store!(state.sr_opt)
        merge_thread_accumulator!(state, c_timer, local_acc)
    end

    # Debug: Check final energy accumulation
    if abs(state.energy.etot) < 1e-15 && real(state.energy.wc) > 0
        @warn "VMCMainCal_fsz: etot is zero at end of sampling. wc=$(state.energy.wc), n_vmc_sample=$n_vmc_sample"
    end
end

"""
    vmc_bf_main_cal!(data::ExpertModeData, state::VMCOptimizationState)

Main calculation with Back Flow.
Equivalent to C's `VMC_BF_MainCal()`.

BackFlow is not supported in Julia-mVMC v0.1: calling this function raises
an error. Inputs that drive `n_proj_bf > 0` should not be passed to v0.1.
"""
function vmc_bf_main_cal!(::ExpertModeData, ::VMCOptimizationState)
    error("BackFlow is not supported in Julia-mVMC v0.1. Remove BackFlow keywords from namelist.def, " *
          "or fall back to the C reference at https://github.com/issp-center-dev/mVMC.")
end
