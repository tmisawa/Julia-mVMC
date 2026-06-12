"""
Stochastic Optimization Functions

Solve S*x = g and update variational parameters using Stochastic Reconfiguration.
"""

"""
    get_opt_flag_for_parameter(data::ExpertModeData, pi::Int) -> Int

Get OptFlag value for a parameter index.
C: OptFlag[pi] (0-based)
Julia: optimization_flags[pi+1] (1-based)

# Arguments
- `data::ExpertModeData`: Expert Mode data
- `pi::Int`: Parameter index (0-based, 0 to 2*NPara-1)

# Returns
- `Int`: 1 if optimized, 0 if fixed
"""
function get_opt_flag_for_parameter(data::ExpertModeData, pi::Int)::Int
    if pi + 1 > length(data.optimization_flags)
        return 0
    end
    return data.optimization_flags[pi+1] ? 1 : 0
end

function _rbm_parameter_sections(data::ExpertModeData)
    return (
        data.charge_rbm_phys_layer_terms,
        data.spin_rbm_phys_layer_terms,
        data.general_rbm_phys_layer_terms,
        data.charge_rbm_hidden_layer_terms,
        data.spin_rbm_hidden_layer_terms,
        data.general_rbm_hidden_layer_terms,
        data.charge_rbm_phys_hidden_terms,
        data.spin_rbm_phys_hidden_terms,
        data.general_rbm_phys_hidden_terms,
    )
end

@inline _parameter_section_width(terms) =
    isempty(terms) ? 0 : maximum(t.idx for t in terms) + 1

function _parameter_count_breakdown(data::ExpertModeData)
    layout = MVMCExpertModeParsers.projection_layout(data)
    n_rbm = MVMCExpertModeParsers.count_rbm_parameters(data)
    n_orbital_idx = MVMCExpertModeParsers.count_orbital_parameters(data)
    n_opt_trans = MVMCExpertModeParsers.count_opt_trans_parameters(data)
    n_para = MVMCExpertModeParsers.count_variational_parameters(data)
    return (;
        layout,
        n_proj = layout.n_proj,
        n_rbm,
        n_orbital_idx,
        n_opt_trans,
        n_para,
    )
end

@enum _ParameterKind begin
    _PARAM_GUTZWILLER
    _PARAM_JASTROW
    _PARAM_DH2
    _PARAM_DH4
    _PARAM_RBM
    _PARAM_ORBITAL
    _PARAM_OPTTRANS
end

struct _ParameterLocation
    kind::_ParameterKind
    section_idx::Int
    item_idx::Int
end

function _foreach_parameter_location(f, data::ExpertModeData, counts)
    layout = counts.layout
    n_para = counts.n_para

    for i = 1:min(layout.n_gutzwiller, length(data.gutzwiller_terms))
        para_idx = layout.gutzwiller_offset + i
        1 <= para_idx <= n_para &&
            f(para_idx, _ParameterLocation(_PARAM_GUTZWILLER, 0, i))
    end
    for i = 1:min(layout.n_jastrow, length(data.jastrow_terms))
        para_idx = layout.jastrow_offset + i
        1 <= para_idx <= n_para &&
            f(para_idx, _ParameterLocation(_PARAM_JASTROW, 0, i))
    end
    for i = 1:min(6 * layout.n_dh2, length(data.doublon_holon_2site_params))
        para_idx = layout.dh2_offset + i
        1 <= para_idx <= n_para &&
            f(para_idx, _ParameterLocation(_PARAM_DH2, 0, i))
    end
    for i = 1:min(10 * layout.n_dh4, length(data.doublon_holon_4site_params))
        para_idx = layout.dh4_offset + i
        1 <= para_idx <= n_para &&
            f(para_idx, _ParameterLocation(_PARAM_DH4, 0, i))
    end

    rbm_offset = counts.n_proj
    for (section_idx, terms) in enumerate(_rbm_parameter_sections(data))
        n_section = _parameter_section_width(terms)
        for term_idx in eachindex(terms)
            para_idx = rbm_offset + terms[term_idx].idx + 1
            1 <= para_idx <= n_para &&
                f(para_idx, _ParameterLocation(_PARAM_RBM, section_idx, term_idx))
        end
        rbm_offset += n_section
    end

    slater_start = counts.n_proj + counts.n_rbm + 1
    for term_idx in eachindex(data.orbital_terms)
        para_idx = slater_start + data.orbital_terms[term_idx].idx
        slater_start <= para_idx < slater_start + counts.n_orbital_idx &&
            f(para_idx, _ParameterLocation(_PARAM_ORBITAL, 0, term_idx))
    end

    opt_trans_start = counts.n_proj + counts.n_rbm + counts.n_orbital_idx + 1
    for i = 1:min(counts.n_opt_trans, length(data.opt_trans))
        para_idx = opt_trans_start + i - 1
        1 <= para_idx <= n_para &&
            f(para_idx, _ParameterLocation(_PARAM_OPTTRANS, 0, i))
    end
    return nothing
end

function _foreach_parameter_location_at(f, data::ExpertModeData, counts, para_idx::Int)
    1 <= para_idx <= counts.n_para || return nothing

    layout = counts.layout
    if para_idx <= counts.n_proj
        if para_idx <= layout.gutzwiller_offset + layout.n_gutzwiller
            local_idx = para_idx - layout.gutzwiller_offset
            1 <= local_idx <= length(data.gutzwiller_terms) &&
                f(_ParameterLocation(_PARAM_GUTZWILLER, 0, local_idx))
        elseif para_idx <= layout.jastrow_offset + layout.n_jastrow
            local_idx = para_idx - layout.jastrow_offset
            1 <= local_idx <= length(data.jastrow_terms) &&
                f(_ParameterLocation(_PARAM_JASTROW, 0, local_idx))
        elseif layout.dh2_offset < para_idx <= layout.dh2_offset + 6 * layout.n_dh2
            local_idx = para_idx - layout.dh2_offset
            1 <= local_idx <= length(data.doublon_holon_2site_params) &&
                f(_ParameterLocation(_PARAM_DH2, 0, local_idx))
        elseif layout.dh4_offset < para_idx <= layout.dh4_offset + 10 * layout.n_dh4
            local_idx = para_idx - layout.dh4_offset
            1 <= local_idx <= length(data.doublon_holon_4site_params) &&
                f(_ParameterLocation(_PARAM_DH4, 0, local_idx))
        end
        return nothing
    elseif para_idx <= counts.n_proj + counts.n_rbm
        rbm_idx = para_idx - counts.n_proj - 1
        section_offset = 0
        for (section_idx, terms) in enumerate(_rbm_parameter_sections(data))
            n_section = _parameter_section_width(terms)
            if rbm_idx < section_offset + n_section
                local_idx = rbm_idx - section_offset
                for term_idx in eachindex(terms)
                    terms[term_idx].idx == local_idx &&
                        f(_ParameterLocation(_PARAM_RBM, section_idx, term_idx))
                end
                return nothing
            end
            section_offset += n_section
        end
        return nothing
    end

    slater_start = counts.n_proj + counts.n_rbm + 1
    opt_trans_start = counts.n_proj + counts.n_rbm + counts.n_orbital_idx + 1
    if slater_start <= para_idx < opt_trans_start
        orbital_idx = para_idx - slater_start
        for term_idx in eachindex(data.orbital_terms)
            data.orbital_terms[term_idx].idx == orbital_idx &&
                f(_ParameterLocation(_PARAM_ORBITAL, 0, term_idx))
        end
    elseif opt_trans_start <= para_idx < opt_trans_start + counts.n_opt_trans
        opt_idx = para_idx - opt_trans_start + 1
        1 <= opt_idx <= length(data.opt_trans) &&
            f(_ParameterLocation(_PARAM_OPTTRANS, 0, opt_idx))
    end
    return nothing
end

function _parameter_location_value(data::ExpertModeData, loc::_ParameterLocation)::ComplexF64
    if loc.kind == _PARAM_GUTZWILLER
        return ComplexF64(data.gutzwiller_terms[loc.item_idx].value)
    elseif loc.kind == _PARAM_JASTROW
        return ComplexF64(data.jastrow_terms[loc.item_idx].value)
    elseif loc.kind == _PARAM_DH2
        return ComplexF64(data.doublon_holon_2site_params[loc.item_idx])
    elseif loc.kind == _PARAM_DH4
        return ComplexF64(data.doublon_holon_4site_params[loc.item_idx])
    elseif loc.kind == _PARAM_RBM
        return ComplexF64(_rbm_parameter_sections(data)[loc.section_idx][loc.item_idx].value)
    elseif loc.kind == _PARAM_ORBITAL
        return ComplexF64(data.orbital_terms[loc.item_idx].value)
    elseif loc.kind == _PARAM_OPTTRANS
        return ComplexF64(data.opt_trans[loc.item_idx])
    end
end

function _set_parameter_location_value!(
    data::ExpertModeData,
    loc::_ParameterLocation,
    value::ComplexF64,
)
    if loc.kind == _PARAM_GUTZWILLER
        data.gutzwiller_terms[loc.item_idx].value = value
    elseif loc.kind == _PARAM_JASTROW
        data.jastrow_terms[loc.item_idx].value = value
    elseif loc.kind == _PARAM_DH2
        data.doublon_holon_2site_params[loc.item_idx] = value
    elseif loc.kind == _PARAM_DH4
        data.doublon_holon_4site_params[loc.item_idx] = value
    elseif loc.kind == _PARAM_RBM
        _rbm_parameter_sections(data)[loc.section_idx][loc.item_idx].value = value
    elseif loc.kind == _PARAM_ORBITAL
        data.orbital_terms[loc.item_idx].value = value
    elseif loc.kind == _PARAM_OPTTRANS
        data.opt_trans[loc.item_idx] = value
    end
    return nothing
end

"""
    get_parameter_value(data, para_idx) -> ComplexF64

flat parameter index（1-based、Proj → RBM → Slater → OptTrans）の現在値を返す。
範囲外や layout の隙間（spin-jastrow 等）は 0 を返す（update 側が no-op になる
index と対応）。duplicate idx は `pack_parameters` と同じ last-wins 規約で返す。
書き込み側は `set_parameter_value!` が同じ layout を使う。
"""
function get_parameter_value(data::ExpertModeData, para_idx::Int)::ComplexF64
    counts = _parameter_count_breakdown(data)
    value = Ref(ComplexF64(0))
    _foreach_parameter_location_at(data, counts, para_idx) do loc
        value[] = _parameter_location_value(data, loc)
    end
    return value[]
end

function _set_parameter_value_direct!(data::ExpertModeData, para_idx::Int, value::ComplexF64)
    counts = _parameter_count_breakdown(data)
    touched_opt_trans = Ref(false)
    _foreach_parameter_location_at(data, counts, para_idx) do loc
        _set_parameter_location_value!(data, loc, value)
        if loc.kind == _PARAM_OPTTRANS
            touched_opt_trans[] = true
        end
    end
    if touched_opt_trans[] && data.qp_weights !== nothing
        MVMCExpertModeParsers.update_qp_weight!(data.qp_weights, data.opt_trans)
    end
    return nothing
end

function _add_parameter_delta_direct!(data::ExpertModeData, para_idx::Int, delta::ComplexF64)
    counts = _parameter_count_breakdown(data)
    touched_opt_trans = Ref(false)
    _foreach_parameter_location_at(data, counts, para_idx) do loc
        _set_parameter_location_value!(data, loc, _parameter_location_value(data, loc) + delta)
        if loc.kind == _PARAM_OPTTRANS
            touched_opt_trans[] = true
        end
    end
    if touched_opt_trans[] && data.qp_weights !== nothing
        MVMCExpertModeParsers.update_qp_weight!(data.qp_weights, data.opt_trans)
    end
    return nothing
end

"""
    update_parameter_value(data::ExpertModeData, para_idx::Int, delta_real::Float64, delta_imag::Float64)

Update a variational parameter value.
para_idx is the real parameter index (0-based in C, but we use 1-based for Julia arrays).

# Arguments
- `data::ExpertModeData`: Expert Mode data
- `para_idx::Int`: Real parameter index (1-based: 1 to NPara)
- `delta_real::Float64`: Real part of delta
- `delta_imag::Float64`: Imaginary part of delta
"""
function update_parameter_value(
    data::ExpertModeData,
    para_idx::Int,
    delta_real::Float64,
    delta_imag::Float64,
)
    delta = ComplexF64(delta_real, delta_imag)
    _add_parameter_delta_direct!(data, para_idx, delta)
    return nothing
end

"""
    build_s_matrix_and_g_vector!(
        S::Matrix{Float64},
        g::Vector{Float64},
        smat_to_para_idx::Vector{Int},
        sr_opt_oo::Vector{ComplexF64},
        sr_opt_ho::Vector{ComplexF64},
        sr_opt_size::Int,
        dsr_opt_sta_del::Float64,
        dsr_opt_step_dt::Float64
    )

Build S matrix and g vector for stochastic reconfiguration.
Equivalent to C's `stcOptInit()`.

# Arguments
- `S::Matrix{Float64}`: Output S matrix (n_smat x n_smat, column-major)
- `g::Vector{Float64}`: Output g vector (n_smat)
- `smat_to_para_idx::Vector{Int}`: Mapping from S matrix index to parameter index (0-based)
- `sr_opt_oo::Vector{ComplexF64}`: Overlap matrix OO
- `sr_opt_ho::Vector{ComplexF64}`: Energy gradient HO
- `sr_opt_size::Int`: Size of SR optimization arrays
- `dsr_opt_sta_del::Float64`: Diagonal stabilization factor
- `dsr_opt_step_dt::Float64`: Step width
"""
function build_s_matrix_and_g_vector!(
    S::Matrix{Float64},
    g::Vector{Float64},
    smat_to_para_idx::Vector{Int},
    sr_opt_oo::Vector{ComplexF64},
    sr_opt_ho::Vector{ComplexF64},
    sr_opt_size::Int,
    dsr_opt_sta_del::Float64,
    dsr_opt_step_dt::Float64,
)
    n_smat = length(smat_to_para_idx)
    ratio_diag = 1.0 + dsr_opt_sta_del

    # Build S matrix: S[i][j] = OO[i+1][j+1] - OO[0][i+1] * OO[0][j+1]
    # C: SROptOO is stored with lda = 2*SROptSize (logical stride for rows)
    # C: calculateOO stores OO[i][j] at index i*(2*SROptSize) + j (0-based)
    # C: So OO[0][j] is at index j (0-based) = j+1 (1-based in Julia)
    # C: OO[pi+2][pj+2] is at index (pi+2)*(2*SROptSize) + (pj+2) (0-based)
    # C: SROptOO[pi+2] = OO[0][pi+1] (pi is 0-based parameter index)
    # NOTE: All C indices need +1 for Julia's 1-based indexing
    lda_oo = 2 * sr_opt_size  # Logical leading dimension used by C code
    for (si, pi) in enumerate(smat_to_para_idx)
        # C: tmp = creal(SROptOO[pi+2])
        # OO[0][pi+1] is at index pi+2 (0-based) = pi+3 (1-based in Julia)
        tmp = real(sr_opt_oo[pi+3])  # OO[0][pi+1], C: SROptOO[pi+2] -> Julia: [pi+2+1]

        for (sj, pj) in enumerate(smat_to_para_idx)
            # C: idx = si + nSmat*sj (column major)
            # Julia: column major, but 1-based indexing
            idx = si + n_smat * (sj - 1)
            # C: SROptOO[offset+(pj+2)] where offset = (pi+2)*(2*SROptSize)
            # OO[pi+1][pj+1] is at index (pi+2)*lda_oo + (pj+2) (0-based) = +1 for Julia
            oo_idx = (pi + 2) * lda_oo + (pj + 2) + 1  # 1-based
            # OO[0][pj+1] is at index pj+2 (0-based) = pj+3 (1-based)
            S[idx] = real(sr_opt_oo[oo_idx]) - tmp * real(sr_opt_oo[pj+3])
        end

        # Modify diagonal elements
        diag_idx = si + n_smat * (si - 1)
        S[diag_idx] *= ratio_diag
    end

    # Build g vector: g[i] = -DSROptStepDt * 2.0 * (HO[i+1] - HO[0] * OO[i+1])
    # C: SROptHO[0] = HO[0], SROptHO[pi+2] = HO[pi+1] (pi is 0-based)
    # Julia: sr_opt_ho[1] = HO[0], sr_opt_ho[pi+3] = HO[pi+1] (pi is 0-based, +1 for Julia 1-based)
    ho_0 = real(sr_opt_ho[1])  # HO[0]
    for (si, pi) in enumerate(smat_to_para_idx)
        ho_idx = pi + 3  # HO[pi+1], C: SROptHO[pi+2] -> Julia: sr_opt_ho[pi+2+1]
        oo_idx = pi + 3  # OO[0][pi+1], C: SROptOO[pi+2] -> Julia: sr_opt_oo[pi+2+1]
        g[si] =
            -dsr_opt_step_dt *
            2.0 *
            (real(sr_opt_ho[ho_idx]) - ho_0 * real(sr_opt_oo[oo_idx]))
    end
end

"""
    stochastic_opt!(data::ExpertModeData, state::VMCOptimizationState) -> Int

Stochastic optimization using LAPACK (direct solver).
Equivalent to C's `StochasticOpt()`.

Solves S*x = g where:
- S is the overlap matrix: S[i][j] = OO[i+1][j+1] - OO[i+1][0] * OO[0][j+1]
- g is the energy gradient: g[i] = 2.0 * (HO[i+1] - H * O[i+1])

# Returns
- `info::Int`: 0 = success, non-zero = error
"""
function stochastic_opt!(data::ExpertModeData, state::VMCOptimizationState, c_timer::CTimer = CTIMER_DISABLED)::Int
    # Get parameters
    n_para = _parameter_count_breakdown(data).n_para
    sr_opt_size = state.sr_opt.sr_opt_size
    sr_opt_oo = state.sr_opt.sr_opt_oo
    sr_opt_ho = state.sr_opt.sr_opt_ho
    all_complex = get_all_complex_flag(data)

    # Initialize optimization flags if empty (same as vmc_para_opt!)
    if isempty(data.optimization_flags)
        data.optimization_flags = fill(true, 2 * n_para)
    end

    # [50] preprocess: real->complex conversion, diag/cut, smat mapping
    ctimer_start!(c_timer, 50)

    # Phase 1: Handle real mode conversion (if needed)
    # In C: AllComplexFlag==0 case converts SROptOO_real to SROptOO and SROptHO_real to SROptHO
    # The C code loops through the entire contiguous memory (SROptOO + SROptHO + SROptO)

    if !all_complex && !isempty(state.sr_opt.sr_opt_oo_real)
        # Convert SROptOO_real to SROptOO
        # C: SROptOO_real is SROptSize*(SROptSize+2) in total size, but DGER uses SROptSize as stride
        # C: SROptOO[i] = SROptOO_real[j] for even indices, 0 for odd indices
        # C: j = int_x/2 + (int_y/2) * SROptSize (NOT SROptSize+2)
        for i = 1:length(sr_opt_oo)
            int_x = (i - 1) % (2 * sr_opt_size)  # Column index (0-based)
            int_y = div(i - 1, 2 * sr_opt_size)  # Row index (0-based)
            if int_x % 2 == 0 && int_y % 2 == 0
                # Convert to real array indices: j = int_x/2 + (int_y/2) * sr_opt_size
                j = div(int_x, 2) + div(int_y, 2) * sr_opt_size + 1  # 1-based
                if j <= length(state.sr_opt.sr_opt_oo_real)
                    sr_opt_oo[i] = ComplexF64(state.sr_opt.sr_opt_oo_real[j], 0.0)
                else
                    sr_opt_oo[i] = 0.0 + 0.0im
                end
            else
                sr_opt_oo[i] = 0.0 + 0.0im
            end
        end

        # Convert SROptHO_real to SROptHO
        # C: SROptHO[i] = SROptHO_real[i/2] for even i, 0 for odd i
        for i = 1:length(sr_opt_ho)
            if (i - 1) % 2 == 0  # even index (0-based)
                j = div(i - 1, 2) + 1  # real index (1-based)
                if j <= length(state.sr_opt.sr_opt_ho_real)
                    sr_opt_ho[i] = ComplexF64(state.sr_opt.sr_opt_ho_real[j], 0.0)
                else
                    sr_opt_ho[i] = 0.0 + 0.0im
                end
            else
                sr_opt_ho[i] = 0.0 + 0.0im
            end
        end
    end

    # Phase 2: Calculate diagonal elements and apply redundant direction cut
    # C: r[pi] = creal(srOptOO[(pi+2)*(2*srOptSize)+(pi+2)]) - creal(srOptOO[pi+2])^2
    # C: SROptOO[(pi+2)*(2*srOptSize)+(pi+2)] = OO[pi+1][pi+1]
    # C: SROptOO[pi+2] = OO[0][pi+1]
    s_diag = zeros(Float64, 2 * n_para)
    for pi = 0:(2*n_para-1)
        # C index: (pi+2)*(2*srOptSize) + (pi+2) (0-based)
        # Julia index: (pi+2) * (2*sr_opt_size) + (pi+2) + 1 (1-based)
        oo_idx_diag = (pi + 2) * (2 * sr_opt_size) + (pi + 2) + 1
        # C: SROptOO[pi+2] (0-based) -> Julia: sr_opt_oo[pi+2+1] (1-based)
        oo_idx_0 = pi + 2 + 1  # +1 for Julia's 1-based indexing

        if oo_idx_diag <= length(sr_opt_oo) && oo_idx_0 <= length(sr_opt_oo)
            s_diag[pi+1] = real(sr_opt_oo[oo_idx_diag]) - real(sr_opt_oo[oo_idx_0])^2
        else
            s_diag[pi+1] = 0.0
        end
    end

    # Find max and min
    s_diag_max = maximum(s_diag)
    s_diag_min = minimum(s_diag)

    # Calculate threshold
    diag_cut_threshold = s_diag_max * data.modpara.dsr_opt_red_cut

    # Build smat_to_para_idx mapping
    smat_to_para_idx = Int[]
    opt_num = 0  # Number of fixed parameters (by OptFlag)
    cut_num = 0  # Number of parameters cut by redundant direction

    for pi = 0:(2*n_para-1)
        opt_flag = get_opt_flag_for_parameter(data, pi)

        if opt_flag != 1  # Fixed parameter (by OptFlag)
            opt_num += 1
            continue
        end

        if s_diag[pi+1] < diag_cut_threshold  # Cut by redundant direction
            cut_num += 1
        else  # Optimize
            push!(smat_to_para_idx, pi)
        end
    end

    n_smat = length(smat_to_para_idx)
    ctimer_stop!(c_timer, 50)

    # If no parameters to optimize, return success
    if n_smat == 0
        @warn "No parameters to optimize (all fixed or cut)"
        return 0
    end

    # [51] stcOptMain: build S/g ([56]) + solve ([57])
    ctimer_start!(c_timer, 51)

    # Phase 3: Build S matrix and g vector ([56] calculate S and g)
    S = zeros(Float64, n_smat, n_smat)
    g = zeros(Float64, n_smat)

    ctimer_start!(c_timer, 56)
    build_s_matrix_and_g_vector!(
        S,
        g,
        smat_to_para_idx,
        sr_opt_oo,
        sr_opt_ho,
        sr_opt_size,
        data.modpara.dsr_opt_sta_del,
        data.modpara.dsr_opt_step_dt,
    )
    ctimer_stop!(c_timer, 56)

    # Phase 4: Solve S*x = g using LAPACK (DPOSV equivalent; [57] DPOSV)
    info = 0
    ctimer_start!(c_timer, 57)
    try
        # Cholesky decomposition (upper triangular)
        potrf!('U', S)

        # Forward/backward substitution (overwrites g with solution x)
        potrs!('U', S, g)
    catch e
        @error "DPOSV failed: $e"
        info = 1
    end
    ctimer_stop!(c_timer, 57)
    ctimer_stop!(c_timer, 51)

    # [52] postprocess: finite check + parameter update
    ctimer_start!(c_timer, 52)

    # Phase 5: Check for inf/nan
    if info == 0
        for si = 1:n_smat
            if !isfinite(g[si])
                @error "StcOpt: r[$si] = $(g[si]) is not finite"
                info = 1
                break
            end
        end
    end

    # Phase 6: Update parameters
    # C: para[pi/2] += r[si] for real part, para[(pi-1)/2] += r[si]*I for imag part
    # Note: g[si] contains the solution r[si] after potrs! (g is overwritten)
    # The g vector already contains -DSROptStepDt*2.0*... from build_s_matrix_and_g_vector!
    if info == 0
        for (si, pi) in enumerate(smat_to_para_idx)
            r_val = g[si]  # Solution vector (g was overwritten by potrs!)

            if pi % 2 == 0  # Real part
                para_idx = div(pi, 2) + 1  # 1-based
                update_parameter_value(data, para_idx, r_val, 0.0)
            else  # Imaginary part
                para_idx = div(pi - 1, 2) + 1  # 1-based
                update_parameter_value(data, para_idx, 0.0, r_val)
            end
        end
    end
    ctimer_stop!(c_timer, 52)

    return info
end

"""
    xdot(p::Vector{Float64}, q::Vector{Float64}) -> Float64

Compute dot product of two vectors.
Equivalent to C's `xdot()`.
"""
function xdot(p::Vector{Float64}, q::Vector{Float64})::Float64
    return dot(p, q)
end

"""
    CGWorkspace

Workspace for Conjugate Gradient method.
Equivalent to C's VecCG layout.
"""
mutable struct CGWorkspace
    x::Vector{Float64}           # Solution vector [nSmat]
    g::Vector{Float64}           # Gradient vector [nSmat]
    sdiag::Vector{Float64}       # Diagonal elements of S [nSmat]
    stcO::Vector{Float64}        # <O> mean values [nSmat]
    stcOs_real::Matrix{Float64}  # Sample data O_real [nSmat, NVMCSample]
    stcOs_imag::Matrix{Float64}  # Sample data O_imag [nSmat, NVMCSample] (for complex)
    q::Vector{Float64}           # CG: S*d [nSmat]
    d::Vector{Float64}           # CG: search direction [nSmat]
    r::Vector{Float64}           # CG: residual [nSmat]
    y_real::Vector{Float64}      # Temporary: [NVMCSample]
    y_imag::Vector{Float64}      # Temporary: [NVMCSample] (for complex)

    function CGWorkspace(n_smat::Int, n_vmc_sample::Int, all_complex::Bool)
        new(
            zeros(Float64, n_smat),
            zeros(Float64, n_smat),
            zeros(Float64, n_smat),
            zeros(Float64, n_smat),
            zeros(Float64, n_smat, n_vmc_sample),
            all_complex ? zeros(Float64, n_smat, n_vmc_sample) : zeros(Float64, 0, 0),
            zeros(Float64, n_smat),
            zeros(Float64, n_smat),
            zeros(Float64, n_smat),
            zeros(Float64, n_vmc_sample),
            all_complex ? zeros(Float64, n_vmc_sample) : zeros(Float64, 0),
        )
    end
end

"""
    stochastic_opt_cg_init!(
        ws::CGWorkspace,
        smat_to_para_idx::Vector{Int},
        sr_opt_o::Vector{ComplexF64},
        sr_opt_oo_diag::Vector{ComplexF64},
        sr_opt_o_store::Vector{ComplexF64},
        sr_opt_ho::Vector{ComplexF64},
        sr_opt_size::Int,
        n_vmc_sample::Int,
        dsr_opt_step_dt::Float64,
        all_complex::Bool
    )

Initialize CG workspace.
Equivalent to C's `StochasticOptCG_Init()`.
"""
function stochastic_opt_cg_init!(
    ws::CGWorkspace,
    smat_to_para_idx::Vector{Int},
    sr_opt_o::Vector{ComplexF64},
    sr_opt_oo_diag::Vector{ComplexF64},
    sr_opt_o_store::Vector{ComplexF64},
    sr_opt_ho::Vector{ComplexF64},
    sr_opt_size::Int,
    n_vmc_sample::Int,
    dsr_opt_step_dt::Float64,
    all_complex::Bool,
)
    n_smat = length(smat_to_para_idx)
    dt = 2.0 * dsr_opt_step_dt
    offset_factor = all_complex ? 2 : 1

    # Get HO[0] (energy)
    sr_opt_ho_0 = real(sr_opt_ho[1])

    # Initialize sample data: stcOs[si, sample] = O[pi+offset][sample]
    for sample = 1:n_vmc_sample
        # C: offset = (sample-1)*OFFSET*SROptSize
        offset = (sample - 1) * offset_factor * sr_opt_size

        for (si, pi) in enumerate(smat_to_para_idx)
            # C: idx = si + (sample-1)*nSmat (column major)
            # C: stcOs_real[idx] = CREAL(srOptO_Store[offset+pi+OFFSET])
            store_idx = offset + pi + offset_factor + 1  # 1-based

            if store_idx <= length(sr_opt_o_store)
                ws.stcOs_real[si, sample] = real(sr_opt_o_store[store_idx])
                if all_complex && !isempty(ws.stcOs_imag)
                    ws.stcOs_imag[si, sample] = imag(sr_opt_o_store[store_idx])
                end
            end
        end
    end

    # Initialize g vector and diagonal elements
    for (si, pi) in enumerate(smat_to_para_idx)
        # C: stcO[si] = CREAL(srOptO[pi+OFFSET])
        # C: g[si] = -dt*(CREAL(srOptHO[pi+OFFSET]) - srOptHO_0 * CREAL(srOptO[pi+OFFSET]))
        # C: sdiag[si] = CREAL(srOptOOdiag[pi+OFFSET]) - CREAL(srOptO[pi+OFFSET])^2
        # C: srOptO = SROptOO (or SROptOO_real), so srOptO[pi+OFFSET] = SROptOO[pi+OFFSET]
        # Julia: sr_opt_o is extracted from sr_opt_oo[1:...], so sr_opt_o[pi+OFFSET] corresponds to sr_opt_oo[pi+OFFSET+1]
        # But since sr_opt_o starts at index 1, we use pi+offset_factor (not +1)
        o_idx = pi + offset_factor  # C: pi+OFFSET -> Julia: pi+offset_factor (sr_opt_o is 1-based but extracted from sr_opt_oo[1:...])

        if o_idx <= length(sr_opt_o)
            o_val = real(sr_opt_o[o_idx])
            ws.stcO[si] = o_val

            # C: srOptHO[pi+OFFSET] -> Julia: sr_opt_ho[pi+offset_factor+1] (since sr_opt_ho[1] = HO[0])
            ho_idx = pi + offset_factor + 1  # C: pi+OFFSET (0-based) -> Julia: pi+offset_factor+1 (1-based)
            ho_val = ho_idx <= length(sr_opt_ho) ? real(sr_opt_ho[ho_idx]) : 0.0
            ws.g[si] = -dt * (ho_val - sr_opt_ho_0 * o_val)

            # C: srOptOOdiag[pi+OFFSET] -> Julia: sr_opt_oo_diag[pi+offset_factor] (since sr_opt_oo_diag is extracted from sr_opt_oo[offset+1:...])
            oo_diag_idx = pi + offset_factor  # C: pi+OFFSET -> Julia: pi+offset_factor (sr_opt_oo_diag is 1-based but extracted)
            oo_diag_val =
                oo_diag_idx <= length(sr_opt_oo_diag) ? real(sr_opt_oo_diag[oo_diag_idx]) :
                0.0
            ws.sdiag[si] = oo_diag_val - o_val * o_val
        end
    end

    # Initialize x to zero
    fill!(ws.x, 0.0)
end

"""
    operate_by_s!(
        z::Vector{Float64},
        x::Vector{Float64},
        ws::CGWorkspace,
        n_smat::Int,
        n_vmc_sample::Int,
        inv_w::Float64,
        dsr_opt_sta_del::Float64,
        all_complex::Bool
    )

Compute z = S*x without explicitly building S matrix.
Equivalent to C's `operate_by_S()`.

S[i][j] = (1/W) * sum{sample} O[i][sample] * O[j][sample] - <O[i]> * <O[j]>
"""
function operate_by_s!(
    z::Vector{Float64},
    x::Vector{Float64},
    ws::CGWorkspace,
    n_smat::Int,
    n_vmc_sample::Int,
    inv_w::Float64,
    dsr_opt_sta_del::Float64,
    all_complex::Bool,
)
    # y_real[sample] = sum{si} x[si] * O_real[si, sample]
    # = stcOs_real^T * x
    mul!(ws.y_real, ws.stcOs_real', x)

    if all_complex && !isempty(ws.stcOs_imag)
        # y_imag[sample] = sum{si} x[si] * O_imag[si, sample]
        mul!(ws.y_imag, ws.stcOs_imag', x)
    end

    # z[si] = sum{sample} O_real[si, sample] * y_real[sample]
    # = stcOs_real * y_real
    mul!(z, ws.stcOs_real, ws.y_real)

    if all_complex && !isempty(ws.stcOs_imag)
        # z[si] += sum{sample} O_imag[si, sample] * y_imag[sample]
        # z += stcOs_imag * y_imag
        mul!(z, ws.stcOs_imag, ws.y_imag, 1.0, 1.0)
    end

    # Compute <O>^T * x
    coef = xdot(ws.stcO, x)

    # z = invW * z - coef * <O> + DSROptStaDel * sdiag * x
    for si = 1:n_smat
        z[si] = inv_w * z[si] - coef * ws.stcO[si] + dsr_opt_sta_del * ws.sdiag[si] * x[si]
    end
end

"""
    stochastic_opt_cg_main!(
        ws::CGWorkspace,
        n_smat::Int,
        n_vmc_sample::Int,
        inv_w::Float64,
        dsr_opt_sta_del::Float64,
        dsr_opt_cg_tol::Float64,
        max_iter::Int,
        all_complex::Bool
    ) -> Int

Main CG iteration loop.
Equivalent to C's `StochasticOptCG_Main()`.

# Returns
- Number of iterations performed
"""
function stochastic_opt_cg_main!(
    ws::CGWorkspace,
    n_smat::Int,
    n_vmc_sample::Int,
    inv_w::Float64,
    dsr_opt_sta_del::Float64,
    dsr_opt_cg_tol::Float64,
    max_iter::Int,
    all_complex::Bool,
)::Int
    # Convergence threshold
    cg_thresh = dsr_opt_cg_tol^2 * Float64(n_smat)^2

    # Initialize: d = r = g
    copyto!(ws.d, ws.g)
    copyto!(ws.r, ws.g)

    # delta = r^T * r
    delta = xdot(ws.r, ws.r)

    iter = 0
    for iter_idx = 1:max_iter
        iter = iter_idx

        # Check convergence
        if delta < cg_thresh
            iter = iter_idx - 1
            break
        end

        # Compute q = S * d
        operate_by_s!(
            ws.q,
            ws.d,
            ws,
            n_smat,
            n_vmc_sample,
            inv_w,
            dsr_opt_sta_del,
            all_complex,
        )

        # alpha = delta / (d^T * q)
        dq = xdot(ws.d, ws.q)
        if abs(dq) < 1e-30
            break
        end
        alpha = delta / dq

        # Update solution: x = x + alpha * d
        for si = 1:n_smat
            ws.x[si] += alpha * ws.d[si]
        end

        # Update residual: r = r - alpha * q
        # (Every 20 iterations, recompute r = g - S*x for numerical stability)
        if iter_idx % 20 == 0
            operate_by_s!(
                ws.r,
                ws.x,
                ws,
                n_smat,
                n_vmc_sample,
                inv_w,
                dsr_opt_sta_del,
                all_complex,
            )
            for si = 1:n_smat
                ws.r[si] = ws.g[si] - ws.r[si]
            end
        else
            for si = 1:n_smat
                ws.r[si] -= alpha * ws.q[si]
            end
        end

        # beta = (r_new^T * r_new) / delta
        delta_new = xdot(ws.r, ws.r)
        beta = delta_new / delta

        # Update delta
        delta = delta_new

        # Update search direction: d = r + beta * d
        for si = 1:n_smat
            ws.d[si] = ws.r[si] + beta * ws.d[si]
        end
    end

    return iter
end

"""
    stochastic_opt_cg!(data::ExpertModeData, state::VMCOptimizationState) -> Int

Stochastic optimization using Conjugate Gradient method.
Equivalent to C's `StochasticOptCG()`.

Uses the CG method to solve S*x = g without explicitly building the S matrix,
which is more memory-efficient for large parameter spaces.

# Returns
- `info::Int`: 0 = success, non-zero = error (number of CG iterations if successful)
"""
function stochastic_opt_cg!(data::ExpertModeData, state::VMCOptimizationState, c_timer::CTimer = CTIMER_DISABLED)::Int
    # NOTE: CG-solver SR internals ([50]-[58]) are not broken down yet; the
    # total CG time is captured by [5] StochasticOpt at the vmc_para_opt! level.
    # The direct solver (stochastic_opt!) carries the [50]/[51]/[56]/[57]/[52]
    # breakdown. c_timer is accepted here for a uniform call from vmc_para_opt!.
    # Get parameters
    n_para = _parameter_count_breakdown(data).n_para
    sr_opt_size = state.sr_opt.sr_opt_size
    n_vmc_sample = data.modpara.nvmc_sample
    all_complex = get_all_complex_flag(data)
    offset_factor = all_complex ? 2 : 1
    n_para_full = offset_factor * n_para

    # Initialize optimization flags if empty (same as vmc_para_opt!)
    if isempty(data.optimization_flags)
        data.optimization_flags = fill(true, 2 * n_para)
    end

    # Get SR arrays
    sr_opt_oo = state.sr_opt.sr_opt_oo
    sr_opt_ho = state.sr_opt.sr_opt_ho
    sr_opt_o_store = state.sr_opt.sr_opt_o_store

    # For CG, we need sr_opt_o (first row of sr_opt_oo) and sr_opt_oo_diag
    # C: srOptO = SROptOO (or SROptOO_real)
    # C: srOptOOdiag = SROptOO + 2*SROptSize (or SROptOO_real + SROptSize)

    # Extract <O> from first row of SROptOO
    sr_opt_o = sr_opt_oo[1:min(offset_factor*sr_opt_size, length(sr_opt_oo))]

    # Extract diagonal: SROptOO + offset*SROptSize
    diag_offset = offset_factor * sr_opt_size
    sr_opt_oo_diag = if diag_offset + offset_factor * sr_opt_size <= length(sr_opt_oo)
        sr_opt_oo[(diag_offset+1):(diag_offset+offset_factor*sr_opt_size)]
    else
        zeros(ComplexF64, offset_factor * sr_opt_size)
    end

    # Phase 1: Calculate diagonal elements and apply redundant direction cut
    # C: sDiagElm[pi] = CREAL(srOptOOdiag[pi+OFFSET]) - CREAL(srOptO[pi+OFFSET])^2
    # C: srOptO = SROptOO, srOptOOdiag = SROptOO + 2*SROptSize (or + SROptSize for real)
    s_diag_elm = zeros(Float64, n_para_full)
    for pi = 0:(n_para_full-1)
        # C: pi+OFFSET (0-based) -> Julia: pi+offset_factor (since arrays are extracted from sr_opt_oo[offset+1:...])
        idx = pi + offset_factor  # C: pi+OFFSET -> Julia: pi+offset_factor
        if idx <= length(sr_opt_oo_diag) && idx <= length(sr_opt_o)
            s_diag_elm[pi+1] = real(sr_opt_oo_diag[idx]) - real(sr_opt_o[idx])^2
        end
    end

    # Find max and min
    s_diag_max = maximum(s_diag_elm)
    s_diag_min = minimum(s_diag_elm)

    # Calculate threshold
    diag_cut_threshold = s_diag_max * data.modpara.dsr_opt_red_cut

    # Build smat_to_para_idx mapping
    smat_to_para_idx = Int[]
    opt_num = 0
    cut_num = 0

    for pi = 0:(n_para_full-1)
        # For real mode, check OptFlag[2*pi], for complex mode check OptFlag[pi]
        opt_flag_idx = all_complex ? pi : 2 * pi
        opt_flag = get_opt_flag_for_parameter(data, opt_flag_idx)

        if opt_flag != 1  # Fixed parameter
            opt_num += 1
            continue
        end

        if s_diag_elm[pi+1] < diag_cut_threshold  # Cut by redundant direction
            cut_num += 1
        else  # Optimize
            push!(smat_to_para_idx, pi)
        end
    end

    n_smat = length(smat_to_para_idx)

    # If no parameters to optimize, return success
    if n_smat == 0
        @warn "StochasticOptCG: No parameters to optimize (all fixed or cut)"
        return 0
    end

    # Phase 2: Create CG workspace and initialize
    ws = CGWorkspace(n_smat, n_vmc_sample, all_complex)

    stochastic_opt_cg_init!(
        ws,
        smat_to_para_idx,
        sr_opt_o,
        sr_opt_oo_diag,
        sr_opt_o_store,
        sr_opt_ho,
        sr_opt_size,
        n_vmc_sample,
        data.modpara.dsr_opt_step_dt,
        all_complex,
    )

    # Phase 3: Run CG iteration
    # inv_w = 1/Wc where Wc is the total weight (for unweighted samples, inv_w = 1/n_vmc_sample)
    inv_w = 1.0 / Float64(n_vmc_sample)

    max_iter =
        data.modpara.nsr_opt_cg_max_iter > 0 ? data.modpara.nsr_opt_cg_max_iter : n_smat

    iter = stochastic_opt_cg_main!(
        ws,
        n_smat,
        n_vmc_sample,
        inv_w,
        data.modpara.dsr_opt_sta_del,
        data.modpara.dsr_opt_cg_tol,
        max_iter,
        all_complex,
    )

    # Phase 4: Check for inf/nan in solution
    info = 0
    for si = 1:n_smat
        if !isfinite(ws.x[si])
            @error "StochasticOptCG: x[$si] = $(ws.x[si]) is not finite"
            info = 1
            break
        end
    end

    # Phase 5: Update parameters
    if info == 0
        for (si, pi) in enumerate(smat_to_para_idx)
            r_val = ws.x[si]  # Solution vector

            if all_complex
                # Complex mode: pi % 2 == 0 -> real, pi % 2 == 1 -> imag
                if pi % 2 == 0
                    para_idx = div(pi, 2) + 1
                    update_parameter_value(data, para_idx, r_val, 0.0)
                else
                    para_idx = div(pi - 1, 2) + 1
                    update_parameter_value(data, para_idx, 0.0, r_val)
                end
            else
                # Real mode: directly update para[pi]
                para_idx = pi + 1
                update_parameter_value(data, para_idx, r_val, 0.0)
            end
        end
    end

    # Log output
    if info == 0 && n_smat > 0
        r_max = maximum(abs, ws.x[1:n_smat])
        # @info "SR-CG Info: NPara=$n_para, nSmat=$n_smat, iter=$iter, " *
        #       "sDiagMax=$s_diag_max, sDiagMin=$s_diag_min, rmax=$r_max"
    end

    return info
end
