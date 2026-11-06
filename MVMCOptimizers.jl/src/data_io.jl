"""
Data I/O Functions

Store and output optimization data.
"""

using Printf

"""
    store_opt_data!(data::ExpertModeData, state::VMCOptimizationState, sample_idx::Int)

Store optimization data for averaging.
Equivalent to C's `StoreOptData()`.

Stores measured pre-SR energy/E² and synchronized post-SR parameters.
"""
function store_opt_data!(data::ExpertModeData, state::VMCOptimizationState, sample_idx::Int)
    sample_idx == length(state.opt_data) ||
        throw(ArgumentError("optimization window must be stored chronologically"))
    # C avevar.c:82–89: unique contiguous Para, not expanded site-pair rows.
    # Own the snapshot: subsequent synchronization must not mutate history.
    push!(state.opt_data, OptDataPoint(state.energy.etot, state.energy.etot2,
        copy(pack_parameters(data))))
end

# Helper: return path, creating output_dir if given
function _output_path(filename::String, output_dir::Union{String,Nothing})
    if output_dir !== nothing && !isempty(output_dir)
        mkpath(output_dir)
        return joinpath(output_dir, filename)
    end
    return filename
end

"""
    output_data!(data::ExpertModeData, state::VMCOptimizationState, step::Int;
                 output_dir=nothing, file_index=nothing)

Output data to files.
Equivalent to C's `outputData()`.

With `file_index === nothing` (optimization) this writes the single
`zvo_out.dat` / `zvo_var.dat` pair; step 0 overwrites and step >= 1 appends.
With `file_index::Int` (PhysCal) it writes the C-indexed per-sample pair
`zvo_out_<file_index:03d>.dat` / `zvo_var_<file_index:03d>.dat` and always
truncates, matching C `InitFilePhysCal`, which opens each sample's file with
`"w"` (mVMC 1.4.0 `src/mVMC/initfile.c:77-95`).
If `output_dir` is set, files are written under that directory (directory is created if needed).
"""
function output_data!(data::ExpertModeData, state::VMCOptimizationState, step::Int; output_dir::Union{String,Nothing}=nothing, file_index::Union{Nothing,Int}=nothing)
    # Get file head from parameters
    data_file_head = data.modpara.c_data_file_head
    if isempty(data_file_head)
        data_file_head = "zvo"
    end

    # Get energy values
    etot = state.energy.etot
    etot2 = state.energy.etot2

    # Calculate variance: (E2 - E^2) / E^2
    variance = if abs(etot) > 1e-14
        real((etot2 - etot * etot) / (etot * etot))
    else
        0.0
    end

    # Sz / Sz^2: weighted-averaged in vmc_main_cal_fsz!. For non-fsz code paths
    # they remain 0, matching the C reference (vmccal.c:576 initialises to 0
    # and never accumulates outside the fsz path).
    sztot = real(state.energy.sztot)
    sztot2 = real(state.energy.sztot2)

    # C PhysCal writes one C-indexed file per sample and truncates it each time
    # (InitFilePhysCal opens "%s_out_%03d.dat" with "w"). Optimization keeps the
    # single non-indexed file with truncate-on-step-0/append semantics.
    suffix = file_index === nothing ? "" : @sprintf("_%03d", file_index)
    write_mode = file_index === nothing ? (step == 0 ? "w" : "a") : "w"

    # Output to zvo_out[_NNN].dat
    out_file = _output_path(data_file_head * "_out" * suffix * ".dat", output_dir)
    open(out_file, write_mode) do f
        # C format: "% .18e % .18e  % .18e % .18e %.18e %.18e\n"
        @printf(
            f,
            "% .18e % .18e  % .18e % .18e %.18e %.18e\n",
            real(etot),
            imag(etot),
            real(etot2),
            variance,
            sztot,
            sztot2
        )
    end

    # Output to zvo_var[_NNN].dat
    var_file = _output_path(data_file_head * "_var" * suffix * ".dat", output_dir)
    open(var_file, write_mode) do f
        # C format: "% .18e % .18e 0.0 % .18e % .18e 0.0 " + parameters
        @printf(
            f,
            "% .18e % .18e 0.0 % .18e % .18e 0.0 ",
            real(etot),
            imag(etot),
            real(etot2),
            imag(etot2)
        )

        # C vmcmain.c:655–657 writes every Para[0:NPara-1] slot:
        # projection (including DH2/DH4), RBM, declared Slater, OptTrans.
        # Mapped orbital rows may repeat an index; they are not extra Para.
        for parameter in pack_parameters(data)
            @printf(f, "% .18e % .18e 0.0 ", real(parameter), imag(parameter))
        end
        @printf(f, "\n")
    end
end

"""
    output_opt_data!(data, state; output_dir=nothing)

C avevar.c OutputOptData: measured E/E² plus unique post-SR Para history.
A one-sample window writes real/zero pairs; larger windows write complex
means and sample deviations. Validate the entire window before opening files.
"""
function output_opt_data!(data::ExpertModeData, state::VMCOptimizationState;
    output_dir::Union{String,Nothing}=nothing)
    counts = _parameter_count_breakdown(data)
    history = state.opt_data
    length(history) == data.modpara.nsr_opt_itr_smp > 0 ||
        throw(ArgumentError("incomplete optimization window"))
    all(point -> length(point.parameters) == counts.n_para, history) ||
        throw(ArgumentError("optimization window parameter shape mismatch"))
    rows = [vcat(ComplexF64[point.energy, point.energy_squared], point.parameters)
            for point in history]
    head = isempty(data.modpara.c_para_file_head) ? "zqp" : data.modpara.c_para_file_head
    # Family ranges follow the same C layout as pack_parameters.
    layout = counts.layout
    families = Tuple{String,String,Int,Int}[
        ("gutzwiller", "NGutzwillerIdx", layout.n_gutzwiller, layout.n_gutzwiller),
        ("jastrow", "NJastrowIdx", layout.n_jastrow, layout.n_jastrow),
        ("doublonHolon2site", "NDoublonHolon2siteIdx", layout.n_dh2, 6 * layout.n_dh2),
        ("doublonHolon4site", "NDoublonHolon4siteIdx", layout.n_dh4, 10 * layout.n_dh4),
    ]
    names = ("chargeRBM_physlayer", "spinRBM_physlayer", "generalRBM_physlayer",
        "chargeRBM_hiddenlayer", "spinRBM_hiddenlayer", "generalRBM_hiddenlayer",
        "chargeRBM_physhidden", "spinRBM_physhidden", "generalRBM_physhidden")
    headers = ("NChargeRBM_PhysLayerIdx", "NSpinRBM_PhysLayerIdx", "NGeneralRBM_PhysLayerIdx",
        "NChargeRBM_HiddenLayerIdx", "NSpinRBM_HiddenLayerIdx", "NGeneralRBM_HiddenLayerIdx",
        "NChargeRBM_PhysHiddenIdx", "NSpinRBM_PhysHiddenIdx", "NGeneralRBM_PhysHiddenIdx")
    for (name, header, terms) in zip(names, headers, _rbm_parameter_sections(data))
        width = _parameter_section_width(terms)
        push!(families, (name, header, width, width))
    end
    if data.i_flg_orbital_general == 0
        push!(families, ("orbital", "NOrbitalIdx", counts.n_orbital_idx, counts.n_orbital_idx))
    elseif data.i_flg_orbital_parallel != 0
        anti = data.n_orbital_anti_parallel
        parallel = counts.n_orbital_idx - anti
        0 <= anti <= counts.n_orbital_idx && parallel > 0 && iseven(parallel) ||
            throw(ArgumentError("invalid AntiParallel + 2*Parallel parameter layout"))
        push!(families, ("orbitalAntiParallel", "NOrbitalAntiParallelIdx", anti, anti))
        push!(families, ("orbitalParallel", "NOrbitalParallelIdx", parallel, parallel))
    else
        push!(families, ("orbital_general", "NOrbitalIdx", counts.n_orbital_idx, counts.n_orbital_idx))
    end
    push!(families, ("trans", "NQPOptTrans", counts.n_opt_trans, counts.n_opt_trans))
    sum(family[4] for family in families) == counts.n_para ||
        throw(ArgumentError("optimization output family layout mismatch"))
    open(_output_path(head * "_opt.dat", output_dir), "w") do root
        if length(rows) == 1
            for value in rows[1]
                @printf(root, "% .18e % .18e ", real(value), 0.0)
            end
        else
            # C CalcAveVar: chronological sums; denominator is n-1, not n.
            statistics = map(eachindex(rows[1])) do index
                mean = 0.0 + 0.0im
                for row in rows
                    mean += row[index]
                end
                mean /= length(rows)
                variance = 0.0
                for row in rows
                    delta = row[index] - mean
                    variance += real(delta * conj(delta))
                end
                (mean, sqrt(variance / (length(rows) - 1)))
            end
            for (mean, deviation) in statistics[1:2]
                @printf(root, "% .18e % .18e % .18e ", real(mean), imag(mean), deviation)
            end
            offset = 2
            for (name, header, header_count, width) in families
                width == 0 && continue
                open(_output_path(head * "_" * name * "_opt.dat", output_dir), "w") do child
                    println(child, "======================")
                    println(child, header, "  ", header_count)
                    for _ in 1:3
                        println(child, "======================")
                    end
                    for index in 1:width
                        mean, deviation = statistics[offset + index]
                        @printf(child, "%d % .18e % .18e \n", index - 1, real(mean), imag(mean))
                        @printf(root, "% .18e % .18e % .18e ", real(mean), imag(mean), deviation)
                    end
                end
                offset += width
            end
        end
        println(root)
    end
end

function _lanczos_energy_by_alpha(
    h1::Float64,
    h2_1::Float64,
    h2_2::Float64,
    h3::Float64,
    h4::Float64,
    alpha::Float64,
)
    tmp_ene = h1 + alpha * (h2_1 + h2_2) + alpha^2 * h3
    dnorm = 1.0 + 2.0 * alpha * h1 + alpha^2 * h2_1
    tmp_ene_v = h2_1 + 2.0 * alpha * h3 + alpha^2 * h4
    (!isfinite(h1) || abs(h1) < eps(Float64)) && return nothing
    norm_ratio = dnorm / h1
    (!isfinite(norm_ratio) || abs(norm_ratio) < 1.0e-12) && return nothing
    ene = tmp_ene / dnorm
    ene_v = ((tmp_ene_v / dnorm) - ene^2) / ene^2
    (!isfinite(ene) || !isfinite(ene_v)) && return nothing
    return ene, ene_v
end

function _invalid_lanczos_energy(reason::AbstractString)
    @warn "Lanczos energy could not be calculated; writing NaN values" reason
    return NaN, NaN, NaN
end

function _lanczos_energy(qqqq::AbstractVector{ComplexF64})
    h1 = real(qqqq[3])
    h2_1 = real(qqqq[4])
    h2_2 = real(qqqq[11])
    h3 = real(qqqq[12])
    h4 = real(qqqq[16])

    tmp_aa = h2_1 * (h2_1 + h2_2) - 2.0 * h1 * h3
    tmp_bb = -h1 * h2_1 + h3
    tmp_cc =
        h2_1 * (h2_1 + h2_2)^2 -
        h1^2 * h2_1 * (h2_1 + 2.0 * h2_2) +
        4.0 * h1^3 * h3 -
        2.0 * h1 * (2.0 * h2_1 + h2_2) * h3 +
        h3^2
    if !(isfinite(tmp_aa) && isfinite(tmp_bb) && isfinite(tmp_cc))
        return _invalid_lanczos_energy("non-finite alpha equation")
    end
    tmp_cc < 0.0 && return _invalid_lanczos_energy("negative alpha discriminant")
    abs(tmp_aa) < eps(Float64) && return _invalid_lanczos_energy("singular alpha equation")

    root = sqrt(tmp_cc)
    alpha_p = (tmp_bb + root) / tmp_aa
    alpha_m = (tmp_bb - root) / tmp_aa
    if !(isfinite(alpha_p) && isfinite(alpha_m))
        return _invalid_lanczos_energy("non-finite alpha")
    end
    energy_p = _lanczos_energy_by_alpha(h1, h2_1, h2_2, h3, h4, alpha_p)
    energy_m = _lanczos_energy_by_alpha(h1, h2_1, h2_2, h3, h4, alpha_m)
    (energy_p === nothing || energy_m === nothing) &&
        return _invalid_lanczos_energy("illegal norm or non-finite energy")
    ene_p, ene_vp = energy_p
    ene_m, ene_vm = energy_m

    if ene_p > ene_m
        return ene_m, ene_vm, alpha_m
    end
    return ene_p, ene_vp, alpha_p
end

function _lanczos_phys_values(
    qqqq::AbstractVector{ComplexF64},
    qphysq::AbstractVector{ComplexF64},
    nphys::Int,
    alpha::Float64,
)
    length(qphysq) == 4 * nphys ||
        throw(
            ArgumentError(
                "Lanczos QPhysQ length $(length(qphysq)) does not match 4 * nphys = " *
                "$(4 * nphys).",
            ),
        )
    values = Vector{ComplexF64}(undef, nphys)
    h1 = qqqq[3]
    h2_1 = qqqq[4]
    dnorm = real(1.0 + 2.0 * alpha * h1 + alpha^2 * h2_1)
    @inbounds for i in 1:nphys
        a0 = qphysq[i]
        a1_01 = qphysq[nphys+i]
        a1_10 = qphysq[2*nphys+i]
        a2_11 = qphysq[3*nphys+i]
        values[i] = (a0 + alpha * (a1_01 + a1_10) + alpha^2 * a2_11) / dnorm
    end
    return values
end

function output_lanczos_func!(
    data::ExpertModeData,
    state::VMCOptimizationState,
    ismp::Int;
    output_dir::Union{String,Nothing} = nothing,
)
    data.modpara.lanczos_mode > 0 || return
    state.phys_quantities === nothing && return

    phys = state.phys_quantities
    data_file_head = data.modpara.c_data_file_head
    isempty(data_file_head) && (data_file_head = "zvo")

    ene, ene_v, alpha = _lanczos_energy(phys.phys_lanczos_qqqq)

    ls_filename = _output_path(@sprintf("%s_ls_out_%03d.dat", data_file_head, ismp), output_dir)
    open(ls_filename, "w") do f
        @printf(f, "% .18e  ", ene)
        @printf(f, "% .18e  ", ene_v)
        @printf(f, "% .18e  ", alpha)
    end

    qqqq_filename = _output_path(@sprintf("%s_ls_qqqq_%03d.dat", data_file_head, ismp), output_dir)
    open(qqqq_filename, "w") do f
        for val in phys.phys_lanczos_qqqq
            @printf(f, "% .18e  ", real(val))
        end
        println(f)
    end

    data.modpara.lanczos_mode > 1 || return

    ls_cis_ajs = _lanczos_phys_values(
        phys.phys_lanczos_qqqq,
        phys.phys_lanczos_qcisajsq,
        length(phys.cis_ajs_idx),
        alpha,
    )
    cisajs_filename =
        _output_path(@sprintf("%s_ls_cisajs_%03d.dat", data_file_head, ismp), output_dir)
    open(cisajs_filename, "w") do f
        for (idx, (ri, si, rj, sj)) in enumerate(phys.cis_ajs_idx)
            val = ls_cis_ajs[idx]
            @printf(
                f,
                "%d %d %d %d % .18e % .18e \n",
                ri,
                si,
                rj,
                sj,
                real(val),
                imag(val)
            )
        end
        println(f)
    end

    ls_cis_ajs_ckt_alt_dc = _lanczos_phys_values(
        phys.phys_lanczos_qqqq,
        phys.phys_lanczos_qcisajscktaltq_dc,
        length(data.green_two_terms),
        alpha,
    )
    cktalt_filename =
        _output_path(@sprintf("%s_ls_cisajscktalt_%03d.dat", data_file_head, ismp), output_dir)
    open(cktalt_filename, "w") do f
        for (idx, term) in enumerate(data.green_two_terms)
            si = term.spin1 == :up ? 0 : 1
            sj = term.spin2 == :up ? 0 : 1
            sk = term.spin3 == :up ? 0 : 1
            sl = term.spin4 == :up ? 0 : 1
            val = ls_cis_ajs_ckt_alt_dc[idx]
            @printf(
                f,
                "%d %d %d %d %d %d %d %d % .18e % .18e\n",
                term.site1,
                si,
                term.site2,
                sj,
                term.site3,
                sk,
                term.site4,
                sl,
                real(val),
                imag(val)
            )
        end
        println(f)
    end

    ls_cis_ajs_ckt_alt = _lanczos_phys_values(
        phys.phys_lanczos_qqqq,
        phys.phys_lanczos_qcisajscktaltq,
        length(phys.cis_ajs_ckt_alt_idx),
        alpha,
    )
    cktaltex_filename =
        _output_path(@sprintf("%s_ls_cisajscktaltex_%03d.dat", data_file_head, ismp), output_dir)
    open(cktaltex_filename, "w") do f
        for val in ls_cis_ajs_ckt_alt
            @printf(f, "% .18e % .18e ", real(val), imag(val))
        end
        println(f)
    end
end

"""
    output_green_func!(data::ExpertModeData, state::VMCOptimizationState, ismp::Int; output_dir=nothing)

Output Green's functions to files.
Equivalent to C's `outputData()` Green's function output section.

Output files (XXX = ismp + NDataIdxStart):
- zvo_cisajs_XXX.dat: 1-body Green's function <c†_i c_j>
- zvo_cisajscktaltex_XXX.dat: factored two-body Green (product / `TwoBodyGEx`)
- zvo_cisajscktalt_XXX.dat: direct two-body Green (`TwoBodyG`)
If `output_dir` is set, files are written under that directory.
"""
function output_green_func!(data::ExpertModeData, state::VMCOptimizationState, ismp::Int; output_dir::Union{String,Nothing}=nothing)
    if state.phys_quantities === nothing
        return  # No physical quantities to output
    end

    phys = state.phys_quantities
    data_file_head = data.modpara.c_data_file_head
    if isempty(data_file_head)
        data_file_head = "zvo"
    end

    # zvo_cisajs_XXX.dat (1-body Green's function) — written from the canonical
    # one-body list (which, with TwoBodyGEx, includes appended factored
    # constituents in C order), not raw data.green_one_terms.
    if !isempty(phys.cis_ajs_idx)
        filename = _output_path(@sprintf("%s_cisajs_%03d.dat", data_file_head, ismp), output_dir)
        open(filename, "w") do f
            for (idx, (ri, si, rj, sj)) in enumerate(phys.cis_ajs_idx)
                val = phys.phys_cis_ajs[idx]
                # C format: "%d %d %d %d % .18e  % .18e \n" (ri, si, rj, sj, real, imag)
                @printf(
                    f,
                    "%d %d %d %d % .18e  % .18e \n",
                    ri,
                    si,
                    rj,
                    sj,
                    real(val),
                    imag(val)
                )
            end
            println(f)  # Empty line at end (C format)
        end
    end

    # zvo_cisajscktaltex_XXX.dat (2-body correlation, product)
    if !isempty(phys.phys_cis_ajs_ckt_alt)
        filename = _output_path(@sprintf("%s_cisajscktaltex_%03d.dat", data_file_head, ismp), output_dir)
        open(filename, "w") do f
            for val in phys.phys_cis_ajs_ckt_alt
                @printf(f, "% .18e  % .18e ", real(val), imag(val))
            end
            println(f)  # Newline at end
        end
    end

    # zvo_cisajscktalt_XXX.dat (2-body correlation, direct)
    if !isempty(data.green_two_terms)
        filename = _output_path(@sprintf("%s_cisajscktalt_%03d.dat", data_file_head, ismp), output_dir)
        open(filename, "w") do f
            for (idx, term) in enumerate(data.green_two_terms)
                val = phys.phys_cis_ajs_ckt_alt_dc[idx]
                # C format: "%d %d %d %d %d %d %d %d % .18e % .18e\n"
                # Format: ri, si, rj, sj, rk, sk, rl, sl, real, imag
                si = term.spin1 == :up ? 0 : 1
                sj = term.spin2 == :up ? 0 : 1
                sk = term.spin3 == :up ? 0 : 1
                sl = term.spin4 == :up ? 0 : 1
                @printf(
                    f,
                    "%d %d %d %d %d %d %d %d % .18e % .18e\n",
                    term.site1,
                    si,
                    term.site2,
                    sj,
                    term.site3,
                    sk,
                    term.site4,
                    sl,
                    real(val),
                    imag(val)
                )
            end
            println(f)  # Empty line at end (C format)
        end
    end

    output_lanczos_func!(data, state, ismp; output_dir = output_dir)
end
