"""
Data I/O Functions

Store and output optimization data.
"""

using Printf

"""
    store_opt_data!(data::ExpertModeData, state::VMCOptimizationState, sample_idx::Int)

Store optimization data for averaging.
Equivalent to C's `StoreOptData()`.

Stores energy and parameters for later averaging.
"""
function store_opt_data!(data::ExpertModeData, state::VMCOptimizationState, sample_idx::Int)
    # Collect current parameters
    parameters = ComplexF64[]

    # Add Gutzwiller parameters
    for term in data.gutzwiller_terms
        push!(parameters, term.value)
    end

    # Add Jastrow parameters
    for term in data.jastrow_terms
        push!(parameters, term.value)
    end

    # Add Orbital parameters
    for term in data.orbital_terms
        push!(parameters, term.value)
    end

    # Store data point
    opt_point = OptDataPoint(state.energy.etot, parameters)

    # Ensure we have enough space
    while length(state.opt_data) <= sample_idx
        push!(state.opt_data, OptDataPoint(0.0 + 0.0im, ComplexF64[]))
    end

    state.opt_data[sample_idx+1] = opt_point
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
    output_data!(data::ExpertModeData, state::VMCOptimizationState, step::Int; output_dir=nothing)

Output data to files.
Equivalent to C's `outputData()`.

Outputs to zvo_out.dat and zvo_var.dat.
Step 0 overwrites the files (new run); step >= 1 appends.
If `output_dir` is set, files are written under that directory (directory is created if needed).
"""
function output_data!(data::ExpertModeData, state::VMCOptimizationState, step::Int; output_dir::Union{String,Nothing}=nothing)
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

    # Step 0: overwrite (new run). Step >= 1: append.
    write_mode = step == 0 ? "w" : "a"

    # Output to zvo_out.dat
    out_file = _output_path(data_file_head * "_out.dat", output_dir)
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

    # Output to zvo_var.dat
    var_file = _output_path(data_file_head * "_var.dat", output_dir)
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

        # Output parameters (Gutzwiller, Jastrow, Orbital)
        for term in data.gutzwiller_terms
            @printf(f, "% .18e % .18e 0.0 ", real(term.value), imag(term.value))
        end
        for term in data.jastrow_terms
            @printf(f, "% .18e % .18e 0.0 ", real(term.value), imag(term.value))
        end
        for term in data.orbital_terms
            @printf(f, "% .18e % .18e 0.0 ", real(term.value), imag(term.value))
        end
        @printf(f, "\n")
    end
end

"""
    output_opt_data!(data::ExpertModeData; output_dir=nothing)

Output optimized parameters.
Equivalent to C's `OutputOptData()`.

Outputs to zqp_opt.dat and individual parameter files.
If `output_dir` is set, files are written under that directory (directory is created if needed).
"""
function output_opt_data!(data::ExpertModeData; output_dir::Union{String,Nothing}=nothing)
    # Get file head from parameters
    para_file_head = data.modpara.c_para_file_head
    if isempty(para_file_head)
        para_file_head = "zqp"
    end

    # Output to zqp_opt.dat
    opt_file = _output_path(para_file_head * "_opt.dat", output_dir)
    open(opt_file, "w") do f
        # Output Gutzwiller parameters
        if !isempty(data.gutzwiller_terms)
            for (i, term) in enumerate(data.gutzwiller_terms)
                @printf(f, "% .18e % .18e \n", real(term.value), imag(term.value))
            end
        end

        # Output Jastrow parameters
        if !isempty(data.jastrow_terms)
            for (i, term) in enumerate(data.jastrow_terms)
                @printf(f, "% .18e % .18e \n", real(term.value), imag(term.value))
            end
        end

        # Output Orbital (Slater) parameters
        if !isempty(data.orbital_terms)
            for (i, term) in enumerate(data.orbital_terms)
                @printf(f, "% .18e % .18e \n", real(term.value), imag(term.value))
            end
        end
    end

    # Output individual parameter files (optional, for detailed analysis)
    # Gutzwiller
    if !isempty(data.gutzwiller_terms)
        gutz_file = _output_path(para_file_head * "_gutzwiller_opt.dat", output_dir)
        open(gutz_file, "w") do f
            println(f, "===============================")
            println(f, "NGutzwillerIdx $(length(data.gutzwiller_terms))")
            println(f, "===============================")
            println(f, "===============================")
            for (i, term) in enumerate(data.gutzwiller_terms)
                @printf(f, "%d % .18e % .18e \n", i-1, real(term.value), imag(term.value))
            end
        end
    end

    # Jastrow
    if !isempty(data.jastrow_terms)
        jast_file = _output_path(para_file_head * "_jastrow_opt.dat", output_dir)
        open(jast_file, "w") do f
            println(f, "===============================")
            println(f, "NJastrowIdx $(length(data.jastrow_terms))")
            println(f, "===============================")
            println(f, "===============================")
            for (i, term) in enumerate(data.jastrow_terms)
                @printf(f, "%d % .18e % .18e \n", i-1, real(term.value), imag(term.value))
            end
        end
    end

    # Orbital (Slater)
    if !isempty(data.orbital_terms)
        orb_file = _output_path(para_file_head * "_orbital_opt.dat", output_dir)
        open(orb_file, "w") do f
            println(f, "===============================")
            println(f, "NOrbitalIdx $(length(data.orbital_terms))")
            println(f, "===============================")
            println(f, "===============================")
            for (i, term) in enumerate(data.orbital_terms)
                @printf(f, "%d % .18e % .18e \n", i-1, real(term.value), imag(term.value))
            end
        end
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
