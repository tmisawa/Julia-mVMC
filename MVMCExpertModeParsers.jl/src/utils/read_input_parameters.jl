"""
Read Input Parameters Utilities

Functions for reading parameter values from InGutzwiller.def, InJastrow.def, etc.
Equivalent to C's ReadInputParameters().
Based on mVMC/src/mVMC/readdef.c
"""

"""
    parse_input_parameter_file(filepath::String) -> Dict{Int, ComplexF64}

Parse an input parameter file (InGutzwiller.def, InJastrow.def, etc.).
Returns a dictionary mapping index to parameter value.

File format:
- First 5 lines are header (skipped)
- Line 2 contains count: "NGutzwillerIdx N"
- After header, lines are: "idx real imag"
"""
function parse_input_parameter_file(filepath::String)::Dict{Int,ComplexF64}
    if !isfile(filepath)
        return Dict{Int,ComplexF64}()
    end

    content = read_def_file(filepath)
    lines = split(content, '\n')

    # Skip first 5 lines (header)
    if length(lines) < 6
        return Dict{Int,ComplexF64}()
    end

    # Read count from line 2 (index 1, 0-based)
    expected_count = 0
    if length(lines) > 1
        header_line = clean_line(lines[2])
        tokens = split_def_line(header_line)
        if length(tokens) >= 2
            expected_count = safe_parse_int(tokens[2], 0)
        end
    end

    # Parse parameter values (starting from line 6, index 5)
    params = Dict{Int,ComplexF64}()
    data_start = 6  # After 5 header lines

    for i = data_start:length(lines)
        line = clean_line(lines[i])
        if isempty(line)
            continue
        end

        tokens = split_def_line(line)
        if length(tokens) >= 3
            idx = safe_parse_int(tokens[1], -1)
            real_val = safe_parse_float(tokens[2], 0.0)
            imag_val = safe_parse_float(tokens[3], 0.0)

            if idx >= 0
                params[idx] = ComplexF64(real_val, imag_val)
            end
        end
    end

    # Validate count
    if expected_count > 0 && length(params) != expected_count
        @warn "Parameter count mismatch in $filepath: expected $expected_count, got $(length(params))"
    end

    return params
end

function _rbm_section_nparam(terms)::Int
    if isempty(terms)
        return 0
    end
    return max(0, maximum(t.idx for t in terms) + 1)
end

function count_rbm_parameters(data::ExpertModeData)::Int
    return _rbm_section_nparam(data.charge_rbm_phys_layer_terms) +
           _rbm_section_nparam(data.spin_rbm_phys_layer_terms) +
           _rbm_section_nparam(data.general_rbm_phys_layer_terms) +
           _rbm_section_nparam(data.charge_rbm_hidden_layer_terms) +
           _rbm_section_nparam(data.spin_rbm_hidden_layer_terms) +
           _rbm_section_nparam(data.general_rbm_hidden_layer_terms) +
           _rbm_section_nparam(data.charge_rbm_phys_hidden_terms) +
           _rbm_section_nparam(data.spin_rbm_phys_hidden_terms) +
           _rbm_section_nparam(data.general_rbm_phys_hidden_terms)
end

function _set_rbm_terms_from_params!(terms, params::Dict{Int, ComplexF64})
    isempty(params) && return
    for term in terms
        if haskey(params, term.idx)
            term.value = params[term.idx]
        end
    end
end

"""
    read_input_parameters!(data::ExpertModeData, namelist_path::String) -> ExpertModeData

Apply optional `In*.def` parameter overlays to `data`. Equivalent
to C's `ReadInputParameters()` in `vmcmain.c:268`.

The base directory for resolving file paths in `namelist.def` is
inferred from `dirname(abspath(namelist_path))` and is **not** a
keyword argument.

# Consumed (writes back into `data`)

| `namelist.def` keyword | Updated field |
|---|---|
| `InGutzwiller` | `data.gutzwiller_terms[i].value` |
| `InJastrow` | `data.jastrow_terms[i].value` |
| `InOrbital` / `InOrbitalAntiParallel` | `data.orbital_terms[i].value` |
| `InOrbitalGeneral` | `data.orbital_terms[i].value` (fsz layout) |
| `InChargeRBM_PhysLayer` / `_HiddenLayer` / `_PhysHidden` | `charge_rbm_*_terms` |
| `InSpinRBM_PhysLayer` / `_HiddenLayer` / `_PhysHidden` | `spin_rbm_*_terms` |
| `InGeneralRBM_PhysLayer` / `_HiddenLayer` / `_PhysHidden` | `general_rbm_*_terms` |

# Recognised but not yet wired (warn-only stubs)

`InDH2`, `InDH4`, `InOrbitalParallel`, `InOptTrans` — these emit an
`@warn` and skip. The corresponding C-side data structures
(doublon-holon, parallel-orbital offset path, OptTrans + FlagOptTrans)
are not yet implemented in Julia.

# Phase ordering

Per `vmcmain.c:264-281`, the canonical order is

    init_parameter! → read_initial_def! → read_input_parameters!
    → sync_modified_parameter! → init_qp_weight!

`read_input_parameters!` runs **after** any `initial.def` overlay so
that selective `In*.def` values take precedence over wholesale
`initial.def` snapshots, matching C's behaviour. Callers (e.g.
`MVMCOptimizers.run_para_opt_from_namelist`) should respect this
order; calling out of sequence is undefined.

# Arguments
- `data::ExpertModeData`: data structure to update in place.
- `namelist_path::String`: path to `namelist.def`. Sibling `In*.def`
  paths listed inside are resolved relative to this file's directory.

Files are optional: missing `In*.def` entries (or files referenced by
`namelist.def` but absent on disk) are silently skipped.
"""
function read_input_parameters!(data::ExpertModeData, namelist_path::String)
    base_dir = dirname(abspath(namelist_path))

    # Read namelist.def to get file list
    namelist_content = read_def_file(namelist_path)
    file_list = parse_namelist_content(namelist_content)

    # Process InGutzwiller.def
    for (file_type, file_path) in file_list
        if file_type == "InGutzwiller"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                params = parse_input_parameter_file(full_path)
                # Update Gutzwiller terms
                # C uses 0-based indexing, Julia uses 1-based
                for (idx_c, value) in params
                    idx_julia = idx_c + 1  # Convert from 0-based to 1-based
                    if 1 <= idx_julia <= length(data.gutzwiller_terms)
                        data.gutzwiller_terms[idx_julia].value = value
                    end
                end
            end
        elseif file_type == "InJastrow"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                params = parse_input_parameter_file(full_path)
                # Update Jastrow terms
                # C uses 0-based indexing, Julia uses 1-based
                for (idx_c, value) in params
                    idx_julia = idx_c + 1  # Convert from 0-based to 1-based
                    if 1 <= idx_julia <= length(data.jastrow_terms)
                        data.jastrow_terms[idx_julia].value = value
                    end
                end
            end
        elseif file_type == "InDH2"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                params = parse_input_parameter_file(full_path)
                # Update DH2 terms (stored in doublon_holon_2site_terms)
                # Note: C implementation uses Proj array with offset
                # We need to map to the appropriate terms
                # For now, we'll skip this as it requires understanding the DH2 structure
                @warn "InDH2.def parsing not yet fully implemented"
            end
        elseif file_type == "InDH4"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                params = parse_input_parameter_file(full_path)
                # Update DH4 terms (stored in doublon_holon_4site_terms)
                @warn "InDH4.def parsing not yet fully implemented"
            end
        elseif file_type == "InOrbital" || file_type == "InOrbitalAntiParallel"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                params = parse_input_parameter_file(full_path)
                # Update Orbital terms (Slater parameters)
                # C uses 0-based indexing, Julia uses 1-based
                for (idx_c, value) in params
                    idx_julia = idx_c + 1  # Convert from 0-based to 1-based
                    if 1 <= idx_julia <= length(data.orbital_terms)
                        data.orbital_terms[idx_julia].value = value
                    end
                end
            end
        elseif file_type == "InOrbitalParallel"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                params = parse_input_parameter_file(full_path)
                # Update OrbitalParallel terms
                # Note: C implementation uses offset iNOrbitalAntiParallel
                # For now, we'll update orbital_terms directly
                for (idx, value) in params
                    # C: Slater[iNOrbitalAntiParallel + idx] = value
                    # We need to find the appropriate orbital term
                    # This is complex and depends on the orbital structure
                    @warn "InOrbitalParallel.def parsing not yet fully implemented"
                end
            end
        elseif file_type == "InOrbitalGeneral"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                params = parse_input_parameter_file(full_path)
                # Update OrbitalGeneral terms
                # C uses 0-based indexing, Julia uses 1-based
                for (idx_c, value) in params
                    idx_julia = idx_c + 1  # Convert from 0-based to 1-based
                    if 1 <= idx_julia <= length(data.orbital_terms)
                        data.orbital_terms[idx_julia].value = value
                    end
                end
            end
        elseif file_type == "InOptTrans"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                params = parse_input_parameter_file(full_path)
                # Update OptTrans parameters
                # Note: OptTrans is stored in data.para_qp_opt_trans or similar
                # We need to check the data structure
                @warn "InOptTrans.def parsing not yet fully implemented"
            end
        elseif startswith(file_type, "InChargeRBM_") ||
               startswith(file_type, "InSpinRBM_") ||
               startswith(file_type, "InGeneralRBM_")
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                params = parse_input_parameter_file(full_path)
                if file_type == "InChargeRBM_PhysLayer"
                    _set_rbm_terms_from_params!(data.charge_rbm_phys_layer_terms, params)
                elseif file_type == "InSpinRBM_PhysLayer"
                    _set_rbm_terms_from_params!(data.spin_rbm_phys_layer_terms, params)
                elseif file_type == "InGeneralRBM_PhysLayer"
                    _set_rbm_terms_from_params!(data.general_rbm_phys_layer_terms, params)
                elseif file_type == "InChargeRBM_HiddenLayer"
                    _set_rbm_terms_from_params!(data.charge_rbm_hidden_layer_terms, params)
                elseif file_type == "InSpinRBM_HiddenLayer"
                    _set_rbm_terms_from_params!(data.spin_rbm_hidden_layer_terms, params)
                elseif file_type == "InGeneralRBM_HiddenLayer"
                    _set_rbm_terms_from_params!(data.general_rbm_hidden_layer_terms, params)
                elseif file_type == "InChargeRBM_PhysHidden"
                    _set_rbm_terms_from_params!(data.charge_rbm_phys_hidden_terms, params)
                elseif file_type == "InSpinRBM_PhysHidden"
                    _set_rbm_terms_from_params!(data.spin_rbm_phys_hidden_terms, params)
                elseif file_type == "InGeneralRBM_PhysHidden"
                    _set_rbm_terms_from_params!(data.general_rbm_phys_hidden_terms, params)
                end
            end
        end
    end

    return data
end

"""
    ensure_optimization_flags_size!(data::ExpertModeData, n_para_full::Int)

Ensure optimization_flags has at least n_para_full elements, defaulting new entries to true.
"""
function ensure_optimization_flags_size!(data::ExpertModeData, n_para_full::Int)
    if isempty(data.optimization_flags)
        data.optimization_flags = fill(true, n_para_full)
    elseif length(data.optimization_flags) < n_para_full
        old_len = length(data.optimization_flags)
        resize!(data.optimization_flags, n_para_full)
        for i in (old_len + 1):n_para_full
            data.optimization_flags[i] = true
        end
    end
end

"""
    set_projection_opt_flags!(data::ExpertModeData, gutzwiller_flags::Dict{Int, Int},
                              jastrow_flags::Dict{Int, Int};
                              gutzwiller_is_complex::Bool=false,
                              jastrow_is_complex::Bool=false)

Set OptFlag for Gutzwiller and Jastrow parameters in data.optimization_flags.
C implementation:
- Gutzwiller: OptFlag[2*i] (real) and OptFlag[2*i+1] (imag, if complex)
- Jastrow: OptFlag[2*(NGutzwillerIdx+i)] (real) and OptFlag[2*(NGutzwillerIdx+i)+1] (imag, if complex)
"""
function set_projection_opt_flags!(
    data::ExpertModeData,
    gutzwiller_flags::Dict{Int, Int},
    jastrow_flags::Dict{Int, Int};
    gutzwiller_is_complex::Bool = false,
    jastrow_is_complex::Bool = false,
)
    if isempty(gutzwiller_flags) && isempty(jastrow_flags)
        return
    end

    n_gutz = max(data.n_gutzwiller_idx, length(data.gutzwiller_terms))
    n_jast = max(data.n_jastrow_idx, length(data.jastrow_terms))
    n_proj = n_gutz + n_jast

    # Ensure size for projection parameters (orbital parameters may be added later)
    ensure_optimization_flags_size!(data, 2 * n_proj)

    for (idx, opt_flag) in gutzwiller_flags
        fidx = idx
        if 2 * fidx + 1 <= length(data.optimization_flags)
            data.optimization_flags[2 * fidx + 1] = (opt_flag != 0)
            if gutzwiller_is_complex
                data.optimization_flags[2 * fidx + 2] = (opt_flag != 0)
            else
                data.optimization_flags[2 * fidx + 2] = false
            end
        end
    end

    for (idx, opt_flag) in jastrow_flags
        fidx = n_gutz + idx
        if 2 * fidx + 1 <= length(data.optimization_flags)
            data.optimization_flags[2 * fidx + 1] = (opt_flag != 0)
            if jastrow_is_complex
                data.optimization_flags[2 * fidx + 2] = (opt_flag != 0)
            else
                data.optimization_flags[2 * fidx + 2] = false
            end
        end
    end
end

"""
    set_orbital_opt_flags!(data::ExpertModeData, opt_flags::Dict{Int, Int})

Set OptFlag for orbital parameters in data.optimization_flags.
C implementation: OptFlag[2*fidx] (real) and OptFlag[2*fidx+1] (imag, if complex)
where fidx = NProj + FlagRBM * NRBM + orbital_idx.
"""
function set_orbital_opt_flags!(data::ExpertModeData, opt_flags::Dict{Int,Int})
    if isempty(opt_flags)
        return
    end

    # Calculate fidx offset: NProj + FlagRBM * NRBM
    n_proj = length(data.gutzwiller_terms) + length(data.jastrow_terms)
    n_rbm = count_rbm_parameters(data)
    flag_rbm = n_rbm > 0 ? 1 : 0
    fidx_offset = n_proj + flag_rbm * n_rbm

    # Get complex flag for orbital parameters
    is_complex = _get_all_complex_flag_local(data)

    # Initialize or resize optimization_flags
    n_para = n_proj + flag_rbm * n_rbm + data.modpara.n_orbital_idx
    n_para_full = 2 * n_para  # Always 2*NPara (real and imaginary)
    ensure_optimization_flags_size!(data, n_para_full)

    # Set OptFlag for each orbital parameter
    for (idx, opt_flag) in opt_flags
        # C: OptFlag[2*fidx] = opt_flag (real part)
        fidx = fidx_offset + idx
        if 2 * fidx + 1 <= length(data.optimization_flags)
            data.optimization_flags[2*fidx+1] = (opt_flag != 0)  # Real part (1-based)
            if is_complex
                # C: OptFlag[2*fidx+1] = opt_flag (imaginary part)
                data.optimization_flags[2*fidx+2] = (opt_flag != 0)  # Imaginary part (1-based)
            end
        end
    end
end

"""
    set_rbm_opt_flags!(data::ExpertModeData, opt_flags::Dict{Int, Int}, fidx_offset::Int; is_complex::Bool=false)

Set OptFlag for RBM parameters in data.optimization_flags.
C implementation uses contiguous RBM blocks in OptFlag with section-specific offsets.
"""
function set_rbm_opt_flags!(
    data::ExpertModeData,
    opt_flags::Dict{Int, Int},
    fidx_offset::Int;
    is_complex::Bool = false,
)
    isempty(opt_flags) && return

    n_proj = length(data.gutzwiller_terms) + length(data.jastrow_terms)
    n_rbm = count_rbm_parameters(data)
    flag_rbm = n_rbm > 0 ? 1 : 0
    n_para = n_proj + flag_rbm * n_rbm + data.modpara.n_orbital_idx
    ensure_optimization_flags_size!(data, 2 * n_para)

    for (idx, opt_flag) in opt_flags
        fidx = fidx_offset + idx
        if 2 * fidx + 1 <= length(data.optimization_flags)
            data.optimization_flags[2 * fidx + 1] = (opt_flag != 0)
            if is_complex
                data.optimization_flags[2 * fidx + 2] = (opt_flag != 0)
            else
                data.optimization_flags[2 * fidx + 2] = false
            end
        end
    end
end

# Helper function to get all_complex flag (local to this module)
function _get_all_complex_flag_local(data::ExpertModeData)::Bool
    # Check if orbital parameters are complex
    if !isempty(data.orbital_terms)
        return any(t -> t.is_complex, data.orbital_terms)
    end
    return false
end
