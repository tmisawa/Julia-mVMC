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

function _parse_float_strict(tok::AbstractString, filepath::String, line_num::Int, field::String)::Float64
    value = tryparse(Float64, tok)
    value === nothing && error("$filepath:$line_num: invalid $field '$tok'")
    isfinite(value) || error("$filepath:$line_num: non-finite $field '$tok'")
    return value
end

function _parse_int_strict(tok::AbstractString, filepath::String, line_num::Int, field::String)::Int
    value = tryparse(Int, tok)
    value === nothing && error("$filepath:$line_num: invalid $field '$tok'")
    return value
end

function parse_indexed_input_parameter_file_strict(
    filepath::String,
    expected_header_count::Int,
    expected_param_count::Int,
    label::AbstractString,
)::Vector{ComplexF64}
    content = read_def_file(filepath)
    lines = split(content, '\n')
    length(lines) >= 5 || error("$label: $filepath must include 5 header lines")

    header_tokens = split_def_line(clean_line(lines[2]))
    length(header_tokens) >= 2 || error("$label: $filepath:2 missing count header")
    header_count = _parse_int_strict(header_tokens[2], filepath, 2, "count")
    header_count == expected_header_count ||
        error("$label: header count mismatch in $filepath: got $header_count, expected $expected_header_count")

    params = fill(0.0 + 0.0im, expected_param_count)
    seen = falses(expected_param_count)
    row_count = 0

    for line_num = 6:length(lines)
        line = clean_line(lines[line_num])
        isempty(line) && continue
        tokens = split_def_line(line)
        length(tokens) == 3 || error("$label: $filepath:$line_num expected 'idx real imag'")
        row_count += 1

        idx = _parse_int_strict(tokens[1], filepath, line_num, "index")
        0 <= idx < expected_param_count ||
            error("$label: index $idx out of range [0, $(expected_param_count - 1)] in $filepath:$line_num")
        !seen[idx+1] || error("$label: duplicated index $idx in $filepath:$line_num")
        real_val = _parse_float_strict(tokens[2], filepath, line_num, "real value")
        imag_val = _parse_float_strict(tokens[3], filepath, line_num, "imag value")

        seen[idx+1] = true
        params[idx+1] = ComplexF64(real_val, imag_val)
    end

    row_count == expected_param_count ||
        error("$label: row count mismatch in $filepath: got $row_count, expected $expected_param_count")
    missing = findfirst(x -> !x, seen)
    missing === nothing || error("$label: missing index $(missing - 1) in $filepath")

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

function count_orbital_parameters(data::ExpertModeData)::Int
    if data.modpara.n_orbital_idx > 0
        return data.modpara.n_orbital_idx
    elseif !isempty(data.orbital_terms)
        return maximum(t.idx for t in data.orbital_terms) + 1
    else
        return 0
    end
end

count_opt_trans_parameters(data::ExpertModeData)::Int = length(data.opt_trans)

function count_variational_parameters(data::ExpertModeData)::Int
    n_rbm = count_rbm_parameters(data)
    return projection_layout(data).n_proj +
           n_rbm +
           count_orbital_parameters(data) +
           count_opt_trans_parameters(data)
end

function _set_rbm_terms_from_params!(terms, params::Dict{Int, ComplexF64})
    isempty(params) && return
    for term in terms
        if haskey(params, term.idx)
            term.value = params[term.idx]
        end
    end
end

function _set_orbital_terms_from_params!(
    data::ExpertModeData,
    params::Dict{Int,ComplexF64},
)
    isempty(params) && return
    for term in data.orbital_terms
        if haskey(params, term.idx)
            term.value = params[term.idx]
        end
    end
    return nothing
end

function _orbital_parallel_offset(data::ExpertModeData)::Int
    # No OrbitalParallel block: InOrbital/InOrbitalAntiParallel cover the whole
    # Slater array, so the overlay window spans all orbital parameters.
    data.i_flg_orbital_parallel == 1 || return count_orbital_parameters(data)
    # Parallel-only (no anti-parallel): the parallel block starts at index 0.
    data.i_flg_orbital_anti_parallel == 1 || return 0
    # Both present: the parallel block begins exactly at NArrayAP, the number of
    # anti-parallel parameters recorded at parse time. This mirrors C readdef.c,
    # `Slater[iNOrbitalAntiParallel + idx]`, instead of inferring the boundary
    # heuristically from consecutive index pairs.
    return data.n_orbital_anti_parallel
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
| `InOrbital` / `InOrbitalAntiParallel` | `data.orbital_terms` matching `.idx` |
| `InOrbitalParallel` | `data.orbital_terms` matching `.idx - offset` |
| `InOrbitalGeneral` | `data.orbital_terms` matching `.idx` (fsz layout) |
| `InOptTrans` | `data.opt_trans` |
| `InDH2` | `data.doublon_holon_2site_params` |
| `InDH4` | `data.doublon_holon_4site_params` |
| `InChargeRBM_PhysLayer` / `_HiddenLayer` / `_PhysHidden` | `charge_rbm_*_terms` |
| `InSpinRBM_PhysLayer` / `_HiddenLayer` / `_PhysHidden` | `spin_rbm_*_terms` |
| `InGeneralRBM_PhysLayer` / `_HiddenLayer` / `_PhysHidden` | `general_rbm_*_terms` |

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
                layout = projection_layout(data)
                expected = 6 * layout.n_dh2
                length(data.doublon_holon_2site_params) == expected ||
                    error("InDH2 target parameter length mismatch: got $(length(data.doublon_holon_2site_params)), expected $expected")
                params = parse_indexed_input_parameter_file_strict(
                    full_path,
                    layout.n_dh2,
                    expected,
                    "InDH2",
                )
                copyto!(data.doublon_holon_2site_params, params)
            end
        elseif file_type == "InDH4"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                layout = projection_layout(data)
                expected = 10 * layout.n_dh4
                length(data.doublon_holon_4site_params) == expected ||
                    error("InDH4 target parameter length mismatch: got $(length(data.doublon_holon_4site_params)), expected $expected")
                params = parse_indexed_input_parameter_file_strict(
                    full_path,
                    layout.n_dh4,
                    expected,
                    "InDH4",
                )
                copyto!(data.doublon_holon_4site_params, params)
            end
        elseif file_type == "InOrbital" || file_type == "InOrbitalAntiParallel"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                params = parse_input_parameter_file(full_path)
                _set_orbital_terms_from_params!(data, params)
            end
        elseif file_type == "InOrbitalParallel"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                offset = _orbital_parallel_offset(data)
                n_orbital = count_orbital_parameters(data)
                expected = n_orbital - offset
                expected > 0 ||
                    error("InOrbitalParallel target parameter length mismatch: got $expected from NOrbitalIdx=$n_orbital and offset=$offset")
                params = parse_indexed_input_parameter_file_strict(
                    full_path,
                    expected,
                    expected,
                    "InOrbitalParallel",
                )
                for term in data.orbital_terms
                    rel_idx = term.idx - offset
                    if 0 <= rel_idx < expected
                        term.value = params[rel_idx+1]
                    end
                end
            end
        elseif file_type == "InOrbitalGeneral"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                params = parse_input_parameter_file(full_path)
                _set_orbital_terms_from_params!(data, params)
            end
        elseif file_type == "InOptTrans"
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                expected = count_opt_trans_parameters(data)
                expected > 0 ||
                    error("InOptTrans target parameter length mismatch: OptTrans is not active")
                params = parse_indexed_input_parameter_file_strict(
                    full_path,
                    expected,
                    expected,
                    "InOptTrans",
                )
                copyto!(data.opt_trans, params)
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

    layout = projection_layout(data)
    n_gutz = layout.n_gutzwiller
    n_proj = layout.n_proj

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
    n_proj = projection_layout(data).n_proj
    n_rbm = count_rbm_parameters(data)
    flag_rbm = n_rbm > 0 ? 1 : 0
    fidx_offset = n_proj + flag_rbm * n_rbm

    # Get complex flag for orbital parameters
    is_complex = _get_all_complex_flag_local(data)

    # Initialize or resize optimization_flags
    n_para =
        n_proj +
        flag_rbm * n_rbm +
        count_orbital_parameters(data) +
        count_opt_trans_parameters(data)
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

function set_opt_trans_opt_flags!(data::ExpertModeData)
    n_opt_trans = count_opt_trans_parameters(data)
    n_opt_trans == 0 && return

    n_proj = projection_layout(data).n_proj
    n_rbm = count_rbm_parameters(data)
    n_orbital = count_orbital_parameters(data)
    n_para = n_proj + n_rbm + n_orbital + n_opt_trans
    ensure_optimization_flags_size!(data, 2 * n_para)

    fidx_offset = n_proj + n_rbm + n_orbital
    for i = 0:(n_opt_trans-1)
        real_idx = 2 * (fidx_offset + i) + 1
        imag_idx = real_idx + 1
        data.optimization_flags[real_idx] = true
        if imag_idx <= length(data.optimization_flags)
            data.optimization_flags[imag_idx] = false
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

    n_proj = projection_layout(data).n_proj
    n_rbm = count_rbm_parameters(data)
    flag_rbm = n_rbm > 0 ? 1 : 0
    n_para =
        n_proj +
        flag_rbm * n_rbm +
        count_orbital_parameters(data) +
        count_opt_trans_parameters(data)
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
    return (
        data.modpara.complex_flag != 0 ||
        any(t -> t.is_complex, data.gutzwiller_terms) ||
        any(t -> t.is_complex, data.jastrow_terms) ||
        data.doublon_holon_2site_complex ||
        data.doublon_holon_4site_complex ||
        any(t -> t.is_complex, data.orbital_terms)
    )
end
