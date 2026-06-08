"""
MVMCExpertModeParsers.jl

A Julia package for parsing and generating mVMC Expert Mode definition files.

This package provides comprehensive functionality for:
- Parsing .def files from mVMC Expert Mode
- Generating .def files for mVMC Expert Mode
- Validating parameter consistency
- Converting between different data formats

Based on the C reference implementation in mVMC.
"""

module MVMCExpertModeParsers

# Standard library imports
using Printf
using Dates
using Random
using Printf
using SFMT

# Include type definitions
include("types/expert_types.jl")

# Include utilities
include("utils/constants.jl")
include("utils/file_utils.jl")
include("utils/validation.jl")
include("utils/parameter_init.jl")
include("utils/read_input_parameters.jl")
include("utils/qp_weight.jl")
include("utils/opt_flag_utils.jl")
include("utils/orbital_qptrans_utils.jl")

# Include parsers
include("parsers/modpara_parser.jl")
include("parsers/locspin_parser.jl")
include("parsers/trans_parser.jl")
include("parsers/coulomb_parser.jl")
include("parsers/hund_parser.jl")
include("parsers/exchange_parser.jl")
include("parsers/pairhop_parser.jl")
include("parsers/interall_parser.jl")
include("parsers/gutzwiller_parser.jl")
include("parsers/jastrow_parser.jl")
include("parsers/orbital_parser.jl")
include("parsers/green_parser.jl")
include("parsers/qptrans_parser.jl")
include("parsers/rbm_parser.jl")
include("parsers/doublon_holon_parser.jl")

# Export main parsing function
export parse_expert_mode_files

# Export quantum projection weight functions
export init_qp_weight!
export update_qp_weight!

function _is_required_if_present_file_type(file_type::String)::Bool
    return (
        file_type == "TwoBodyGEx" ||
        file_type == "DH2" ||
        file_type == "DH4" ||
        file_type == "DoublonHolon2Site" ||
        file_type == "DoublonHolon4Site"
    )
end

function _is_unsupported_projection_file_type(file_type::String)::Bool
    return file_type == "SpinJastrow"
end

"""
    parse_expert_mode_files(namelist_path::String) -> ExpertModeData

Parse all Expert Mode files referenced by a `namelist.def` file. The
base directory for relative `.def` paths is derived from
`dirname(abspath(namelist_path))` automatically.

# Arguments
- `namelist_path`: Path to the `namelist.def` file (relative or absolute).

# Returns
- `ExpertModeData`: Complete Expert Mode data structure.

# Example
```julia
data = parse_expert_mode_files("inputs/namelist.def")
```
"""
function parse_expert_mode_files(namelist_path::String)::ExpertModeData
    base_dir = dirname(abspath(namelist_path))

    data = ExpertModeData()
    errors = String[]
    warnings = String[]
    # Some inputs are required-if-present: record a fatal failure and rethrow
    # AFTER the outer try/catch below, which would otherwise swallow an in-loop
    # error and return a data object with a missing required block.
    required_file_fatal = nothing

    try
        # Read namelist.def to get file list
        namelist_content = read_def_file(namelist_path)
        file_list = parse_namelist_content(namelist_content)

        # C reads every orbital-block count in a first pass, so OrbitalParallel may
        # appear in any namelist position. This single-pass parser instead builds
        # the parallel block on top of the already-parsed anti-parallel block (whose
        # NArrayAP offsets the parallel indices), so the anti-parallel file must be
        # listed first. Enforce that explicitly rather than emit a wrong layout.
        required_file_fatal = _orbital_file_order_error(file_list)

        # Parse each file (continue on errors, like C implementation)
        for (file_type, file_path) in file_list
            required_file_fatal === nothing || break
            full_path = joinpath(base_dir, file_path)

            if _is_unsupported_projection_file_type(file_type)
                required_file_fatal =
                    "$file_type inputs are not supported yet; SpinJastrow must hard-fail because projection layout would otherwise be wrong"
                break
            end

            if !validate_file_exists(full_path)
                if _is_required_if_present_file_type(file_type)
                    required_file_fatal = "Required $file_type file not found: $full_path"
                    break   # stop parsing; post-loop work is guarded below
                end
                warning_msg = "File not found: $full_path"
                push!(warnings, warning_msg)
                @warn warning_msg
                continue
            end

            try
                parse_file_by_type!(data, file_type, full_path)
            catch e
                if _is_required_if_present_file_type(file_type)
                    required_file_fatal = "Error parsing required $file_type file $full_path: $e"
                    break   # stop parsing; post-loop work is guarded below
                end
                error_msg = "Error parsing $file_type file $full_path: $e"
                push!(warnings, error_msg)  # Treat as warning, not error
                @warn error_msg
                # Continue processing other files
            end
        end

        # Skip all post-loop processing when a required-if-present file failed:
        # the run is about to error, so avoid spurious warnings / side effects.
        if required_file_fatal === nothing
            # RBM OptFlag is defined in RBM idx files and uses C-specific block offsets.
            # Apply projection-family-dependent flags after all files are parsed,
            # independent of namelist order.
            set_dh_opt_flags!(data)
            apply_rbm_opt_flags_from_files!(data, file_list, base_dir)
            apply_orbital_opt_flags_from_files!(data, file_list, base_dir)
            set_opt_trans_opt_flags!(data)

            # Judge orbital mode (equivalent to C's JudgeOrbitalMode)
            judge_orbital_mode!(data)

            # Build orbital matrices AFTER judge_orbital_mode! so i_flg_orbital_general is correctly set
            # This ensures FSZ mode uses 2*Nsite x 2*Nsite matrices
            if !isempty(data.orbital_terms)
                build_orbital_sgn_matrix!(data)
            end

            # Calculate Ne from NLocSpin and NCond if NCond != -1 (C code: readdef.c:593)
            # This matches C implementation's behavior in ReadDefFileNInt
            if data.modpara.ncond != -1
                if data.modpara.ncond % 2 != 0
                    push!(warnings, "NCond must be even, got $(data.modpara.ncond)")
                elseif data.modpara.nelec == 0
                    # Calculate Ne = (NLocSpin + NCond) / 2
                    data.modpara.nelec = (data.modpara.nlocspin + data.modpara.ncond) ÷ 2
                end
            end

            # Read input parameter files (InGutzwiller.def, InJastrow.def, InOrbital.def, etc.)
            # This updates parameter values in gutzwiller_terms, jastrow_terms, orbital_terms
            # Equivalent to C's ReadInputParameters()
            read_input_parameters!(data, namelist_path)

            # Print summary like C implementation
            if !isempty(warnings)
                @info "Parsing completed with $(length(warnings)) warnings"
            else
                @info "Parsing completed successfully"
            end
        end  # if required_file_fatal === nothing

    catch e
        error_msg = "Critical error in parse_expert_mode_files: $e"
        push!(errors, error_msg)
        @error error_msg
    end

    # Required-if-present failures are fatal on the public path. Throw here,
    # outside the outer try/catch, so the error is not swallowed.
    if required_file_fatal !== nothing
        error(required_file_fatal)
    end

    return data
end

"""
    _orbital_file_order_error(file_list) -> Union{String, Nothing}

Return an error message if `OrbitalParallel` is listed before `Orbital` /
`OrbitalAntiParallel` in the namelist, otherwise `nothing`. The single-pass
parser offsets the parallel orbital block by the anti-parallel `NArrayAP`, so the
anti-parallel block must be parsed first. A pure-parallel system (no anti-parallel
file) is valid (`NArrayAP == 0`).
"""
function _orbital_file_order_error(file_list::Vector{Tuple{String,String}})
    parallel_pos = nothing
    anti_pos = nothing
    for (i, (file_type, _)) in enumerate(file_list)
        if file_type == "OrbitalParallel" && parallel_pos === nothing
            parallel_pos = i
        elseif (file_type == "Orbital" || file_type == "OrbitalAntiParallel") &&
               anti_pos === nothing
            anti_pos = i
        end
    end
    if parallel_pos !== nothing && anti_pos !== nothing && parallel_pos < anti_pos
        return "OrbitalParallel must be listed after Orbital/OrbitalAntiParallel in " *
               "namelist.def: the anti-parallel block defines the NArrayAP offset " *
               "for the parallel orbital block"
    end
    return nothing
end

function apply_orbital_opt_flags_from_files!(
    data::ExpertModeData,
    file_list::Vector{Tuple{String,String}},
    base_dir::String,
)
    orbital_opt_flags = Dict{Int,Int}()
    n_orbital_anti_parallel = 0

    for (file_type, file_path) in file_list
        if !(
            file_type == "Orbital" ||
            file_type == "OrbitalAntiParallel" ||
            file_type == "OrbitalParallel" ||
            file_type == "OrbitalGeneral"
        )
            continue
        end

        full_path = joinpath(base_dir, file_path)
        validate_file_exists(full_path) || continue

        result, opt_flags, header_count = parse_orbital_def(full_path)
        result.success || continue

        # Record NArrayAP from the header (C's iNOrbitalAntiParallel) before the
        # opt-flag emptiness check, so the parallel opt-flag offset below is
        # correct even when the anti-parallel file carries no opt flags.
        if file_type == "Orbital" || file_type == "OrbitalAntiParallel"
            n_orbital_anti_parallel = header_count
        end

        isempty(opt_flags) && continue

        if file_type == "Orbital" || file_type == "OrbitalAntiParallel"
            for (idx, opt) in opt_flags
                orbital_opt_flags[idx] = opt
            end
        elseif file_type == "OrbitalParallel"
            for (idx, opt) in opt_flags
                orbital_opt_flags[n_orbital_anti_parallel+2*idx] = opt
                orbital_opt_flags[n_orbital_anti_parallel+2*idx+1] = opt
            end
        elseif file_type == "OrbitalGeneral"
            for (idx, opt) in opt_flags
                orbital_opt_flags[idx] = opt
            end
        end
    end

    set_orbital_opt_flags!(data, orbital_opt_flags)
end

function apply_rbm_opt_flags_from_files!(
    data::ExpertModeData,
    file_list::Vector{Tuple{String, String}},
    base_dir::String,
)
    rbm_files = Dict{String, String}()
    for (file_type, file_path) in file_list
        if startswith(file_type, "ChargeRBM_") || startswith(file_type, "SpinRBM_") || startswith(file_type, "GeneralRBM_")
            full_path = joinpath(base_dir, file_path)
            if validate_file_exists(full_path)
                rbm_files[file_type] = full_path
            end
        end
    end
    isempty(rbm_files) && return

    rbm_opt = Dict{String, Tuple{Dict{Int, Int}, Bool}}()
    for (file_type, full_path) in rbm_files
        content = read_def_file(full_path)
        if file_type == "ChargeRBM_PhysLayer"
            r = parse_charge_rbm_phys_layer_content_extended(content)
            rbm_opt[file_type] = (r.opt_flags, r.is_complex_flag)
        elseif file_type == "SpinRBM_PhysLayer"
            r = parse_spin_rbm_phys_layer_content_extended(content)
            rbm_opt[file_type] = (r.opt_flags, r.is_complex_flag)
        elseif file_type == "GeneralRBM_PhysLayer"
            r = parse_general_rbm_phys_layer_content_extended(content)
            rbm_opt[file_type] = (r.opt_flags, r.is_complex_flag)
        elseif file_type == "ChargeRBM_HiddenLayer"
            r = parse_charge_rbm_hidden_layer_content_extended(content)
            rbm_opt[file_type] = (r.opt_flags, r.is_complex_flag)
        elseif file_type == "SpinRBM_HiddenLayer"
            r = parse_spin_rbm_hidden_layer_content_extended(content)
            rbm_opt[file_type] = (r.opt_flags, r.is_complex_flag)
        elseif file_type == "GeneralRBM_HiddenLayer"
            r = parse_general_rbm_hidden_layer_content_extended(content)
            rbm_opt[file_type] = (r.opt_flags, r.is_complex_flag)
        elseif file_type == "ChargeRBM_PhysHidden"
            r = parse_charge_rbm_phys_hidden_content_extended(content)
            rbm_opt[file_type] = (r.opt_flags, r.is_complex_flag)
        elseif file_type == "SpinRBM_PhysHidden"
            r = parse_spin_rbm_phys_hidden_content_extended(content)
            rbm_opt[file_type] = (r.opt_flags, r.is_complex_flag)
        elseif file_type == "GeneralRBM_PhysHidden"
            r = parse_general_rbm_phys_hidden_content_extended(content)
            rbm_opt[file_type] = (r.opt_flags, r.is_complex_flag)
        end
    end

    n_proj = projection_layout(data).n_proj
    n_charge_phys = _rbm_section_nparam(data.charge_rbm_phys_layer_terms)
    n_spin_phys = _rbm_section_nparam(data.spin_rbm_phys_layer_terms)
    n_general_phys = _rbm_section_nparam(data.general_rbm_phys_layer_terms)
    n_charge_hidden = _rbm_section_nparam(data.charge_rbm_hidden_layer_terms)
    n_spin_hidden = _rbm_section_nparam(data.spin_rbm_hidden_layer_terms)
    n_general_hidden = _rbm_section_nparam(data.general_rbm_hidden_layer_terms)
    n_charge_physhidden = _rbm_section_nparam(data.charge_rbm_phys_hidden_terms)
    n_spin_physhidden = _rbm_section_nparam(data.spin_rbm_phys_hidden_terms)
    n_general_physhidden = _rbm_section_nparam(data.general_rbm_phys_hidden_terms)

    phys_offset = n_proj
    hidden_offset = n_proj + n_charge_phys + n_spin_phys + n_general_phys
    physhidden_offset = hidden_offset + n_charge_hidden + n_spin_hidden + n_general_hidden

    sections = (
        ("ChargeRBM_PhysLayer", phys_offset, n_charge_phys),
        ("SpinRBM_PhysLayer", phys_offset + n_charge_phys, n_spin_phys),
        ("GeneralRBM_PhysLayer", phys_offset + n_charge_phys + n_spin_phys, n_general_phys),
        ("ChargeRBM_HiddenLayer", hidden_offset, n_charge_hidden),
        ("SpinRBM_HiddenLayer", hidden_offset + n_charge_hidden, n_spin_hidden),
        ("GeneralRBM_HiddenLayer", hidden_offset + n_charge_hidden + n_spin_hidden, n_general_hidden),
        ("ChargeRBM_PhysHidden", physhidden_offset, n_charge_physhidden),
        ("SpinRBM_PhysHidden", physhidden_offset + n_charge_physhidden, n_spin_physhidden),
        ("GeneralRBM_PhysHidden", physhidden_offset + n_charge_physhidden + n_spin_physhidden, n_general_physhidden),
    )

    for (section_name, offset, n_section) in sections
        n_section == 0 && continue
        opt_and_flag = get(rbm_opt, section_name, (Dict{Int, Int}(), false))
        set_rbm_opt_flags!(data, opt_and_flag[1], offset; is_complex = opt_and_flag[2])
    end
end

"""
    judge_orbital_mode!(data::ExpertModeData)

Judge orbital mode based on which orbital files are present.
Equivalent to C's JudgeOrbitalMode() function.

Sets i_flg_orbital_general based on:
- If OrbitalGeneral file exists: i_flg_orbital_general = 1
- If both OrbitalAntiParallel and OrbitalParallel exist: i_flg_orbital_general = 1
- Otherwise: i_flg_orbital_general = 0 (sz conserved)
"""
function judge_orbital_mode!(data::ExpertModeData)
    # C implementation logic:
    # if (iFlgOrbitalGeneral == 1) {
    #   if (iFlgOrbitalAP == 0 && iFlgOrbitalP == 0) { //(1, 0, 0)
    #     iret = 0;
    #   } else { //(1, 1, 0) or (1, 0, 1) or (1, 1, 1)
    #     iret = -1;  // Error: multiple definition
    #   }
    # } else {
    #   if (iFlgOrbitalAP == 1) {
    #     if (iFlgOrbitalP == 1) { //(0, 1, 1)
    #       iFlgOrbitalGeneral = 1;
    #       iret = 0;
    #     } else { //(0, 1, 0)
    #       iret = 0;
    #     }
    #   } else { // (0, 0, 0) or (0, 0, 1)
    #     iret = -2;  // Error: no orbital file
    #   }
    # }

    if data.i_flg_orbital_general == 1
        # OrbitalGeneral is already set
        if data.i_flg_orbital_anti_parallel == 1 || data.i_flg_orbital_parallel == 1
            @warn "Multiple definition of Orbital files: OrbitalGeneral conflicts with OrbitalAntiParallel/OrbitalParallel"
        end
    else
        # Check if both OrbitalAntiParallel and OrbitalParallel exist
        if data.i_flg_orbital_anti_parallel == 1 && data.i_flg_orbital_parallel == 1
            # Both exist: set to general mode
            data.i_flg_orbital_general = 1
        elseif data.i_flg_orbital_anti_parallel == 0 && data.i_flg_orbital_parallel == 0
            # No orbital files: this is an error in C, but we'll just keep it as 0
            # (sz conserved mode, but no orbital terms)
        end
    end
end

"""
    parse_file_by_type!(data::ExpertModeData, file_type::String, file_path::String)

Parse a specific file type and update the ExpertModeData.
"""
function parse_file_by_type!(data::ExpertModeData, file_type::String, file_path::String)
    if file_type == "ModPara"
        result = parse_modpara_def(file_path)
        if result.success
            data.modpara = result.data
        end
    elseif file_type == "LocSpin"
        # Read NLocSpin from the first line of locspn.def (C code: ReadBuffInt)
        # Format: first line after separator contains "NlocalSpin N"
        content = read_def_file(file_path)
        lines = split(content, '\n')
        # Look for "NlocalSpin" or "NLocSpin" in the first few lines
        for line in lines[1:min(5, length(lines))]
            clean_line_str = clean_line(line)
            if isempty(clean_line_str)
                continue
            end
            tokens = split_def_line(clean_line_str)
            if length(tokens) >= 2 && (tokens[1] == "NlocalSpin" || tokens[1] == "NLocSpin")
                n_locspin = safe_parse_int(tokens[2], 0)
                data.modpara.nlocspin = n_locspin
                break
            end
        end

        result = parse_locspin_def(file_path)
        if result.success
            data.locspin_terms = result.data
        end
    elseif file_type == "Trans"
        result = parse_trans_def(file_path)
        if result.success
            data.transfer_terms = result.data
        end
    elseif file_type == "CoulombIntra"
        result = parse_coulomb_intra_def(file_path)
        if result.success
            data.coulomb_intra_terms = result.data
        end
    elseif file_type == "CoulombInter"
        result = parse_coulomb_inter_def(file_path)
        if result.success
            data.coulomb_inter_terms = result.data
        end
    elseif file_type == "Hund"
        result = parse_hund_def(file_path)
        if result.success
            data.hund_terms = result.data
        end
    elseif file_type == "Exchange"
        result = parse_exchange_def(file_path)
        if result.success
            data.exchange_terms = result.data
        end
    elseif file_type == "PairHop"
        result = parse_pairhop_def(file_path)
        if result.success
            data.pair_hop_terms = result.data
        end
    elseif file_type == "InterAll"
        result = parse_interall_def(file_path)
        if result.success
            data.inter_all_terms = result.data
        end
    elseif file_type == "Gutzwiller"
        # Use extended parser to get site->idx mapping
        content = read_def_file(file_path)
        extended_result = parse_gutzwiller_content_extended(content)
        if extended_result.success
            data.gutzwiller_terms = extended_result.terms
            data.n_gutzwiller_idx = extended_result.n_gutzwiller_idx
            # Build gutzwiller_idx array from site_idx_map
            if !isempty(extended_result.site_idx_map)
                max_site = maximum(keys(extended_result.site_idx_map))
                data.gutzwiller_idx = zeros(Int, max_site + 1)
                for (site, idx) in extended_result.site_idx_map
                    data.gutzwiller_idx[site+1] = idx  # 1-based indexing for Julia array
                end
            end
            set_projection_opt_flags!(
                data,
                extended_result.opt_flags,
                Dict{Int, Int}();
                gutzwiller_is_complex = extended_result.is_complex_flag,
            )
        end
    elseif file_type == "Jastrow"
        content = read_def_file(file_path)
        extended_result = parse_jastrow_content_extended(content)
        if extended_result.success
            data.jastrow_terms = extended_result.terms
            data.n_jastrow_idx = extended_result.n_jastrow_idx
            # Build JastrowIdx matrix for C-compatible projection calculations
            nsite = data.modpara.nsite
            if nsite > 0
                jastrow_idx_matrix, n_jastrow_idx =
                    build_jastrow_idx_matrix(file_path, nsite)
                data.jastrow_idx = jastrow_idx_matrix
                data.n_jastrow_idx = n_jastrow_idx
            end
            set_projection_opt_flags!(
                data,
                Dict{Int, Int}(),
                extended_result.opt_flags;
                jastrow_is_complex = extended_result.is_complex_flag,
            )
        end
    elseif file_type == "Orbital"
        result, opt_flags, anti_count = parse_orbital_def(file_path)
        if result.success
            data.orbital_terms = result.data
            # C implementation: KWOrbital sets iFlgOrbitalAntiParallel = 1
            data.i_flg_orbital_anti_parallel = 1
            # NOTE: build_orbital_sgn_matrix! is called after judge_orbital_mode!
            # to ensure correct matrix size for FSZ mode.
            # Use the header declared count (C's iNOrbitalAntiParallel), not
            # max(idx)+1, so unreferenced trailing parameters still reserve slots.
            if anti_count > 0
                data.modpara.n_orbital_idx = anti_count
                data.n_orbital_anti_parallel = anti_count  # NArrayAP
            end
        end
    elseif file_type == "OrbitalAntiParallel"
        result, opt_flags, anti_count = parse_orbital_def(file_path)
        if result.success
            data.orbital_terms = result.data
            data.i_flg_orbital_anti_parallel = 1
            # NOTE: build_orbital_sgn_matrix! is called after judge_orbital_mode!
            # Use the header declared count (C's iNOrbitalAntiParallel).
            if anti_count > 0
                data.modpara.n_orbital_idx = anti_count
                data.n_orbital_anti_parallel = anti_count  # NArrayAP
            end
        end
    elseif file_type == "OrbitalParallel"
        result, opt_flags, parallel_count = parse_orbital_def(file_path)
        if result.success
            # C implementation: OrbitalParallel indices are interleaved for up-up and down-down
            # See readdef.c GetInfoOrbitalParallel (lines 2608-2620):
            #   for (spn_i = 0; spn_i < 2; spn_i++) {
            #     fij = NArrayAP + 2 * fij_org + spn_i;
            #     Array[all_i][all_j] = fij;
            #   }
            # So indices are: [NArrayAP + 2*0 = up0, NArrayAP + 2*0 + 1 = down0,
            #                  NArrayAP + 2*1 = up1, NArrayAP + 2*1 + 1 = down1, ...]

            # NArrayAP is the anti-parallel parameter count recorded by the
            # preceding Orbital/OrbitalAntiParallel parse (0 for a pure-parallel
            # system). C reads OrbitalParallel after the anti-parallel block;
            # the namelist order is enforced by parse_expert_mode_files.
            n_orbital_anti_parallel = data.n_orbital_anti_parallel

            # Create interleaved parallel orbital terms matching C implementation
            for term in result.data
                fij_org = term.idx  # Original index from the file (0-based)

                # up-up orbital term: fij = NArrayAP + 2 * fij_org
                offset_idx_up = n_orbital_anti_parallel + 2 * fij_org
                term_up = OrbitalTerm(
                    term.site1,
                    term.site2,
                    offset_idx_up,
                    term.value,
                    term.is_complex,
                    term.sign,
                )
                push!(data.orbital_terms, term_up)

                # down-down orbital term: fij = NArrayAP + 2 * fij_org + 1
                offset_idx_down = n_orbital_anti_parallel + 2 * fij_org + 1
                term_down = OrbitalTerm(
                    term.site1,
                    term.site2,
                    offset_idx_down,
                    term.value,
                    term.is_complex,
                    term.sign,
                )
                push!(data.orbital_terms, term_down)
            end

            data.i_flg_orbital_parallel = 1
            # Total parameter count = NArrayAP + 2*iNOrbitalParallel (header count),
            # matching C readdef.c (bufInt[IdxNOrbit] += 2*iNOrbitalParallel).
            data.modpara.n_orbital_idx = n_orbital_anti_parallel + 2 * parallel_count
        end
    elseif file_type == "OrbitalGeneral"
        result, opt_flags, general_count = parse_orbital_def(file_path)
        if result.success
            data.orbital_terms = result.data
            data.i_flg_orbital_general = 1
            # NOTE: build_orbital_sgn_matrix! is called after judge_orbital_mode!
            # Use the header declared count (C's NOrbitalIdx for OrbitalGeneral).
            if general_count > 0
                data.modpara.n_orbital_idx = general_count
            end
        end
    elseif file_type == "OneBodyG"
        result = parse_green_one_def(file_path)
        if result.success
            data.green_one_terms = result.data
        end
    elseif file_type == "TwoBodyG"
        result = parse_green_two_def(file_path)
        if result.success
            data.green_two_terms = result.data
        end
    elseif file_type == "TwoBodyGEx"
        result = parse_green_two_ex_def(file_path)
        if result.success
            data.green_two_ex_terms = result.data
        else
            # Fatal: a present-but-unparseable factored definition must never
            # degrade silently into an empty list (spec Finding 4).
            error("Failed to parse TwoBodyGEx file '$file_path': $(result.error_message)")
        end
    elseif file_type == "QPTrans"
        result = parse_qptrans_def(file_path)
        if result.success
            data.qptrans_terms = result.data
        end
    elseif file_type == "TransSym"
        # qptransidx.def format: first NQPTrans lines are idx value (ParaQPTrans),
        # then Nsite * NQPTrans lines are i j itmp itmpsgn (QPTrans indices)
        result = parse_qptransidx_def(file_path)
        if result.success
            data.qptrans_terms = result.data
            # Also read ParaQPTrans and NQPTrans from qptransidx.def
            content = read_def_file(file_path)
            lines = split(content, '\n')
            IGNORE_LINES_IN_DEF = 5
            if length(lines) > IGNORE_LINES_IN_DEF
                # Read NQPTrans from header (line 2: "NQPTrans N")
                n_qp_trans = 0
                if length(lines) > 1
                    header_line = clean_line(lines[2])
                    tokens = split_def_line(header_line)
                    if length(tokens) >= 2 && tokens[1] == "NQPTrans"
                        n_qp_trans = safe_parse_int(tokens[2], 0)
                    end
                end
                data.n_qp_trans = n_qp_trans

                # Read ParaQPTrans values (first NQPTrans lines after header)
                para_qp_trans = ComplexF64[]
                line_idx = IGNORE_LINES_IN_DEF + 1
                for i = 1:n_qp_trans
                    while line_idx <= length(lines)
                        line = clean_line(lines[line_idx])
                        if isempty(line)
                            line_idx += 1
                            continue
                        end
                        tokens = split_def_line(line)
                        if length(tokens) >= 2
                            idx = safe_parse_int(tokens[1], -1)
                            value = safe_parse_float(tokens[2])
                            if idx >= 0
                                while length(para_qp_trans) <= idx
                                    push!(para_qp_trans, ComplexF64(0.0))
                                end
                                para_qp_trans[idx+1] = ComplexF64(value)
                            end
                            line_idx += 1
                            break
                        end
                        line_idx += 1
                    end
                end
                data.para_qp_trans = para_qp_trans

                # Build QPTrans mappings
                build_qp_trans_mappings!(data, file_path)
            end
        end
    elseif file_type == "OptTrans"
        parse_opttrans_def!(data, file_path)
    elseif file_type == "ChargeRBM_PhysLayer"
        result = parse_charge_rbm_phys_layer_def(file_path)
        if result.success
            data.charge_rbm_phys_layer_terms = result.data
        end
    elseif file_type == "SpinRBM_PhysLayer"
        result = parse_spin_rbm_phys_layer_def(file_path)
        if result.success
            data.spin_rbm_phys_layer_terms = result.data
        end
    elseif file_type == "GeneralRBM_PhysLayer"
        result = parse_general_rbm_phys_layer_def(file_path)
        if result.success
            data.general_rbm_phys_layer_terms = result.data
        end
    elseif file_type == "ChargeRBM_HiddenLayer"
        result = parse_charge_rbm_hidden_layer_def(file_path)
        if result.success
            data.charge_rbm_hidden_layer_terms = result.data
        end
    elseif file_type == "SpinRBM_HiddenLayer"
        result = parse_spin_rbm_hidden_layer_def(file_path)
        if result.success
            data.spin_rbm_hidden_layer_terms = result.data
        end
    elseif file_type == "GeneralRBM_HiddenLayer"
        result = parse_general_rbm_hidden_layer_def(file_path)
        if result.success
            data.general_rbm_hidden_layer_terms = result.data
        end
    elseif file_type == "ChargeRBM_PhysHidden"
        result = parse_charge_rbm_phys_hidden_def(file_path)
        if result.success
            data.charge_rbm_phys_hidden_terms = result.data
        end
    elseif file_type == "SpinRBM_PhysHidden"
        result = parse_spin_rbm_phys_hidden_def(file_path)
        if result.success
            data.spin_rbm_phys_hidden_terms = result.data
        end
    elseif file_type == "GeneralRBM_PhysHidden"
        result = parse_general_rbm_phys_hidden_def(file_path)
        if result.success
            data.general_rbm_phys_hidden_terms = result.data
        end
    elseif file_type == "DH2" || file_type == "DoublonHolon2Site"
        result = parse_doublon_holon_2site_def(file_path, data.modpara.nsite)
        if result.success
            data.doublon_holon_2site_indices = result.data.indices
            data.doublon_holon_2site_opt_flags = result.data.opt_flags
            data.doublon_holon_2site_complex = result.data.is_complex
            data.doublon_holon_2site_params = fill(0.0 + 0.0im, length(result.data.opt_flags))
        else
            error("Failed to parse DH2 file '$file_path': $(result.error_message)")
        end
    elseif file_type == "DH4" || file_type == "DoublonHolon4Site"
        result = parse_doublon_holon_4site_def(file_path, data.modpara.nsite)
        if result.success
            data.doublon_holon_4site_indices = result.data.indices
            data.doublon_holon_4site_opt_flags = result.data.opt_flags
            data.doublon_holon_4site_complex = result.data.is_complex
            data.doublon_holon_4site_params = fill(0.0 + 0.0im, length(result.data.opt_flags))
        else
            error("Failed to parse DH4 file '$file_path': $(result.error_message)")
        end
    end
end

end # module
