"""
RBM Parser

Parser for RBM (Restricted Boltzmann Machine) related definition files.
RBM def files are index mapping tables + optimization flags (C-mVMC format).
"""

"""
    RBMParseResult

Extended parse result for RBM index files.
"""
struct RBMParseResult{T<:RBMTerm}
    success::Bool
    terms::Union{Vector{T}, Nothing}
    n_rbm_idx::Int
    opt_flags::Dict{Int, Int}  # idx -> opt_flag
    is_complex_flag::Bool
    error_message::String
    line_number::Int
end

function _to_simple_result(result::RBMParseResult{T})::ParseResult{Vector{T}} where {T<:RBMTerm}
    return ParseResult{Vector{T}}(
        result.success,
        result.terms,
        result.error_message,
        result.line_number,
    )
end

function _parse_rbm_header(lines::AbstractVector{<:AbstractString})::Tuple{Int, Bool, Int}
    # Default: no explicit header
    n_rbm_idx = 0
    is_complex_flag = false
    start_line = 1

    if length(lines) >= 3
        line2 = split_def_line(clean_line(lines[2]))
        line3 = split_def_line(clean_line(lines[3]))
        if length(line2) >= 2 && length(line3) >= 2 && line3[1] == "ComplexType"
            n_rbm_idx = safe_parse_int(line2[2], 0)
            is_complex_flag = safe_parse_int(line3[2], 0) != 0
            # mVMC idx.def format: first 5 lines are header/separators
            start_line = min(6, length(lines) + 1)
        end
    end

    return n_rbm_idx, is_complex_flag, start_line
end

function _collect_rbm_sections(
    lines::AbstractVector{<:AbstractString},
    start_line::Int,
    map_cols::Int,
)::Tuple{Vector{Tuple{Int, Vector{String}}}, Vector{Tuple{Int, Vector{String}}}}
    map_entries = Tuple{Int, Vector{String}}[]
    opt_entries = Tuple{Int, Vector{String}}[]

    in_opt_section = false
    seen_first_col = Set{Int}()

    for line_num in start_line:length(lines)
        tokens = split_def_line(clean_line(lines[line_num]))
        isempty(tokens) && continue

        if map_cols == 2
            # 2-column formats are ambiguous (map and opt both have 2 cols).
            # Split at the first repeated first-column index.
            if !in_opt_section && length(tokens) == 2
                first_col = safe_parse_int(tokens[1], typemin(Int))
                if first_col in seen_first_col
                    in_opt_section = true
                else
                    push!(seen_first_col, first_col)
                end
            end

            if in_opt_section
                if length(tokens) >= 2
                    push!(opt_entries, (line_num, tokens))
                end
            else
                if length(tokens) == 2
                    push!(map_entries, (line_num, tokens))
                end
            end
        else
            # For 3/4-column maps, opt section is 2-column.
            if in_opt_section
                if length(tokens) >= 2
                    push!(opt_entries, (line_num, tokens))
                end
            else
                if length(tokens) == map_cols
                    push!(map_entries, (line_num, tokens))
                elseif length(tokens) == 2
                    in_opt_section = true
                    push!(opt_entries, (line_num, tokens))
                end
            end
        end
    end

    return map_entries, opt_entries
end

function _parse_opt_flags(opt_entries::Vector{Tuple{Int, Vector{String}}})::Dict{Int, Int}
    opt_flags = Dict{Int, Int}()
    for (_line_num, tokens) in opt_entries
        if length(tokens) < 2
            continue
        end
        idx = safe_parse_int(tokens[1], -1)
        opt_flag = safe_parse_int(tokens[2], -1)
        if idx >= 0 && opt_flag >= 0
            opt_flags[idx] = opt_flag
        end
    end
    return opt_flags
end

function _parse_rbm_content(
    make_term::Function,
    content::String,
    context_name::String,
    map_cols::Int,
    term_type::Type{T},
)::RBMParseResult{T} where {T<:RBMTerm}
    context = ParsingContext(context_name)
    terms = T[]

    lines = split(content, '\n')
    n_rbm_idx, is_complex_flag, start_line = _parse_rbm_header(lines)
    map_entries, opt_entries = _collect_rbm_sections(lines, start_line, map_cols)

    for (line_num, tokens) in map_entries
        context.line_number = line_num
        try
            term = make_term(tokens, is_complex_flag)
            push!(terms, term)
        catch e
            push!(context.errors, "Line $line_num: Error parsing $context_name term: $e")
        end
    end

    if n_rbm_idx <= 0 && !isempty(terms)
        n_rbm_idx = maximum(t.idx for t in terms) + 1
    end

    opt_flags = _parse_opt_flags(opt_entries)

    success = isempty(context.errors)
    return RBMParseResult{T}(
        success,
        success ? terms : nothing,
        n_rbm_idx,
        opt_flags,
        is_complex_flag,
        join(context.errors, "; "),
        context.line_number,
    )
end

function _read_and_parse_rbm(filepath::String, parser::Function, term_type::Type{T})::ParseResult{Vector{T}} where {T<:RBMTerm}
    try
        content = read_def_file(filepath)
        return _to_simple_result(parser(content))
    catch e
        return ParseResult{Vector{T}}(false, nothing, "Error reading file: $e", 0)
    end
end

"""
    parse_charge_rbm_phys_layer_content_extended(content::String) -> RBMParseResult{ChargeRBMPhysLayerTerm}
"""
function parse_charge_rbm_phys_layer_content_extended(content::String)::RBMParseResult{ChargeRBMPhysLayerTerm}
    return _parse_rbm_content(
        content,
        "ChargeRBM_PhysLayer",
        2,
        ChargeRBMPhysLayerTerm,
    ) do tokens, is_complex_flag
        site = safe_parse_int(tokens[1], -1)
        idx = safe_parse_int(tokens[2], -1)
        if site < 0 || idx < 0
            error("invalid site/idx")
        end
        ChargeRBMPhysLayerTerm(site, ComplexF64(0.0), is_complex_flag, idx)
    end
end

"""
    parse_charge_rbm_phys_layer_content(content::String) -> ParseResult{Vector{ChargeRBMPhysLayerTerm}}
"""
function parse_charge_rbm_phys_layer_content(content::String)::ParseResult{Vector{ChargeRBMPhysLayerTerm}}
    return _to_simple_result(parse_charge_rbm_phys_layer_content_extended(content))
end

"""
    parse_charge_rbm_phys_layer_def(filepath::String) -> ParseResult{Vector{ChargeRBMPhysLayerTerm}}
"""
function parse_charge_rbm_phys_layer_def(filepath::String)::ParseResult{Vector{ChargeRBMPhysLayerTerm}}
    return _read_and_parse_rbm(filepath, parse_charge_rbm_phys_layer_content_extended, ChargeRBMPhysLayerTerm)
end

function parse_spin_rbm_phys_layer_content_extended(content::String)::RBMParseResult{SpinRBMPhysLayerTerm}
    return _parse_rbm_content(
        content,
        "SpinRBM_PhysLayer",
        2,
        SpinRBMPhysLayerTerm,
    ) do tokens, is_complex_flag
        site = safe_parse_int(tokens[1], -1)
        idx = safe_parse_int(tokens[2], -1)
        if site < 0 || idx < 0
            error("invalid site/idx")
        end
        SpinRBMPhysLayerTerm(site, ComplexF64(0.0), is_complex_flag, idx)
    end
end

function parse_spin_rbm_phys_layer_content(content::String)::ParseResult{Vector{SpinRBMPhysLayerTerm}}
    return _to_simple_result(parse_spin_rbm_phys_layer_content_extended(content))
end

function parse_spin_rbm_phys_layer_def(filepath::String)::ParseResult{Vector{SpinRBMPhysLayerTerm}}
    return _read_and_parse_rbm(filepath, parse_spin_rbm_phys_layer_content_extended, SpinRBMPhysLayerTerm)
end

function parse_general_rbm_phys_layer_content_extended(content::String)::RBMParseResult{GeneralRBMPhysLayerTerm}
    return _parse_rbm_content(
        content,
        "GeneralRBM_PhysLayer",
        3,
        GeneralRBMPhysLayerTerm,
    ) do tokens, is_complex_flag
        site = safe_parse_int(tokens[1], -1)
        spin = safe_parse_int(tokens[2], -1)
        idx = safe_parse_int(tokens[3], -1)
        if site < 0 || spin < 0 || idx < 0
            error("invalid site/spin/idx")
        end
        GeneralRBMPhysLayerTerm(site, spin, ComplexF64(0.0), is_complex_flag, idx)
    end
end

function parse_general_rbm_phys_layer_content(content::String)::ParseResult{Vector{GeneralRBMPhysLayerTerm}}
    return _to_simple_result(parse_general_rbm_phys_layer_content_extended(content))
end

function parse_general_rbm_phys_layer_def(filepath::String)::ParseResult{Vector{GeneralRBMPhysLayerTerm}}
    return _read_and_parse_rbm(filepath, parse_general_rbm_phys_layer_content_extended, GeneralRBMPhysLayerTerm)
end

function parse_charge_rbm_hidden_layer_content_extended(content::String)::RBMParseResult{ChargeRBMHiddenLayerTerm}
    return _parse_rbm_content(
        content,
        "ChargeRBM_HiddenLayer",
        2,
        ChargeRBMHiddenLayerTerm,
    ) do tokens, is_complex_flag
        hidden = safe_parse_int(tokens[1], -1)
        idx = safe_parse_int(tokens[2], -1)
        if hidden < 0 || idx < 0
            error("invalid hidden/idx")
        end
        ChargeRBMHiddenLayerTerm(hidden, ComplexF64(0.0), is_complex_flag, idx)
    end
end

function parse_charge_rbm_hidden_layer_content(content::String)::ParseResult{Vector{ChargeRBMHiddenLayerTerm}}
    return _to_simple_result(parse_charge_rbm_hidden_layer_content_extended(content))
end

function parse_charge_rbm_hidden_layer_def(filepath::String)::ParseResult{Vector{ChargeRBMHiddenLayerTerm}}
    return _read_and_parse_rbm(filepath, parse_charge_rbm_hidden_layer_content_extended, ChargeRBMHiddenLayerTerm)
end

function parse_spin_rbm_hidden_layer_content_extended(content::String)::RBMParseResult{SpinRBMHiddenLayerTerm}
    return _parse_rbm_content(
        content,
        "SpinRBM_HiddenLayer",
        2,
        SpinRBMHiddenLayerTerm,
    ) do tokens, is_complex_flag
        hidden = safe_parse_int(tokens[1], -1)
        idx = safe_parse_int(tokens[2], -1)
        if hidden < 0 || idx < 0
            error("invalid hidden/idx")
        end
        SpinRBMHiddenLayerTerm(hidden, ComplexF64(0.0), is_complex_flag, idx)
    end
end

function parse_spin_rbm_hidden_layer_content(content::String)::ParseResult{Vector{SpinRBMHiddenLayerTerm}}
    return _to_simple_result(parse_spin_rbm_hidden_layer_content_extended(content))
end

function parse_spin_rbm_hidden_layer_def(filepath::String)::ParseResult{Vector{SpinRBMHiddenLayerTerm}}
    return _read_and_parse_rbm(filepath, parse_spin_rbm_hidden_layer_content_extended, SpinRBMHiddenLayerTerm)
end

function parse_general_rbm_hidden_layer_content_extended(content::String)::RBMParseResult{GeneralRBMHiddenLayerTerm}
    return _parse_rbm_content(
        content,
        "GeneralRBM_HiddenLayer",
        2,
        GeneralRBMHiddenLayerTerm,
    ) do tokens, is_complex_flag
        hidden = safe_parse_int(tokens[1], -1)
        idx = safe_parse_int(tokens[2], -1)
        if hidden < 0 || idx < 0
            error("invalid hidden/idx")
        end
        GeneralRBMHiddenLayerTerm(hidden, ComplexF64(0.0), is_complex_flag, idx)
    end
end

function parse_general_rbm_hidden_layer_content(content::String)::ParseResult{Vector{GeneralRBMHiddenLayerTerm}}
    return _to_simple_result(parse_general_rbm_hidden_layer_content_extended(content))
end

function parse_general_rbm_hidden_layer_def(filepath::String)::ParseResult{Vector{GeneralRBMHiddenLayerTerm}}
    return _read_and_parse_rbm(filepath, parse_general_rbm_hidden_layer_content_extended, GeneralRBMHiddenLayerTerm)
end

function parse_charge_rbm_phys_hidden_content_extended(content::String)::RBMParseResult{ChargeRBMPhysHiddenTerm}
    return _parse_rbm_content(
        content,
        "ChargeRBM_PhysHidden",
        3,
        ChargeRBMPhysHiddenTerm,
    ) do tokens, is_complex_flag
        site = safe_parse_int(tokens[1], -1)
        hidden = safe_parse_int(tokens[2], -1)
        idx = safe_parse_int(tokens[3], -1)
        if site < 0 || hidden < 0 || idx < 0
            error("invalid site/hidden/idx")
        end
        ChargeRBMPhysHiddenTerm(site, hidden, ComplexF64(0.0), is_complex_flag, idx)
    end
end

function parse_charge_rbm_phys_hidden_content(content::String)::ParseResult{Vector{ChargeRBMPhysHiddenTerm}}
    return _to_simple_result(parse_charge_rbm_phys_hidden_content_extended(content))
end

function parse_charge_rbm_phys_hidden_def(filepath::String)::ParseResult{Vector{ChargeRBMPhysHiddenTerm}}
    return _read_and_parse_rbm(filepath, parse_charge_rbm_phys_hidden_content_extended, ChargeRBMPhysHiddenTerm)
end

function parse_spin_rbm_phys_hidden_content_extended(content::String)::RBMParseResult{SpinRBMPhysHiddenTerm}
    return _parse_rbm_content(
        content,
        "SpinRBM_PhysHidden",
        3,
        SpinRBMPhysHiddenTerm,
    ) do tokens, is_complex_flag
        site = safe_parse_int(tokens[1], -1)
        hidden = safe_parse_int(tokens[2], -1)
        idx = safe_parse_int(tokens[3], -1)
        if site < 0 || hidden < 0 || idx < 0
            error("invalid site/hidden/idx")
        end
        SpinRBMPhysHiddenTerm(site, hidden, ComplexF64(0.0), is_complex_flag, idx)
    end
end

function parse_spin_rbm_phys_hidden_content(content::String)::ParseResult{Vector{SpinRBMPhysHiddenTerm}}
    return _to_simple_result(parse_spin_rbm_phys_hidden_content_extended(content))
end

function parse_spin_rbm_phys_hidden_def(filepath::String)::ParseResult{Vector{SpinRBMPhysHiddenTerm}}
    return _read_and_parse_rbm(filepath, parse_spin_rbm_phys_hidden_content_extended, SpinRBMPhysHiddenTerm)
end

function parse_general_rbm_phys_hidden_content_extended(content::String)::RBMParseResult{GeneralRBMPhysHiddenTerm}
    return _parse_rbm_content(
        content,
        "GeneralRBM_PhysHidden",
        4,
        GeneralRBMPhysHiddenTerm,
    ) do tokens, is_complex_flag
        site = safe_parse_int(tokens[1], -1)
        spin = safe_parse_int(tokens[2], -1)
        hidden = safe_parse_int(tokens[3], -1)
        idx = safe_parse_int(tokens[4], -1)
        if site < 0 || spin < 0 || hidden < 0 || idx < 0
            error("invalid site/spin/hidden/idx")
        end
        GeneralRBMPhysHiddenTerm(site, spin, hidden, ComplexF64(0.0), is_complex_flag, idx)
    end
end

function parse_general_rbm_phys_hidden_content(content::String)::ParseResult{Vector{GeneralRBMPhysHiddenTerm}}
    return _to_simple_result(parse_general_rbm_phys_hidden_content_extended(content))
end

function parse_general_rbm_phys_hidden_def(filepath::String)::ParseResult{Vector{GeneralRBMPhysHiddenTerm}}
    return _read_and_parse_rbm(filepath, parse_general_rbm_phys_hidden_content_extended, GeneralRBMPhysHiddenTerm)
end
