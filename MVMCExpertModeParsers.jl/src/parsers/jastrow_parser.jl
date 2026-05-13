"""
Jastrow Parser

Parser for jastrow.def files containing Jastrow factor terms.
"""

"""
    parse_jastrow_def(filepath::String) -> ParseResult{Vector{JastrowTerm}}

Parse jastrow.def file from file path.
"""
function parse_jastrow_def(filepath::String)::ParseResult{Vector{JastrowTerm}}
    try
        content = read_def_file(filepath)
        return parse_jastrow_content(content)
    catch e
        return ParseResult{Vector{JastrowTerm}}(false, nothing, "Error reading file: $e", 0)
    end
end

"""
    parse_jastrow_content(content::String) -> ParseResult{Vector{JastrowTerm}}

Parse jastrow.def content from string.
"""
function parse_jastrow_content(content::String)::ParseResult{Vector{JastrowTerm}}
    result = parse_jastrow_content_extended(content)
    return ParseResult{Vector{JastrowTerm}}(
        result.success,
        result.terms,
        result.error_message,
        result.line_number
    )
end

"""
    JastrowParseResult

Extended result for Jastrow parsing that includes opt flags.
"""
struct JastrowParseResult
    success::Bool
    terms::Union{Vector{JastrowTerm}, Nothing}
    n_jastrow_idx::Int
    opt_flags::Dict{Int, Int}  # idx -> opt_flag mapping
    is_complex_flag::Bool
    error_message::String
    line_number::Int
end

"""
    parse_jastrow_content_extended(content::String) -> JastrowParseResult

Parse jastrow.def content from string, returning extended result with opt flags.
"""
function parse_jastrow_content_extended(content::String)::JastrowParseResult
    context = ParsingContext("jastrow.def")
    terms = JastrowTerm[]
    opt_flags = Dict{Int, Int}()  # idx -> opt_flag mapping

    lines = split(content, '\n')

    # Skip header lines (first 5 lines) for jastrowidx.def format
    # Format: =============================================
    #         NJastrowIdx          N
    #         ComplexType          flag
    #         =============================================
    #         =============================================
    # Then data lines: site1 site2 idx (Nsite*(Nsite-1) lines), then idx opt_flag (NPara lines)
    IGNORE_LINES_IN_DEF = 5
    start_line = 1

    # Read NJastrowIdx from header (line 2: "NJastrowIdx N")
    n_jastrow_idx = 0
    has_header = false
    if length(lines) > 1
        header_line = clean_line(lines[2])
        tokens = split_def_line(header_line)
        if length(tokens) >= 2 && tokens[1] == "NJastrowIdx"
            n_jastrow_idx = safe_parse_int(tokens[2], 0)
            has_header = true
        end
    end

    # Read ComplexType from header (line 3: "ComplexType flag")
    # This corresponds to C implementation's iComplexFlgJastrow
    complex_type = 0
    if has_header && length(lines) > 2
        complex_type_line = clean_line(lines[3])
        tokens = split_def_line(complex_type_line)
        if length(tokens) >= 2 && tokens[1] == "ComplexType"
            complex_type = safe_parse_int(tokens[2], 0)
        end
    end
    # Convert to boolean: 0 = false (real), != 0 = true (complex)
    is_complex_flag = (complex_type != 0)

    # Check if this looks like jastrowidx.def format (has header)
    if length(lines) > IGNORE_LINES_IN_DEF && has_header
        first_data_line = clean_line(lines[IGNORE_LINES_IN_DEF+1])
        if !isempty(first_data_line)
            tokens = split_def_line(first_data_line)
            # If first data line has at least 3 tokens (site1 site2 idx), skip header
            if length(tokens) >= 3
                start_line = IGNORE_LINES_IN_DEF + 1
            end
        end
    end

    # If no header, process all lines normally (jastrow.def format)
    if !has_header
        start_line = 1
    end

    in_opt_section = false
    for line_num = start_line:length(lines)
        line = lines[line_num]
        context.line_number = line_num
        clean_line_str = clean_line(line)

        if isempty(clean_line_str)
            continue
        end

        tokens = split_def_line(clean_line_str)

        if length(tokens) == 2
            if !isempty(terms)
                in_opt_section = true
            else
                continue
            end
        end

        if in_opt_section
            if length(tokens) >= 2
                idx = safe_parse_int(tokens[1], -1)
                opt_flag = safe_parse_int(tokens[2], -1)
                if idx >= 0 && opt_flag >= 0
                    opt_flags[idx] = opt_flag
                end
            end
            continue
        end

        if length(tokens) < 3
            continue
        end

        try
            # Pass is_complex_flag to parse_jastrow_term
            term = parse_jastrow_term(tokens, context, is_complex_flag)
            if term !== nothing
                push!(terms, term)
            end
        catch e
            push!(context.errors, "Line $line_num: Error parsing Jastrow term: $e")
        end
    end

    # If NJastrowIdx is specified and we have more terms, limit to unique idx values
    if n_jastrow_idx > 0 && length(terms) > n_jastrow_idx
        # For now, just take first n_jastrow_idx terms
        terms = terms[1:min(n_jastrow_idx, length(terms))]
    end

    success = length(context.errors) == 0
    return JastrowParseResult(
        success,
        success ? terms : nothing,
        n_jastrow_idx,
        opt_flags,
        is_complex_flag,
        join(context.errors, "; "),
        context.line_number,
    )
end

"""
    parse_jastrow_term(tokens::Vector{String}, context::ParsingContext, is_complex_flag::Bool) -> Union{JastrowTerm, Nothing}

Parse a single Jastrow term from tokens.
For jastrowidx.def format, tokens[3] is the idx, not a value.
The is_complex_flag comes from the ComplexType header in the file.
"""
function parse_jastrow_term(
    tokens::Vector{String},
    context::ParsingContext,
    is_complex_flag::Bool = false,
)::Union{JastrowTerm,Nothing}
    if length(tokens) < 3
        push!(context.warnings, "Insufficient tokens for Jastrow term")
        return nothing
    end

    # Parse site indices
    site1 = safe_parse_int(tokens[1], -1)
    site2 = safe_parse_int(tokens[2], -1)

    if site1 < 0 || site2 < 0
        push!(context.errors, "Invalid site indices: $site1, $site2")
        return nothing
    end

    # Parse third token
    # For jastrowidx.def: this is the idx, not a value
    # For jastrow.def (legacy format): this might be a value
    # Try to parse as integer first (jastrowidx.def format)
    idx = safe_parse_int(tokens[3], -1)

    if idx >= 0
        # This is jastrowidx.def format: site1 site2 idx
        # Value is not specified here, initialize to 0.0
        # The actual value will be set from InJastrow.def or during initialization
        value = ComplexF64(0.0, 0.0)
    else
        # Try to parse as a value (legacy jastrow.def format)
        value = safe_parse_complex(tokens[3])
        # For legacy format, check if value has imaginary part
        # But prefer the is_complex_flag if provided
        if !is_complex_flag && imag(value) != 0.0
            # Fallback to checking value if flag not provided
            is_complex_flag = true
        end
    end

    # Use is_complex_flag from ComplexType header (C implementation's iComplexFlgJastrow)
    # This matches C implementation: ComplexType 0 = real, ComplexType != 0 = complex
    return JastrowTerm(site1, site2, value, is_complex_flag)
end

"""
    build_jastrow_idx_matrix(filepath::String, nsite::Int) -> Tuple{Matrix{Int}, Int}

Build the JastrowIdx matrix from jastrowidx.def file.
Returns (jastrow_idx_matrix, n_jastrow_idx).
The matrix has size (nsite, nsite) where jastrow_idx[site1+1, site2+1] = idx.
"""
function build_jastrow_idx_matrix(filepath::String, nsite::Int)::Tuple{Matrix{Int},Int}
    if nsite <= 0
        return (Matrix{Int}(undef, 0, 0), 0)
    end

    # Initialize with -1 (invalid index)
    jastrow_idx = fill(-1, nsite, nsite)
    n_jastrow_idx = 0

    try
        content = read_def_file(filepath)
        lines = split(content, '\n')

        # Read NJastrowIdx from header (line 2: "NJastrowIdx N")
        if length(lines) > 1
            header_line = clean_line(lines[2])
            tokens = split_def_line(header_line)
            if length(tokens) >= 2 && tokens[1] == "NJastrowIdx"
                n_jastrow_idx = safe_parse_int(tokens[2], 0)
            end
        end

        # Skip header lines (first 5 lines)
        IGNORE_LINES_IN_DEF = 5

        for line_num = (IGNORE_LINES_IN_DEF+1):length(lines)
            line = lines[line_num]
            clean_line_str = clean_line(line)

            if isempty(clean_line_str)
                continue
            end

            tokens = split_def_line(clean_line_str)

            # Skip 2-column lines (idx opt_flag format)
            if length(tokens) == 2
                break
            end

            if length(tokens) < 3
                continue
            end

            # Parse site1 site2 idx
            site1 = safe_parse_int(tokens[1], -1)
            site2 = safe_parse_int(tokens[2], -1)
            idx = safe_parse_int(tokens[3], -1)

            if site1 >= 0 && site1 < nsite && site2 >= 0 && site2 < nsite && idx >= 0
                # Store both directions for symmetric access
                jastrow_idx[site1+1, site2+1] = idx
                jastrow_idx[site2+1, site1+1] = idx  # Symmetric
            end
        end

    catch e
        @warn "Error building JastrowIdx matrix: $e"
    end

    return (jastrow_idx, n_jastrow_idx)
end
