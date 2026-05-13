"""
Gutzwiller Parser

Parser for gutzwiller.def files containing Gutzwiller projection terms.
"""

"""
    parse_gutzwiller_def(filepath::String) -> ParseResult{Vector{GutzwillerTerm}}

Parse gutzwiller.def file from file path.
"""
function parse_gutzwiller_def(filepath::String)::ParseResult{Vector{GutzwillerTerm}}
    try
        content = read_def_file(filepath)
        return parse_gutzwiller_content(content)
    catch e
        return ParseResult{Vector{GutzwillerTerm}}(
            false,
            nothing,
            "Error reading file: $e",
            0,
        )
    end
end

"""
    GutzwillerParseResult

Extended result for Gutzwiller parsing that includes site->idx mapping.
"""
struct GutzwillerParseResult
    success::Bool
    terms::Union{Vector{GutzwillerTerm},Nothing}
    site_idx_map::Dict{Int,Int}  # site -> idx mapping
    n_gutzwiller_idx::Int
    opt_flags::Dict{Int, Int}  # idx -> opt_flag mapping
    is_complex_flag::Bool
    error_message::String
    line_number::Int
end

"""
    parse_gutzwiller_content(content::String) -> ParseResult{Vector{GutzwillerTerm}}

Parse gutzwiller.def content from string.
"""
function parse_gutzwiller_content(content::String)::ParseResult{Vector{GutzwillerTerm}}
    result = parse_gutzwiller_content_extended(content)
    return ParseResult{Vector{GutzwillerTerm}}(
        result.success,
        result.terms,
        result.error_message,
        result.line_number,
    )
end

"""
    parse_gutzwiller_content_extended(content::String) -> GutzwillerParseResult

Parse gutzwiller.def content from string, returning extended result with site->idx mapping.
"""
function parse_gutzwiller_content_extended(content::String)::GutzwillerParseResult
    context = ParsingContext("gutzwiller.def")
    terms = GutzwillerTerm[]
    site_idx_map = Dict{Int,Int}()  # site -> idx mapping
    opt_flags = Dict{Int,Int}()  # idx -> opt_flag mapping

    lines = split(content, '\n')

    # Skip header lines (first 5 lines) for gutzwilleridx.def format
    # Format: =============================================
    #         NGutzwillerIdx          N
    #         ComplexType          flag
    #         =============================================
    #         =============================================
    # Then data lines: site idx (Nsite lines), then idx opt_flag (NPara lines)
    IGNORE_LINES_IN_DEF = 5
    start_line = 1

    # Read NGutzwillerIdx from header (line 2: "NGutzwillerIdx N")
    n_gutzwiller_idx = 0
    has_header = false
    if length(lines) > 1
        header_line = clean_line(lines[2])
        tokens = split_def_line(header_line)
        if length(tokens) >= 2 && tokens[1] == "NGutzwillerIdx"
            n_gutzwiller_idx = safe_parse_int(tokens[2], 0)
            has_header = true
        end
    end

    # Read ComplexType from header (line 3: "ComplexType flag")
    # This corresponds to C implementation's iComplexFlgGutzwiller
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

    # Check if this looks like gutzwilleridx.def format (has header)
    if length(lines) > IGNORE_LINES_IN_DEF && has_header
        first_data_line = clean_line(lines[IGNORE_LINES_IN_DEF+1])
        if !isempty(first_data_line)
            tokens = split_def_line(first_data_line)
            # If first data line has exactly 2 tokens (site idx), skip header
            if length(tokens) == 2
                start_line = IGNORE_LINES_IN_DEF + 1
            end
        end
    end

    # If no header, process all lines normally (gutzwiller.def format)
    if !has_header
        # Process all lines as site value format
        for line_num = 1:length(lines)
            line = lines[line_num]
            context.line_number = line_num
            clean_line_str = clean_line(line)

            if isempty(clean_line_str)
                continue
            end

            tokens = split_def_line(clean_line_str)
            if length(tokens) < 2
                continue
            end

            try
                # Pass is_complex_flag to parse_gutzwiller_term
                term = parse_gutzwiller_term(tokens, context, is_complex_flag)
                if term !== nothing
                    push!(terms, term)
                end
            catch e
                push!(context.errors, "Line $line_num: Error parsing Gutzwiller term: $e")
            end
        end

        success = length(context.errors) == 0
        return GutzwillerParseResult(
            success,
            success ? terms : nothing,
            site_idx_map,
            n_gutzwiller_idx,
            opt_flags,
            is_complex_flag,
            join(context.errors, "; "),
            context.line_number,
        )
    end

    # Read site idx format lines (Nsite lines)
    # Then create NGutzwillerIdx terms based on unique idx values
    seen_idx_values = Set{Int}()
    seen_sites = Set{Int}()
    in_opt_section = false

    for line_num = start_line:length(lines)
        line = lines[line_num]
        context.line_number = line_num
        clean_line_str = clean_line(line)

        if isempty(clean_line_str)
            continue
        end

        tokens = split_def_line(clean_line_str)
        if length(tokens) < 2
            continue
        end

        site = safe_parse_int(tokens[1], -1)
        idx = safe_parse_int(tokens[2], -1)

        if site >= 0 && idx >= 0
            if !in_opt_section && (site in seen_sites)
                in_opt_section = true
            end

            if in_opt_section
                # idx opt_flag format
                opt_flags[site] = idx
            else
                # site idx format
                site_idx_map[site] = idx
                push!(seen_idx_values, idx)
                push!(seen_sites, site)
            end
        end
    end

    # Create NGutzwillerIdx terms based on unique idx values
    # If NGutzwillerIdx is specified, create that many terms
    # Otherwise, create one term per unique idx value
    if n_gutzwiller_idx > 0
        # Create NGutzwillerIdx terms
        for i = 1:n_gutzwiller_idx
            # Find a site with idx = i-1 (0-indexed)
            site = -1
            for (s, idx_val) in site_idx_map
                if idx_val == i - 1
                    site = s
                    break
                end
            end
            if site < 0
                # Use default site 0 if not found
                site = 0
            end
            # Use idx value as the term value
            value = ComplexF64(i - 1)
            term = GutzwillerTerm(site, value, is_complex_flag)
            push!(terms, term)
        end
    else
        # Create one term per unique idx value
        for idx_val in sort(collect(seen_idx_values))
            # Find a site with this idx value
            site = -1
            for (s, idx_val_map) in site_idx_map
                if idx_val_map == idx_val
                    site = s
                    break
                end
            end
            if site < 0
                site = 0
            end
            value = ComplexF64(idx_val)
            term = GutzwillerTerm(site, value, is_complex_flag)
            push!(terms, term)
        end
    end

    success = length(context.errors) == 0
    return GutzwillerParseResult(
        success,
        success ? terms : nothing,
        site_idx_map,
        n_gutzwiller_idx,
        opt_flags,
        is_complex_flag,
        join(context.errors, "; "),
        context.line_number,
    )
end

"""
    parse_gutzwiller_term(tokens::Vector{String}, context::ParsingContext, is_complex_flag::Bool) -> Union{GutzwillerTerm, Nothing}

Parse a single Gutzwiller term from tokens.
For gutzwilleridx.def format, tokens[2] is the idx, not a value.
The is_complex_flag comes from the ComplexType header in the file.
"""
function parse_gutzwiller_term(
    tokens::Vector{String},
    context::ParsingContext,
    is_complex_flag::Bool = false,
)::Union{GutzwillerTerm,Nothing}
    if length(tokens) < 2
        push!(context.warnings, "Insufficient tokens for Gutzwiller term")
        return nothing
    end

    # Parse site index
    site = safe_parse_int(tokens[1], -1)

    if site < 0
        push!(context.errors, "Invalid site index: $site")
        return nothing
    end

    # Parse second token
    # For gutzwilleridx.def: this is the idx, not a value
    # For gutzwiller.def (legacy format): this might be a value
    # Try to parse as integer first (gutzwilleridx.def format)
    idx = safe_parse_int(tokens[2], -1)

    if idx >= 0
        # This is gutzwilleridx.def format: site idx
        # Value is not specified here, initialize to 0.0
        # The actual value will be set from InGutzwiller.def or during initialization
        value = ComplexF64(0.0, 0.0)
    else
        # Try to parse as a value (legacy gutzwiller.def format)
        value = safe_parse_complex(tokens[2])
        # For legacy format, check if value has imaginary part
        # But prefer the is_complex_flag if provided
        if !is_complex_flag && imag(value) != 0.0
            # Fallback to checking value if flag not provided
            is_complex_flag = true
        end
    end

    # Use is_complex_flag from ComplexType header (C implementation's iComplexFlgGutzwiller)
    # This matches C implementation: ComplexType 0 = real, ComplexType != 0 = complex
    return GutzwillerTerm(site, value, is_complex_flag)
end
