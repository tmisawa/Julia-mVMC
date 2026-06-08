"""
Orbital Parser

Parser for orbital.def files containing orbital parameter terms.
"""

"""
    parse_orbital_def(filepath::String) -> Tuple{ParseResult{Vector{OrbitalTerm}}, Dict{Int, Int}, Int}

Parse orbital.def file from file path.
Returns orbital terms, OptFlag dictionary (idx -> opt_flag), and the header
declared parameter count (`NOrbitalIdx` / `NOrbitalParallel`). The declared count
is what C reads via `ReadBuffIntCmpFlg` (`iNOrbitalAntiParallel`/`iNOrbitalParallel`),
and is the authoritative parameter count even when some indices are unreferenced.
"""
function parse_orbital_def(
    filepath::String,
)::Tuple{ParseResult{Vector{OrbitalTerm}},Dict{Int,Int},Int}
    try
        content = read_def_file(filepath)
        return parse_orbital_content(content)
    catch e
        return (
            ParseResult{Vector{OrbitalTerm}}(false, nothing, "Error reading file: $e", 0),
            Dict{Int,Int}(),
            0,
        )
    end
end

"""
    parse_orbital_content(content::String) -> Tuple{ParseResult{Vector{OrbitalTerm}}, Dict{Int, Int}, Int}

Parse orbital.def content from string.
Returns orbital terms, OptFlag dictionary (idx -> opt_flag), and the header
declared parameter count (0 for the headerless legacy format, where it falls
back to `max(idx)+1`).
"""
function parse_orbital_content(
    content::String,
)::Tuple{ParseResult{Vector{OrbitalTerm}},Dict{Int,Int},Int}
    context = ParsingContext("orbital.def")
    terms = OrbitalTerm[]

    lines = split(content, '\n')

    # Skip header lines (first 5 lines) for orbitalidx.def format
    # Format: =============================================
    #         NOrbitalIdx          N
    #         ComplexType          flag
    #         =============================================
    #         =============================================
    # Then data lines: site1 site2 idx (NOrbitalIdx lines), then idx opt_flag (NPara lines)
    IGNORE_LINES_IN_DEF = 5
    start_line = 1

    # Read NOrbitalIdx from header (line 2: "NOrbitalIdx N")
    n_orbital_idx = 0
    has_header = false
    if length(lines) > 1
        header_line = clean_line(lines[2])
        tokens = split_def_line(header_line)
        # Accept any NOrbital* header keyword (NOrbitalIdx, NOrbitalParallel,
        # NOrbitalAntiParallel, NOrbitalGeneral). C reads the count positionally
        # via ReadBuffIntCmpFlg regardless of the keyword spelling.
        if length(tokens) >= 2 && startswith(tokens[1], "NOrbital")
            n_orbital_idx = safe_parse_int(tokens[2], 0)
            has_header = true
        end
    end

    # Read ComplexType from header (line 3: "ComplexType flag")
    # This corresponds to C implementation's iComplexFlgOrbital
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

    # Check if this looks like orbitalidx.def format (has header)
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

    # If no header, process all lines normally (orbital.def format)
    if !has_header
        start_line = 1
    end

    # Parse all OrbitalIdx entries (3-column lines: site1 site2 idx)
    # Note: orbitalidx.def has N_site * N_site OrbitalIdx entries (256 for 16 sites),
    # followed by NOrbitalIdx OptFlag entries (idx opt_flag format, 2 columns).
    # The NOrbitalIdx value is the number of unique orbital parameters, NOT the number of OrbitalIdx entries.
    processing_orbital_idx = true  # Flag to track if we're still in OrbitalIdx section
    opt_flags = Dict{Int,Int}()  # Store OptFlag: idx -> opt_flag

    for line_num = start_line:length(lines)
        line = lines[line_num]
        context.line_number = line_num
        clean_line_str = clean_line(line)

        if isempty(clean_line_str)
            continue
        end

        tokens = split_def_line(clean_line_str)

        # 2-column lines are OptFlag (idx opt_flag format)
        # Once we see 2-column lines after 3-column lines, we're done with OrbitalIdx
        if length(tokens) == 2
            if processing_orbital_idx && !isempty(terms)
                # We've finished OrbitalIdx section, now in OptFlag section
                processing_orbital_idx = false
            end

            # Parse OptFlag: idx opt_flag
            if !processing_orbital_idx
                try
                    idx = safe_parse_int(tokens[1], -1)
                    opt_flag = safe_parse_int(tokens[2], 0)
                    if idx >= 0
                        opt_flags[idx] = opt_flag
                    end
                catch e
                    push!(context.errors, "Line $line_num: Error parsing OptFlag: $e")
                end
            end
            continue
        end

        if length(tokens) < 3
            continue
        end

        # Only process 3-column lines in OrbitalIdx section
        if !processing_orbital_idx
            continue
        end

        try
            # Pass is_complex_flag to parse_orbital_term
            term = parse_orbital_term(tokens, context, is_complex_flag)
            if term !== nothing
                push!(terms, term)
            end
        catch e
            push!(context.errors, "Line $line_num: Error parsing orbital term: $e")
        end
    end

    # Note: NOrbitalIdx is the number of unique orbital parameters, not the number of OrbitalIdx entries.
    # All OrbitalIdx entries (site pairs) should be parsed, even if there are more than NOrbitalIdx.

    success = length(context.errors) == 0
    result = ParseResult{Vector{OrbitalTerm}}(
        success,
        success ? terms : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
    # Declared parameter count from the header. C reads this directly as
    # iNOrbitalAntiParallel / iNOrbitalParallel (ReadBuffIntCmpFlg), so it is the
    # authoritative count even when some indices are unreferenced by site pairs.
    # Fall back to max(idx)+1 only for the headerless legacy orbital.def format.
    declared_count = if has_header
        n_orbital_idx
    elseif !isempty(terms)
        maximum(t.idx for t in terms) + 1
    else
        0
    end
    return (result, opt_flags, declared_count)
end

"""
    parse_orbital_term(tokens::Vector{String}, context::ParsingContext, is_complex_flag::Bool) -> Union{OrbitalTerm, Nothing}

Parse a single orbital term from tokens.
For orbitalidx.def format, tokens[3] is the orbital index (idx), not a value.
The is_complex_flag comes from the ComplexType header in the file.
"""
function parse_orbital_term(
    tokens::Vector{String},
    context::ParsingContext,
    is_complex_flag::Bool = false,
)::Union{OrbitalTerm,Nothing}
    if length(tokens) < 3
        push!(context.warnings, "Insufficient tokens for orbital term")
        return nothing
    end

    # Parse site indices
    site1 = safe_parse_int(tokens[1], -1)
    site2 = safe_parse_int(tokens[2], -1)

    if site1 < 0 || site2 < 0
        push!(context.errors, "Invalid site indices: $site1, $site2")
        return nothing
    end

    # Parse third token and optional fourth token (sign)
    # For orbitalidx.def: this is the orbital index (idx), not a value
    # For orbital.def (legacy format): this might be a value
    # Try to parse as integer first (orbitalidx.def format)
    idx = safe_parse_int(tokens[3], -1)

    # Parse sign (4th token, optional, default: 1)
    sign = 1
    if length(tokens) >= 4
        sign_parsed = safe_parse_int(tokens[4], 1)
        if sign_parsed == 1 || sign_parsed == -1
            sign = sign_parsed
        end
    end

    if idx >= 0
        # This is orbitalidx.def format: site1 site2 idx [sgn]
        # Value is not specified here, initialize to 0.0
        # The actual value will be set from InOrbital.def or during initialization
        value = ComplexF64(0.0, 0.0)
        # Store idx in OrbitalTerm
        return OrbitalTerm(site1, site2, idx, value, is_complex_flag, sign)
    else
        # Try to parse as a value (legacy orbital.def format)
        value = safe_parse_complex(tokens[3])
        # For legacy format, check if value has imaginary part
        # But prefer the is_complex_flag if provided
        if is_complex_flag
            # Use flag from header
        elseif imag(value) != 0.0
            # Fallback to checking value if flag not provided
            is_complex_flag = true
        end
        # Use is_complex_flag from ComplexType header (C implementation's iComplexFlgOrbital)
        # This matches C implementation: ComplexType 0 = real, ComplexType != 0 = complex
        return OrbitalTerm(site1, site2, 0, value, is_complex_flag, sign)
    end
end
