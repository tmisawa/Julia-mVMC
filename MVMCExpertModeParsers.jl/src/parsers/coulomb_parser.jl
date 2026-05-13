"""
Coulomb Parser

Parser for coulombintra.def and coulombinter.def files.
"""

"""
    parse_coulomb_intra_def(filepath::String) -> ParseResult{Vector{CoulombIntraTerm}}

Parse coulombintra.def file from file path.
"""
function parse_coulomb_intra_def(filepath::String)::ParseResult{Vector{CoulombIntraTerm}}
    try
        content = read_def_file(filepath)
        return parse_coulomb_intra_content(content)
    catch e
        return ParseResult{Vector{CoulombIntraTerm}}(
            false,
            nothing,
            "Error reading file: $e",
            0,
        )
    end
end

"""
    parse_coulomb_intra_content(content::String) -> ParseResult{Vector{CoulombIntraTerm}}

Parse coulombintra.def content from string.
"""
function parse_coulomb_intra_content(content::String)::ParseResult{Vector{CoulombIntraTerm}}
    context = ParsingContext("coulombintra.def")
    terms = CoulombIntraTerm[]

    lines = split(content, '\n')

    for (line_num, line) in enumerate(lines)
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
            term = parse_coulomb_intra_term(tokens, context)
            if term !== nothing
                push!(terms, term)
            end
        catch e
            push!(context.errors, "Line $line_num: Error parsing Coulomb intra term: $e")
        end
    end

    success = length(context.errors) == 0
    return ParseResult{Vector{CoulombIntraTerm}}(
        success,
        success ? terms : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
end

"""
    parse_coulomb_intra_term(tokens::Vector{String}, context::ParsingContext) -> Union{CoulombIntraTerm, Nothing}

Parse a single Coulomb intra term from tokens.
Format: site value
"""
function parse_coulomb_intra_term(
    tokens::Vector{String},
    context::ParsingContext,
)::Union{CoulombIntraTerm,Nothing}
    if length(tokens) < 2
        # Skip lines with insufficient tokens (headers, etc.)
        return nothing
    end

    # Skip lines that don't start with a number (header lines)
    # This handles lines like "NCoulombIntra 8" or "=================="
    site = tryparse(Int, tokens[1])
    if site === nothing
        # Not a valid integer - likely a header line, skip silently
        return nothing
    end

    if site < 0
        push!(context.warnings, "Negative site index: $site")
        return nothing
    end

    # Parse value
    value = safe_parse_float(tokens[2])

    return CoulombIntraTerm(site, value)
end

"""
    parse_coulomb_inter_def(filepath::String) -> ParseResult{Vector{CoulombInterTerm}}

Parse coulombinter.def file from file path.
"""
function parse_coulomb_inter_def(filepath::String)::ParseResult{Vector{CoulombInterTerm}}
    try
        content = read_def_file(filepath)
        return parse_coulomb_inter_content(content)
    catch e
        return ParseResult{Vector{CoulombInterTerm}}(
            false,
            nothing,
            "Error reading file: $e",
            0,
        )
    end
end

"""
    parse_coulomb_inter_content(content::String) -> ParseResult{Vector{CoulombInterTerm}}

Parse coulombinter.def content from string.
"""
function parse_coulomb_inter_content(content::String)::ParseResult{Vector{CoulombInterTerm}}
    context = ParsingContext("coulombinter.def")
    terms = CoulombInterTerm[]

    lines = split(content, '\n')

    # Skip header lines (first 5 lines) for coulombinter.def format
    IGNORE_LINES_IN_DEF = 5
    start_line = 1

    # Check if this looks like coulombinter.def format (has header)
    if length(lines) > IGNORE_LINES_IN_DEF
        first_data_line = clean_line(lines[IGNORE_LINES_IN_DEF+1])
        if !isempty(first_data_line)
            tokens = split_def_line(first_data_line)
            # If first data line has at least 3 tokens (site1 site2 value), skip header
            if length(tokens) >= 3
                start_line = IGNORE_LINES_IN_DEF + 1
            end
        end
    end

    for line_num = start_line:length(lines)
        line = lines[line_num]
        context.line_number = line_num
        clean_line_str = clean_line(line)

        if isempty(clean_line_str)
            continue
        end

        tokens = split_def_line(clean_line_str)
        if length(tokens) < 3
            continue
        end

        try
            term = parse_coulomb_inter_term(tokens, context)
            if term !== nothing
                push!(terms, term)
            end
        catch e
            push!(context.errors, "Line $line_num: Error parsing Coulomb inter term: $e")
        end
    end

    success = length(context.errors) == 0
    return ParseResult{Vector{CoulombInterTerm}}(
        success,
        success ? terms : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
end

"""
    parse_coulomb_inter_term(tokens::Vector{String}, context::ParsingContext) -> Union{CoulombInterTerm, Nothing}

Parse a single Coulomb inter term from tokens.
"""
function parse_coulomb_inter_term(
    tokens::Vector{String},
    context::ParsingContext,
)::Union{CoulombInterTerm,Nothing}
    if length(tokens) < 3
        push!(context.warnings, "Insufficient tokens for Coulomb inter term")
        return nothing
    end

    # Parse site indices
    site1 = safe_parse_int(tokens[1], -1)
    site2 = safe_parse_int(tokens[2], -1)

    if site1 < 0 || site2 < 0
        push!(context.errors, "Invalid site indices: $site1, $site2")
        return nothing
    end

    # Parse value
    value = safe_parse_float(tokens[3])

    return CoulombInterTerm(site1, site2, value)
end
