"""
Green's Function Parser

Parser for greenone.def and greentwo.def files containing Green's function measurement terms.
"""

"""
    parse_green_one_def(filepath::String) -> ParseResult{Vector{GreenOneTerm}}

Parse greenone.def file from file path.
"""
function parse_green_one_def(filepath::String)::ParseResult{Vector{GreenOneTerm}}
    try
        content = read_def_file(filepath)
        return parse_green_one_content(content)
    catch e
        return ParseResult{Vector{GreenOneTerm}}(
            false,
            nothing,
            "Error reading file: $e",
            0,
        )
    end
end

"""
    parse_green_one_content(content::String) -> ParseResult{Vector{GreenOneTerm}}

Parse greenone.def content from string.
"""
function parse_green_one_content(content::String)::ParseResult{Vector{GreenOneTerm}}
    context = ParsingContext("greenone.def")
    terms = GreenOneTerm[]

    lines = split(content, '\n')

    # Skip header lines (first 5 lines) for greenone.def format
    IGNORE_LINES_IN_DEF = 5
    start_line = 1

    # Check if this looks like greenone.def format (has header)
    if length(lines) > IGNORE_LINES_IN_DEF
        first_data_line = clean_line(lines[IGNORE_LINES_IN_DEF+1])
        if !isempty(first_data_line)
            tokens = split_def_line(first_data_line)
            # If first data line has at least 2 tokens (site1 site2), skip header
            if length(tokens) >= 2
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
        if length(tokens) < 2
            continue
        end

        try
            term = parse_green_one_term(tokens, context)
            if term !== nothing
                push!(terms, term)
            end
        catch e
            push!(context.errors, "Line $line_num: Error parsing Green's function term: $e")
        end
    end

    success = length(context.errors) == 0
    return ParseResult{Vector{GreenOneTerm}}(
        success,
        success ? terms : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
end

"""
    parse_green_one_term(tokens::Vector{String}, context::ParsingContext) -> Union{GreenOneTerm, Nothing}

Parse a single Green's function term from tokens.

C format (cisajs.def / greenone.def): ri si rj sj
- Column 1: ri (site i)
- Column 2: si (spin i, 0=up, 1=down)
- Column 3: rj (site j)
- Column 4: sj (spin j, 0=up, 1=down)
"""
function parse_green_one_term(
    tokens::Vector{String},
    context::ParsingContext,
)::Union{GreenOneTerm,Nothing}
    if length(tokens) < 4
        push!(
            context.warnings,
            "Insufficient tokens for Green's function term (need ri si rj sj)",
        )
        return nothing
    end

    # Parse: ri si rj sj (C format)
    site1 = safe_parse_int(tokens[1], -1)      # ri
    spin1_int = safe_parse_int(tokens[2], -1)  # si (0=up, 1=down)
    site2 = safe_parse_int(tokens[3], -1)      # rj
    spin2_int = safe_parse_int(tokens[4], -1)  # sj (0=up, 1=down)

    if site1 < 0 || site2 < 0
        push!(context.errors, "Invalid site indices: $site1, $site2")
        return nothing
    end

    if spin1_int < 0 || spin1_int > 1 || spin2_int < 0 || spin2_int > 1
        push!(
            context.errors,
            "Invalid spin indices (must be 0 or 1): $spin1_int, $spin2_int",
        )
        return nothing
    end

    # Convert spin integer to symbol (0=up, 1=down)
    spin1 = spin1_int == 0 ? :up : :down
    spin2 = spin2_int == 0 ? :up : :down

    return GreenOneTerm(site1, site2, spin1, spin2)
end

"""
    parse_green_two_def(filepath::String) -> ParseResult{Vector{GreenTwoTerm}}

Parse greentwo.def file from file path.
"""
function parse_green_two_def(filepath::String)::ParseResult{Vector{GreenTwoTerm}}
    try
        content = read_def_file(filepath)
        return parse_green_two_content(content)
    catch e
        return ParseResult{Vector{GreenTwoTerm}}(
            false,
            nothing,
            "Error reading file: $e",
            0,
        )
    end
end

"""
    parse_green_two_content(content::String) -> ParseResult{Vector{GreenTwoTerm}}

Parse greentwo.def content from string.
"""
function parse_green_two_content(content::String)::ParseResult{Vector{GreenTwoTerm}}
    context = ParsingContext("greentwo.def")
    terms = GreenTwoTerm[]

    lines = split(content, '\n')

    # Skip header lines (first 5 lines) for greentwo.def format
    IGNORE_LINES_IN_DEF = 5
    start_line = 1

    # Check if this looks like greentwo.def format (has header)
    if length(lines) > IGNORE_LINES_IN_DEF
        first_data_line = clean_line(lines[IGNORE_LINES_IN_DEF+1])
        if !isempty(first_data_line)
            tokens = split_def_line(first_data_line)
            # If first data line has at least 4 tokens (site1 site2 site3 site4), skip header
            if length(tokens) >= 4
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
        if length(tokens) < 4
            continue
        end

        try
            term = parse_green_two_term(tokens, context)
            if term !== nothing
                push!(terms, term)
            end
        catch e
            push!(
                context.errors,
                "Line $line_num: Error parsing two-body Green's function term: $e",
            )
        end
    end

    success = length(context.errors) == 0
    return ParseResult{Vector{GreenTwoTerm}}(
        success,
        success ? terms : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
end

"""
    parse_green_two_term(tokens::Vector{String}, context::ParsingContext) -> Union{GreenTwoTerm, Nothing}

Parse a single two-body Green's function term from tokens.

C format (greentwo.def / cisajscktalt.def): ri si rj sj rk sk rl sl
- Columns 1-2: ri si (site i, spin i)
- Columns 3-4: rj sj (site j, spin j)
- Columns 5-6: rk sk (site k, spin k)
- Columns 7-8: rl sl (site l, spin l)

This represents <c†_{ri,si} c_{rj,sj} c†_{rk,sk} c_{rl,sl}>
"""
function parse_green_two_term(
    tokens::Vector{String},
    context::ParsingContext,
)::Union{GreenTwoTerm,Nothing}
    if length(tokens) < 8
        push!(
            context.warnings,
            "Insufficient tokens for two-body Green's function term (need ri si rj sj rk sk rl sl)",
        )
        return nothing
    end

    # Parse: ri si rj sj rk sk rl sl (C format)
    site1 = safe_parse_int(tokens[1], -1)      # ri
    spin1_int = safe_parse_int(tokens[2], -1)  # si
    site2 = safe_parse_int(tokens[3], -1)      # rj
    spin2_int = safe_parse_int(tokens[4], -1)  # sj
    site3 = safe_parse_int(tokens[5], -1)      # rk
    spin3_int = safe_parse_int(tokens[6], -1)  # sk
    site4 = safe_parse_int(tokens[7], -1)      # rl
    spin4_int = safe_parse_int(tokens[8], -1)  # sl

    if site1 < 0 || site2 < 0 || site3 < 0 || site4 < 0
        push!(context.errors, "Invalid site indices: $site1, $site2, $site3, $site4")
        return nothing
    end

    # Validate spin indices (must be 0 or 1)
    for (i, spin_int) in enumerate([spin1_int, spin2_int, spin3_int, spin4_int])
        if spin_int < 0 || spin_int > 1
            push!(context.errors, "Invalid spin$i index (must be 0 or 1): $spin_int")
            return nothing
        end
    end

    # Convert spin integer to symbol (0=up, 1=down)
    spin1 = spin1_int == 0 ? :up : :down
    spin2 = spin2_int == 0 ? :up : :down
    spin3 = spin3_int == 0 ? :up : :down
    spin4 = spin4_int == 0 ? :up : :down

    return GreenTwoTerm(site1, site2, site3, site4, spin1, spin2, spin3, spin4)
end
