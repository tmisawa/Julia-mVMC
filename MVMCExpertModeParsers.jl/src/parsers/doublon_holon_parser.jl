"""
Doublon-Holon Parser

Parser for Doublon-Holon related definition files.
Handles 2-site and 4-site Doublon-Holon terms.
"""

"""
    parse_doublon_holon_2site_def(filepath::String) -> ParseResult{Vector{DoublonHolon2SiteTerm}}

Parse Doublon-Holon 2-site definition file.
"""
function parse_doublon_holon_2site_def(
    filepath::String,
)::ParseResult{Vector{DoublonHolon2SiteTerm}}
    try
        content = read_def_file(filepath)
        return parse_doublon_holon_2site_content(content)
    catch e
        return ParseResult{Vector{DoublonHolon2SiteTerm}}(
            false,
            nothing,
            "Error reading file: $e",
            0,
        )
    end
end

"""
    parse_doublon_holon_2site_content(content::String) -> ParseResult{Vector{DoublonHolon2SiteTerm}}

Parse Doublon-Holon 2-site content from string.
"""
function parse_doublon_holon_2site_content(
    content::String,
)::ParseResult{Vector{DoublonHolon2SiteTerm}}
    context = ParsingContext("DoublonHolon2Site")
    terms = DoublonHolon2SiteTerm[]

    lines = split(content, '\n')

    for (line_num, line) in enumerate(lines)
        context.line_number = line_num
        clean_line_str = clean_line(line)

        if isempty(clean_line_str) || startswith(clean_line_str, "#")
            continue
        end

        tokens = split_def_line(clean_line_str)
        if length(tokens) < 4
            continue
        end

        try
            site1 = safe_parse_int(tokens[1], 0)
            site2 = safe_parse_int(tokens[2], 0)
            value = safe_parse_complex(tokens[3], 0.0 + 0.0im)
            is_complex = length(tokens) >= 5 ? safe_parse_int(tokens[5], 0) != 0 : false

            push!(terms, DoublonHolon2SiteTerm(site1, site2, value, is_complex))
        catch e
            push!(
                context.errors,
                "Line $line_num: Error parsing DoublonHolon2Site term: $e",
            )
        end
    end

    success = length(context.errors) == 0
    return ParseResult{Vector{DoublonHolon2SiteTerm}}(
        success,
        success ? terms : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
end

"""
    parse_doublon_holon_4site_def(filepath::String) -> ParseResult{Vector{DoublonHolon4SiteTerm}}

Parse Doublon-Holon 4-site definition file.
"""
function parse_doublon_holon_4site_def(
    filepath::String,
)::ParseResult{Vector{DoublonHolon4SiteTerm}}
    try
        content = read_def_file(filepath)
        return parse_doublon_holon_4site_content(content)
    catch e
        return ParseResult{Vector{DoublonHolon4SiteTerm}}(
            false,
            nothing,
            "Error reading file: $e",
            0,
        )
    end
end

"""
    parse_doublon_holon_4site_content(content::String) -> ParseResult{Vector{DoublonHolon4SiteTerm}}

Parse Doublon-Holon 4-site content from string.
"""
function parse_doublon_holon_4site_content(
    content::String,
)::ParseResult{Vector{DoublonHolon4SiteTerm}}
    context = ParsingContext("DoublonHolon4Site")
    terms = DoublonHolon4SiteTerm[]

    lines = split(content, '\n')

    for (line_num, line) in enumerate(lines)
        context.line_number = line_num
        clean_line_str = clean_line(line)

        if isempty(clean_line_str) || startswith(clean_line_str, "#")
            continue
        end

        tokens = split_def_line(clean_line_str)
        if length(tokens) < 5
            continue
        end

        try
            site1 = safe_parse_int(tokens[1], 0)
            site2 = safe_parse_int(tokens[2], 0)
            site3 = safe_parse_int(tokens[3], 0)
            site4 = safe_parse_int(tokens[4], 0)
            value = safe_parse_complex(tokens[5], 0.0 + 0.0im)
            is_complex = length(tokens) >= 6 ? safe_parse_int(tokens[6], 0) != 0 : false

            push!(
                terms,
                DoublonHolon4SiteTerm(site1, site2, site3, site4, value, is_complex),
            )
        catch e
            push!(
                context.errors,
                "Line $line_num: Error parsing DoublonHolon4Site term: $e",
            )
        end
    end

    success = length(context.errors) == 0
    return ParseResult{Vector{DoublonHolon4SiteTerm}}(
        success,
        success ? terms : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
end
