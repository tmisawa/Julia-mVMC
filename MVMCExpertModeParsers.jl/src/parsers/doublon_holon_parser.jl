"""
Doublon-Holon Parser

C-compatible parsers for `DH2` / `DH4` (`doublonholon2siteidx.def` /
`doublonholon4siteidx.def`). These files define neighbor index tables; the
projection parameters live in the C `Proj` slice and are stored separately.
"""

function _dh_parse_int(tok::AbstractString, line_num::Int, field::AbstractString)::Int
    value = tryparse(Int, tok)
    value === nothing && error("line $line_num: invalid integer for $field: '$tok'")
    return value
end

function _dh_header_value(
    lines::Vector{SubString{String}},
    line_num::Int,
    fallback_name::AbstractString,
)::Int
    length(lines) >= line_num || error("missing header line $line_num")
    tokens = split_def_line(lines[line_num])
    length(tokens) >= 2 || error("line $line_num: expected '$fallback_name <int>'")
    return _dh_parse_int(tokens[2], line_num, fallback_name)
end

function _dh_body_rows(lines::Vector{SubString{String}})
    rows = Tuple{Int,Vector{String}}[]
    for line_num = 6:length(lines)
        tokens = split_def_line(lines[line_num])
        isempty(tokens) && continue
        push!(rows, (line_num, tokens))
    end
    return rows
end

function _dh_check_site(site::Int, nsite::Int, line_num::Int, field::AbstractString)
    0 <= site < nsite || error("line $line_num: $field site $site out of range [0, $(nsite - 1)]")
end

function _dh_check_opt_flag(flag::Int, line_num::Int)
    (flag == 0 || flag == 1) || error("line $line_num: opt flag must be 0 or 1, got $flag")
end

"""
    parse_doublon_holon_2site_def(filepath::String, nsite::Int)

Parse a C-compatible DH2 definition file. The optimization-table first column
is intentionally ignored, matching C `GetInfoOpt`; file row order controls the
local DH parameter order.
"""
function parse_doublon_holon_2site_def(
    filepath::String,
    nsite::Int,
)::ParseResult{DoublonHolon2SiteDefinition}
    try
        content = read_def_file(filepath)
        return parse_doublon_holon_2site_content(content, nsite)
    catch e
        return ParseResult{DoublonHolon2SiteDefinition}(
            false,
            nothing,
            "Error reading/parsing file: $e",
            0,
        )
    end
end

function parse_doublon_holon_2site_content(
    content::String,
    nsite::Int,
)::ParseResult{DoublonHolon2SiteDefinition}
    context = ParsingContext("DH2")
    try
        nsite > 0 || error("Nsite must be positive before parsing DH2")
        lines = split(content, '\n')
        length(lines) >= 5 || error("DH2 file must include 5 header lines")

        n_dh2 = _dh_header_value(lines, 2, "NDoublonHolon2siteIdx")
        complex_type = _dh_header_value(lines, 3, "ComplexType")
        n_dh2 >= 0 || error("NDoublonHolon2siteIdx must be non-negative")
        is_complex = complex_type != 0

        rows = _dh_body_rows(lines)
        expected_main = nsite * n_dh2
        expected_opt = 6 * n_dh2
        expected_total = expected_main + expected_opt
        length(rows) == expected_total || error(
            "DH2 row count mismatch: got $(length(rows)), expected $expected_total " *
            "($expected_main neighbor rows + $expected_opt opt rows)",
        )

        indices = [DoublonHolon2SiteIndex(fill(-1, nsite, 2)) for _ = 1:n_dh2]
        seen = falses(n_dh2, nsite)

        for row_i = 1:expected_main
            line_num, tokens = rows[row_i]
            context.line_number = line_num
            length(tokens) == 4 || error("line $line_num: expected 'i x0 x1 n'")
            site = _dh_parse_int(tokens[1], line_num, "center site")
            x0 = _dh_parse_int(tokens[2], line_num, "neighbor 0")
            x1 = _dh_parse_int(tokens[3], line_num, "neighbor 1")
            idx = _dh_parse_int(tokens[4], line_num, "DH2 index")
            _dh_check_site(site, nsite, line_num, "center")
            _dh_check_site(x0, nsite, line_num, "neighbor 0")
            _dh_check_site(x1, nsite, line_num, "neighbor 1")
            0 <= idx < n_dh2 || error("line $line_num: DH2 index $idx out of range [0, $(n_dh2 - 1)]")
            !seen[idx+1, site+1] || error("line $line_num: duplicate DH2 row for index $idx site $site")
            seen[idx+1, site+1] = true
            indices[idx+1].neighbors[site+1, 1] = x0
            indices[idx+1].neighbors[site+1, 2] = x1
        end

        if any(x -> !x, seen)
            missing = findfirst(x -> !x, seen)
            error("missing DH2 neighbor row for index $(missing[1] - 1) site $(missing[2] - 1)")
        end

        opt_flags = Vector{Bool}(undef, expected_opt)
        for local_i = 1:expected_opt
            line_num, tokens = rows[expected_main+local_i]
            context.line_number = line_num
            length(tokens) == 2 || error("line $line_num: expected 'local_param_index opt_flag'")
            _dh_parse_int(tokens[1], line_num, "ignored local parameter index")
            flag = _dh_parse_int(tokens[2], line_num, "opt flag")
            _dh_check_opt_flag(flag, line_num)
            opt_flags[local_i] = flag != 0
        end

        data = DoublonHolon2SiteDefinition(indices, opt_flags, is_complex)
        return ParseResult{DoublonHolon2SiteDefinition}(true, data, "", context.line_number)
    catch e
        return ParseResult{DoublonHolon2SiteDefinition}(
            false,
            nothing,
            sprint(showerror, e),
            context.line_number,
        )
    end
end

"""
    parse_doublon_holon_4site_def(filepath::String, nsite::Int)

Parse a C-compatible DH4 definition file. The optimization-table first column
is intentionally ignored, matching C `GetInfoOpt`; file row order controls the
local DH parameter order.
"""
function parse_doublon_holon_4site_def(
    filepath::String,
    nsite::Int,
)::ParseResult{DoublonHolon4SiteDefinition}
    try
        content = read_def_file(filepath)
        return parse_doublon_holon_4site_content(content, nsite)
    catch e
        return ParseResult{DoublonHolon4SiteDefinition}(
            false,
            nothing,
            "Error reading/parsing file: $e",
            0,
        )
    end
end

function parse_doublon_holon_4site_content(
    content::String,
    nsite::Int,
)::ParseResult{DoublonHolon4SiteDefinition}
    context = ParsingContext("DH4")
    try
        nsite > 0 || error("Nsite must be positive before parsing DH4")
        lines = split(content, '\n')
        length(lines) >= 5 || error("DH4 file must include 5 header lines")

        n_dh4 = _dh_header_value(lines, 2, "NDoublonHolon4siteIdx")
        complex_type = _dh_header_value(lines, 3, "ComplexType")
        n_dh4 >= 0 || error("NDoublonHolon4siteIdx must be non-negative")
        is_complex = complex_type != 0

        rows = _dh_body_rows(lines)
        expected_main = nsite * n_dh4
        expected_opt = 10 * n_dh4
        expected_total = expected_main + expected_opt
        length(rows) == expected_total || error(
            "DH4 row count mismatch: got $(length(rows)), expected $expected_total " *
            "($expected_main neighbor rows + $expected_opt opt rows)",
        )

        indices = [DoublonHolon4SiteIndex(fill(-1, nsite, 4)) for _ = 1:n_dh4]
        seen = falses(n_dh4, nsite)

        for row_i = 1:expected_main
            line_num, tokens = rows[row_i]
            context.line_number = line_num
            length(tokens) == 6 || error("line $line_num: expected 'i x0 x1 x2 x3 n'")
            site = _dh_parse_int(tokens[1], line_num, "center site")
            x0 = _dh_parse_int(tokens[2], line_num, "neighbor 0")
            x1 = _dh_parse_int(tokens[3], line_num, "neighbor 1")
            x2 = _dh_parse_int(tokens[4], line_num, "neighbor 2")
            x3 = _dh_parse_int(tokens[5], line_num, "neighbor 3")
            idx = _dh_parse_int(tokens[6], line_num, "DH4 index")
            _dh_check_site(site, nsite, line_num, "center")
            _dh_check_site(x0, nsite, line_num, "neighbor 0")
            _dh_check_site(x1, nsite, line_num, "neighbor 1")
            _dh_check_site(x2, nsite, line_num, "neighbor 2")
            _dh_check_site(x3, nsite, line_num, "neighbor 3")
            0 <= idx < n_dh4 || error("line $line_num: DH4 index $idx out of range [0, $(n_dh4 - 1)]")
            !seen[idx+1, site+1] || error("line $line_num: duplicate DH4 row for index $idx site $site")
            seen[idx+1, site+1] = true
            indices[idx+1].neighbors[site+1, 1] = x0
            indices[idx+1].neighbors[site+1, 2] = x1
            indices[idx+1].neighbors[site+1, 3] = x2
            indices[idx+1].neighbors[site+1, 4] = x3
        end

        if any(x -> !x, seen)
            missing = findfirst(x -> !x, seen)
            error("missing DH4 neighbor row for index $(missing[1] - 1) site $(missing[2] - 1)")
        end

        opt_flags = Vector{Bool}(undef, expected_opt)
        for local_i = 1:expected_opt
            line_num, tokens = rows[expected_main+local_i]
            context.line_number = line_num
            length(tokens) == 2 || error("line $line_num: expected 'local_param_index opt_flag'")
            _dh_parse_int(tokens[1], line_num, "ignored local parameter index")
            flag = _dh_parse_int(tokens[2], line_num, "opt flag")
            _dh_check_opt_flag(flag, line_num)
            opt_flags[local_i] = flag != 0
        end

        data = DoublonHolon4SiteDefinition(indices, opt_flags, is_complex)
        return ParseResult{DoublonHolon4SiteDefinition}(true, data, "", context.line_number)
    catch e
        return ParseResult{DoublonHolon4SiteDefinition}(
            false,
            nothing,
            sprint(showerror, e),
            context.line_number,
        )
    end
end

"""
    parse_doublon_holon_2site_def(filepath::String)

Deprecated compatibility shim for the pre-DH-1 value-bearing DH2 parser. The
runtime parser is `parse_doublon_holon_2site_def(filepath, nsite)`.
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
        length(tokens) < 4 && continue

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

    success = isempty(context.errors)
    return ParseResult{Vector{DoublonHolon2SiteTerm}}(
        success,
        success ? terms : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
end

"""
    parse_doublon_holon_4site_def(filepath::String)

Deprecated compatibility shim for the pre-DH-1 value-bearing DH4 parser. The
runtime parser is `parse_doublon_holon_4site_def(filepath, nsite)`.
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
        length(tokens) < 5 && continue

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

    success = isempty(context.errors)
    return ParseResult{Vector{DoublonHolon4SiteTerm}}(
        success,
        success ? terms : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
end
