"""
Trans.def Parser

Parser for trans.def files containing transfer (hopping) terms.
"""

"""
    parse_trans_def(filepath::String) -> ParseResult{Vector{TransferTerm}}

Parse trans.def file from file path.
"""
function parse_trans_def(filepath::String)::ParseResult{Vector{TransferTerm}}
    try
        content = read_def_file(filepath)
        return parse_trans_content(content)
    catch e
        return ParseResult{Vector{TransferTerm}}(
            false,
            nothing,
            "Error reading file: $e",
            0,
        )
    end
end

"""
    parse_trans_content(content::String) -> ParseResult{Vector{TransferTerm}}

Parse trans.def content from string.
"""
function parse_trans_content(content::String)::ParseResult{Vector{TransferTerm}}
    context = ParsingContext("trans.def")
    terms = TransferTerm[]

    lines = split(content, '\n')

    for (line_num, line) in enumerate(lines)
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
            term = parse_transfer_term(tokens, context)
            if term !== nothing
                push!(terms, term)
            end
        catch e
            push!(context.errors, "Line $line_num: Error parsing transfer term: $e")
        end
    end

    success = length(context.errors) == 0
    return ParseResult{Vector{TransferTerm}}(
        success,
        success ? terms : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
end

"""
    parse_transfer_term(tokens::Vector{String}, context::ParsingContext) -> Union{TransferTerm, Nothing}

Parse a single transfer term from tokens.
Format: site1 spin1 site2 spin2 real_value imag_value
Where spin1 and spin2 should be equal (same-spin hopping).
"""
function parse_transfer_term(
    tokens::Vector{String},
    context::ParsingContext,
)::Union{TransferTerm,Nothing}
    if length(tokens) < 6
        push!(
            context.warnings,
            "Insufficient tokens for transfer term (need 6: site1 spin1 site2 spin2 real imag)",
        )
        return nothing
    end

    # Parse site indices and spin indices
    # Format: site1 spin1 site2 spin2 real imag
    site1 = safe_parse_int(tokens[1], -1)
    spin1 = safe_parse_int(tokens[2], -1)
    site2 = safe_parse_int(tokens[3], -1)
    spin2 = safe_parse_int(tokens[4], -1)

    if site1 < 0 || site2 < 0
        push!(context.errors, "Invalid site indices: $site1, $site2")
        return nothing
    end

    if spin1 < 0 || spin1 > 1 || spin2 < 0 || spin2 > 1
        push!(context.errors, "Invalid spin indices: $spin1, $spin2 (must be 0 or 1)")
        return nothing
    end

    # Validate that spin1 == spin2 (same-spin hopping)
    if spin1 != spin2
        push!(
            context.warnings,
            "Spin indices differ ($spin1 != $spin2), using spin1=$spin1",
        )
    end

    # Parse complex value (real and imaginary parts)
    real_val = safe_parse_float(tokens[5], 0.0)
    imag_val = safe_parse_float(tokens[6], 0.0)
    value = ComplexF64(real_val, imag_val)

    return TransferTerm(site1, spin1, site2, spin2, value)
end
