"""
QPTrans Parser

Parser for qptrans.def files containing quantum projection translation terms.
"""

"""
    parse_qptrans_def(filepath::String) -> ParseResult{Vector{QPTransTerm}}

Parse qptrans.def file from file path.
"""
function parse_qptrans_def(filepath::String)::ParseResult{Vector{QPTransTerm}}
    try
        content = read_def_file(filepath)
        return parse_qptrans_content(content)
    catch e
        return ParseResult{Vector{QPTransTerm}}(false, nothing, "Error reading file: $e", 0)
    end
end

"""
    parse_qptrans_content(content::String) -> ParseResult{Vector{QPTransTerm}}

Parse qptrans.def content from string.
"""
function parse_qptrans_content(content::String)::ParseResult{Vector{QPTransTerm}}
    context = ParsingContext("qptrans.def")
    terms = QPTransTerm[]

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
            term = parse_qptrans_term(tokens, context)
            if term !== nothing
                push!(terms, term)
            end
        catch e
            push!(context.errors, "Line $line_num: Error parsing QPTrans term: $e")
        end
    end

    success = length(context.errors) == 0
    return ParseResult{Vector{QPTransTerm}}(
        success,
        success ? terms : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
end

"""
    parse_qptrans_term(tokens::Vector{String}, context::ParsingContext) -> Union{QPTransTerm, Nothing}

Parse a single QPTrans term from tokens.
"""
function parse_qptrans_term(
    tokens::Vector{String},
    context::ParsingContext,
)::Union{QPTransTerm,Nothing}
    if length(tokens) < 2
        push!(context.warnings, "Insufficient tokens for QPTrans term")
        return nothing
    end

    # Parse site index
    site = safe_parse_int(tokens[1], -1)

    if site < 0
        push!(context.errors, "Invalid site index: $site")
        return nothing
    end

    # Parse momentum vector (can be 1D, 2D, or 3D)
    momentum = Float64[]
    if length(tokens) >= 2
        for i = 2:length(tokens)
            push!(momentum, safe_parse_float(tokens[i]))
        end
    end

    # Parse phase (optional, default to 0.0)
    phase = 0.0
    if length(tokens) >= 3
        phase = safe_parse_float(tokens[end])
    end

    return QPTransTerm(site, momentum, phase)
end

"""
    parse_qptransidx_def(filepath::String) -> ParseResult{Vector{QPTransTerm}}

Parse qptransidx.def file from file path.
Format: first 5 lines are header, then NQPTrans lines are idx value (ParaQPTrans),
then Nsite * NQPTrans lines are i j itmp itmpsgn (QPTrans indices).
"""
function parse_qptransidx_def(filepath::String)::ParseResult{Vector{QPTransTerm}}
    try
        content = read_def_file(filepath)
        return parse_qptransidx_content(content)
    catch e
        return ParseResult{Vector{QPTransTerm}}(false, nothing, "Error reading file: $e", 0)
    end
end

"""
    parse_qptransidx_content(content::String) -> ParseResult{Vector{QPTransTerm}}

Parse qptransidx.def content from string.
"""
function parse_qptransidx_content(content::String)::ParseResult{Vector{QPTransTerm}}
    context = ParsingContext("qptransidx.def")
    terms = QPTransTerm[]

    lines = split(content, '\n')

    # Skip header lines (first 5 lines)
    IGNORE_LINES_IN_DEF = 5
    if length(lines) <= IGNORE_LINES_IN_DEF
        return ParseResult{Vector{QPTransTerm}}(false, nothing, "File has too few lines", 0)
    end

    # Read NQPTrans from header (line 2: "NQPTrans N")
    n_qp_trans = 0
    if length(lines) > 1
        header_line = clean_line(lines[2])
        tokens = split_def_line(header_line)
        if length(tokens) >= 2 && tokens[1] == "NQPTrans"
            n_qp_trans = safe_parse_int(tokens[2], 0)
        end
    end

    if n_qp_trans == 0
        return ParseResult{Vector{QPTransTerm}}(
            false,
            nothing,
            "Could not parse NQPTrans from header",
            0,
        )
    end

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
                    # Store ParaQPTrans value (imaginary part is 0)
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

    # Create QPTransTerm for each ParaQPTrans value
    # Each term represents one quantum projection translation
    for i = 1:n_qp_trans
        if i <= length(para_qp_trans)
            momentum = Float64[]  # Empty momentum vector
            phase = real(para_qp_trans[i])  # Use ParaQPTrans value as phase
            # Use site=0 as default (actual site mapping is in the index part)
            term = QPTransTerm(0, momentum, phase)
            push!(terms, term)
        end
    end

    success = length(context.errors) == 0
    return ParseResult{Vector{QPTransTerm}}(
        success,
        success ? terms : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
end
