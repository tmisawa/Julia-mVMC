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

"""
    parse_opttrans_def!(data::ExpertModeData, filepath::String) -> ExpertModeData

Parse C-mVMC's `opttrans.def` / `QPOptTrans` file and initialise both
`data.para_qp_opt_trans` and the mutable runtime `data.opt_trans`.
"""
function parse_opttrans_def!(data::ExpertModeData, filepath::String)
    n_site = data.modpara.nsite
    n_site > 0 || error("OptTrans requires ModPara/Nsite to be parsed before $filepath")

    content = read_def_file(filepath)
    lines = split(content, '\n')
    ignore_lines = 5
    length(lines) > ignore_lines || error("$filepath must include 5 header lines")

    header_tokens = split_def_line(clean_line(lines[2]))
    length(header_tokens) >= 2 || error("$filepath:2 missing NQPOptTrans header")
    n_qp_opt_trans = _parse_int_strict(header_tokens[2], filepath, 2, "NQPOptTrans")
    n_qp_opt_trans >= 1 || error("NQPOptTrans should be larger than 0 in $filepath")

    para = fill(0.0 + 0.0im, n_qp_opt_trans)
    seen_para = falses(n_qp_opt_trans)
    line_idx = ignore_lines + 1

    for _ = 1:n_qp_opt_trans
        while line_idx <= length(lines) && isempty(clean_line(lines[line_idx]))
            line_idx += 1
        end
        line_idx <= length(lines) || error("$filepath ended before ParaQPOptTrans block was complete")

        tokens = split_def_line(clean_line(lines[line_idx]))
        length(tokens) >= 2 || error("$filepath:$line_idx expected 'idx value'")
        idx = _parse_int_strict(tokens[1], filepath, line_idx, "ParaQPOptTrans index")
        0 <= idx < n_qp_opt_trans ||
            error("$filepath:$line_idx ParaQPOptTrans index $idx out of range [0, $(n_qp_opt_trans - 1)]")
        !seen_para[idx+1] || error("$filepath:$line_idx duplicated ParaQPOptTrans index $idx")
        value = _parse_float_strict(tokens[2], filepath, line_idx, "ParaQPOptTrans value")
        seen_para[idx+1] = true
        para[idx+1] = ComplexF64(value, 0.0)
        line_idx += 1
    end

    missing_para = findfirst(x -> !x, seen_para)
    missing_para === nothing ||
        error("$filepath missing ParaQPOptTrans index $(missing_para - 1)")

    qp_opt_trans = [fill(-1, n_site) for _ = 1:n_qp_opt_trans]
    qp_opt_trans_sgn = [ones(Int, n_site) for _ = 1:n_qp_opt_trans]
    seen_map = falses(n_qp_opt_trans, n_site)

    while line_idx <= length(lines)
        line = clean_line(lines[line_idx])
        if isempty(line)
            line_idx += 1
            continue
        end

        tokens = split_def_line(line)
        length(tokens) >= 4 || error("$filepath:$line_idx expected 'optidx site mapped sign'")
        optidx = _parse_int_strict(tokens[1], filepath, line_idx, "QPOptTrans index")
        site = _parse_int_strict(tokens[2], filepath, line_idx, "site index")
        mapped = _parse_int_strict(tokens[3], filepath, line_idx, "mapped site index")
        sign = _parse_int_strict(tokens[4], filepath, line_idx, "QPOptTrans sign")

        0 <= optidx < n_qp_opt_trans ||
            error("$filepath:$line_idx QPOptTrans index $optidx out of range [0, $(n_qp_opt_trans - 1)]")
        0 <= site < n_site || error("$filepath:$line_idx site index $site out of range [0, $(n_site - 1)]")
        0 <= mapped < n_site || error("$filepath:$line_idx mapped site index $mapped out of range [0, $(n_site - 1)]")
        sign == 1 || sign == -1 || error("$filepath:$line_idx sign must be +1 or -1")
        !seen_map[optidx+1, site+1] || error("$filepath:$line_idx duplicated QPOptTrans entry optidx=$optidx site=$site")

        seen_map[optidx+1, site+1] = true
        qp_opt_trans[optidx+1][site+1] = mapped
        qp_opt_trans_sgn[optidx+1][site+1] = sign
        line_idx += 1
    end

    for optidx = 1:n_qp_opt_trans
        for site = 1:n_site
            if !seen_map[optidx, site]
                error("$filepath missing QPOptTrans entry optidx=$(optidx - 1) site=$(site - 1)")
            end
        end
    end

    if data.modpara.nmp_trans >= 0
        for signs in qp_opt_trans_sgn
            fill!(signs, 1)
        end
    end

    data.n_qp_opt_trans = n_qp_opt_trans
    data.para_qp_opt_trans = para
    data.opt_trans = copy(para)
    data.qp_opt_trans = qp_opt_trans
    data.qp_opt_trans_sgn = qp_opt_trans_sgn

    return data
end
