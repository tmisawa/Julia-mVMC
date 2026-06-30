"""
Pair Hopping Parser

Parses pairhop.def files for pair hopping terms.
"""

"""
    parse_pairhop_def(file_path::String) -> ParseResult{Vector{PairHopTerm}}

Parse a pairhop.def file and return the pair hopping terms.

# Arguments
- `file_path`: Path to the pairhop.def file

# Returns
- `ParseResult{Vector{PairHopTerm}}`: Parsed pair hopping terms

# Example
```julia
result = parse_pairhop_def("pairhop.def")
if result.success
    pairhop_terms = result.data
end
```
"""
function parse_pairhop_def(file_path::String)::ParseResult{Vector{PairHopTerm}}
    try
        content = read_def_file(file_path)
        return parse_pairhop_content(content)
    catch e
        return ParseResult{Vector{PairHopTerm}}(false, nothing, "Error reading file: $e", 0)
    end
end

"""
    parse_pairhop_content(content::String) -> ParseResult{Vector{PairHopTerm}}

Parse pairhop.def content and return the pair hopping terms.

C-mVMC expands each input PairHop line `(i, j, value)` into two internal terms,
`(i, j, value)` and `(j, i, value)`. Return the expanded internal list so the
runtime sees the same term ordering as C.

# Arguments
- `content`: Content of the pairhop.def file

# Returns
- `ParseResult{Vector{PairHopTerm}}`: Parsed pair hopping terms
"""
function parse_pairhop_content(content::String)::ParseResult{Vector{PairHopTerm}}
    pairhop_terms = PairHopTerm[]
    errors = String[]

    lines = split(content, '\n')

    # Skip header lines (first 5 lines) for StdFace-generated pairhopp.def.
    IGNORE_LINES_IN_DEF = 5
    start_line = 1
    if length(lines) > IGNORE_LINES_IN_DEF
        first_line = clean_line(lines[1])
        second_line = clean_line(lines[2])
        if startswith(first_line, "=") || startswith(second_line, "NPairHopp")
            start_line = IGNORE_LINES_IN_DEF + 1
        end
    end

    line_number = 0
    for line_num = start_line:length(lines)
        line = lines[line_num]
        line_number = line_num
        clean_line_str = clean_line(line)

        # Skip empty lines and comments
        if isempty(clean_line_str) || startswith(clean_line_str, "#")
            continue
        end

        try
            # Parse line: site1 site2 value
            parts = split(clean_line_str)
            if length(parts) < 3
                push!(
                    errors,
                    "Line $line_number: Invalid format, expected 'site1 site2 value'",
                )
                continue
            end

            # Parse site1, site2, and value
            site1 = try
                parse(Int, parts[1])
            catch
                push!(errors, "Line $line_number: Invalid site1 number '$(parts[1])'")
                continue
            end

            site2 = try
                parse(Int, parts[2])
            catch
                push!(errors, "Line $line_number: Invalid site2 number '$(parts[2])'")
                continue
            end

            value = try
                parse(Float64, parts[3])
            catch
                push!(errors, "Line $line_number: Invalid value '$(parts[3])'")
                continue
            end

            # Validate site numbers (0-based indexing)
            if site1 < 0
                push!(
                    errors,
                    "Line $line_number: Site1 number must be non-negative, got $site1",
                )
                continue
            end

            if site2 < 0
                push!(
                    errors,
                    "Line $line_number: Site2 number must be non-negative, got $site2",
                )
                continue
            end

            push!(pairhop_terms, PairHopTerm(site1, site2, value))
            push!(pairhop_terms, PairHopTerm(site2, site1, value))

        catch e
            push!(errors, "Line $line_number: Error parsing line '$clean_line_str': $e")
        end
    end

    success = isempty(errors)
    return ParseResult{Vector{PairHopTerm}}(
        success,
        pairhop_terms,
        isempty(errors) ? "" : join(errors, "; "),
        0,
    )
end
