"""
InterAll Parser

Parses interall.def files for general interaction terms.

The format of interall.def is:
    site0 spin0 site1 spin1 site2 spin2 site3 spin3 real_value imag_value

where site0-site3 are site indices (0-based) and spin0-spin3 are spin indices (0 or 1).
The sites are stored as combined indices: site + spin * Nsite.
"""

"""
    parse_interall_def(file_path::String) -> ParseResult{Vector{InterAllTerm}}

Parse an interall.def file and return the general interaction terms.

# Arguments
- `file_path`: Path to the interall.def file

# Returns
- `ParseResult{Vector{InterAllTerm}}`: Parsed general interaction terms

# Example
```julia
result = parse_interall_def("interall.def")
if result.success
    interall_terms = result.data
end
```
"""
function parse_interall_def(file_path::String)::ParseResult{Vector{InterAllTerm}}
    try
        content = read_def_file(file_path)
        return parse_interall_content(content)
    catch e
        return ParseResult{Vector{InterAllTerm}}(
            false,
            nothing,
            "Error reading file: $e",
            0,
        )
    end
end

"""
    parse_interall_content(content::String) -> ParseResult{Vector{InterAllTerm}}

Parse interall.def content and return the general interaction terms.

# Arguments
- `content`: Content of the interall.def file

# Returns
- `ParseResult{Vector{InterAllTerm}}`: Parsed general interaction terms
"""
function parse_interall_content(content::String)::ParseResult{Vector{InterAllTerm}}
    interall_terms = InterAllTerm[]
    errors = String[]

    lines = split(content, '\n')
    line_number = 0

    for line in lines
        line_number += 1
        clean_line_str = clean_line(line)

        # Skip empty lines, comments, and header lines
        if isempty(clean_line_str) ||
           startswith(clean_line_str, "#") ||
           startswith(clean_line_str, "=")
            continue
        end

        # Skip lines that look like headers (contain letters but not just numbers)
        if occursin(r"[A-Za-z]", clean_line_str)
            continue
        end

        try
            # Parse line: site0 spin0 site1 spin1 site2 spin2 site3 spin3 real_value imag_value
            parts = split(clean_line_str)
            if length(parts) < 10
                # Skip lines with fewer than 10 parts (not data lines)
                continue
            end

            # Parse all 10 values (matching C implementation exactly)
            site0 = parse(Int, parts[1])
            spin0 = parse(Int, parts[2])
            site1 = parse(Int, parts[3])
            spin1 = parse(Int, parts[4])
            site2 = parse(Int, parts[5])
            spin2 = parse(Int, parts[6])
            site3 = parse(Int, parts[7])
            spin3 = parse(Int, parts[8])
            real_val = parse(Float64, parts[9])
            imag_val = parse(Float64, parts[10])

            value = complex(real_val, imag_val)
            is_complex = abs(imag_val) > 1e-14

            push!(
                interall_terms,
                InterAllTerm(
                    site0,
                    spin0,
                    site1,
                    spin1,
                    site2,
                    spin2,
                    site3,
                    spin3,
                    value,
                    is_complex,
                ),
            )

        catch e
            # Skip non-data lines silently
            continue
        end
    end

    return ParseResult{Vector{InterAllTerm}}(true, interall_terms, "", 0)
end
