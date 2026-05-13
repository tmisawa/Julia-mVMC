"""
Local Spin Parser

Parses locspn.def files for local spin information.
"""

"""
    parse_locspin_def(file_path::String) -> ParseResult{Vector{LocSpinTerm}}

Parse a locspn.def file and return the local spin terms.

# Arguments
- `file_path`: Path to the locspn.def file

# Returns
- `ParseResult{Vector{LocSpinTerm}}`: Parsed local spin terms

# Example
```julia
result = parse_locspin_def("locspn.def")
if result.success
    locspin_terms = result.data
end
```
"""
function parse_locspin_def(file_path::String)::ParseResult{Vector{LocSpinTerm}}
    try
        content = read_def_file(file_path)
        return parse_locspin_content(content)
    catch e
        return ParseResult{Vector{LocSpinTerm}}(false, nothing, "Error reading file: $e", 0)
    end
end

"""
    parse_locspin_content(content::String) -> ParseResult{Vector{LocSpinTerm}}

Parse locspn.def content and return the local spin terms.

# Arguments
- `content`: Content of the locspn.def file

# Returns
- `ParseResult{Vector{LocSpinTerm}}`: Parsed local spin terms
"""
function parse_locspin_content(content::String)::ParseResult{Vector{LocSpinTerm}}
    locspin_terms = LocSpinTerm[]
    errors = String[]

    lines = split(content, '\n')
    line_number = 0

    # Skip header lines (equivalent to C's IgnoreLinesInDef, typically 5 lines)
    # C implementation: for (i = 0; i < IgnoreLinesInDef; i++) fgets(...)
    IGNORE_LINES_IN_DEF = 5
    start_line = IGNORE_LINES_IN_DEF + 1

    for line_idx = start_line:length(lines)
        line = lines[line_idx]
        line_number = line_idx
        clean_line_str = clean_line(line)

        # Skip empty lines and comments
        if isempty(clean_line_str) || startswith(clean_line_str, "#")
            continue
        end

        try
            # Parse line: site spin_value
            parts = split(clean_line_str)
            if length(parts) < 2
                push!(
                    errors,
                    "Line $line_number: Invalid format, expected 'site spin_value'",
                )
                continue
            end

            # Parse site and spin value
            site = try
                parse(Int, parts[1])
            catch
                push!(errors, "Line $line_number: Invalid site number '$(parts[1])'")
                continue
            end

            spin_value = try
                parse(Int, parts[2])
            catch
                push!(errors, "Line $line_number: Invalid spin value '$(parts[2])'")
                continue
            end

            # Validate site number (0-based indexing)
            if site < 0
                push!(
                    errors,
                    "Line $line_number: Site number must be non-negative, got $site",
                )
                continue
            end

            # Validate spin value
            if spin_value < 0
                push!(
                    errors,
                    "Line $line_number: Spin value must be non-negative, got $spin_value",
                )
                continue
            end

            push!(locspin_terms, LocSpinTerm(site, spin_value))

        catch e
            push!(errors, "Line $line_number: Error parsing line '$clean_line_str': $e")
        end
    end

    success = isempty(errors)
    return ParseResult{Vector{LocSpinTerm}}(
        success,
        locspin_terms,
        isempty(errors) ? "" : join(errors, "; "),
        0,
    )
end
