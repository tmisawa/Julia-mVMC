"""
Exchange Coupling Parser

Parses exchange.def files for exchange coupling terms.
"""

"""
    parse_exchange_def(file_path::String) -> ParseResult{Vector{ExchangeTerm}}

Parse an exchange.def file and return the exchange coupling terms.

# Arguments
- `file_path`: Path to the exchange.def file

# Returns
- `ParseResult{Vector{ExchangeTerm}}`: Parsed exchange coupling terms

# Example
```julia
result = parse_exchange_def("exchange.def")
if result.success
    exchange_terms = result.data
end
```
"""
function parse_exchange_def(file_path::String)::ParseResult{Vector{ExchangeTerm}}
    try
        content = read_def_file(file_path)
        return parse_exchange_content(content)
    catch e
        return ParseResult{Vector{ExchangeTerm}}(
            false,
            nothing,
            "Error reading file: $e",
            0,
        )
    end
end

"""
    parse_exchange_content(content::String) -> ParseResult{Vector{ExchangeTerm}}

Parse exchange.def content and return the exchange coupling terms.

# Arguments
- `content`: Content of the exchange.def file

# Returns
- `ParseResult{Vector{ExchangeTerm}}`: Parsed exchange coupling terms
"""
function parse_exchange_content(content::String)::ParseResult{Vector{ExchangeTerm}}
    exchange_terms = ExchangeTerm[]
    errors = String[]

    lines = split(content, '\n')

    # Skip header lines (first 5 lines) for exchange.def format
    IGNORE_LINES_IN_DEF = 5
    start_line = 1

    # Check if this looks like exchange.def format (has header)
    if length(lines) > IGNORE_LINES_IN_DEF
        first_data_line = clean_line(lines[IGNORE_LINES_IN_DEF+1])
        if !isempty(first_data_line)
            parts = split(first_data_line)
            # If first data line has at least 3 tokens (site1 site2 value), skip header
            if length(parts) >= 3
                start_line = IGNORE_LINES_IN_DEF + 1
            end
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

            # Check for self-interaction
            if site1 == site2
                push!(
                    errors,
                    "Line $line_number: Self-interaction not allowed for exchange coupling (site1 = site2 = $site1)",
                )
                continue
            end

            push!(exchange_terms, ExchangeTerm(site1, site2, value))

        catch e
            push!(errors, "Line $line_number: Error parsing line '$clean_line_str': $e")
        end
    end

    success = isempty(errors)
    return ParseResult{Vector{ExchangeTerm}}(
        success,
        exchange_terms,
        isempty(errors) ? "" : join(errors, "; "),
        0,
    )
end
