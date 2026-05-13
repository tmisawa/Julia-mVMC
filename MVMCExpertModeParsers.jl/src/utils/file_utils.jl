"""
File Utilities for Expert Mode Parsing

File I/O operations and path management utilities.
"""

"""
    read_def_file(filepath::String) -> String

Read the entire contents of a definition file.
"""
function read_def_file(filepath::String)::String
    if !isfile(filepath)
        error("File not found: $filepath")
    end

    content = ""
    open(filepath, "r") do file
        content = read(file, String)
    end

    return content
end

"""
    get_file_info(filepath::String) -> FileInfo

Get information about a file.
"""
function get_file_info(filepath::String)::FileInfo
    exists = isfile(filepath)
    size_bytes = exists ? stat(filepath).size : 0
    last_modified = exists ? unix2datetime(stat(filepath).mtime) : DateTime(1970)

    return FileInfo(basename(filepath), filepath, exists, size_bytes, last_modified)
end

"""
    validate_file_exists(filepath::String) -> Bool

Check if a file exists and is readable.
"""
function validate_file_exists(filepath::String)::Bool
    return isfile(filepath) && isreadable(filepath)
end

"""
    clean_line(line::String) -> String

Clean a line by removing comments and extra whitespace.
"""
function clean_line(line::AbstractString)::String
    # Remove comments (lines starting with # or //)
    if startswith(strip(line), "#") || startswith(strip(line), "//")
        return ""
    end

    # Remove inline comments
    comment_pos = findfirst("#", line)
    if comment_pos !== nothing
        line = line[1:(comment_pos[1]-1)]
    end

    comment_pos = findfirst("//", line)
    if comment_pos !== nothing
        line = line[1:(comment_pos[1]-1)]
    end

    # Strip whitespace
    return strip(line)
end

"""
    split_def_line(line::String) -> Vector{String}

Split a definition file line into tokens.
"""
function split_def_line(line::AbstractString)::Vector{String}
    # Clean the line first
    clean_line_str = clean_line(line)

    if isempty(clean_line_str)
        return String[]
    end

    # Split by whitespace
    tokens = split(clean_line_str)

    # Filter out empty tokens
    return filter(!isempty, tokens)
end

"""
    parse_complex_value(str::String) -> ComplexF64

Parse a complex number from string format.
"""
function parse_complex_value(str::String)::ComplexF64
    # Handle different formats: "1.0+2.0i", "1.0 2.0", "1.0", etc.
    if occursin("i", str) || occursin("j", str)
        # Complex format with i or j
        return parse(ComplexF64, str)
    elseif occursin(" ", str)
        # Space-separated real and imaginary parts
        parts = split(str)
        if length(parts) == 2
            real_part = parse(Float64, parts[1])
            imag_part = parse(Float64, parts[2])
            return ComplexF64(real_part, imag_part)
        else
            error("Invalid complex format: $str")
        end
    else
        # Real number
        return ComplexF64(parse(Float64, str), 0.0)
    end
end

"""
    safe_parse_int(str::String, default::Int = 0) -> Int

Safely parse an integer with a default value.
"""
function safe_parse_int(str::String, default::Int = 0)::Int
    try
        return parse(Int, str)
    catch
        return default
    end
end

"""
    parse_namelist_content(content::String) -> Vector{Tuple{String, String}}

Parse namelist content and extract file types and file names.
Returns a vector of tuples (file_type, file_name).
"""
function parse_namelist_content(content::String)::Vector{Tuple{String,String}}
    files = Tuple{String,String}[]
    lines = split(content, '\n')

    for line in lines
        clean_line_str = clean_line(line)
        if isempty(clean_line_str)
            continue
        end

        tokens = split_def_line(clean_line_str)
        if length(tokens) >= 2
            # Format: "Type filename.def"
            file_type = tokens[1]
            file_name = tokens[2]
            push!(files, (file_type, file_name))
        end
    end

    return files
end

"""
    safe_parse_float(str::String, default::Float64 = 0.0) -> Float64

Safely parse a float with a default value.
"""
function safe_parse_float(str::String, default::Float64 = 0.0)::Float64
    try
        return parse(Float64, str)
    catch
        return default
    end
end

"""
    safe_parse_complex(str::String, default::ComplexF64 = 0.0+0.0im) -> ComplexF64

Safely parse a complex number with a default value.
"""
function safe_parse_complex(str::String, default::ComplexF64 = 0.0+0.0im)::ComplexF64
    try
        return parse_complex_value(str)
    catch
        return default
    end
end
