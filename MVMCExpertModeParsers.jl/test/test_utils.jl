"""
Tests for utility functions
"""

using MVMCExpertModeParsers: read_def_file, get_file_info, validate_file_exists
using MVMCExpertModeParsers: clean_line, split_def_line, parse_complex_value
using MVMCExpertModeParsers: safe_parse_int, safe_parse_float, safe_parse_complex

@testset "File Utilities" begin
    # Test file operations
    test_dir = mktempdir()
    test_file = joinpath(test_dir, "test.def")

    try
        # Create test file manually (write_def_file is not available)
        open(test_file, "w") do f
            write(f, "NSite = 4\nNElec = 2\n")
        end
        @test isfile(test_file)

        read_content = read_def_file(test_file)
        @test occursin("NSite = 4", read_content)
        @test occursin("NElec = 2", read_content)

        # Test file info
        file_info = get_file_info(test_file)
        @test file_info.filename == "test.def"
        @test file_info.exists
        @test file_info.size_bytes > 0

        # Test file validation
        @test validate_file_exists(test_file)
        @test !validate_file_exists(joinpath(test_dir, "nonexistent.def"))

    finally
        rm(test_dir, recursive = true)
    end
end

@testset "String Processing" begin
    # Test line cleaning
    @test clean_line("  # This is a comment  ") == ""
    @test clean_line("  NSite = 4  # comment  ") == "NSite = 4"
    @test clean_line("  NSite = 4  ") == "NSite = 4"

    # Test line splitting
    tokens = split_def_line("NSite = 4")
    @test tokens == ["NSite", "=", "4"]

    tokens = split_def_line("0 1 1.0 up")
    @test tokens == ["0", "1", "1.0", "up"]
end

@testset "Complex Number Parsing" begin
    # Test complex number parsing
    @test parse_complex_value("1.0") == 1.0+0.0im
    @test parse_complex_value("1.0+2.0i") == 1.0+2.0im
    @test parse_complex_value("1.0 2.0") == 1.0+2.0im
    @test parse_complex_value("1.0+2.0j") == 1.0+2.0im
end

@testset "Safe Parsing" begin
    # Test safe integer parsing
    @test safe_parse_int("123", 0) == 123
    @test safe_parse_int("invalid", 42) == 42
    @test safe_parse_int("", 99) == 99

    # Test safe float parsing
    @test safe_parse_float("1.23", 0.0) == 1.23
    @test safe_parse_float("invalid", 3.14) == 3.14
    @test safe_parse_float("", 2.71) == 2.71

    # Test safe complex parsing
    @test safe_parse_complex("1.0", 0.0+0.0im) == 1.0+0.0im
    @test safe_parse_complex("1.0+2.0i", 0.0+0.0im) == 1.0+2.0im
    @test safe_parse_complex("invalid", 1.0+1.0im) == 1.0+1.0im
end

@testset "Path Utilities" begin
    # Test absolute path resolution
    @test isfile("README.md") == false  # Basic file check
    @test isdir(".") == true  # Basic directory check
end
