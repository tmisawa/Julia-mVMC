"""
Read Input Parameters Tests

Tests for read_input_parameters! function (equivalent to C's ReadInputParameters()).
"""

using Test
using MVMCExpertModeParsers

import MVMCExpertModeParsers: ExpertModeData, GutzwillerTerm, JastrowTerm, OrbitalTerm
import MVMCExpertModeParsers: read_input_parameters!

"""
    test_read_input_parameters_basic()

Test basic functionality of read_input_parameters!().
"""
function test_read_input_parameters_basic()
    @testset "Basic read_input_parameters! Tests" begin
        # Create test data
        data = ExpertModeData()

        # Add Gutzwiller and Jastrow terms
        push!(data.gutzwiller_terms, GutzwillerTerm(0, 0.0+0.0im, false))
        push!(data.gutzwiller_terms, GutzwillerTerm(1, 0.0+0.0im, false))
        push!(data.jastrow_terms, JastrowTerm(0, 1, 0.0+0.0im, false))

        # Add Orbital terms
        for i = 1:3
            push!(data.orbital_terms, OrbitalTerm(0, i, 0.0+0.0im, false))
        end

        # Create a temporary namelist.def file
        test_dir = mktempdir()
        namelist_path = joinpath(test_dir, "namelist.def")

        # Create namelist.def with InGutzwiller entry
        open(namelist_path, "w") do f
            write(f, "ModPara modpara.def\n")
            write(f, "InGutzwiller ingutzwiller.def\n")
        end

        # Create InGutzwiller.def file
        ingutzwiller_path = joinpath(test_dir, "ingutzwiller.def")
        open(ingutzwiller_path, "w") do f
            write(f, "=============================================\n")
            write(f, "NGutzwillerIdx          2\n")
            write(f, "ComplexType         0\n")
            write(f, "=============================================\n")
            write(f, "=============================================\n")
            write(f, "0 0.5 0.0\n")  # idx=0, real=0.5, imag=0.0
            write(f, "1 0.3 0.0\n")  # idx=1, real=0.3, imag=0.0
        end

        # Read input parameters
        read_input_parameters!(data, namelist_path)

        # Test: Gutzwiller terms should be updated
        @test data.gutzwiller_terms[1].value == ComplexF64(0.5, 0.0)
        @test data.gutzwiller_terms[2].value == ComplexF64(0.3, 0.0)

        # Cleanup
        rm(test_dir, recursive = true)
    end
end

"""
    test_read_input_parameters_orbital()

Test read_input_parameters!() with InOrbital file.
"""
function test_read_input_parameters_orbital()
    @testset "read_input_parameters! Orbital Tests" begin
        data = ExpertModeData()

        # Add Orbital terms
        for i = 1:3
            push!(data.orbital_terms, OrbitalTerm(0, i, 0.0+0.0im, false))
        end

        # Create a temporary namelist.def file
        test_dir = mktempdir()
        namelist_path = joinpath(test_dir, "namelist.def")

        open(namelist_path, "w") do f
            write(f, "ModPara modpara.def\n")
            write(f, "InOrbital inorbital.def\n")
        end

        # Create InOrbital.def file
        inorbital_path = joinpath(test_dir, "inorbital.def")
        open(inorbital_path, "w") do f
            write(f, "=============================================\n")
            write(f, "NOrbitalIdx          3\n")
            write(f, "ComplexType         0\n")
            write(f, "=============================================\n")
            write(f, "=============================================\n")
            write(f, "0 1.0 0.0\n")
            write(f, "1 2.0 0.0\n")
            write(f, "2 3.0 0.0\n")
        end

        # Read input parameters
        read_input_parameters!(data, namelist_path)

        # Test: Orbital terms should be updated
        @test data.orbital_terms[1].value == ComplexF64(1.0, 0.0)
        @test data.orbital_terms[2].value == ComplexF64(2.0, 0.0)
        @test data.orbital_terms[3].value == ComplexF64(3.0, 0.0)

        # Cleanup
        rm(test_dir, recursive = true)
    end
end

"""
    test_read_input_parameters_complex()

Test read_input_parameters!() with complex values.
"""
function test_read_input_parameters_complex()
    @testset "read_input_parameters! Complex Tests" begin
        data = ExpertModeData()

        push!(data.gutzwiller_terms, GutzwillerTerm(0, 0.0+0.0im, true))

        test_dir = mktempdir()
        namelist_path = joinpath(test_dir, "namelist.def")

        open(namelist_path, "w") do f
            write(f, "ModPara modpara.def\n")
            write(f, "InGutzwiller ingutzwiller.def\n")
        end

        ingutzwiller_path = joinpath(test_dir, "ingutzwiller.def")
        open(ingutzwiller_path, "w") do f
            write(f, "=============================================\n")
            write(f, "NGutzwillerIdx          1\n")
            write(f, "ComplexType         1\n")
            write(f, "=============================================\n")
            write(f, "=============================================\n")
            write(f, "0 0.5 0.3\n")  # idx=0, real=0.5, imag=0.3
        end

        read_input_parameters!(data, namelist_path)

        @test data.gutzwiller_terms[1].value == ComplexF64(0.5, 0.3)

        rm(test_dir, recursive = true)
    end
end

"""
    test_read_input_parameters_missing_file()

Test read_input_parameters!() with missing file (should not error).
"""
function test_read_input_parameters_missing_file()
    @testset "read_input_parameters! Missing File Tests" begin
        data = ExpertModeData()

        push!(data.gutzwiller_terms, GutzwillerTerm(0, 0.0+0.0im, false))

        test_dir = mktempdir()
        namelist_path = joinpath(test_dir, "namelist.def")

        # Create namelist.def with InGutzwiller entry, but don't create the file
        open(namelist_path, "w") do f
            write(f, "ModPara modpara.def\n")
            write(f, "InGutzwiller ingutzwiller.def\n")
        end

        # Should not error even if file doesn't exist
        read_input_parameters!(data, namelist_path)

        # Value should remain unchanged
        @test data.gutzwiller_terms[1].value == ComplexF64(0.0, 0.0)

        rm(test_dir, recursive = true)
    end
end

# Run all tests
@testset "Read Input Parameters Tests" begin
    test_read_input_parameters_basic()
    test_read_input_parameters_orbital()
    test_read_input_parameters_complex()
    test_read_input_parameters_missing_file()
end
