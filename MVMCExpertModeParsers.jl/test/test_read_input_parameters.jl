"""
Read Input Parameters Tests

Tests for read_input_parameters! function (equivalent to C's ReadInputParameters()).
"""

using Test
using MVMCExpertModeParsers

import MVMCExpertModeParsers:
    ExpertModeData,
    GutzwillerTerm,
    JastrowTerm,
    OrbitalTerm,
    DoublonHolon2SiteIndex,
    DoublonHolon4SiteIndex
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

function _write_indexed_input_parameter_file(path, header_name, header_count, rows)
    open(path, "w") do f
        write(f, "=============================================\n")
        write(f, "$header_name          $header_count\n")
        write(f, "ComplexType         1\n")
        write(f, "=============================================\n")
        write(f, "=============================================\n")
        for row in rows
            write(f, row * "\n")
        end
    end
end

function test_read_input_parameters_dh_overlays()
    @testset "read_input_parameters! DH overlays" begin
        mktempdir() do test_dir
            namelist_path = joinpath(test_dir, "namelist.def")
            open(namelist_path, "w") do f
                write(f, "ModPara modpara.def\n")
                write(f, "InDH2 indh2.def\n")
                write(f, "InDH4 indh4.def\n")
            end

            data = ExpertModeData()
            data.doublon_holon_2site_indices = [DoublonHolon2SiteIndex([1 0; 0 1])]
            data.doublon_holon_2site_params = fill(0.0 + 0.0im, 6)
            data.doublon_holon_4site_indices = [DoublonHolon4SiteIndex([1 0 1 0; 0 1 0 1])]
            data.doublon_holon_4site_params = fill(0.0 + 0.0im, 10)

            _write_indexed_input_parameter_file(
                joinpath(test_dir, "indh2.def"),
                "NDoublonHolon2siteIdx",
                1,
                [
                    "5 5.0 -5.0",
                    "3 3.0 -3.0",
                    "1 1.0 -1.0",
                    "0 0.0 0.5",
                    "2 2.0 -2.0",
                    "4 4.0 -4.0",
                ],
            )
            _write_indexed_input_parameter_file(
                joinpath(test_dir, "indh4.def"),
                "NDoublonHolon4siteIdx",
                1,
                ["$(i - 1) $(10 + i).0 -$(10 + i).0" for i = 1:10],
            )

            read_input_parameters!(data, namelist_path)

            @test data.doublon_holon_2site_params == ComplexF64[
                0.0 + 0.5im,
                1.0 - 1.0im,
                2.0 - 2.0im,
                3.0 - 3.0im,
                4.0 - 4.0im,
                5.0 - 5.0im,
            ]
            @test data.doublon_holon_4site_params == ComplexF64[
                11.0 - 11.0im,
                12.0 - 12.0im,
                13.0 - 13.0im,
                14.0 - 14.0im,
                15.0 - 15.0im,
                16.0 - 16.0im,
                17.0 - 17.0im,
                18.0 - 18.0im,
                19.0 - 19.0im,
                20.0 - 20.0im,
            ]
        end

        mktempdir() do test_dir
            namelist_path = joinpath(test_dir, "namelist.def")
            write(namelist_path, "InDH2 indh2.def\n")

            data = ExpertModeData()
            data.doublon_holon_2site_indices = [DoublonHolon2SiteIndex([1 0; 0 1])]
            data.doublon_holon_2site_params = fill(0.0 + 0.0im, 6)

            bad_cases = Dict(
                "duplicate" => [
                    "0 0.0 0.0",
                    "1 1.0 0.0",
                    "1 2.0 0.0",
                    "3 3.0 0.0",
                    "4 4.0 0.0",
                    "5 5.0 0.0",
                ],
                "range" => [
                    "0 0.0 0.0",
                    "1 1.0 0.0",
                    "2 2.0 0.0",
                    "3 3.0 0.0",
                    "4 4.0 0.0",
                    "6 6.0 0.0",
                ],
                "short" => [
                    "0 0.0 0.0",
                    "1 1.0 0.0",
                    "2 2.0 0.0",
                    "3 3.0 0.0",
                    "4 4.0 0.0",
                ],
            )

            for rows in values(bad_cases)
                _write_indexed_input_parameter_file(
                    joinpath(test_dir, "indh2.def"),
                    "NDoublonHolon2siteIdx",
                    1,
                    rows,
                )
                @test_throws ErrorException read_input_parameters!(data, namelist_path)
                @test data.doublon_holon_2site_params == fill(0.0 + 0.0im, 6)
            end
        end

        mktempdir() do test_dir
            dh_specs = (
                (
                    file_type = "InDH2",
                    file_name = "indh2.def",
                    header_name = "NDoublonHolon2siteIdx",
                    n_index = 1,
                    n_param = 6,
                    setup! = data -> begin
                        data.doublon_holon_2site_indices = [DoublonHolon2SiteIndex([1 0; 0 1])]
                        data.doublon_holon_2site_params = fill(0.0 + 0.0im, 6)
                    end,
                    params = data -> data.doublon_holon_2site_params,
                ),
                (
                    file_type = "InDH4",
                    file_name = "indh4.def",
                    header_name = "NDoublonHolon4siteIdx",
                    n_index = 1,
                    n_param = 10,
                    setup! = data -> begin
                        data.doublon_holon_4site_indices = [DoublonHolon4SiteIndex([1 0 1 0; 0 1 0 1])]
                        data.doublon_holon_4site_params = fill(0.0 + 0.0im, 10)
                    end,
                    params = data -> data.doublon_holon_4site_params,
                ),
            )

            for spec in dh_specs
                namelist_path = joinpath(test_dir, "namelist.def")
                write(namelist_path, "$(spec.file_type) $(spec.file_name)\n")
                def_path = joinpath(test_dir, spec.file_name)
                good_rows = ["$(i - 1) $(i).0 -$(i).0" for i = 1:spec.n_param]
                bad_cases = (
                    (header_count = 2, rows = good_rows),
                    (header_count = spec.n_index, rows = vcat(["0 0.0"], good_rows[2:end])),
                    (header_count = spec.n_index, rows = vcat(["0 NaN 0.0"], good_rows[2:end])),
                    (header_count = spec.n_index, rows = vcat(["bad-index 0.0 0.0"], good_rows[2:end])),
                )

                for case in bad_cases
                    data = ExpertModeData()
                    spec.setup!(data)
                    _write_indexed_input_parameter_file(
                        def_path,
                        spec.header_name,
                        case.header_count,
                        case.rows,
                    )
                    @test_throws ErrorException read_input_parameters!(data, namelist_path)
                    @test spec.params(data) == fill(0.0 + 0.0im, spec.n_param)
                end
            end
        end
    end
end

# Run all tests
@testset "Read Input Parameters Tests" begin
    test_read_input_parameters_basic()
    test_read_input_parameters_orbital()
    test_read_input_parameters_complex()
    test_read_input_parameters_missing_file()
    test_read_input_parameters_dh_overlays()
end
