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
    ModParaParameters,
    DoublonHolon2SiteIndex,
    DoublonHolon4SiteIndex
import MVMCExpertModeParsers: read_input_parameters!
import MVMCExpertModeParsers:
    _orbital_parallel_offset,
    count_orbital_parameters,
    parse_file_by_type!,
    _orbital_file_order_error,
    parse_expert_mode_files

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
            push!(data.orbital_terms, OrbitalTerm(0, i, i - 1, 0.0+0.0im, false))
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
            push!(data.orbital_terms, OrbitalTerm(0, i, i - 1, 0.0+0.0im, false))
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

function test_read_input_parameters_orbital_idx_overlay()
    @testset "read_input_parameters! Orbital overlays by parameter idx" begin
        mktempdir() do test_dir
            namelist_path = joinpath(test_dir, "namelist.def")
            write(namelist_path, "InOrbital inorbital.def\n")

            data = ExpertModeData()
            data.modpara.n_orbital_idx = 2
            data.orbital_terms = [
                OrbitalTerm(0, 0, 1, 0.0 + 0.0im, true, 1),
                OrbitalTerm(0, 1, 0, 0.0 + 0.0im, true, 1),
                OrbitalTerm(1, 0, 1, 0.0 + 0.0im, true, -1),
            ]

            _write_indexed_input_parameter_file(
                joinpath(test_dir, "inorbital.def"),
                "NOrbitalIdx",
                2,
                ["1 20.0 -0.2", "0 10.0 -0.1"],
            )

            read_input_parameters!(data, namelist_path)

            @test data.orbital_terms[1].value == 20.0 - 0.2im
            @test data.orbital_terms[2].value == 10.0 - 0.1im
            @test data.orbital_terms[3].value == 20.0 - 0.2im
        end
    end
end

function test_read_input_parameters_orbital_general_idx_overlay()
    @testset "read_input_parameters! OrbitalGeneral overlays by parameter idx" begin
        mktempdir() do test_dir
            namelist_path = joinpath(test_dir, "namelist.def")
            write(namelist_path, "InOrbitalGeneral inorbitalgeneral.def\n")

            data = ExpertModeData()
            data.modpara.n_orbital_idx = 2
            data.orbital_terms = [
                OrbitalTerm(0, 0, 1, 0.0 + 0.0im, true, 1),
                OrbitalTerm(0, 1, 0, 0.0 + 0.0im, true, 1),
                OrbitalTerm(1, 0, 1, 0.0 + 0.0im, true, -1),
            ]

            _write_indexed_input_parameter_file(
                joinpath(test_dir, "inorbitalgeneral.def"),
                "NOrbitalIdx",
                2,
                ["1 22.0 -0.4", "0 11.0 -0.3"],
            )

            read_input_parameters!(data, namelist_path)

            @test data.orbital_terms[1].value == 22.0 - 0.4im
            @test data.orbital_terms[2].value == 11.0 - 0.3im
            @test data.orbital_terms[3].value == 22.0 - 0.4im
        end
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

function test_read_input_parameters_orbital_parallel_and_opttrans()
    @testset "read_input_parameters! OrbitalParallel and OptTrans overlays" begin
        mktempdir() do test_dir
            namelist_path = joinpath(test_dir, "namelist.def")
            write(
                namelist_path,
                "InOrbitalParallel inorbitalparallel.def\nInOptTrans inopttrans.def\n",
            )

            data = ExpertModeData()
            data.i_flg_orbital_anti_parallel = 1
            data.i_flg_orbital_parallel = 1
            data.modpara.n_orbital_idx = 5
            # NArrayAP: one anti-parallel parameter (idx 0); parallel block starts at 1.
            data.n_orbital_anti_parallel = 1
            data.orbital_terms = [
                OrbitalTerm(0, 1, 0, 0.0 + 0.0im, false, 1),
                OrbitalTerm(0, 0, 1, 0.0 + 0.0im, false, 1),
                OrbitalTerm(0, 0, 2, 0.0 + 0.0im, false, 1),
                OrbitalTerm(1, 1, 3, 0.0 + 0.0im, false, 1),
                OrbitalTerm(1, 1, 4, 0.0 + 0.0im, false, 1),
            ]
            data.opt_trans = [1.0 + 0.0im, 2.0 + 0.0im]

            _write_indexed_input_parameter_file(
                joinpath(test_dir, "inorbitalparallel.def"),
                "NOrbitalParallel",
                4,
                [
                    "0 10.0 0.1",
                    "1 20.0 0.2",
                    "2 30.0 0.3",
                    "3 40.0 0.4",
                ],
            )
            _write_indexed_input_parameter_file(
                joinpath(test_dir, "inopttrans.def"),
                "NQPOptTrans",
                2,
                [
                    "1 0.25 -0.25",
                    "0 0.50 -0.50",
                ],
            )

            read_input_parameters!(data, namelist_path)

            @test data.orbital_terms[1].value == 0.0 + 0.0im
            @test [t.value for t in data.orbital_terms[2:end]] == ComplexF64[
                10.0 + 0.1im,
                20.0 + 0.2im,
                30.0 + 0.3im,
                40.0 + 0.4im,
            ]
            @test data.opt_trans == ComplexF64[
                0.50 - 0.50im,
                0.25 - 0.25im,
            ]
        end
    end
end

function test_parse_opttrans_def()
    @testset "OptTrans definition parser" begin
        mktempdir() do test_dir
            data = ExpertModeData()
            data.modpara = ModParaParameters(nsite = 2, nmp_trans = 1)
            path = joinpath(test_dir, "opttrans.def")
            open(path, "w") do f
                write(f, "=============================================\n")
                write(f, "NQPOptTrans          2\n")
                write(f, "=============================================\n")
                write(f, "=============================================\n")
                write(f, "=============================================\n")
                write(f, "0 0.25\n")
                write(f, "1 0.75\n")
                write(f, "0 0 1 -1\n")
                write(f, "0 1 0 -1\n")
                write(f, "1 0 0 1\n")
                write(f, "1 1 1 -1\n")
            end

            MVMCExpertModeParsers.parse_file_by_type!(data, "OptTrans", path)

            @test data.n_qp_opt_trans == 2
            @test data.para_qp_opt_trans == ComplexF64[0.25 + 0.0im, 0.75 + 0.0im]
            @test data.opt_trans == data.para_qp_opt_trans
            @test data.qp_opt_trans == [[1, 0], [0, 1]]
            # APFlag is off when NMPTrans >= 0, so C forces all signs to +1.
            @test data.qp_opt_trans_sgn == [[1, 1], [1, 1]]
        end
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

function test_orbital_parallel_offset_cases()
    @testset "_orbital_parallel_offset branches" begin
        # No OrbitalParallel block: the overlay window spans all orbital params,
        # so the offset equals the total orbital parameter count.
        data = ExpertModeData()
        data.i_flg_orbital_anti_parallel = 1
        data.i_flg_orbital_parallel = 0
        data.modpara.n_orbital_idx = 4
        data.n_orbital_anti_parallel = 4
        @test _orbital_parallel_offset(data) == count_orbital_parameters(data) == 4

        # Parallel-only (no anti-parallel): the parallel block starts at index 0.
        data = ExpertModeData()
        data.i_flg_orbital_anti_parallel = 0
        data.i_flg_orbital_parallel = 1
        data.modpara.n_orbital_idx = 6
        data.n_orbital_anti_parallel = 0
        @test _orbital_parallel_offset(data) == 0

        # Both present: the offset is exactly NArrayAP, independent of how the
        # orbital_terms indices happen to be grouped (no heuristic inference).
        for narrayap in (1, 2, 3, 5)
            data = ExpertModeData()
            data.i_flg_orbital_anti_parallel = 1
            data.i_flg_orbital_parallel = 1
            data.modpara.n_orbital_idx = narrayap + 4
            data.n_orbital_anti_parallel = narrayap
            @test _orbital_parallel_offset(data) == narrayap
        end
    end
end

function test_read_input_parameters_orbital_parallel_offset_regression()
    @testset "InOrbitalParallel exact offset (heuristic-fragile config)" begin
        mktempdir() do test_dir
            namelist_path = joinpath(test_dir, "namelist.def")
            write(namelist_path, "InOrbitalParallel inorbitalparallel.def\n")

            # Anti-parallel block (NArrayAP = 3) deliberately contains two terms
            # that share (site1, site2, sign) with consecutive indices (idx 0, 1).
            # The previous heuristic inferred the parallel-block start from the
            # smallest such consecutive pair and would have returned offset 0,
            # corrupting the anti-parallel parameters. The exact NArrayAP-based
            # offset (3) confines the overlay to the parallel block (idx 3, 4).
            data = ExpertModeData()
            data.i_flg_orbital_anti_parallel = 1
            data.i_flg_orbital_parallel = 1
            data.modpara.n_orbital_idx = 5
            data.n_orbital_anti_parallel = 3
            data.orbital_terms = [
                OrbitalTerm(0, 1, 0, 7.0 + 0.0im, false, 1),  # anti
                OrbitalTerm(0, 1, 1, 8.0 + 0.0im, false, 1),  # anti (same key as idx 0)
                OrbitalTerm(0, 2, 2, 9.0 + 0.0im, false, 1),  # anti
                OrbitalTerm(1, 0, 3, 0.0 + 0.0im, false, 1),  # parallel up
                OrbitalTerm(1, 0, 4, 0.0 + 0.0im, false, 1),  # parallel down
            ]

            @test _orbital_parallel_offset(data) == 3

            _write_indexed_input_parameter_file(
                joinpath(test_dir, "inorbitalparallel.def"),
                "NOrbitalParallel",
                2,
                ["0 10.0 0.1", "1 20.0 0.2"],
            )

            read_input_parameters!(data, namelist_path)

            # Anti-parallel parameters are untouched.
            @test [t.value for t in data.orbital_terms[1:3]] ==
                  ComplexF64[7.0, 8.0, 9.0]
            # Parallel parameters receive the overlay values.
            @test data.orbital_terms[4].value == 10.0 + 0.1im
            @test data.orbital_terms[5].value == 20.0 + 0.2im
        end
    end
end

function test_orbital_parallel_offset_from_parse()
    @testset "parse records NArrayAP for OrbitalParallel offset" begin
        mktempdir() do test_dir
            anti_path = joinpath(test_dir, "orbitalidx.def")
            write(
                anti_path,
                """
                =============================================
                NOrbitalIdx          2
                ComplexType          0
                =============================================
                =============================================
                    0      0      0
                    0      1      1
                    1      0      1
                    1      1      0
                    0      1
                    1      1
                """,
            )

            data = ExpertModeData()
            parse_file_by_type!(data, "OrbitalAntiParallel", anti_path)
            # OrbitalAntiParallel alone records NArrayAP = max idx + 1 = 2.
            @test data.n_orbital_anti_parallel == 2
            @test data.i_flg_orbital_anti_parallel == 1

            par_path = joinpath(test_dir, "orbitalidxpara.def")
            write(
                par_path,
                """
                =============================================
                NOrbitalParallel     1
                ComplexType          0
                =============================================
                =============================================
                    0      1      0
                    1      0      0
                    0      1
                """,
            )
            parse_file_by_type!(data, "OrbitalParallel", par_path)

            # Parsing OrbitalParallel preserves the exact NArrayAP recorded from
            # the anti-parallel block, which is the offset used by the overlay.
            @test data.i_flg_orbital_parallel == 1
            @test data.n_orbital_anti_parallel == 2
            @test _orbital_parallel_offset(data) == 2
        end
    end
end

function test_orbital_parallel_offset_sparse_header()
    @testset "NArrayAP from header count, not max(idx)+1" begin
        mktempdir() do test_dir
            # Anti-parallel file declares NOrbitalIdx=3 but only idx 0,1 are
            # referenced by site pairs (idx 2 is a declared-but-unreferenced
            # parameter, listed only in the OptFlag section). C uses the header
            # count (iNOrbitalAntiParallel=3); max(idx)+1 would give 2.
            anti_path = joinpath(test_dir, "orbitalidx.def")
            write(
                anti_path,
                """
                =============================================
                NOrbitalIdx          3
                ComplexType          0
                =============================================
                =============================================
                    0      1      0
                    1      0      1
                    0      1
                    1      1
                    2      1
                """,
            )

            data = ExpertModeData()
            parse_file_by_type!(data, "OrbitalAntiParallel", anti_path)
            @test data.n_orbital_anti_parallel == 3
            @test data.modpara.n_orbital_idx == 3

            # Parallel file declares 1 parallel parameter.
            par_path = joinpath(test_dir, "orbitalidxpara.def")
            write(
                par_path,
                """
                =============================================
                NOrbitalIdx          1
                ComplexType          0
                =============================================
                =============================================
                    0      1      0
                    0      1
                """,
            )
            parse_file_by_type!(data, "OrbitalParallel", par_path)

            # NArrayAP stays 3 (header), so parallel slots are 3,4 and the
            # reserved-but-unused anti slot 2 is not stolen by the parallel block.
            @test data.n_orbital_anti_parallel == 3
            @test _orbital_parallel_offset(data) == 3
            @test data.modpara.n_orbital_idx == 3 + 2 * 1
            @test any(t -> t.idx == 3, data.orbital_terms)
            @test any(t -> t.idx == 4, data.orbital_terms)
            @test all(t -> t.idx != 2, data.orbital_terms)

            # InOrbitalParallel overlay must land on slots 3,4 only.
            namelist_path = joinpath(test_dir, "namelist.def")
            write(namelist_path, "InOrbitalParallel inorbitalparallel.def\n")
            _write_indexed_input_parameter_file(
                joinpath(test_dir, "inorbitalparallel.def"),
                "NOrbitalParallel",
                2,
                ["0 11.0 0.0", "1 22.0 0.0"],
            )
            read_input_parameters!(data, namelist_path)
            @test all(
                t.value == 11.0 + 0.0im for t in data.orbital_terms if t.idx == 3
            )
            @test all(
                t.value == 22.0 + 0.0im for t in data.orbital_terms if t.idx == 4
            )
        end
    end
end

function test_orbital_file_order_error()
    @testset "OrbitalParallel-before-AntiParallel is rejected" begin
        # Helper: parallel before anti is an error; anti-first or pure-parallel ok.
        @test _orbital_file_order_error([
            ("OrbitalParallel", "p.def"),
            ("OrbitalAntiParallel", "a.def"),
        ]) !== nothing
        @test _orbital_file_order_error([
            ("OrbitalParallel", "p.def"),
            ("Orbital", "a.def"),
        ]) !== nothing
        @test _orbital_file_order_error([
            ("OrbitalAntiParallel", "a.def"),
            ("OrbitalParallel", "p.def"),
        ]) === nothing
        @test _orbital_file_order_error([("OrbitalParallel", "p.def")]) === nothing
        @test _orbital_file_order_error([("Gutzwiller", "g.def")]) === nothing

        # End-to-end: parse_expert_mode_files hard-fails on the bad order, before
        # any referenced file is even opened.
        mktempdir() do test_dir
            namelist_path = joinpath(test_dir, "namelist.def")
            write(
                namelist_path,
                "OrbitalParallel orbitalidxpara.def\nOrbitalAntiParallel orbitalidx.def\n",
            )
            err = nothing
            try
                parse_expert_mode_files(namelist_path)
            catch e
                err = e
            end
            @test err isa ErrorException
            @test occursin("OrbitalParallel must be listed after", sprint(showerror, err))
        end
    end
end

# Run all tests
@testset "Read Input Parameters Tests" begin
    test_read_input_parameters_basic()
    test_read_input_parameters_orbital()
    test_read_input_parameters_orbital_idx_overlay()
    test_read_input_parameters_orbital_general_idx_overlay()
    test_read_input_parameters_complex()
    test_read_input_parameters_missing_file()
    test_read_input_parameters_orbital_parallel_and_opttrans()
    test_orbital_parallel_offset_cases()
    test_read_input_parameters_orbital_parallel_offset_regression()
    test_orbital_parallel_offset_from_parse()
    test_orbital_parallel_offset_sparse_header()
    test_orbital_file_order_error()
    test_parse_opttrans_def()
    test_read_input_parameters_dh_overlays()
end
