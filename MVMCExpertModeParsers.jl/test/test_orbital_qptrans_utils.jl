"""
Orbital and QPTrans Utilities Tests

Tests for build_orbital_sgn_matrix!() and build_qp_trans_mappings!() functions.
"""

using Test
using MVMCExpertModeParsers

import MVMCExpertModeParsers: ExpertModeData, OrbitalTerm, ModParaParameters
import MVMCExpertModeParsers: build_orbital_sgn_matrix!, build_qp_trans_mappings!

# Test data directory
const TEST_DATA_DIR = joinpath(@__DIR__, "samples", "HeisenbergChain")

"""
    test_build_orbital_sgn_matrix_sz_conserved()

Test build_orbital_sgn_matrix!() for sz-conserved case (iFlgOrbitalGeneral == 0).
"""
function test_build_orbital_sgn_matrix_sz_conserved()
    @testset "build_orbital_sgn_matrix! sz-conserved Tests" begin
        # Create test data
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
        data.i_flg_orbital_general = 0  # sz-conserved
        data.modpara.nmp_trans = 4  # APFlag = 0 (periodic)

        # Add orbital terms with signs
        # OrbitalTerm(site1, site2, idx, value, is_complex, sign)
        push!(data.orbital_terms, OrbitalTerm(0, 0, 0, 0.5+0.0im, false, 1))  # sign = 1
        push!(data.orbital_terms, OrbitalTerm(0, 1, 1, 0.3+0.0im, false, -1))  # sign = -1
        push!(data.orbital_terms, OrbitalTerm(1, 2, 2, 0.2+0.0im, false, 1))  # sign = 1

        # Build orbital sign matrix
        build_orbital_sgn_matrix!(data)

        # Test: orbital_sgn should be set
        @test data.orbital_sgn !== nothing
        @test size(data.orbital_sgn) == (4, 4)  # n_site x n_site

        # Test: APFlag == 0, so all signs should be +1 (C implementation behavior)
        @test all(data.orbital_sgn .== 1)

        # Test with APFlag == 1 (anti-periodic)
        data2 = ExpertModeData()
        data2.modpara = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
        data2.i_flg_orbital_general = 0
        data2.modpara.nmp_trans = -4  # APFlag = 1 (anti-periodic)

        push!(data2.orbital_terms, OrbitalTerm(0, 0, 0, 0.5+0.0im, false, 1))
        push!(data2.orbital_terms, OrbitalTerm(0, 1, 1, 0.3+0.0im, false, -1))
        push!(data2.orbital_terms, OrbitalTerm(1, 2, 2, 0.2+0.0im, false, 1))

        build_orbital_sgn_matrix!(data2)

        # Test: APFlag == 1, so signs from orbital_terms should be preserved
        @test data2.orbital_sgn !== nothing
        @test data2.orbital_sgn[1, 1] == 1  # site 0, site 0
        @test data2.orbital_sgn[1, 2] == -1  # site 0, site 1
        @test data2.orbital_sgn[2, 3] == 1  # site 1, site 2
    end
end

"""
    test_build_orbital_sgn_matrix_fsz()

Test build_orbital_sgn_matrix!() for fsz case (iFlgOrbitalGeneral == 1).
"""
function test_build_orbital_sgn_matrix_fsz()
    @testset "build_orbital_sgn_matrix! fsz Tests" begin
        # Create test data for fsz case
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite = 2, nelec = 2, nlocspin = 0)
        data.i_flg_orbital_general = 1  # fsz
        data.modpara.nmp_trans = 2  # APFlag = 0

        # Add orbital terms (for fsz, site indices can be 0..2*Nsite-1)
        push!(data.orbital_terms, OrbitalTerm(0, 1, 0, 0.5+0.0im, false, 1))
        push!(data.orbital_terms, OrbitalTerm(2, 3, 1, 0.3+0.0im, false, -1))

        build_orbital_sgn_matrix!(data)

        # Test: orbital_sgn should be 2*Nsite x 2*Nsite
        @test data.orbital_sgn !== nothing
        @test size(data.orbital_sgn) == (4, 4)  # 2*2 x 2*2

        # Test: APFlag == 0, so signs should follow C implementation
        # For fsz with APFlag == 0: F_{IJ} = 1, F_{JI} = -1 for I < J
        @test data.orbital_sgn[1, 2] == 1  # I=0, J=1
        @test data.orbital_sgn[2, 1] == -1  # J=1, I=0
    end
end

"""
    test_build_orbital_sgn_matrix_with_real_data()

Test build_orbital_sgn_matrix!() with real HeisenbergChain data.
"""
function test_build_orbital_sgn_matrix_with_real_data()
    @testset "build_orbital_sgn_matrix! with real data Tests" begin
        # Skip if test data directory doesn't exist
        if !isdir(TEST_DATA_DIR)
            @test_skip "Test data directory not found: $TEST_DATA_DIR"
            return
        end

        orbitalidx_file = joinpath(TEST_DATA_DIR, "orbitalidx.def")
        if !isfile(orbitalidx_file)
            @test_skip "orbitalidx.def not found in test directory"
            return
        end

        # Parse orbitalidx.def
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite = 16, nelec = 8, nlocspin = 16)
        data.i_flg_orbital_general = 0  # sz-conserved
        data.modpara.nmp_trans = 4  # APFlag = 0

        # Parse orbitalidx.def (same format as orbital.def)
        # parse_orbital_content returns (ParseResult, Dict{Int,Int})
        parse_result, _opt_flags = MVMCExpertModeParsers.parse_orbital_content(
            read(orbitalidx_file, String),
        )

        if parse_result.success
            data.orbital_terms = parse_result.data

            # Build orbital sign matrix
            build_orbital_sgn_matrix!(data)

            # Test: orbital_sgn should be set
            @test data.orbital_sgn !== nothing
            @test size(data.orbital_sgn) == (16, 16)  # n_site x n_site

            # Test: APFlag == 0, so all signs should be +1
            @test all(data.orbital_sgn .== 1)
        else
            @warn "Failed to parse orbitalidx.def: $(parse_result.error_message)"
        end
    end
end

"""
    test_build_qp_trans_mappings_basic()

Test build_qp_trans_mappings!() with basic data.
"""
function test_build_qp_trans_mappings_basic()
    @testset "build_qp_trans_mappings! basic Tests" begin
        # Create test data
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
        data.n_qp_trans = 2
        data.n_qp_opt_trans = 1
        data.modpara.nmp_trans = 2  # APFlag = 0

        # Create temporary qptransidx.def file
        test_dir = mktempdir()
        qptransidx_file = joinpath(test_dir, "qptransidx.def")

        open(qptransidx_file, "w") do f
            write(f, "=============================================\n")
            write(f, "NQPTrans          2\n")
            write(f, "=============================================\n")
            write(f, "======== TrIdx_TrWeight_and_TrIdx_i_xi ======\n")
            write(f, "=============================================\n")
            write(f, "0    1.00000\n")
            write(f, "1    1.00000\n")
            # QPTrans mappings: mpidx original_site translated_site sign
            write(f, "    0      0      0      1\n")  # mpidx=0, site 0 -> 0, sign=1
            write(f, "    0      1      1      1\n")  # mpidx=0, site 1 -> 1, sign=1
            write(f, "    0      2      2      1\n")  # mpidx=0, site 2 -> 2, sign=1
            write(f, "    0      3      3      1\n")  # mpidx=0, site 3 -> 3, sign=1
            write(f, "    1      0      1      1\n")  # mpidx=1, site 0 -> 1, sign=1
            write(f, "    1      1      2      1\n")  # mpidx=1, site 1 -> 2, sign=1
            write(f, "    1      2      3      1\n")  # mpidx=1, site 2 -> 3, sign=1
            write(f, "    1      3      0      1\n")  # mpidx=1, site 3 -> 0, sign=1
        end

        # Build QPTrans mappings
        build_qp_trans_mappings!(data, qptransidx_file)

        # Test: QPTrans should be set
        @test !isempty(data.qp_trans)
        @test length(data.qp_trans) == 2  # NQPTrans = 2
        @test length(data.qp_trans[1]) == 4  # Nsite = 4
        @test length(data.qp_trans[2]) == 4

        # Test: QPTrans mappings
        # mpidx=0: identity mapping
        @test data.qp_trans[1][1] == 0  # site 0 -> 0
        @test data.qp_trans[1][2] == 1  # site 1 -> 1
        @test data.qp_trans[1][3] == 2  # site 2 -> 2
        @test data.qp_trans[1][4] == 3  # site 3 -> 3

        # mpidx=1: shift by 1
        @test data.qp_trans[2][1] == 1  # site 0 -> 1
        @test data.qp_trans[2][2] == 2  # site 1 -> 2
        @test data.qp_trans[2][3] == 3  # site 2 -> 3
        @test data.qp_trans[2][4] == 0  # site 3 -> 0

        # Test: QPTransInv (inverse mapping)
        @test !isempty(data.qp_trans_inv)
        @test length(data.qp_trans_inv) == 2
        # For mpidx=1: inverse of site 0 -> 1 is site 1 -> 0
        @test data.qp_trans_inv[2][2] == 0  # translated site 1 -> original site 0

        # Test: QPTransSgn
        @test !isempty(data.qp_trans_sgn)
        @test length(data.qp_trans_sgn) == 2
        # APFlag == 0, so all signs should be +1
        @test all(data.qp_trans_sgn[1] .== 1)
        @test all(data.qp_trans_sgn[2] .== 1)

        # Test: QPOptTrans (identity mapping by default)
        @test !isempty(data.qp_opt_trans)
        @test length(data.qp_opt_trans) == 1  # NQPOptTrans = 1
        @test data.qp_opt_trans[1] == [0, 1, 2, 3]  # Identity mapping

        # Cleanup
        rm(test_dir, recursive = true)
    end
end

"""
    test_build_qp_trans_mappings_with_real_data()

Test build_qp_trans_mappings!() with real HeisenbergChain data.
"""
function test_build_qp_trans_mappings_with_real_data()
    @testset "build_qp_trans_mappings! with real data Tests" begin
        # Skip if test data directory doesn't exist
        if !isdir(TEST_DATA_DIR)
            @test_skip "Test data directory not found: $TEST_DATA_DIR"
            return
        end

        qptransidx_file = joinpath(TEST_DATA_DIR, "qptransidx.def")
        if !isfile(qptransidx_file)
            @test_skip "qptransidx.def not found in test directory"
            return
        end

        # Create test data
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite = 16, nelec = 8, nlocspin = 16)
        data.n_qp_trans = 4  # From qptransidx.def
        data.n_qp_opt_trans = 1
        data.modpara.nmp_trans = 4  # APFlag = 0

        # Build QPTrans mappings
        build_qp_trans_mappings!(data, qptransidx_file)

        # Test: QPTrans should be set
        @test !isempty(data.qp_trans)
        @test length(data.qp_trans) == 4  # NQPTrans = 4
        @test length(data.qp_trans[1]) == 16  # Nsite = 16

        # Test: QPTransInv should be set
        @test !isempty(data.qp_trans_inv)
        @test length(data.qp_trans_inv) == 4

        # Test: QPTransSgn should be set
        @test !isempty(data.qp_trans_sgn)
        @test length(data.qp_trans_sgn) == 4

        # Test: APFlag == 0, so all signs should be +1
        for mpidx = 1:4
            @test all(data.qp_trans_sgn[mpidx] .== 1)
        end

        # Test: QPTrans and QPTransInv are inverse mappings
        # For each mpidx, QPTransInv[QPTrans[site]] should equal site
        for mpidx = 1:4
            for site = 1:16
                translated = data.qp_trans[mpidx][site]
                original = data.qp_trans_inv[mpidx][translated+1]  # +1 for 1-based indexing
                @test original == (site - 1)  # Convert to 0-based for comparison
            end
        end

        # Test: QPOptTrans (identity mapping by default)
        @test !isempty(data.qp_opt_trans)
        @test length(data.qp_opt_trans) == 1
        @test data.qp_opt_trans[1] == collect(0:15)  # Identity mapping: [0, 1, 2, ..., 15]
    end
end

"""
    test_build_qp_trans_mappings_apflag()

Test build_qp_trans_mappings!() with APFlag == 1 (anti-periodic).
"""
function test_build_qp_trans_mappings_apflag()
    @testset "build_qp_trans_mappings! APFlag Tests" begin
        # Create test data with APFlag == 1
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
        data.n_qp_trans = 1
        data.n_qp_opt_trans = 1
        data.modpara.nmp_trans = -4  # APFlag = 1 (anti-periodic)

        # Create temporary qptransidx.def file with negative signs
        test_dir = mktempdir()
        qptransidx_file = joinpath(test_dir, "qptransidx.def")

        open(qptransidx_file, "w") do f
            write(f, "=============================================\n")
            write(f, "NQPTrans          1\n")
            write(f, "=============================================\n")
            write(f, "======== TrIdx_TrWeight_and_TrIdx_i_xi ======\n")
            write(f, "=============================================\n")
            write(f, "0    1.00000\n")
            # QPTrans mappings with negative signs
            write(f, "    0      0      1     -1\n")  # mpidx=0, site 0 -> 1, sign=-1
            write(f, "    0      1      2      1\n")   # mpidx=0, site 1 -> 2, sign=1
            write(f, "    0      2      3     -1\n")   # mpidx=0, site 2 -> 3, sign=-1
            write(f, "    0      3      0      1\n")   # mpidx=0, site 3 -> 0, sign=1
        end

        # Build QPTrans mappings
        build_qp_trans_mappings!(data, qptransidx_file)

        # Test: APFlag == 1, so signs from file should be preserved
        @test data.qp_trans_sgn[1][1] == -1  # site 0 -> sign -1
        @test data.qp_trans_sgn[1][2] == 1   # site 1 -> sign 1
        @test data.qp_trans_sgn[1][3] == -1  # site 2 -> sign -1
        @test data.qp_trans_sgn[1][4] == 1   # site 3 -> sign 1

        # Cleanup
        rm(test_dir, recursive = true)
    end
end

# Run all tests
@testset "Orbital and QPTrans Utilities Tests" begin
    test_build_orbital_sgn_matrix_sz_conserved()
    test_build_orbital_sgn_matrix_fsz()
    test_build_orbital_sgn_matrix_with_real_data()
    test_build_qp_trans_mappings_basic()
    test_build_qp_trans_mappings_with_real_data()
    test_build_qp_trans_mappings_apflag()
end
