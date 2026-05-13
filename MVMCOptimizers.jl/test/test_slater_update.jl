"""
Slater Update Tests

Tests for build_qp_trans_matrices() and update_slater_elm_fcmp!() functions.
"""

using Test
using Random
using Logging
using MVMCExpertModeParsers: initialize_parameters!, read_input_parameters!
using SFMT

using MVMCExpertModeParsers: ExpertModeData, OrbitalTerm, ModParaParameters
using MVMCExpertModeParsers: build_orbital_sgn_matrix!, build_qp_trans_mappings!
using MVMCExpertModeParsers: init_qp_weight!

using MVMCOptimizers: build_qp_trans_matrices, update_slater_elm_fcmp!
using MVMCOptimizers: VMCOptimizationState

# Test data directory
const TEST_DATA_DIR = joinpath(
    @__DIR__,
    "..",
    "..",
    "MVMCExpertModeParsers.jl",
    "test",
    "samples",
    "HeisenbergChain",
)

"""
    test_build_qp_trans_matrices_success()

Test build_qp_trans_matrices() when QPTrans mappings exist.
"""
function test_build_qp_trans_matrices_success()
    @testset "build_qp_trans_matrices success Tests" begin
        # Create test data with QPTrans mappings
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
        data.n_qp_trans = 2
        data.n_qp_opt_trans = 1

        # Set up QPTrans mappings
        data.qp_trans = [
            [0, 1, 2, 3],  # mpidx=0: identity
            [1, 2, 3, 0],   # mpidx=1: shift by 1
        ]
        data.qp_trans_inv = [
            [0, 1, 2, 3],  # mpidx=0: identity
            [3, 0, 1, 2],   # mpidx=1: inverse of shift
        ]
        data.qp_trans_sgn = [
            [1, 1, 1, 1],  # mpidx=0: all +1
            [1, 1, 1, 1],   # mpidx=1: all +1
        ]

        # Build QPTrans matrices
        qp_trans, qp_trans_inv, qp_trans_sgn, qp_opt_trans, qp_opt_trans_sgn =
            build_qp_trans_matrices(data)

        # Test: QPTrans should match input
        @test length(qp_trans) == 2
        @test qp_trans[1] == [0, 1, 2, 3]
        @test qp_trans[2] == [1, 2, 3, 0]

        # Test: QPTransInv should match input
        @test length(qp_trans_inv) == 2
        @test qp_trans_inv[1] == [0, 1, 2, 3]
        @test qp_trans_inv[2] == [3, 0, 1, 2]

        # Test: QPTransSgn should match input
        @test length(qp_trans_sgn) == 2
        @test qp_trans_sgn[1] == [1, 1, 1, 1]
        @test qp_trans_sgn[2] == [1, 1, 1, 1]

        # Test: QPOptTrans should be identity mapping (default)
        @test length(qp_opt_trans) == 1
        @test qp_opt_trans[1] == [0, 1, 2, 3]  # Identity mapping

        # Test: QPOptTransSgn should be all +1
        @test length(qp_opt_trans_sgn) == 1
        @test qp_opt_trans_sgn[1] == [1, 1, 1, 1]
    end
end

"""
    test_build_qp_trans_matrices_error()

Test build_qp_trans_matrices() when QPTrans mappings are missing (should throw error).
"""
function test_build_qp_trans_matrices_error()
    @testset "build_qp_trans_matrices error Tests" begin
        # Create test data without QPTrans mappings
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
        data.n_qp_trans = 2
        data.n_qp_opt_trans = 1

        # QPTrans mappings are empty (not set)
        data.qp_trans = Vector{Vector{Int}}()
        data.qp_trans_inv = Vector{Vector{Int}}()
        data.qp_trans_sgn = Vector{Vector{Int}}()

        # Test: Should throw ArgumentError
        # Suppress expected @error log during this test
        with_logger(NullLogger()) do
            @test_throws ArgumentError build_qp_trans_matrices(data)
        end
    end
end

"""
    test_build_qp_trans_matrices_qp_opt_trans()

Test build_qp_trans_matrices() with QPOptTrans mappings.
"""
function test_build_qp_trans_matrices_qp_opt_trans()
    @testset "build_qp_trans_matrices QPOptTrans Tests" begin
        # Create test data with both QPTrans and QPOptTrans mappings
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
        data.n_qp_trans = 1
        data.n_qp_opt_trans = 2

        # Set up QPTrans mappings
        data.qp_trans = [[0, 1, 2, 3]]
        data.qp_trans_inv = [[0, 1, 2, 3]]
        data.qp_trans_sgn = [[1, 1, 1, 1]]

        # Set up QPOptTrans mappings
        data.qp_opt_trans = [
            [0, 1, 2, 3],  # optidx=0: identity
            [1, 2, 3, 0],   # optidx=1: shift by 1
        ]
        data.qp_opt_trans_sgn = [
            [1, 1, 1, 1],  # optidx=0: all +1
            [1, 1, 1, 1],   # optidx=1: all +1
        ]

        # Build QPTrans matrices
        qp_trans, qp_trans_inv, qp_trans_sgn, qp_opt_trans, qp_opt_trans_sgn =
            build_qp_trans_matrices(data)

        # Test: QPOptTrans should use parsed data
        @test length(qp_opt_trans) == 2
        @test qp_opt_trans[1] == [0, 1, 2, 3]
        @test qp_opt_trans[2] == [1, 2, 3, 0]

        # Test: QPOptTransSgn should use parsed data
        @test length(qp_opt_trans_sgn) == 2
        @test qp_opt_trans_sgn[1] == [1, 1, 1, 1]
        @test qp_opt_trans_sgn[2] == [1, 1, 1, 1]
    end
end

"""
    test_update_slater_elm_fcmp_basic()

Test update_slater_elm_fcmp!() with basic data.
"""
function test_update_slater_elm_fcmp_basic()
    @testset "update_slater_elm_fcmp! basic Tests" begin
        # Create test data
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
        data.modpara.nsp_gauss_leg = 1
        data.modpara.nmp_trans = 2
        data.i_flg_orbital_general = 0  # sz-conserved
        data.n_qp_trans = 2
        data.n_qp_opt_trans = 1

        # Add orbital terms (use 4-arg constructor: site1, site2, value, is_complex)
        for i = 0:3
            for j = 0:3
                push!(data.orbital_terms, OrbitalTerm(i, j, 0.1+0.0im, false))
            end
        end

        # Build orbital sign matrix
        build_orbital_sgn_matrix!(data)

        # Set up QPTrans mappings
        data.qp_trans = [
            [0, 1, 2, 3],  # mpidx=0: identity
            [1, 2, 3, 0],   # mpidx=1: shift by 1
        ]
        data.qp_trans_inv = [
            [0, 1, 2, 3],  # mpidx=0: identity
            [3, 0, 1, 2],   # mpidx=1: inverse
        ]
        data.qp_trans_sgn = [
            [1, 1, 1, 1],  # mpidx=0: all +1
            [1, 1, 1, 1],   # mpidx=1: all +1
        ]

        # Initialize quantum projection weights
        data.para_qp_trans = [1.0+0.0im, 1.0+0.0im]
        init_qp_weight!(data)

        # Create VMCOptimizationState
        n_qp_full =
            data.modpara.nsp_gauss_leg * data.modpara.nmp_trans * data.n_qp_opt_trans
        state = VMCOptimizationState(
            data.modpara.nsite,
            data.modpara.nelec,
            0,  # n_proj
            length(data.orbital_terms),  # n_para
            n_qp_full,
            1,  # n_vmc_sample
            false,  # all_complex
            false,  # use_fsz
        )

        # Test: Should not throw error
        try
            update_slater_elm_fcmp!(data, state)
            @test true  # Success
        catch e
            @test false  # Should not throw
            @error "update_slater_elm_fcmp! failed: $e"
        end

        # Test: SlaterElm should be updated (not all zeros)
        n_site2 = 2 * data.modpara.nsite
        slater_elm_size = n_qp_full * n_site2 * n_site2
        @test length(state.slater_matrix.slater_elm) == slater_elm_size
        # Note: We don't check specific values here, just that the function runs without error
    end
end

"""
    test_update_slater_elm_fcmp_missing_qp_weights()

Test update_slater_elm_fcmp!() when quantum projection weights are missing.
"""
function test_update_slater_elm_fcmp_missing_qp_weights()
    @testset "update_slater_elm_fcmp! missing qp_weights Tests" begin
        # Create test data without qp_weights
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
        data.modpara.nsp_gauss_leg = 1
        data.modpara.nmp_trans = 2
        data.i_flg_orbital_general = 0
        data.n_qp_trans = 2
        data.n_qp_opt_trans = 1

        # Add orbital terms (use 4-arg constructor: site1, site2, value, is_complex)
        for i = 0:3
            for j = 0:3
                push!(data.orbital_terms, OrbitalTerm(i, j, 0.1+0.0im, false))
            end
        end

        build_orbital_sgn_matrix!(data)

        # Set up QPTrans mappings
        data.qp_trans = [[0, 1, 2, 3], [1, 2, 3, 0]]
        data.qp_trans_inv = [[0, 1, 2, 3], [3, 0, 1, 2]]
        data.qp_trans_sgn = [[1, 1, 1, 1], [1, 1, 1, 1]]

        # qp_weights is not initialized (should be nothing)
        data.qp_weights = nothing

        # Create VMCOptimizationState
        n_qp_full = 1 * 2 * 1
        state = VMCOptimizationState(
            data.modpara.nsite,
            data.modpara.nelec,
            0,
            length(data.orbital_terms),
            n_qp_full,
            1,
            false,
            false,
        )

        # Test: Should return early (no error, but function returns)
        # The function logs an error and returns, so we just check it doesn't throw
        # Suppress expected @error log during this test
        with_logger(NullLogger()) do
            try
                update_slater_elm_fcmp!(data, state)
                @test true  # Function returns without throwing
            catch e
                @test false  # Should not throw
            end
        end
    end
end

"""
    test_update_slater_elm_fcmp_with_real_data()

Test update_slater_elm_fcmp!() with real HeisenbergChain data.
"""
function test_update_slater_elm_fcmp_with_real_data()
    @testset "update_slater_elm_fcmp! with real data Tests" begin
        # Skip if test data directory doesn't exist
        if !isdir(TEST_DATA_DIR)
            @test_skip "Test data directory not found: $TEST_DATA_DIR"
            return
        end

        namelist_file = joinpath(TEST_DATA_DIR, "namelist.def")
        if !isfile(namelist_file)
            @test_skip "namelist.def not found in test directory"
            return
        end

        # Parse Expert Mode files
        data = parse_expert_mode_files(namelist_file)

        # Initialize parameters with seeded RNG
        rng = SFMT19937RNG()
        Random.seed!(rng, data.modpara.rnd_seed)
        initialize_parameters!(data; rng = rng)

        # Initialize quantum projection weights
        init_qp_weight!(data)

        # Verify that QPTrans mappings are set (from build_qp_trans_mappings! in parse_expert_mode_files)
        if isempty(data.qp_trans)
            @test_skip "QPTrans mappings not found. qptransidx.def may not be parsed correctly."
            return
        end

        # Create VMCOptimizationState
        n_qp_full =
            data.modpara.nsp_gauss_leg * abs(data.modpara.nmp_trans) * data.n_qp_opt_trans
        state = VMCOptimizationState(
            data.modpara.nsite,
            data.modpara.nelec,
            length(data.gutzwiller_terms) + length(data.jastrow_terms),
            length(data.orbital_terms),
            n_qp_full,
            1,
            data.modpara.complex_flag != 0,
            data.i_flg_orbital_general != 0,
        )

        # Test: Should not throw error
        try
            update_slater_elm_fcmp!(data, state)
            @test true  # Success
        catch e
            @test false  # Should not throw
            @error "update_slater_elm_fcmp! failed with real data: $e"
        end

        # Test: SlaterElm should be updated
        n_site2 = 2 * data.modpara.nsite
        slater_elm_size = n_qp_full * n_site2 * n_site2
        @test length(state.slater_matrix.slater_elm) == slater_elm_size
    end
end

# Run all tests
@testset "Slater Update Tests" begin
    test_build_qp_trans_matrices_success()
    test_build_qp_trans_matrices_error()
    test_build_qp_trans_matrices_qp_opt_trans()
    test_update_slater_elm_fcmp_basic()
    test_update_slater_elm_fcmp_missing_qp_weights()
    test_update_slater_elm_fcmp_with_real_data()
end
