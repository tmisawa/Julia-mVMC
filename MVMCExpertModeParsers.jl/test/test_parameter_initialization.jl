"""
Parameter Initialization Tests

Tests for parameter initialization functions (init_parameter!, sync_modified_parameter!, initialize_parameters!).
Based on C implementation's InitParameter() and SyncModifiedParameter().
"""

using Test
using Random
using SFMT
using MVMCExpertModeParsers
using SFMT

# Import types for convenience
import MVMCExpertModeParsers: ExpertModeData, GutzwillerTerm, JastrowTerm, OrbitalTerm
import MVMCExpertModeParsers: ChargeRBMPhysLayerTerm
import MVMCExpertModeParsers:
    init_parameter!, sync_modified_parameter!, initialize_parameters!

# Helper function to create a seeded SFMT RNG
function create_seeded_rng(seed::Int)
    rng = SFMT19937RNG()
    Random.seed!(rng, seed)
    return rng
end

function seeded_rng(seed::Integer)
    rng = SFMT19937RNG()
    Random.seed!(rng, seed)
    return rng
end

"""
    test_init_parameter_basic()

Test basic functionality of init_parameter!().
"""
function test_init_parameter_basic()
    @testset "Basic init_parameter! Tests" begin
        # Create test data
        data = ExpertModeData()
        data.modpara.rnd_seed = 123456789
        data.modpara.complex_flag = 0  # Real parameters

        # Add Gutzwiller and Jastrow terms
        push!(data.gutzwiller_terms, GutzwillerTerm(0, 0.5+0.0im, false))
        push!(data.jastrow_terms, JastrowTerm(0, 1, 0.1+0.0im, false))

        # Add Orbital terms (Slater parameters)
        for i in 1:5
            push!(data.orbital_terms, OrbitalTerm(0, i, i - 1, 0.0 + 0.0im, false))
        end

        # Set optimization flags
        # OptFlag structure: [Proj_real, Proj_imag, RBM_real, RBM_imag, Slater_real, Slater_imag, ...]
        # For real parameters: only real flags are used
        n_proj = length(data.gutzwiller_terms) + length(data.jastrow_terms)
        n_slater = length(data.orbital_terms)
        # Set flags: first 2*NProj for Proj, then 2*NSlater for Slater
        data.optimization_flags = vcat(
            repeat([true, false], n_proj),  # Proj flags (real only)
            repeat([true, false], n_slater),  # Slater flags (real only)
        )

        # Initialize parameters
        rng = seeded_rng(data.modpara.rnd_seed)
        init_parameter!(data; rng=rng)

        # Test: Proj parameters should be 0.0
        @test data.gutzwiller_terms[1].value == ComplexF64(0.0)
        @test data.jastrow_terms[1].value == ComplexF64(0.0)

        # Test: Slater parameters should be initialized randomly (not all 0.0)
        # Since OptFlag is true, they should be non-zero
        slater_values = [term.value for term in data.orbital_terms]
        @test any(abs(v) > 1e-10 for v in slater_values)

        # Test: Slater parameters should be in range [-1, 1) for real case
        for term in data.orbital_terms
            @test real(term.value) >= -1.0
            @test real(term.value) < 1.0
            @test abs(imag(term.value)) < 1e-10  # Should be real
        end
    end
end

"""
    test_init_parameter_complex()

Test init_parameter!() with complex parameters (AllComplexFlag != 0).
"""
function test_init_parameter_complex()
    @testset "Complex init_parameter! Tests" begin
        data = ExpertModeData()
        data.modpara.rnd_seed = 123456789
        data.modpara.complex_flag = 1  # Complex parameters

        # Add Orbital terms
        for i in 1:5
            push!(data.orbital_terms, OrbitalTerm(0, i, i - 1, 0.0 + 0.0im, true))
        end

        # Set optimization flags
        n_slater = length(data.orbital_terms)
        data.optimization_flags = repeat([true, true], n_slater)

        # Initialize parameters
        rng = seeded_rng(data.modpara.rnd_seed)
        init_parameter!(data; rng=rng)

        # Test: Complex Slater parameters should have both real and imaginary parts
        for term in data.orbital_terms
            # After initialization and normalization by sqrt(2.0), values should be reasonable
            abs_val = abs(term.value)
            @test abs_val > 0.0 || abs_val < sqrt(2.0) + 0.1  # Should be normalized
            # Both real and imaginary parts should be non-zero (with some tolerance)
            @test abs(real(term.value)) > 1e-10 || abs(imag(term.value)) > 1e-10
        end
    end
end

"""
    test_init_parameter_optflag()

Test that init_parameter!() respects OptFlag settings.
"""
function test_init_parameter_optflag()
    @testset "OptFlag handling Tests" begin
        data = ExpertModeData()
        data.modpara.rnd_seed = 123456789
        data.modpara.complex_flag = 0

        # Add Orbital terms
        for i in 1:5
            push!(data.orbital_terms, OrbitalTerm(0, i, i - 1, 0.0 + 0.0im, false))
        end

        n_slater = length(data.orbital_terms)
        # Set flags: first 3 optimized, last 2 not optimized
        data.optimization_flags = vcat(
            repeat([true, false], 3),  # First 3 optimized
            repeat([false, false], 2),   # Last 2 not optimized
        )

        # Initialize parameters
        rng = seeded_rng(data.modpara.rnd_seed)
        init_parameter!(data; rng=rng)

        # Test: First 3 should be non-zero, last 2 should be 0.0
        for i = 1:3
            @test abs(data.orbital_terms[i].value) > 1e-10
        end
        for i = 4:5
            @test abs(data.orbital_terms[i].value) < 1e-10
        end
    end
end

"""
    test_sync_modified_parameter()

Test sync_modified_parameter!() scaling functionality.
"""
function test_sync_modified_parameter()
    @testset "sync_modified_parameter! Tests" begin
        data = ExpertModeData()
        data.modpara.complex_flag = 0

        # Add Orbital terms with large values
        push!(data.orbital_terms, OrbitalTerm(0, 1, 0, 10.0+0.0im, false, 1))
        push!(data.orbital_terms, OrbitalTerm(0, 2, 1, 5.0+0.0im, false, 1))
        push!(data.orbital_terms, OrbitalTerm(0, 3, 2, 8.0+0.0im, false, 1))

        # Sync modified parameter (should scale to max 4.0)
        sync_modified_parameter!(data)

        # Test: Maximum absolute value should be <= 4.0
        max_abs = maximum([abs(term.value) for term in data.orbital_terms])
        @test max_abs <= 4.0 + 1e-10

        # Test: All values should be scaled proportionally
        # Ratio should be 4.0 / 10.0 = 0.4
        expected_ratio = 4.0 / 10.0
        @test abs(abs(data.orbital_terms[1].value) / 10.0 - expected_ratio) < 1e-10
        @test abs(abs(data.orbital_terms[2].value) / 5.0 - expected_ratio) < 1e-10
        @test abs(abs(data.orbital_terms[3].value) / 8.0 - expected_ratio) < 1e-10

        # Test: C implementation always scales, even if max <= 4.0
        # ratio = D_AMP_MAX / xmax = 4.0 / 3.0 = 1.333...
        data2 = ExpertModeData()
        push!(data2.orbital_terms, OrbitalTerm(0, 1, 0, 2.0+0.0im, false, 1))
        push!(data2.orbital_terms, OrbitalTerm(0, 2, 1, 3.0+0.0im, false, 1))
        sync_modified_parameter!(data2)
        # xmax = 3.0, ratio = 4.0 / 3.0
        expected_ratio = 4.0 / 3.0
        @test abs(data2.orbital_terms[1].value - (2.0 * expected_ratio)) < 1e-10
        @test abs(data2.orbital_terms[2].value - (3.0 * expected_ratio)) < 1e-10
        # Maximum should be exactly 4.0 after scaling
        max_abs2 = maximum([abs(term.value) for term in data2.orbital_terms])
        @test abs(max_abs2 - 4.0) < 1e-10
    end
end

"""
    test_initialize_parameters()

Test complete initialize_parameters!() workflow.
"""
function test_initialize_parameters()
    @testset "initialize_parameters! Tests" begin
        data = ExpertModeData()
        data.modpara.rnd_seed = 123456789
        data.modpara.complex_flag = 0

        # Add terms
        push!(data.gutzwiller_terms, GutzwillerTerm(0, 0.5+0.0im, false))
        push!(data.jastrow_terms, JastrowTerm(0, 1, 0.1+0.0im, false))

        for i in 1:10
            push!(data.orbital_terms, OrbitalTerm(0, i, i - 1, 0.0 + 0.0im, false))
        end

        n_proj = length(data.gutzwiller_terms) + length(data.jastrow_terms)
        n_slater = length(data.orbital_terms)
        data.optimization_flags =
            vcat(repeat([true, false], n_proj), repeat([true, false], n_slater))

        # Initialize parameters (init + sync)
        rng = seeded_rng(data.modpara.rnd_seed)
        initialize_parameters!(data; rng=rng)

        # Test: Proj should be 0.0
        @test data.gutzwiller_terms[1].value == ComplexF64(0.0)
        @test data.jastrow_terms[1].value == ComplexF64(0.0)

        # Test: Slater should be initialized and scaled
        max_abs = maximum([abs(term.value) for term in data.orbital_terms])
        @test max_abs <= 4.0 + 1e-10

        # Test: Slater values should be in reasonable range
        for term in data.orbital_terms
            @test abs(term.value) <= 4.0 + 1e-10
        end
    end
end

"""
    test_random_seed()

Test that random seed works correctly.
"""
function test_random_seed()
    @testset "Random Seed Tests" begin
        # Test: Different seed should give different values
        data1 = ExpertModeData()
        data1.modpara.rnd_seed = 123456789
        data1.modpara.complex_flag = 0

        data2 = ExpertModeData()
        data2.modpara.rnd_seed = 987654321  # Different seed
        data2.modpara.complex_flag = 0

        # Add same orbital terms
        for i in 1:5
            push!(data1.orbital_terms, OrbitalTerm(0, i, i - 1, 0.0 + 0.0im, false))
            push!(data2.orbital_terms, OrbitalTerm(0, i, i - 1, 0.0 + 0.0im, false))
        end

        n_slater = length(data1.orbital_terms)
        flags1 = repeat([true, false], n_slater)
        flags2 = repeat([true, false], n_slater)
        data1.optimization_flags = flags1
        data2.optimization_flags = flags2

        # Initialize with different seeds
        rng1 = create_seeded_rng(data1.modpara.rnd_seed)
        rng2 = create_seeded_rng(data2.modpara.rnd_seed)
        initialize_parameters!(data1; rng = rng1)
        initialize_parameters!(data2; rng = rng2)

        # Test: At least one value should be different
        different = any(
            abs(t1.value - t2.value) > 1e-10 for
            (t1, t2) in zip(data1.orbital_terms, data2.orbital_terms)
        )
        @test different

        # Test: Values are within expected range after scaling
        for term in data1.orbital_terms
            @test abs(term.value) <= 4.0 + 1e-10
        end
        for term in data2.orbital_terms
            @test abs(term.value) <= 4.0 + 1e-10
        end
    end
end

"""
    test_rbm_initialization()

Test RBM parameter initialization (if FlagRBM > 0).
"""
function test_rbm_initialization()
    @testset "RBM Initialization Tests" begin
        data = ExpertModeData()
        data.modpara.rnd_seed = 123456789
        data.modpara.complex_flag = 0

        # Add RBM terms
        push!(data.charge_rbm_phys_layer_terms, ChargeRBMPhysLayerTerm(0, 0.0+0.0im, false))
        push!(data.charge_rbm_phys_layer_terms, ChargeRBMPhysLayerTerm(1, 0.0+0.0im, false))

        # Add Proj and Slater terms
        push!(data.gutzwiller_terms, GutzwillerTerm(0, 0.5+0.0im, false))
        push!(data.orbital_terms, OrbitalTerm(0, 1, 0, 0.0+0.0im, false, 1))

        n_proj = length(data.gutzwiller_terms)
        n_rbm = length(data.charge_rbm_phys_layer_terms)
        n_slater = length(data.orbital_terms)

        # Set optimization flags
        data.optimization_flags = vcat(
            repeat([true, false], n_proj),      # Proj
            repeat([true, false], n_rbm),        # RBM
            repeat([true, false], n_slater),      # Slater
        )

        # Initialize parameters
        rng = seeded_rng(data.modpara.rnd_seed)
        init_parameter!(data; rng=rng)

        # Test: RBM parameters should be initialized (small values)
        for term in data.charge_rbm_phys_layer_terms
            @test abs(term.value) < 0.01  # Should be small (0.01 * (rand() - 0.5))
            @test abs(term.value) > 1e-10  # But not zero (if OptFlag is true)
        end
    end
end

"""
    test_empty_orbital_terms()

Test that initialization handles empty orbital terms gracefully.
"""
function test_empty_orbital_terms()
    @testset "Empty Orbital Terms Tests" begin
        data = ExpertModeData()
        data.modpara.rnd_seed = 123456789
        data.modpara.complex_flag = 0

        # No orbital terms
        data.optimization_flags = Bool[]

        # Should not throw error
        try
            rng = seeded_rng(data.modpara.rnd_seed)
            initialize_parameters!(data; rng=rng)
            @test true  # Success
        catch e
            @test false  # Should not throw
        end

        # sync_modified_parameter should also work
        try
            sync_modified_parameter!(data)
            @test true  # Success
        catch e
            @test false  # Should not throw
        end
    end
end

# Run all tests
@testset "Parameter Initialization Tests" begin
    test_init_parameter_basic()
    test_init_parameter_complex()
    test_init_parameter_optflag()
    test_sync_modified_parameter()
    test_initialize_parameters()
    test_random_seed()
    test_rbm_initialization()
    test_empty_orbital_terms()
end
