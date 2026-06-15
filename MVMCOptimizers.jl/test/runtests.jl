using Test
using MVMCOptimizers
using MVMCExpertModeParsers  # For legendre_poly and gauss_legendre
using LinearAlgebra

# Import functions from MVMCExpertModeParsers
using MVMCExpertModeParsers: legendre_poly, gauss_legendre

@testset "MVMCOptimizers.jl" begin
    @testset "Legendre Polynomial" begin
        # Test base cases
        @test legendre_poly(0.5, 0) == 1.0
        @test legendre_poly(0.5, 1) == 0.5

        # Test known values of Legendre polynomials
        # P_2(x) = (3x^2 - 1)/2
        x = 0.5
        @test legendre_poly(x, 2) ≈ (3*x^2 - 1)/2

        # P_3(x) = (5x^3 - 3x)/2
        @test legendre_poly(x, 3) ≈ (5*x^3 - 3*x)/2

        # Test at x = 1 (all Legendre polynomials equal 1)
        for n = 0:5
            @test legendre_poly(1.0, n) ≈ 1.0
        end
    end

    @testset "Gauss-Legendre Quadrature" begin
        # Test n=1 case
        # gauss_legendre returns a tuple (x, w), not in-place
        x, w = gauss_legendre(0.0, Float64(π), 1)
        @test x[1] ≈ π/2
        @test w[1] ≈ π

        # Test standard interval [-1, 1]
        n = 5
        x, w = gauss_legendre(-1.0, 1.0, n)

        # Check that weights sum to interval length
        @test sum(w) ≈ 2.0

        # Check that points are in interval
        @test all(-1.0 .<= x .<= 1.0)

        # Test exact integration of polynomials up to degree 2n-1
        # For n=5, should integrate polynomials up to degree 9 exactly
        # Test with x^4 (degree 4)
        integral = sum(w .* (x .^ 4))
        exact = 2/5  # ∫_{-1}^{1} x^4 dx = 2/5
        @test integral ≈ exact rtol=1e-10
    end

end

@testset "Slater Update Tests" begin
    include("test_slater_update.jl")
end

@testset "Unit Tests" begin
    include("../test_unit/helpers/mock_state.jl")
    include("../test_unit/helpers/mock_data.jl")
    include("../test_unit/test_unit_stochastic_opt.jl")
    include("../test_unit/test_unit_vmc_sampling_rbm.jl")
    include("../test_unit/test_unit_vmc_sampling_proj.jl")
    include("../test_unit/test_unit_vmc_sampling_misc.jl")
    include("../test_unit/test_unit_slater_update.jl")
    include("../test_unit/test_unit_vmc_main_cal_sr.jl")
    include("../test_unit/test_unit_parameter_sync.jl")
    include("../test_unit/test_unit_types.jl")
    include("../test_unit/test_unit_threading.jl")
    include("../test_unit/test_unit_parallel.jl")
    include("../test_unit/test_unit_weight_average.jl")
    include("../test_unit/test_unit_unsupported_inputs.jl")
    include("../test_unit/test_unit_physcal_factored_green.jl")
    include("../test_unit/test_unit_read_opt_para.jl")
    include("../test_unit/test_unit_run_phys_cal_runner.jl")
end

# Integration tests against the C reference live at the workspace root
# (test/integration/runtests.jl) so that input fixtures and reference data
# can be bundled alongside the public release. They are not part of the
# subpackage test suite.
