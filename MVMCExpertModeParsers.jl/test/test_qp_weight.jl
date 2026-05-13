"""
Quantum Projection Weight Tests

Tests for init_qp_weight! and update_qp_weight! functions (equivalent to C's InitQPWeight() and UpdateQPWeight()).
"""

using Test
using MVMCExpertModeParsers

import MVMCExpertModeParsers: ExpertModeData, QuantumProjectionWeights
import MVMCExpertModeParsers:
    init_qp_weight!, update_qp_weight!, gauss_legendre, legendre_poly

"""
    test_gauss_legendre()

Test Gauss-Legendre quadrature function.
"""
function test_gauss_legendre()
    @testset "Gauss-Legendre Quadrature Tests" begin
        # Test with n=1
        x, w = gauss_legendre(0.0, 1.0, 1)
        @test length(x) == 1
        @test length(w) == 1
        @test abs(x[1] - 0.5) < 1e-10
        @test abs(w[1] - 1.0) < 1e-10

        # Test with n=4
        x, w = gauss_legendre(0.0, Float64(π), 4)
        @test length(x) == 4
        @test length(w) == 4
        # Check symmetry: x[i] + x[n-i+1] = π for integration range [0, π]
        @test abs(x[1] + x[4] - Float64(π)) < 1e-10
        @test abs(x[2] + x[3] - Float64(π)) < 1e-10
        @test abs(w[1] - w[4]) < 1e-10
        @test abs(w[2] - w[3]) < 1e-10

        # Test integration: ∫₀^π sin(x) dx = 2
        x, w = gauss_legendre(0.0, Float64(π), 8)
        integral = sum(w .* sin.(x))
        @test abs(integral - 2.0) < 1e-10
    end
end

"""
    test_legendre_poly()

Test Legendre polynomial function.
"""
function test_legendre_poly()
    @testset "Legendre Polynomial Tests" begin
        # P_0(x) = 1
        @test abs(legendre_poly(0.5, 0) - 1.0) < 1e-10

        # P_1(x) = x
        @test abs(legendre_poly(0.5, 1) - 0.5) < 1e-10

        # P_2(x) = (3x^2 - 1)/2
        x_val = 0.5
        p2_expected = (3.0 * x_val^2 - 1.0) / 2.0
        @test abs(legendre_poly(x_val, 2) - p2_expected) < 1e-10

        # P_3(x) = (5x^3 - 3x)/2
        p3_expected = (5.0 * x_val^3 - 3.0 * x_val) / 2.0
        @test abs(legendre_poly(x_val, 3) - p3_expected) < 1e-10
    end
end

"""
    test_init_qp_weight_basic()

Test basic init_qp_weight! functionality.
"""
function test_init_qp_weight_basic()
    @testset "Basic init_qp_weight! Tests" begin
        data = ExpertModeData()
        data.modpara.nsp_gauss_leg = 4
        data.modpara.nsp_stot = 1
        data.modpara.nmp_trans = 1
        data.para_qp_trans = [1.0 + 0.0im]

        init_qp_weight!(data)

        @test data.qp_weights !== nothing
        @test length(data.qp_weights.qp_full_weight) == 4
        @test length(data.qp_weights.qp_fix_weight) == 4
        @test length(data.qp_weights.spgl_cos) == 4
        @test length(data.qp_weights.spgl_sin) == 4

        # Check that weights are non-zero
        @test any(abs(w) > 1e-10 for w in data.qp_weights.qp_fix_weight)
    end
end

"""
    test_init_qp_weight_nsp_gauss_leg_1()

Test init_qp_weight! with NSPGaussLeg=1 (special case).
"""
function test_init_qp_weight_nsp_gauss_leg_1()
    @testset "init_qp_weight! NSPGaussLeg=1 Tests" begin
        data = ExpertModeData()
        data.modpara.nsp_gauss_leg = 1
        data.modpara.nsp_stot = 0
        data.modpara.nmp_trans = 2
        data.para_qp_trans = [1.0 + 0.0im, 0.5 + 0.0im]

        init_qp_weight!(data)

        @test data.qp_weights !== nothing
        @test length(data.qp_weights.qp_fix_weight) == 2
        @test data.qp_weights.spgl_cos[1] == ComplexF64(1.0, 0.0)
        @test data.qp_weights.spgl_sin[1] == ComplexF64(0.0, 0.0)
        @test data.qp_weights.qp_fix_weight[1] == ComplexF64(1.0, 0.0)
        @test data.qp_weights.qp_fix_weight[2] == ComplexF64(0.5, 0.0)
    end
end

"""
    test_update_qp_weight()

Test update_qp_weight! function.
"""
function test_update_qp_weight()
    @testset "update_qp_weight! Tests" begin
        weights = QuantumProjectionWeights()
        weights.qp_fix_weight = [1.0 + 0.0im, 2.0 + 0.0im, 3.0 + 0.0im]

        # Test without OptTrans (FlagOptTrans == 0)
        update_qp_weight!(weights, ComplexF64[])
        @test length(weights.qp_full_weight) == 3
        @test weights.qp_full_weight == weights.qp_fix_weight

        # Test with OptTrans (FlagOptTrans > 0)
        opt_trans = [0.5 + 0.0im, 1.5 + 0.0im]
        update_qp_weight!(weights, opt_trans)
        @test length(weights.qp_full_weight) == 6  # 3 * 2
        @test weights.qp_full_weight[1] == ComplexF64(0.5, 0.0)  # 0.5 * 1.0
        @test weights.qp_full_weight[2] == ComplexF64(1.0, 0.0)  # 0.5 * 2.0
        @test weights.qp_full_weight[3] == ComplexF64(1.5, 0.0)  # 0.5 * 3.0
        @test weights.qp_full_weight[4] == ComplexF64(1.5, 0.0)  # 1.5 * 1.0
        @test weights.qp_full_weight[5] == ComplexF64(3.0, 0.0)  # 1.5 * 2.0
        @test weights.qp_full_weight[6] == ComplexF64(4.5, 0.0)  # 1.5 * 3.0
    end
end

"""
    test_init_qp_weight_trigonometric()

Test trigonometric values calculation.
"""
function test_init_qp_weight_trigonometric()
    @testset "Trigonometric Values Tests" begin
        data = ExpertModeData()
        data.modpara.nsp_gauss_leg = 4
        data.modpara.nsp_stot = 1
        data.modpara.nmp_trans = 1
        data.para_qp_trans = [1.0 + 0.0im]

        init_qp_weight!(data)

        weights = data.qp_weights

        # Check trigonometric identities
        for i = 1:length(weights.spgl_cos)
            cos_val = weights.spgl_cos[i]
            sin_val = weights.spgl_sin[i]

            # cos^2 + sin^2 = 1 (for beta/2)
            cos_sq = abs(cos_val)^2
            sin_sq = abs(sin_val)^2
            @test abs(cos_sq + sin_sq - 1.0) < 1e-10

            # Check computed products
            @test abs(weights.spgl_cos_sin[i] - cos_val * sin_val) < 1e-10
            @test abs(weights.spgl_cos_cos[i] - cos_val * cos_val) < 1e-10
            @test abs(weights.spgl_sin_sin[i] - sin_val * sin_val) < 1e-10
        end
    end
end

# Run all tests
@testset "Quantum Projection Weight Tests" begin
    test_gauss_legendre()
    test_legendre_poly()
    test_init_qp_weight_basic()
    test_init_qp_weight_nsp_gauss_leg_1()
    test_update_qp_weight()
    test_init_qp_weight_trigonometric()
end
