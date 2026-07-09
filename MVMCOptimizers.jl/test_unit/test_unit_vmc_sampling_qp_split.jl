using Test
using MVMCOptimizers
using MVMCOptimizers: calculate_ip_fcmp, calculate_log_ip_fcmp,
                      calculate_ip_real, calculate_log_ip_real
using MVMCExpertModeParsers

function qp_weight_data(weights::Vector{ComplexF64})
    data = MVMCExpertModeParsers.ExpertModeData()
    data.qp_weights = MVMCExpertModeParsers.QuantumProjectionWeights()
    data.qp_weights.qp_full_weight = copy(weights)
    return data
end

@testset "IP helpers default to no reduction and allow empty ranges" begin
    data = qp_weight_data(ComplexF64[2.0 + 0.0im, 3.0 + 0.0im, 5.0 + 0.0im])
    pf = ComplexF64[7.0 + 0.0im, 11.0 + 0.0im, 13.0 + 0.0im]
    @test calculate_ip_fcmp(pf, 1, 4, data) == 112.0 + 0.0im
    @test calculate_ip_fcmp(pf, 2, 3, data) == 33.0 + 0.0im
    @test calculate_ip_fcmp(pf, 2, 2, data) == 0.0 + 0.0im
    @test calculate_log_ip_fcmp(pf, 1, 4, data) ≈ log(112.0 + 0.0im)
end

@testset "real IP helpers default to no reduction and allow empty ranges" begin
    data = qp_weight_data(ComplexF64[2.0 + 9.0im, 3.0 + 0.0im])
    pf = [7.0, 11.0]
    @test calculate_ip_real(pf, 1, 3, data) == 47.0
    @test calculate_ip_real(pf, 2, 2, data) == 0.0
    @test calculate_log_ip_real(pf, 1, 3, data) ≈ log(abs(47.0) + 1e-100)
end

@testset "comm1 reduction requires a context" begin
    data = qp_weight_data(ComplexF64[1.0 + 0.0im])
    pf = ComplexF64[2.0 + 0.0im]
    @test_throws ArgumentError calculate_ip_fcmp(pf, 1, 2, data; reduce = :comm1)
    @test calculate_ip_fcmp(pf, 1, 2, data; ctx = serial_context(), reduce = :comm1) ==
          2.0 + 0.0im
end
