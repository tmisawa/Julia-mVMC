using Test
using MVMCOptimizers
using MVMCExpertModeParsers
using MVMCExpertModeParsers: ExpertModeData, ModParaParameters, OrbitalTerm

@testset "unit/slater_update: build_orbital_idx_sgn_matrices" begin
    data = ExpertModeData()
    data.modpara = ModParaParameters(nsite = 3)
    data.orbital_terms = [
        OrbitalTerm(0, 1, 0, 1.0 + 0.0im, true, 1),
        OrbitalTerm(1, 2, 2, 2.5 + 0.0im, true, -1),
        OrbitalTerm(2, 0, 2, 0.0 + 0.0im, true, 1), # duplicate idx with zero value (should not overwrite)
    ]

    orbital_idx, orbital_sgn, slater = MVMCOptimizers.build_orbital_idx_sgn_matrices(data)

    @test size(orbital_idx) == (3, 3)
    @test size(orbital_sgn) == (3, 3)
    @test length(slater) == 3  # max idx=2 => 3 entries

    @test orbital_idx[0 + 1, 1 + 1] == 0
    @test orbital_idx[1 + 1, 2 + 1] == 2

    # Fallback sign matrix built from orbital_terms
    @test orbital_sgn[0 + 1, 1 + 1] == 1
    @test orbital_sgn[1 + 1, 2 + 1] == -1

    # Slater values placed by idx; zero-valued duplicate does not overwrite.
    @test slater[0 + 1] == 1.0 + 0.0im
    @test slater[2 + 1] == 2.5 + 0.0im
end

@testset "unit/slater_update: build_orbital_idx_sgn_matrices uses data.orbital_sgn when provided" begin
    data = ExpertModeData()
    data.modpara = ModParaParameters(nsite = 2)
    data.orbital_terms = [
        OrbitalTerm(0, 1, 0, 1.0 + 0.0im, true, 1),
    ]

    provided = [1 -1; -1 1]
    data.orbital_sgn = provided

    _, orbital_sgn, _ = MVMCOptimizers.build_orbital_idx_sgn_matrices(data)
    @test orbital_sgn === provided
end

@testset "unit/slater_update: build_qp_trans_matrices requires qp_trans mappings" begin
    data = ExpertModeData()
    data.modpara = ModParaParameters(nsite = 2)
    data.n_qp_trans = 1
    data.n_qp_opt_trans = 1
    data.qp_trans = Vector{Vector{Int}}()
    data.qp_trans_inv = Vector{Vector{Int}}()
    data.qp_trans_sgn = Vector{Vector{Int}}()

    @test_logs (:error,) @test_throws ArgumentError MVMCOptimizers.build_qp_trans_matrices(data)
end

@testset "unit/slater_update: build_qp_trans_matrices opt fallback is identity" begin
    data = ExpertModeData()
    data.modpara = ModParaParameters(nsite = 3)
    data.n_qp_trans = 1
    data.n_qp_opt_trans = 2

    data.qp_trans = [collect(0:2)]
    data.qp_trans_inv = [collect(0:2)]
    data.qp_trans_sgn = [ones(Int, 3)]

    # Leave qp_opt_trans empty => identity fallback
    data.qp_opt_trans = Vector{Vector{Int}}()
    data.qp_opt_trans_sgn = Vector{Vector{Int}}()

    _, _, _, qp_opt_trans, qp_opt_trans_sgn = MVMCOptimizers.build_qp_trans_matrices(data)

    @test length(qp_opt_trans) == 2
    @test length(qp_opt_trans_sgn) == 2
    for k = 1:2
        @test qp_opt_trans[k] == [0, 1, 2]
        @test qp_opt_trans_sgn[k] == [1, 1, 1]
    end
end
