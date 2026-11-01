using Test, Random
using MVMCOptimizers
using MVMCExpertModeParsers: ExpertModeData

@testset "Lanczos intermediate Pfaffian updates" begin
    M = MVMCOptimizers
    data = ExpertModeData()
    data.modpara.nsite = 6
    data.modpara.nelec = 3
    state = M.VMCOptimizationState(6, 3, 0, 0, 2, 1, false, false)
    rng = MersenneTwister(89)
    slater = state.slater_matrix.slater_elm_real
    # Full antisymmetric Slater matrices, two independent QP sectors.
    for q in 0:1
        a = randn(rng, 12, 12)
        a = a - transpose(a)
        for i in 0:11, j in 0:11
            slater[q * 144 + i * 12 + j + 1] = a[i+1, j+1]
        end
    end
    original = [0, 1, 2, 0, 1, 2]
    @test M.calculate_m_all_real!(original, 1, 3, data, state) == 0
    original_pf = copy(state.slater_matrix.pf_m_real)
    original_inv = copy(state.slater_matrix.inv_m_real)
    updated = similar(original_pf)
    direct = similar(original_pf)
    for moved in ([0,1,2,0,1,2], [3,1,2,0,1,2], [3,4,2,0,1,2],
                  [3,1,2,4,1,2], [3,4,5,0,1,2])
        @test M._lanczos_updated_pfaffians_real!(moved, original, data, state, updated) == 0
        @test M._lanczos_pfaffians_real!(moved, data, state, direct) == 0
        @test updated ≈ direct rtol=1e-11 atol=1e-11
        @test state.slater_matrix.pf_m_real == original_pf
        @test state.slater_matrix.inv_m_real == original_inv
    end
    for _ in 1:100
        moved = [randperm(rng, 6)[1:3] .- 1; randperm(rng, 6)[1:3] .- 1]
        @test M._lanczos_updated_pfaffians_real!(moved, original, data, state, updated) == 0
        @test M._lanczos_pfaffians_real!(moved, data, state, direct) == 0
        @test updated ≈ direct rtol=1e-11 atol=1e-11
    end
    # A zero row in one QP sector makes a moved configuration singular. The
    # update must take the same failure path as a fresh factorization.
    for j in 0:11
        slater[3 * 12 + j + 1] = 0.0
        slater[j * 12 + 3 + 1] = 0.0
    end
    nodal = [3,1,2,0,1,2]
    @test M._lanczos_updated_pfaffians_real!(nodal, original, data, state, updated) != 0
    @test M._lanczos_pfaffians_real!(nodal, data, state, direct) != 0
end
