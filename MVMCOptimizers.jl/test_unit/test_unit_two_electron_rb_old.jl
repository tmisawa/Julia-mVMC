using Test
using LinearAlgebra
using Random
using MVMCOptimizers
using MVMCExpertModeParsers: ExpertModeData, ModParaParameters

# C mVMC 1.4.0 pfupdate_two_real.c:230 (upstream #100): the two-electron update
# reads mOld_ab = sltE[rsaOld][rsbOld] with rsbOld = rbOld + t*Nsite. 1.3.0 used
# raOld, i.e. the diagonal element sltE[rsaOld][raOld + t*Nsite] for s == t.
# The element cancels mathematically, so the check poisons the 1.3.0 element
# with NaN (it is read nowhere else for this move) after the old inverse is
# formed: the 1.4.0 update must stay finite and match the explicit new inverse.
@testset "unit/two-electron real update reads the old b-site (C 1.4.0 rbOld)" begin
    nsite, ne = 4, 2
    n2, nsz = 2nsite, 2ne
    pos(idx) = [idx[k] + (k > ne ? nsite : 0) for k = 1:nsz] .+ 1
    rng = MersenneTwister(20261106)
    A = randn(rng, n2, n2)
    F = A - A'
    old = [0, 1, 0, 2]       # up 0,1; down 0,2
    new = [2, 3, 0, 2]       # up 0 -> 2 (ma=0) and 1 -> 3 (mb=1), s = t = 0
    Mold = F[pos(old), pos(old)]
    Mnew = F[pos(new), pos(new)]
    Fpoison = copy(F)
    Fpoison[1, 1] = NaN      # sltE[raOld][raOld]: the element 1.3.0 read
    for (fill_f, label) in ((F, "plain"), (Fpoison, "1.3.0 element poisoned"))
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite = nsite, nelec = ne)
        state = MVMCOptimizers.VMCOptimizationState(nsite, ne, 0, 1, 1, 1, false, false)
        sm = state.slater_matrix
        resize!(sm.slater_elm_real, n2 * n2)
        resize!(sm.inv_m_real, nsz * nsz)
        resize!(sm.pf_m_real, 1)
        for i = 1:n2, j = 1:n2
            sm.slater_elm_real[(i-1)*n2+j] = fill_f[i, j]   # C row-major
        end
        Mi = inv(Mold)
        for i = 1:nsz, j = 1:nsz
            sm.inv_m_real[(i-1)*nsz+j] = Mi[i, j]
        end
        sm.pf_m_real[1] = 1.0
        MVMCOptimizers.update_m_all_two_real!(0, 0, 1, 0, 0, 1, new, 1, 2, data, state)
        updated = permutedims(reshape(copy(sm.inv_m_real), nsz, nsz))
        @test all(isfinite, updated)
        @test isfinite(sm.pf_m_real[1])
        # Rank-2 update of a 4x4 inverse with O(1) entries: 64 eps * cond budget.
        exact = inv(Mnew)
        budget = 64 * eps(Float64) * cond(Mold) * cond(Mnew)
        @test maximum(abs.(updated - exact)) <= budget * maximum(abs.(exact))
        # |Pf(Mnew)/Pf(Mold)| = sqrt(|det(Mnew)/det(Mold)|)
        @test isapprox(abs(sm.pf_m_real[1]), sqrt(abs(det(Mnew) / det(Mold))); rtol = budget)
    end
end
