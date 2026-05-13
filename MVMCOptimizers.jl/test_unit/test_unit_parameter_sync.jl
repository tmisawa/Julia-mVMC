using Test
using MVMCOptimizers
using MVMCExpertModeParsers
using MVMCExpertModeParsers: ExpertModeData, GutzwillerTerm, JastrowTerm, OrbitalTerm

@testset "unit/parameter_sync: sync_modified_parameter!" begin
    @testset "shift gutzwiller/jastrow when all optimized (default)" begin
        data = ExpertModeData()
        data.gutzwiller_terms = [
            GutzwillerTerm(0, 2.0 + 1.5im, true),
            GutzwillerTerm(1, -1.0 + 0.25im, true),
        ]
        data.jastrow_terms = [
            JastrowTerm(0, 1, 5.0 - 0.75im, true),
        ]

        # optimization_flags empty => treat all as optimized => shift should be applied
        data.optimization_flags = Bool[]

        orig = vcat(
            (t.value for t in data.gutzwiller_terms)...,
            (t.value for t in data.jastrow_terms)...,
        )
        shift = sum(real.(orig)) / length(orig)

        ret = MVMCOptimizers.sync_modified_parameter!(data)
        @test ret == 0

        newvals = vcat(
            (t.value for t in data.gutzwiller_terms)...,
            (t.value for t in data.jastrow_terms)...,
        )
        @test (sum(real.(newvals)) / length(newvals)) ≈ 0.0 atol = 1e-14

        # Shift is real-only: imag parts preserved
        @test imag.(newvals) == imag.(orig)
        @test real.(newvals) ≈ (real.(orig) .- shift) atol = 1e-14
    end

    @testset "no shift when not all optimized (explicit flags)" begin
        data = ExpertModeData()
        data.gutzwiller_terms = [
            GutzwillerTerm(0, 2.0 + 0.0im, true),
            GutzwillerTerm(1, -1.0 + 0.0im, true),
        ]
        data.jastrow_terms = [
            JastrowTerm(0, 1, 5.0 + 0.0im, true),
        ]

        # Need length >= max index accessed by is_gutzwiller_optimized/is_jastrow_optimized.
        # gutz idx 1 -> 1, gutz idx 2 -> 3, jast idx 1 (offset by n_gutz=2) -> 5
        flags = falses(5)
        flags[1] = true
        flags[3] = false  # not all gutz optimized => no shift
        flags[5] = true
        data.optimization_flags = flags

        orig_g = [t.value for t in data.gutzwiller_terms]
        orig_j = [t.value for t in data.jastrow_terms]

        MVMCOptimizers.sync_modified_parameter!(data)

        @test [t.value for t in data.gutzwiller_terms] == orig_g
        @test [t.value for t in data.jastrow_terms] == orig_j
    end

    @testset "rescale orbital terms to D_AMP_MAX" begin
        data = ExpertModeData()
        data.orbital_terms = [
            OrbitalTerm(0, 0, 0, 2.0 + 0.0im, true),
            OrbitalTerm(0, 1, 1, -6.0 + 8.0im, true), # abs = 10
            OrbitalTerm(1, 1, 2, 0.5 - 0.5im, true),
        ]

        orig = [t.value for t in data.orbital_terms]
        xmax = maximum(abs.(orig))
        ratio = MVMCOptimizers.D_AMP_MAX / xmax

        MVMCOptimizers.sync_modified_parameter!(data)

        newvals = [t.value for t in data.orbital_terms]
        @test maximum(abs.(newvals)) ≈ MVMCOptimizers.D_AMP_MAX atol = 1e-14
        @test newvals ≈ (orig .* ratio) atol = 1e-14
    end

    @testset "para_qp_trans is NOT normalized (C normalizes OptTrans, not ParaQPTrans)" begin
        # C's SyncModifiedParameter (parameter.c:163-175) normalizes the
        # dedicated OptTrans[] array gated by FlagOptTrans>0. ParaQPTrans
        # (the QPTrans phase factors consumed by qp_weight_update.jl) is a
        # different array and must remain untouched. The Julia loader had
        # historically rescaled para_qp_trans, which silently corrupted
        # quantum projection weights — see B5 review round 2.
        data = ExpertModeData()
        data.n_qp_opt_trans = 1
        data.para_qp_trans = ComplexF64[
            1.0 + 0.0im,
            -2.0 + 2.0im, # abs = sqrt(8); deliberately not normalized
        ]

        orig = copy(data.para_qp_trans)

        MVMCOptimizers.sync_modified_parameter!(data)

        # Must be byte-identical to the input.
        @test data.para_qp_trans == orig
    end
end
