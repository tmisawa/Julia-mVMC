using Test
using MVMCOptimizers
using MVMCExpertModeParsers
using MVMCExpertModeParsers:
    ExpertModeData,
    GutzwillerTerm,
    JastrowTerm,
    OrbitalTerm,
    DoublonHolon2SiteIndex,
    DoublonHolon4SiteIndex,
    projection_layout

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

    @testset "correlation shifts can be disabled while Slater rescale remains enabled" begin
        data = ExpertModeData()
        data.gutzwiller_terms = [GutzwillerTerm(0, 10.0 + 0.25im, true)]
        data.jastrow_terms = [JastrowTerm(0, 1, 20.0 - 0.5im, true)]
        data.doublon_holon_2site_indices = [DoublonHolon2SiteIndex([1 0; 0 1])]
        data.doublon_holon_2site_params = ComplexF64[i + i * im for i = 1:6]
        data.orbital_terms = [
            OrbitalTerm(0, 0, 0, 2.0 + 0.0im, true),
            OrbitalTerm(0, 1, 1, 0.0 + 8.0im, true),
        ]
        data.optimization_flags = fill(true, 2 * projection_layout(data).n_proj)

        orig_g = [t.value for t in data.gutzwiller_terms]
        orig_j = [t.value for t in data.jastrow_terms]
        orig_dh = copy(data.doublon_holon_2site_params)
        orig_orb = [t.value for t in data.orbital_terms]
        ratio = MVMCOptimizers.D_AMP_MAX / maximum(abs.(orig_orb))

        MVMCOptimizers.sync_modified_parameter!(data; shift_correlations = false)

        @test [t.value for t in data.gutzwiller_terms] == orig_g
        @test [t.value for t in data.jastrow_terms] == orig_j
        @test data.doublon_holon_2site_params == orig_dh
        @test [t.value for t in data.orbital_terms] ≈ (orig_orb .* ratio) atol = 1e-14
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

    @testset "DH2 shift subtracts bin averages and shifts Gutzwiller" begin
        data = ExpertModeData()
        data.gutzwiller_terms = [GutzwillerTerm(0, 10.0 + 0.25im, true)]
        data.doublon_holon_2site_indices = [DoublonHolon2SiteIndex([1 0; 0 1])]
        data.doublon_holon_2site_params = ComplexF64[
            1 + 1im,
            2 + 2im,
            3 + 3im,
            4 + 4im,
            5 + 5im,
            6 + 6im,
        ]
        data.optimization_flags = fill(true, 2 * projection_layout(data).n_proj)

        MVMCOptimizers.sync_modified_parameter!(data)

        @test data.doublon_holon_2site_params == ComplexF64[
            -2 + 1im,
            -2 + 2im,
            0 + 3im,
            0 + 4im,
            2 + 5im,
            2 + 6im,
        ]
        @test data.gutzwiller_terms[1].value == 17.0 + 0.25im
    end

    @testset "DH4 shift subtracts five-bin averages" begin
        data = ExpertModeData()
        data.gutzwiller_terms = [GutzwillerTerm(0, 1.0 + 0.0im, true)]
        data.doublon_holon_4site_indices = [DoublonHolon4SiteIndex([1 0 1 0; 0 1 0 1])]
        data.doublon_holon_4site_params = ComplexF64[
            1 + 0im,
            2 + 0im,
            3 + 0im,
            4 + 0im,
            5 + 0im,
            6 + 0im,
            7 + 0im,
            8 + 0im,
            9 + 0im,
            10 + 0im,
        ]
        data.optimization_flags = fill(true, 2 * projection_layout(data).n_proj)

        MVMCOptimizers.sync_modified_parameter!(data)

        @test real.(data.doublon_holon_4site_params) ≈ [-4, -4, -2, -2, 0, 0, 2, 2, 4, 4]
        @test data.gutzwiller_terms[1].value == 12.0 + 0.0im
    end

    @testset "DH shift is disabled when any Gutzwiller real flag is fixed" begin
        data = ExpertModeData()
        data.gutzwiller_terms = [GutzwillerTerm(0, 10.0 + 0.0im, true)]
        data.doublon_holon_2site_indices = [DoublonHolon2SiteIndex([1 0; 0 1])]
        data.doublon_holon_2site_params = ComplexF64[i + 0im for i = 1:6]
        data.optimization_flags = fill(true, 2 * projection_layout(data).n_proj)
        data.optimization_flags[1] = false

        orig_dh = copy(data.doublon_holon_2site_params)
        MVMCOptimizers.sync_modified_parameter!(data)

        @test data.doublon_holon_2site_params == orig_dh
        @test data.gutzwiller_terms[1].value == 10.0 + 0.0im
    end

    @testset "DH shift is disabled when any DH real flag is fixed" begin
        data2 = ExpertModeData()
        data2.gutzwiller_terms = [GutzwillerTerm(0, 10.0 + 0.0im, true)]
        data2.doublon_holon_2site_indices = [DoublonHolon2SiteIndex([1 0; 0 1])]
        data2.doublon_holon_2site_params = ComplexF64[i + 0im for i = 1:6]
        layout2 = projection_layout(data2)
        data2.optimization_flags = fill(true, 2 * layout2.n_proj)
        data2.optimization_flags[2 * layout2.dh2_offset + 1] = false

        orig_dh2 = copy(data2.doublon_holon_2site_params)
        @test !MVMCOptimizers.flag_shift_dh2(data2, layout2)
        MVMCOptimizers.sync_modified_parameter!(data2)
        @test data2.doublon_holon_2site_params == orig_dh2
        @test data2.gutzwiller_terms[1].value == 10.0 + 0.0im

        data4 = ExpertModeData()
        data4.gutzwiller_terms = [GutzwillerTerm(0, 20.0 + 0.0im, true)]
        data4.doublon_holon_4site_indices = [DoublonHolon4SiteIndex([1 0 1 0; 0 1 0 1])]
        data4.doublon_holon_4site_params = ComplexF64[i + 0im for i = 1:10]
        layout4 = projection_layout(data4)
        data4.optimization_flags = fill(true, 2 * layout4.n_proj)
        data4.optimization_flags[2 * layout4.dh4_offset + 1] = false

        orig_dh4 = copy(data4.doublon_holon_4site_params)
        @test !MVMCOptimizers.flag_shift_dh4(data4, layout4)
        MVMCOptimizers.sync_modified_parameter!(data4)
        @test data4.doublon_holon_4site_params == orig_dh4
        @test data4.gutzwiller_terms[1].value == 20.0 + 0.0im
    end

    @testset "DH shift happens before Gutzwiller-Jastrow shift" begin
        data = ExpertModeData()
        data.gutzwiller_terms = [GutzwillerTerm(0, 10.0 + 0.0im, true)]
        data.jastrow_terms = [JastrowTerm(0, 1, 20.0 + 0.0im, true)]
        data.doublon_holon_2site_indices = [DoublonHolon2SiteIndex([1 0; 0 1])]
        data.doublon_holon_2site_params = ComplexF64[i + 0im for i = 1:6]
        data.optimization_flags = fill(true, 2 * projection_layout(data).n_proj)

        MVMCOptimizers.sync_modified_parameter!(data)

        # DH2 contributes gShift = 3 + 4 = 7, then GJ shifts average(17, 20)=18.5.
        @test data.gutzwiller_terms[1].value == -1.5 + 0.0im
        @test data.jastrow_terms[1].value == 1.5 + 0.0im
    end
end
