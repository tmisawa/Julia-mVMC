using Test
using Random: MersenneTwister

using MVMCExpertModeParsers
using MVMCExpertModeParsers:
    ExpertModeData,
    ModParaParameters,
    GutzwillerTerm,
    JastrowTerm,
    ChargeRBMPhysLayerTerm,
    ChargeRBMHiddenLayerTerm,
    ChargeRBMPhysHiddenTerm

section_nparam(terms) = isempty(terms) ? 0 : (maximum(t.idx for t in terms) + 1)

function set_rbm_opt_flags!(flags::Vector{Bool}, n_proj::Int, rbm_opts::AbstractVector{Bool})
    for (i, is_opt) in enumerate(rbm_opts)
        i0 = i - 1
        opt_flag_idx = 2 * i0 + 2 * n_proj + 1
        flags[opt_flag_idx] = is_opt
    end
end

function make_data_for_rbm_init_tests(; modpara_complex_flag::Int = 0, gutz_is_complex::Bool = false)
    data = ExpertModeData()
    data.modpara = ModParaParameters(complex_flag = modpara_complex_flag, nneuron = 1)

    # NProj = 2 (also used to test OptFlag offset)
    data.gutzwiller_terms = [GutzwillerTerm(0, 1.0 + 1.0im, gutz_is_complex)]
    data.jastrow_terms = [JastrowTerm(0, 1, 2.0 - 3.0im, false)]

    # RBM sections: charge phys (2) + charge hidden (3) + charge phys-hidden (1) = 6 params
    data.charge_rbm_phys_layer_terms = [
        ChargeRBMPhysLayerTerm(0, 9.0 + 9.0im, false, 0),
        ChargeRBMPhysLayerTerm(0, 8.0 + 8.0im, false, 1),
    ]
    data.charge_rbm_hidden_layer_terms = [
        ChargeRBMHiddenLayerTerm(0, 7.0 + 7.0im, false, 0),
        ChargeRBMHiddenLayerTerm(0, 6.0 + 6.0im, false, 1),
        ChargeRBMHiddenLayerTerm(0, 5.0 + 5.0im, false, 2),
    ]
    data.charge_rbm_phys_hidden_terms = [
        ChargeRBMPhysHiddenTerm(0, 0, 4.0 + 4.0im, false, 0),
    ]

    n_proj = length(data.gutzwiller_terms) + length(data.jastrow_terms)
    rbm_sections = (
        data.charge_rbm_phys_layer_terms,
        data.spin_rbm_phys_layer_terms,
        data.general_rbm_phys_layer_terms,
        data.charge_rbm_hidden_layer_terms,
        data.spin_rbm_hidden_layer_terms,
        data.general_rbm_hidden_layer_terms,
        data.charge_rbm_phys_hidden_terms,
        data.spin_rbm_phys_hidden_terms,
        data.general_rbm_phys_hidden_terms,
    )
    n_rbm = sum(section_nparam(s) for s in rbm_sections)

    # Allocate optimization flags large enough for RBM OptFlag indexing.
    data.optimization_flags = falses(2 * (n_proj + n_rbm))
    return data, rbm_sections, n_proj, n_rbm
end

function expected_rbm_values_real(rng, n_rbm::Int, n_proj::Int, opt::Vector{Bool}; nneuron_divisor::Float64 = 1.0)
    vals = fill(ComplexF64(0.0), n_rbm)
    for i0 in 0:(n_rbm - 1)
        opt_flag_idx = 2 * i0 + 2 * n_proj + 1
        should_optimize = (opt_flag_idx <= length(opt) && opt[opt_flag_idx])
        if should_optimize
            vals[i0 + 1] = ComplexF64(0.01 * (rand(rng) - 0.5) / nneuron_divisor)
        end
    end
    return vals
end

function expected_rbm_values_complex(rng, n_rbm::Int, n_proj::Int, opt::Vector{Bool})
    vals = fill(ComplexF64(0.0), n_rbm)
    for i0 in 0:(n_rbm - 1)
        opt_flag_idx = 2 * i0 + 2 * n_proj + 1
        should_optimize = (opt_flag_idx <= length(opt) && opt[opt_flag_idx])
        if should_optimize
            r1 = rand(rng)
            r2 = rand(rng)
            vals[i0 + 1] = ComplexF64(1e-2 * r1 * exp(2.0im * π * r2))
        end
    end
    return vals
end

@testset "parameter_init: complex flag + RBM initialization" begin
    @testset "real RBM init uses NProj offset and OptFlag" begin
        data, rbm_sections, n_proj, n_rbm =
            make_data_for_rbm_init_tests(modpara_complex_flag = 0, gutz_is_complex = false)
        @test n_proj == 2
        @test n_rbm == 6

        # RBM opt flags for global indices i0=0..5
        rbm_opts = Bool[true, false, true, false, true, true]
        set_rbm_opt_flags!(data.optimization_flags, n_proj, rbm_opts)

        rng_init = MersenneTwister(1234)
        rng_exp = copy(rng_init)
        MVMCExpertModeParsers.init_parameter!(data; rng = rng_init)

        # Proj params are always zeroed.
        @test all(t -> t.value == 0.0 + 0.0im, data.gutzwiller_terms)
        @test all(t -> t.value == 0.0 + 0.0im, data.jastrow_terms)

        expected = expected_rbm_values_real(rng_exp, n_rbm, n_proj, data.optimization_flags; nneuron_divisor = 1.0)

        # Check scatter by section offsets.
        section_offset = 0
        for terms in rbm_sections
            n_section = section_nparam(terms)
            for term in terms
                @test term.value == expected[section_offset + term.idx + 1]
                @test imag(term.value) == 0.0
            end
            section_offset += n_section
        end

        # OptFlag=false entries should remain exactly zero.
        @test expected[2] == 0.0 + 0.0im
        @test expected[4] == 0.0 + 0.0im
    end

    @testset "complex RBM init is enabled by ModPara.complex_flag" begin
        data, rbm_sections, n_proj, n_rbm =
            make_data_for_rbm_init_tests(modpara_complex_flag = 1, gutz_is_complex = false)
        rbm_opts = trues(n_rbm)
        set_rbm_opt_flags!(data.optimization_flags, n_proj, rbm_opts)

        rng_init = MersenneTwister(2025)
        rng_exp = copy(rng_init)
        MVMCExpertModeParsers.init_parameter!(data; rng = rng_init)

        expected = expected_rbm_values_complex(rng_exp, n_rbm, n_proj, data.optimization_flags)

        assigned = ComplexF64[]
        section_offset = 0
        for terms in rbm_sections
            n_section = section_nparam(terms)
            for term in terms
                push!(assigned, term.value)
                @test term.value == expected[section_offset + term.idx + 1]
                @test abs(term.value) <= 1e-2 + 1e-14
            end
            section_offset += n_section
        end
        @test any(v -> imag(v) != 0.0, assigned)
    end

    @testset "complex RBM init is enabled by complex gutzwiller flag (AllComplexFlag)" begin
        data, rbm_sections, n_proj, n_rbm =
            make_data_for_rbm_init_tests(modpara_complex_flag = 0, gutz_is_complex = true)
        rbm_opts = trues(n_rbm)
        set_rbm_opt_flags!(data.optimization_flags, n_proj, rbm_opts)

        rng_init = MersenneTwister(7)
        rng_exp = copy(rng_init)
        MVMCExpertModeParsers.init_parameter!(data; rng = rng_init)

        expected = expected_rbm_values_complex(rng_exp, n_rbm, n_proj, data.optimization_flags)

        section_offset = 0
        for terms in rbm_sections
            n_section = section_nparam(terms)
            for term in terms
                @test term.value == expected[section_offset + term.idx + 1]
            end
            section_offset += n_section
        end
    end
end
