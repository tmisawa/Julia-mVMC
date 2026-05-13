"""
Tests for validation functions
"""

using MVMCExpertModeParsers:
    ModParaParameters, TransferTerm, CoulombIntraTerm, CoulombInterTerm
using MVMCExpertModeParsers: GutzwillerTerm, JastrowTerm, OrbitalTerm, ExpertModeData
using MVMCExpertModeParsers: validate_modpara_params, validate_transfer_terms
using MVMCExpertModeParsers: validate_coulomb_intra_terms, validate_coulomb_inter_terms
using MVMCExpertModeParsers: validate_gutzwiller_terms, validate_jastrow_terms
using MVMCExpertModeParsers: validate_orbital_terms, validate_expert_mode_data

@testset "ModPara Validation" begin
    # Test valid parameters
    valid_params = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0, ncond = -1)
    result = validate_modpara_params(valid_params)
    @test result.is_valid
    @test isempty(result.errors)

    # Test invalid parameters
    invalid_params = ModParaParameters(nsite = -1, nelec = -1, nlocspin = -1)
    result = validate_modpara_params(invalid_params)
    @test !result.is_valid
    @test !isempty(result.errors)
    @test any(contains(e, "NSite must be positive") for e in result.errors)
    @test any(contains(e, "NElec must be non-negative") for e in result.errors)
    @test any(contains(e, "NLocSpin must be non-negative") for e in result.errors)
end

@testset "Transfer Terms Validation" begin
    nsite = 4

    # Test valid terms
    valid_terms = [
        TransferTerm(0, 1, 1.0+0.0im, :up),
        TransferTerm(1, 2, 1.0+0.0im, :down),
        TransferTerm(2, 3, 1.0+0.0im, :both),
    ]
    result = validate_transfer_terms(valid_terms, nsite)
    @test result.is_valid
    @test isempty(result.errors)

    # Test invalid terms (out of range)
    invalid_terms =
        [TransferTerm(-1, 1, 1.0+0.0im, :up), TransferTerm(0, 4, 1.0+0.0im, :up)]
    result = validate_transfer_terms(invalid_terms, nsite)
    @test !result.is_valid
    @test !isempty(result.errors)
    @test any(contains(e, "out of range") for e in result.errors)
end

@testset "Coulomb Terms Validation" begin
    nsite = 4

    # Test valid Coulomb intra terms
    valid_intra_terms = [CoulombIntraTerm(0, 4.0), CoulombIntraTerm(1, 4.0)]
    result = validate_coulomb_intra_terms(valid_intra_terms, nsite)
    @test result.is_valid
    @test isempty(result.errors)

    # Test invalid Coulomb intra terms
    invalid_intra_terms = [CoulombIntraTerm(-1, 4.0), CoulombIntraTerm(4, 4.0)]
    result = validate_coulomb_intra_terms(invalid_intra_terms, nsite)
    @test !result.is_valid
    @test !isempty(result.errors)

    # Test valid Coulomb inter terms
    valid_inter_terms = [CoulombInterTerm(0, 1, 1.0), CoulombInterTerm(1, 2, 1.0)]
    result = validate_coulomb_inter_terms(valid_inter_terms, nsite)
    @test result.is_valid
    @test isempty(result.errors)

    # Test invalid Coulomb inter terms
    invalid_inter_terms = [CoulombInterTerm(-1, 1, 1.0), CoulombInterTerm(0, 4, 1.0)]
    result = validate_coulomb_inter_terms(invalid_inter_terms, nsite)
    @test !result.is_valid
    @test !isempty(result.errors)
end

@testset "Gutzwiller Terms Validation" begin
    nsite = 4

    # Test valid terms
    valid_terms = [GutzwillerTerm(0, 0.5+0.0im, false), GutzwillerTerm(1, 0.3+0.0im, false)]
    result = validate_gutzwiller_terms(valid_terms, nsite)
    @test result.is_valid
    @test isempty(result.errors)

    # Test invalid terms
    invalid_terms =
        [GutzwillerTerm(-1, 0.5+0.0im, false), GutzwillerTerm(4, 0.5+0.0im, false)]
    result = validate_gutzwiller_terms(invalid_terms, nsite)
    @test !result.is_valid
    @test !isempty(result.errors)
end

@testset "Jastrow Terms Validation" begin
    nsite = 4

    # Test valid terms
    valid_terms = [JastrowTerm(0, 1, 0.1+0.0im, false), JastrowTerm(1, 2, 0.1+0.0im, false)]
    result = validate_jastrow_terms(valid_terms, nsite)
    @test result.is_valid
    @test isempty(result.errors)

    # Test invalid terms
    invalid_terms =
        [JastrowTerm(-1, 1, 0.1+0.0im, false), JastrowTerm(0, 4, 0.1+0.0im, false)]
    result = validate_jastrow_terms(invalid_terms, nsite)
    @test !result.is_valid
    @test !isempty(result.errors)
end

@testset "Orbital Terms Validation" begin
    nsite = 4

    # Test valid terms
    valid_terms = [OrbitalTerm(0, 1, 0.5+0.0im, false), OrbitalTerm(1, 2, 0.5+0.0im, false)]
    result = validate_orbital_terms(valid_terms, nsite)
    @test result.is_valid
    @test isempty(result.errors)

    # Test invalid terms
    invalid_terms =
        [OrbitalTerm(-1, 1, 0.5+0.0im, false), OrbitalTerm(0, 4, 0.5+0.0im, false)]
    result = validate_orbital_terms(invalid_terms, nsite)
    @test !result.is_valid
    @test !isempty(result.errors)
end

@testset "Expert Mode Data Validation" begin
    # Test valid data
    data = ExpertModeData()
    data.modpara = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
    push!(data.transfer_terms, TransferTerm(0, 1, 1.0+0.0im, :up))
    push!(data.coulomb_intra_terms, CoulombIntraTerm(0, 4.0))
    push!(data.gutzwiller_terms, GutzwillerTerm(0, 0.5+0.0im, false))

    result = validate_expert_mode_data(data)
    @test result.is_valid
    @test isempty(result.errors)

    # Test invalid data
    invalid_data = ExpertModeData()
    invalid_data.modpara = ModParaParameters(nsite = -1, nelec = -1, nlocspin = -1)
    push!(invalid_data.transfer_terms, TransferTerm(-1, 1, 1.0+0.0im, :up))

    result = validate_expert_mode_data(invalid_data)
    @test !result.is_valid
    @test !isempty(result.errors)
end
