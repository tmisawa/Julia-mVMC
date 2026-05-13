"""
Tests for integration functions
"""

using MVMCExpertModeParsers: ExpertModeData, ModParaParameters, TransferTerm
using MVMCExpertModeParsers: CoulombIntraTerm, GutzwillerTerm, parse_namelist_content

@testset "Expert Mode Data Integration" begin
    # Test ExpertModeData creation
    data = ExpertModeData()
    @test data.modpara.nsite == 0
    @test isempty(data.transfer_terms)
    @test isempty(data.coulomb_intra_terms)
    @test isempty(data.gutzwiller_terms)

    # Test data population
    data.modpara = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
    push!(data.transfer_terms, TransferTerm(0, 1, 1.0+0.0im, :up))
    push!(data.coulomb_intra_terms, CoulombIntraTerm(0, 4.0))
    push!(data.gutzwiller_terms, GutzwillerTerm(0, 0.5+0.0im, false))

    @test data.modpara.nsite == 4
    @test length(data.transfer_terms) == 1
    @test length(data.coulomb_intra_terms) == 1
    @test length(data.gutzwiller_terms) == 1
end

@testset "Namelist Parsing" begin
    # Test namelist content parsing
    namelist_content = """
    ModPara modpara.def
    Trans trans.def
    CoulombIntra coulombintra.def
    Gutzwiller gutzwiller.def
    """

    file_list = parse_namelist_content(namelist_content)
    # file_list is now Vector{Tuple{String, String}}
    file_names = [fname for (_, fname) in file_list]
    @test "modpara.def" in file_names
    @test "trans.def" in file_names
    @test "coulombintra.def" in file_names
    @test "gutzwiller.def" in file_names
end
