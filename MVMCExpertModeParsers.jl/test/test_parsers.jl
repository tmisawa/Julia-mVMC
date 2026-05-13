"""
Tests for parser functions
"""

using MVMCExpertModeParsers: ModParaParameters, TransferTerm, CoulombIntraTerm
using MVMCExpertModeParsers: GutzwillerTerm, JastrowTerm, OrbitalTerm
using MVMCExpertModeParsers: parse_modpara_content, parse_trans_content
using MVMCExpertModeParsers: parse_coulomb_intra_content, parse_gutzwiller_content
using MVMCExpertModeParsers: parse_jastrow_content, parse_orbital_content

@testset "ModPara Parser" begin
    # Test ModParaParameters creation
    params = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
    @test params.nsite == 4
    @test params.nelec == 2
    @test params.nlocspin == 0

    # Test parsing from content
    content = """
    NSite = 4
    NElec = 2
    NLocSpin = 0
    VMCCalMode = 0
    LanczosMode = 0
    """

    result = parse_modpara_content(content)
    @test result.success
    @test result.data.nsite == 4
    @test result.data.nelec == 2
    @test result.data.nlocspin == 0
end

@testset "Transfer Parser" begin
    # Test TransferTerm creation
    term = TransferTerm(0, 1, 1.0+0.0im, :up)
    @test term.site1 == 0
    @test term.site2 == 1
    @test term.value == 1.0+0.0im
    @test term.spin == :up

    # Test parsing from content
    # Format: site1 spin1 site2 spin2 real imag
    # spin: 0 = up, 1 = down
    content = """
    0 0 1 0  1.0  0.0
    1 1 2 1  2.0  0.0
    2 0 3 1  3.0  0.0
    """

    result = parse_trans_content(content)
    @test result.success
    @test length(result.data) == 3
    @test result.data[1].site1 == 0
    @test result.data[1].spin1 == 0
    @test result.data[1].site2 == 1
    @test result.data[1].spin2 == 0
    @test result.data[1].spin == :up
end

@testset "Coulomb Parser" begin
    # Test CoulombIntraTerm creation
    term = CoulombIntraTerm(0, 4.0)
    @test term.site == 0
    @test term.value == 4.0

    # Test parsing from content
    content = """
    0 4.0
    1 4.0
    """

    result = parse_coulomb_intra_content(content)
    @test result.success
    @test length(result.data) == 2
    @test result.data[1].site == 0
    @test result.data[1].value == 4.0
end

@testset "Gutzwiller Parser" begin
    # Test GutzwillerTerm creation
    term = GutzwillerTerm(0, 0.5+0.0im, false)
    @test term.site == 0
    @test term.value == 0.5+0.0im
    @test term.is_complex == false

    # Test parsing from content
    content = """
    0 0.5
    1 0.3
    """

    result = parse_gutzwiller_content(content)
    @test result.success
    @test length(result.data) == 2
    @test result.data[1].site == 0
    @test result.data[1].value == 0.5+0.0im
end

@testset "Jastrow Parser" begin
    # Test JastrowTerm creation
    term = JastrowTerm(0, 1, 0.1+0.0im, false)
    @test term.site1 == 0
    @test term.site2 == 1
    @test term.value == 0.1+0.0im
    @test term.is_complex == false

    # Test parsing from content
    content = """
    0 1 0.1
    1 2 0.1
    """

    result = parse_jastrow_content(content)
    @test result.success
    @test length(result.data) == 2
    @test result.data[1].site1 == 0
    @test result.data[1].site2 == 1
end

@testset "Orbital Parser" begin
    # Test OrbitalTerm creation
    # Constructor: OrbitalTerm(site1, site2, idx, value, is_complex, sign)
    term = OrbitalTerm(0, 1, 0, 0.5+0.0im, false, 1)
    @test term.site1 == 0
    @test term.site2 == 1
    @test term.value == 0.5+0.0im
    @test term.is_complex == false

    # Test parsing from content
    # Format: site1 site2 idx [sign]
    content = """
    0 1 0
    1 2 1
    """

    # parse_orbital_content returns (ParseResult, Dict{Int,Int})
    result, opt_flags = MVMCExpertModeParsers.parse_orbital_content(content)
    @test result.success
    @test length(result.data) == 2
    @test result.data[1].site1 == 0
    @test result.data[1].site2 == 1
    @test isempty(opt_flags)
end
