"""
Tests for parser functions
"""

using MVMCExpertModeParsers: ModParaParameters, TransferTerm, CoulombIntraTerm
using MVMCExpertModeParsers: GutzwillerTerm, JastrowTerm, OrbitalTerm, PairHopTerm
using MVMCExpertModeParsers: parse_modpara_content, parse_trans_content
using MVMCExpertModeParsers: parse_coulomb_intra_content, parse_gutzwiller_content
using MVMCExpertModeParsers: parse_jastrow_content, parse_orbital_content
using MVMCExpertModeParsers: parse_pairhop_content

@testset "ModPara Parser" begin
    # Test ModParaParameters creation
    params = ModParaParameters(nsite = 4, nelec = 2, nlocspin = 0)
    @test params.nsite == 4
    @test params.nelec == 2
    @test params.nlocspin == 0

    defaults = ModParaParameters()
    @test defaults.dsr_opt_cg_tol == 1e-10
    @test defaults.nsr_opt_cg_max_iter == 0
    @test defaults.nsrcg == 0
    @test defaults.nstore_o == 1
    @test defaults.use_diag_scale == 0
    @test defaults.rescale_smat == 0

    # Test parsing from content
    content = """
    NSite = 4
    NElec = 2
    NLocSpin = 0
    VMCCalMode = 0
    LanczosMode = 0
    NStore = 0
    NSRCG = 1
    useDiagScale = 1
    RescaleSmat = 1
    """

    result = parse_modpara_content(content)
    @test result.success
    @test result.data.nsite == 4
    @test result.data.nelec == 2
    @test result.data.nlocspin == 0
    @test result.data.dsr_opt_cg_tol == 1e-10
    @test result.data.nsr_opt_cg_max_iter == 0
    @test result.data.nstore_o == 0
    @test result.data.nsrcg == 1
    @test result.data.use_diag_scale == 1
    @test result.data.rescale_smat == 1
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

@testset "PairHop Parser" begin
    term = PairHopTerm(0, 1, 0.25)
    @test term.site1 == 0
    @test term.site2 == 1
    @test term.value == 0.25

    result = parse_pairhop_content("""
    0 1 0.25
    2 3 -0.5
    """)

    @test result.success
    @test length(result.data) == 4
    @test [(t.site1, t.site2, t.value) for t in result.data] == [
        (0, 1, 0.25),
        (1, 0, 0.25),
        (2, 3, -0.5),
        (3, 2, -0.5),
    ]

    raw_long = parse_pairhop_content(join(["$i $(i + 1) 0.1" for i = 0:5], "\n"))
    @test raw_long.success
    @test length(raw_long.data) == 12

    with_header = parse_pairhop_content("""
    =============================================
    NPairHopp          1
    =============================================
    ====== Pair-Hopping term ============
    =============================================
        4     5         1.250000000000000
    """)

    @test with_header.success
    @test [(t.site1, t.site2, t.value) for t in with_header.data] == [
        (4, 5, 1.25),
        (5, 4, 1.25),
    ]

    self_term = parse_pairhop_content("2 2 0.75")
    @test self_term.success
    @test [(t.site1, t.site2, t.value) for t in self_term.data] == [
        (2, 2, 0.75),
        (2, 2, 0.75),
    ]
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
