using Test
using MVMCExpertModeParsers

@testset "contract/trans.def: spin1/spin2 are preserved" begin
    # Format per line: site1 spin1 site2 spin2 real imag
    content = """
    # same-spin hopping (up -> up)
    0 0 1 0  1.25  -0.5
    # spin-flip term (up -> down): should not be rejected; spin1/spin2 must be preserved
    2 0 3 1  -0.75  0.125
    """

    result = MVMCExpertModeParsers.parse_trans_content(content)
    @test result.success
    @test result.data !== nothing

    terms = result.data
    @test length(terms) == 2

    t1 = terms[1]
    @test t1.site1 == 0
    @test t1.spin1 == 0
    @test t1.site2 == 1
    @test t1.spin2 == 0
    @test t1.value == ComplexF64(1.25, -0.5)
    @test t1.spin == :up

    t2 = terms[2]
    @test t2.site1 == 2
    @test t2.spin1 == 0
    @test t2.site2 == 3
    @test t2.spin2 == 1
    @test t2.value == ComplexF64(-0.75, 0.125)
    # `spin` is derived from spin1 by TransferTerm constructor.
    @test t2.spin == :up
end

