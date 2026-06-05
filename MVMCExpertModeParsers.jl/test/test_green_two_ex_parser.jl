using Test
using MVMCExpertModeParsers
using MVMCExpertModeParsers: GreenTwoExTerm, ExpertModeData

@testset "GreenTwoExTerm type and ExpertModeData field" begin
    # The struct stores two one-body Green specs (C reorder already absorbed).
    t = GreenTwoExTerm(0, 0, 1, 0, 2, 1, 3, 1)
    @test t.site_i1 == 0
    @test t.spin_i1 == 0
    @test t.site_j1 == 1
    @test t.spin_j1 == 0
    @test t.site_i2 == 2
    @test t.spin_i2 == 1
    @test t.site_j2 == 3
    @test t.spin_j2 == 1

    # A fresh ExpertModeData has an empty factored-term list by default.
    data = ExpertModeData()
    @test isa(data.green_two_ex_terms, Vector{GreenTwoExTerm})
    @test isempty(data.green_two_ex_terms)
end
