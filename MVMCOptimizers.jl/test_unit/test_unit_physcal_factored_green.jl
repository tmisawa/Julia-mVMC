using Test
using MVMCOptimizers
using MVMCExpertModeParsers: ExpertModeData, GreenOneTerm, GreenTwoExTerm

@testset "PhysicalQuantities index fields" begin
    pq = MVMCOptimizers.PhysicalQuantities(2, 1, 3)
    @test pq.cis_ajs_idx == NTuple{4,Int}[]
    @test pq.cis_ajs_ckt_alt_idx == Tuple{Int,Int}[]
    @test length(pq.phys_cis_ajs) == 2
    @test length(pq.phys_cis_ajs_ckt_alt) == 1
    @test length(pq.phys_cis_ajs_ckt_alt_dc) == 3
end
