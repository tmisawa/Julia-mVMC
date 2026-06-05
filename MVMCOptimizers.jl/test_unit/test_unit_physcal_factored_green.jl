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

@testset "canonical one-body list" begin
    # No TwoBodyGEx: canonical == greenone.def order, no dedup.
    g1 = [GreenOneTerm(0, 1, :up, :up), GreenOneTerm(1, 0, :down, :down)]
    canon = MVMCOptimizers.build_canonical_cis_ajs_idx(g1, GreenTwoExTerm[], 2)
    @test canon == NTuple{4,Int}[(0, 0, 1, 0), (1, 1, 0, 1)]

    # TwoBodyGEx present: explicit terms first, then appended factored
    # constituents in file order, de-duplicated. GreenTwoExTerm stores the two
    # one-body Greens directly (reorder already absorbed at parse time): here
    # second = <c†_{2,1} c_{3,1}> = key (2,1,3,1), which is NOT in greenone.def
    # and must be appended.
    g1b = [GreenOneTerm(0, 1, :up, :up)]                  # (0,0,1,0)
    ex = [GreenTwoExTerm(0, 0, 1, 0, 2, 1, 3, 1)]         # A=(0,0,1,0) present; B=(2,1,3,1) new
    canon2 = MVMCOptimizers.build_canonical_cis_ajs_idx(g1b, ex, 4)
    @test canon2 == NTuple{4,Int}[(0, 0, 1, 0), (2, 1, 3, 1)]

    # Out-of-range site is rejected.
    @test_throws ErrorException MVMCOptimizers.build_canonical_cis_ajs_idx(
        [GreenOneTerm(0, 5, :up, :up)], GreenTwoExTerm[], 2)

    # An invalid spin symbol is rejected, not silently treated as down (Finding 4).
    @test_throws ErrorException MVMCOptimizers.build_canonical_cis_ajs_idx(
        [GreenOneTerm(0, 1, :both, :up)], GreenTwoExTerm[], 2)
end

@testset "factored index resolution is 1-based" begin
    canon = NTuple{4,Int}[(0, 0, 1, 0), (2, 1, 3, 1)]
    ex = [GreenTwoExTerm(0, 0, 1, 0, 2, 1, 3, 1)]  # A=(0,0,1,0)->canon[1], B=(2,1,3,1)->canon[2]
    pairs = MVMCOptimizers.resolve_cis_ajs_ckt_alt_idx(canon, ex)
    @test pairs == [(1, 2)]

    # C index 0 (first one-body Green) must resolve to Julia index 1.
    canon2 = NTuple{4,Int}[(0, 0, 0, 0)]
    ex2 = [GreenTwoExTerm(0, 0, 0, 0, 0, 0, 0, 0)]  # both constituents = canon[1]
    @test MVMCOptimizers.resolve_cis_ajs_ckt_alt_idx(canon2, ex2) == [(1, 1)]
end
