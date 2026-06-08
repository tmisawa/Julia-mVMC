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

@testset "factored accumulation: w * local[idx0] * conj(local[idx1])" begin
    pq = MVMCOptimizers.PhysicalQuantities(2, 1, 0)
    pq.cis_ajs_ckt_alt_idx = [(1, 2)]
    pq.local_cis_ajs[1] = 2.0 + 1.0im
    pq.local_cis_ajs[2] = 3.0 - 4.0im
    MVMCOptimizers.accumulate_factored_green!(pq, 0.5)
    # 0.5 * (2+1im) * conj(3-4im) = 0.5 * (2+1im) * (3+4im) = 0.5 * (2+11im)
    @test pq.phys_cis_ajs_ckt_alt[1] ≈ 0.5 * ((2.0 + 1.0im) * conj(3.0 - 4.0im))
    @test pq.phys_cis_ajs_ckt_alt[1] ≈ (1.0 + 5.5im)

    # Accumulates (adds), not overwrites.
    MVMCOptimizers.accumulate_factored_green!(pq, 0.5)
    @test pq.phys_cis_ajs_ckt_alt[1] ≈ (2.0 + 11.0im)
end

@testset "calculate_green_func! legacy no-acc path writes to state phys" begin
    data = ExpertModeData()
    data.modpara.nsite = 2
    data.modpara.nelec = 1
    data.green_one_terms = [
        GreenOneTerm(0, 0, :up, :up),
        GreenOneTerm(1, 1, :down, :down),
    ]
    data.green_two_ex_terms = [
        GreenTwoExTerm(0, 0, 0, 0, 1, 1, 1, 1),
    ]

    state = MVMCOptimizers.VMCOptimizationState(2, 1, 0, 0, 1, 1, true, false)
    MVMCOptimizers.initialize_phys_quantities!(state, data)
    phys = state.phys_quantities

    ele_idx = Int[0, 1]
    ele_cfg = Int[0, -1, -1, 1]
    ele_num = Int[1, 0, 0, 1]

    MVMCOptimizers.calculate_green_func!(
        data,
        state,
        0.5,
        1.0 + 0.0im,
        ele_idx,
        ele_cfg,
        ele_num,
        Int[],
    )

    @test phys.local_cis_ajs == [1.0 + 0.0im, 1.0 + 0.0im]
    @test phys.phys_cis_ajs == [0.5 + 0.0im, 0.5 + 0.0im]
    @test phys.phys_cis_ajs_ckt_alt == [0.5 + 0.0im]
end

using MVMCOptimizers: VMCOptimizationState
using Printf

@testset "output: canonical cisajs + factored ex, with output_dir and numbering" begin
    data = ExpertModeData()
    data.modpara.nsite = 4
    data.modpara.c_data_file_head = "zvo"
    state = VMCOptimizationState(4, 2, 0, 0, 1, 4, true, false)

    pq = MVMCOptimizers.PhysicalQuantities(1, 1, 0)
    pq.cis_ajs_idx = NTuple{4,Int}[(0, 0, 1, 0)]
    pq.cis_ajs_ckt_alt_idx = [(1, 1)]
    pq.phys_cis_ajs[1] = 1.25 + 0.0im
    pq.phys_cis_ajs_ckt_alt[1] = 2.5 - 1.0im
    state.phys_quantities = pq

    mktempdir() do dir
        MVMCOptimizers.output_green_func!(data, state, 1; output_dir = dir)

        cisajs = joinpath(dir, "zvo_cisajs_001.dat")
        ex = joinpath(dir, "zvo_cisajscktaltex_001.dat")
        @test isfile(cisajs)
        @test isfile(ex)

        line = first(filter(!isempty, split(read(cisajs, String), "\n")))
        cols = split(line)
        @test cols[1] == "0" && cols[2] == "0" && cols[3] == "1" && cols[4] == "0"

        exline = first(filter(!isempty, split(read(ex, String), "\n")))
        vals = parse.(Float64, split(exline))
        @test length(vals) == 2
        @test vals[1] ≈ 2.5
        @test vals[2] ≈ -1.0
    end
end

@testset "PhysCal output file index uses NDataIdxStart" begin
    data = ExpertModeData()
    data.modpara.n_data_idx_start = 1
    @test MVMCOptimizers.physcal_output_file_index(data, 0) == 1
    @test MVMCOptimizers.physcal_output_file_index(data, 3) == 4

    data.modpara.n_data_idx_start = 7
    @test MVMCOptimizers.physcal_output_file_index(data, 0) == 7
    @test MVMCOptimizers.physcal_output_file_index(data, 2) == 9
end

@testset "output_data_phys!: out/var truncate-on-first-sample, Green files use NDataIdxStart" begin
    # Regression guard for the fmt-1 fix: output_data_phys! must drive the
    # energy/param write mode from the 0-based sample index (so the first sample
    # truncates and a re-run does not accumulate stale lines), while numbering the
    # Green files with ismp + NDataIdxStart. Previously the C-visible file index
    # (>= 1 for NDataIdxStart >= 1) was forwarded as the write-mode selector, so
    # the first sample appended instead of truncating.
    data = ExpertModeData()
    data.modpara.nsite = 2
    data.modpara.c_data_file_head = "zvo"
    data.modpara.n_data_idx_start = 1

    state = VMCOptimizationState(2, 2, 0, 0, 1, 2, true, false)
    pq = MVMCOptimizers.PhysicalQuantities(1, 0, 0)
    pq.cis_ajs_idx = NTuple{4,Int}[(0, 0, 1, 0)]
    state.phys_quantities = pq

    mktempdir() do dir
        outpath = joinpath(dir, "zvo_out.dat")
        varpath = joinpath(dir, "zvo_var.dat")

        # One run of two samples: ismp = 0 (truncate) then 1 (append).
        MVMCOptimizers.output_data_phys!(data, state, 0; output_dir = dir)
        MVMCOptimizers.output_data_phys!(data, state, 1; output_dir = dir)
        @test count(==('\n'), read(outpath, String)) == 2  # exactly 2 lines, no stale leading line

        # Green files numbered ismp + NDataIdxStart (= 1, 2), not 0-based.
        @test isfile(joinpath(dir, "zvo_cisajs_001.dat"))
        @test isfile(joinpath(dir, "zvo_cisajs_002.dat"))
        @test !isfile(joinpath(dir, "zvo_cisajs_000.dat"))

        # A fresh run (ismp = 0 again) re-truncates out/var — no cross-run pollution.
        MVMCOptimizers.output_data_phys!(data, state, 0; output_dir = dir)
        @test count(==('\n'), read(outpath, String)) == 1
        @test isfile(varpath)
    end
end

@testset "no TwoBodyGEx preserves greenone order and duplicates in output" begin
    data = ExpertModeData()
    data.modpara.nsite = 2
    data.modpara.c_data_file_head = "zvo"
    data.green_one_terms = [
        GreenOneTerm(0, 1, :up, :up),
        GreenOneTerm(0, 1, :up, :up),
        GreenOneTerm(1, 0, :down, :down),
    ]

    state = VMCOptimizationState(2, 2, 0, 0, 1, 2, true, false)
    MVMCOptimizers.initialize_phys_quantities!(state, data)
    pq = state.phys_quantities

    @test pq.cis_ajs_idx == NTuple{4,Int}[
        (0, 0, 1, 0),
        (0, 0, 1, 0),
        (1, 1, 0, 1),
    ]

    pq.phys_cis_ajs .= ComplexF64[1.0 + 0im, 2.0 + 0im, 3.0 + 0im]
    mktempdir() do dir
        MVMCOptimizers.output_green_func!(data, state, 1; output_dir = dir)
        lines = filter(!isempty, split(read(joinpath(dir, "zvo_cisajs_001.dat"), String), "\n"))
        @test length(lines) == 3
        @test split(lines[1])[1:4] == ["0", "0", "1", "0"]
        @test split(lines[2])[1:4] == ["0", "0", "1", "0"]
        @test split(lines[3])[1:4] == ["1", "1", "0", "1"]
    end
end

@testset "initialize_phys_quantities! wires the factored canonical list and pairs" begin
    # Integration of build_canonical_cis_ajs_idx + resolve_cis_ajs_ckt_alt_idx
    # through initialize_phys_quantities! (the helpers are unit-tested above, but
    # the wiring that stores them on PhysicalQuantities and sizes the buffers was
    # previously only exercised with an empty factored list).

    # Single factored term: explicit greenone first, one appended constituent.
    data = ExpertModeData()
    data.modpara.nsite = 4
    data.green_one_terms = [GreenOneTerm(0, 1, :up, :up)]            # (0,0,1,0)
    data.green_two_ex_terms = [GreenTwoExTerm(0, 0, 1, 0, 2, 1, 3, 1)]  # c1=(0,0,1,0); c2=(2,1,3,1) new

    state = VMCOptimizationState(4, 2, 0, 0, 1, 4, true, false)
    MVMCOptimizers.initialize_phys_quantities!(state, data)
    pq = state.phys_quantities

    @test pq.cis_ajs_idx == NTuple{4,Int}[(0, 0, 1, 0), (2, 1, 3, 1)]
    @test pq.cis_ajs_ckt_alt_idx == [(1, 2)]
    @test length(pq.phys_cis_ajs) == 2          # canonical includes the appended constituent
    @test length(pq.local_cis_ajs) == 2         # local buffer sized for accumulate_factored_green!
    @test length(pq.phys_cis_ajs_ckt_alt) == 1  # == number of factored terms
    @test length(pq.phys_cis_ajs_ckt_alt_dc) == 0

    # Two factored terms sharing constituents (swapped): exercises cross-term
    # dedup (canonical stays length 2) and ordered multi-pair resolution.
    data2 = ExpertModeData()
    data2.modpara.nsite = 4
    data2.green_one_terms = [GreenOneTerm(0, 1, :up, :up)]           # (0,0,1,0)
    data2.green_two_ex_terms = [
        GreenTwoExTerm(0, 0, 1, 0, 2, 1, 3, 1),   # c1=(0,0,1,0)=idx1; c2=(2,1,3,1)=idx2
        GreenTwoExTerm(2, 1, 3, 1, 0, 0, 1, 0),   # c1=(2,1,3,1)=idx2; c2=(0,0,1,0)=idx1
    ]

    state2 = VMCOptimizationState(4, 2, 0, 0, 1, 4, true, false)
    MVMCOptimizers.initialize_phys_quantities!(state2, data2)
    pq2 = state2.phys_quantities

    @test pq2.cis_ajs_idx == NTuple{4,Int}[(0, 0, 1, 0), (2, 1, 3, 1)]  # deduped, length 2
    @test pq2.cis_ajs_ckt_alt_idx == [(1, 2), (2, 1)]                   # ordered, shared indices
    @test length(pq2.phys_cis_ajs) == 2
    @test length(pq2.phys_cis_ajs_ckt_alt) == 2
end

@testset "FSZ + factored is rejected before sampling" begin
    data = ExpertModeData()
    data.green_two_ex_terms = [GreenTwoExTerm(0, 0, 1, 0, 2, 1, 3, 1)]

    # sz-conserved (i_flg_orbital_general == 0): allowed.
    data.i_flg_orbital_general = 0
    @test MVMCOptimizers.validate_factored_green_supported(data) === nothing

    # FSZ / general-orbital (== 1): rejected with a clear message.
    data.i_flg_orbital_general = 1
    threw = false
    msg = ""
    try
        MVMCOptimizers.validate_factored_green_supported(data)
    catch err
        threw = true
        msg = sprint(showerror, err)
    end
    @test threw
    @test occursin("TwoBodyGEx", msg)
    @test occursin("FSZ", msg)

    # No factored terms: FSZ is fine for the rest of PhysCal.
    data.green_two_ex_terms = GreenTwoExTerm[]
    @test MVMCOptimizers.validate_factored_green_supported(data) === nothing
end

@testset "vmc_phys_cal! rejects FSZ + factored before sampling" begin
    # The guard runs before RNG / init_parameter!, so a minimal (physically
    # incomplete) input still fails for the intended reason (Finding 2).
    data = ExpertModeData()
    data.green_two_ex_terms = [GreenTwoExTerm(0, 0, 1, 0, 2, 1, 3, 1)]
    data.i_flg_orbital_general = 1

    err = try
        MVMCOptimizers.vmc_phys_cal!(data)
        nothing
    catch e
        e
    end

    @test err isa ErrorException
    msg = sprint(showerror, err)
    @test occursin("TwoBodyGEx", msg)
    @test occursin("FSZ", msg)
end
