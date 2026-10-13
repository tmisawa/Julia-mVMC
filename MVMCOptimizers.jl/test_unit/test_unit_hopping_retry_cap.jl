using Test, Random
using MVMCOptimizers

# Supplied raw-word tape; running past its end is a draw-count defect.
mutable struct HoppingRetryTape <: AbstractRNG
    raw::Vector{UInt32}
    cursor::Int
end
Random.rand(rng::HoppingRetryTape, ::Type{UInt32}) = begin
    rng.cursor < length(rng.raw) || error("tape exhausted after $(rng.cursor) words")
    rng.cursor += 1
    rng.raw[rng.cursor]
end

function hopping_retry_cases(path)
    cases = Dict{String,Any}[]
    case = Dict{String,Any}()
    for line in eachline(path)
        key, rest = occursin(' ', line) ? split(line, ' '; limit = 2) : (line, "")
        ints() = parse.(Int, split(rest))
        if key == "case"
            name, kind = split(rest)
            case = Dict{String,Any}("name" => String(name), "kind" => String(kind))
        elseif key == "words"
            case["words"] = parse.(UInt32, split(rest))
        elseif key == "nsite" || key == "nelec"
            case[key] = parse(Int, rest)
        elseif key == "end"
            push!(cases, case)
        else
            case[String(key)] = ints()
        end
    end
    return cases
end

@testset "C candidate retry-cap order (mvmc-rs #559)" begin
    fixture = joinpath(@__DIR__, "fixtures", "c_hopping_retry_cap", "expected.txt")
    cases = hopping_retry_cases(fixture)
    @test length(cases) == 13
    for c in cases
        nsite, nelec = c["nsite"], c["nelec"]
        idx, cfg, spn, loc = c["ele_idx"], c["ele_cfg"], c["ele_spn"], c["loc_spn"]
        num = Int.(cfg .!= -1)
        # Metropolis marker, drawn only for a non-rejected candidate as in C.
        rng = HoppingRetryTape(vcat(c["words"], UInt32(0x7fffffff)), 0)
        kind = c["kind"]
        if kind == "hopping"
            mi, ri, rj, s, reject = MVMCOptimizers.make_candidate_hopping(
                idx, cfg, nsite, nelec, loc, rng)
            t = s
        elseif kind == "hopping_fsz" || kind == "hopping_csz"
            two_sz = kind == "hopping_fsz" ? -1 : 0
            mi, ri, rj, s, t, reject = MVMCOptimizers.make_candidate_hopping_fsz(
                idx, cfg, num, spn, loc, nsite, nelec, two_sz, rng)
        else
            @test kind == "conduction"
            mi, ri, rj, s, t, reject = MVMCOptimizers.make_candidate_local_spin_flip_conduction(
                idx, cfg, num, spn, loc, nsite, nelec, rng)
        end
        candidate_words = rng.cursor
        Bool(reject) || MVMCOptimizers.rng_real2(rng)
        got = [mi, ri, rj, s, t, Int(reject), candidate_words, rng.cursor]
        @test (c["name"], got) == (c["name"], c["result"])
    end

    # No itinerant electron: C's source loop never terminates (unsupported
    # input); Julia rejects after Nsite^2 + 1 electron draws.
    rng = HoppingRetryTape(zeros(UInt32, 2 * (2 * 2 + 1)), 0)
    out = MVMCOptimizers.make_candidate_hopping([0, 1], [0, -1, -1, 0], 2, 1, [1, 1], rng)
    @test out[5] == 1
    @test rng.cursor == 2 * (2 * 2 + 1)
end
