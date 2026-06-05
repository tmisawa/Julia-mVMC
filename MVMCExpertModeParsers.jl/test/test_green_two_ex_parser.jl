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

using MVMCExpertModeParsers: parse_green_two_ex_content

# Standard 5-line .def header followed by data rows.
function _greentwoex(count::Int, rows::Vector{String})
    header = join([
        "=============================================",
        "NCisAjsCktAlt   $count",
        "=============================================",
        "======== Factored two-body Green ============",
        "=============================================",
    ], "\n")
    return header * "\n" * join(rows, "\n") * "\n"
end

# Build content with an arbitrary line-2 (count) string, for header tests.
function _greentwoex_raw(line2::String, rows::Vector{String})
    header = join([
        "=============================================",
        line2,
        "=============================================",
        "======== Factored two-body Green ============",
        "=============================================",
    ], "\n")
    return header * "\n" * join(rows, "\n") * "\n"
end

@testset "parse_green_two_ex_content" begin
    @testset "valid content maps columns with C reorder" begin
        # Row x0..x7 = 0 0 1 0  2 1 3 1
        #   first  Green ⟨c†_{0,0} c_{1,0}⟩ -> (0,0,1,0)
        #   second Green ⟨c†_{3,1} c_{2,1}⟩ -> (3,1,2,1)   (x6,x7,x4,x5)
        content = _greentwoex(1, ["    0     0     1     0     2     1     3     1"])
        res = parse_green_two_ex_content(content)
        @test res.success
        @test length(res.data) == 1
        t = res.data[1]
        @test (t.site_i1, t.spin_i1, t.site_j1, t.spin_j1) == (0, 0, 1, 0)
        @test (t.site_i2, t.spin_i2, t.site_j2, t.spin_j2) == (3, 1, 2, 1)
    end

    @testset "header count must match parsed rows" begin
        content = _greentwoex(3, [
            "0 0 0 0 0 0 0 0",
            "0 0 0 0 0 1 0 1",
        ])  # header says 3, only 2 rows
        res = parse_green_two_ex_content(content)
        @test !res.success
        @test occursin("count", lowercase(res.error_message))
    end

    @testset "a data row that is not exactly 8 integers is rejected" begin
        content = _greentwoex(1, ["0 0 0 0 0 0 0"])  # 7 fields
        res = parse_green_two_ex_content(content)
        @test !res.success
        @test occursin("8", res.error_message)
    end

    @testset "non-integer field is rejected" begin
        content = _greentwoex(1, ["0 0 0 0 0 0 0 x"])
        res = parse_green_two_ex_content(content)
        @test !res.success
    end

    @testset "spin outside {0,1} is rejected" begin
        content = _greentwoex(1, ["0 2 0 0 0 0 0 0"])
        res = parse_green_two_ex_content(content)
        @test !res.success
    end

    @testset "header count must be a non-negative integer in token 2" begin
        # Missing header entirely.
        @test !parse_green_two_ex_content("no header here\njust one line\n").success
        # Keyword present but the count token is not an integer.
        @test !parse_green_two_ex_content(
            _greentwoex_raw("NCisAjsCktAlt x", ["0 0 0 0 0 0 0 0"])).success
        # First token is an integer but token 2 (the count) is not — guards
        # against the permissive "first integer token" shortcut (Finding 2).
        @test !parse_green_two_ex_content(
            _greentwoex_raw("0 bogus", ["0 0 0 0 0 0 0 0"])).success
        # Negative count.
        @test !parse_green_two_ex_content(
            _greentwoex_raw("NCisAjsCktAlt -1", ["0 0 0 0 0 0 0 0"])).success
    end
end

using MVMCExpertModeParsers: parse_file_by_type!

@testset "TwoBodyGEx dispatch" begin
    mktempdir() do dir
        good = joinpath(dir, "greentwoex.def")
        write(good, _greentwoex(1, ["0 0 1 0 2 1 3 1"]))
        data = ExpertModeData()
        parse_file_by_type!(data, "TwoBodyGEx", good)
        @test length(data.green_two_ex_terms) == 1
        @test data.green_two_ex_terms[1].site_i2 == 3  # reorder preserved through dispatch

        # A malformed factored file must be fatal, not silently empty.
        bad = joinpath(dir, "bad_greentwoex.def")
        write(bad, _greentwoex(2, ["0 0 1 0 2 1 3 1"]))  # header 2, one row
        data2 = ExpertModeData()
        @test_throws ErrorException parse_file_by_type!(data2, "TwoBodyGEx", bad)
    end
end

# parse_expert_mode_files is exported; `using MVMCExpertModeParsers` (top of file)
# brings it into scope.
@testset "TwoBodyGEx is fatal on the public parse_expert_mode_files path" begin
    # Valid: factored terms reach the public ExpertModeData.
    mktempdir() do dir
        write(joinpath(dir, "greentwoex.def"), _greentwoex(1, ["0 0 1 0 2 1 3 1"]))
        write(joinpath(dir, "namelist.def"), "      TwoBodyGEx  greentwoex.def\n")
        d = parse_expert_mode_files(joinpath(dir, "namelist.def"))
        @test length(d.green_two_ex_terms) == 1
    end
    # Malformed greentwoex.def is fatal (not a swallowed warning).
    mktempdir() do dir
        write(joinpath(dir, "greentwoex.def"), _greentwoex(2, ["0 0 1 0 2 1 3 1"]))
        write(joinpath(dir, "namelist.def"), "      TwoBodyGEx  greentwoex.def\n")
        @test_throws ErrorException parse_expert_mode_files(joinpath(dir, "namelist.def"))
    end
    # Missing greentwoex.def is fatal.
    mktempdir() do dir
        write(joinpath(dir, "namelist.def"), "      TwoBodyGEx  greentwoex.def\n")
        @test_throws ErrorException parse_expert_mode_files(joinpath(dir, "namelist.def"))
    end
end
