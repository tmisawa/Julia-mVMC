"""
Tests for C-compatible DH2/DH4 index-table parsing and projection layout.
"""

using Test
using MVMCExpertModeParsers
using MVMCExpertModeParsers: ExpertModeData, GutzwillerTerm, JastrowTerm
using MVMCExpertModeParsers: DoublonHolon2SiteTerm, DoublonHolon4SiteTerm
using MVMCExpertModeParsers: DoublonHolon2SiteIndex, DoublonHolon4SiteIndex
using MVMCExpertModeParsers: parse_doublon_holon_2site_content
using MVMCExpertModeParsers: parse_doublon_holon_4site_content
using MVMCExpertModeParsers: parse_file_by_type!, parse_expert_mode_files
using MVMCExpertModeParsers: projection_layout, projection_parameters, set_dh_opt_flags!
using MVMCExpertModeParsers: validate_doublon_holon_2site_terms
using MVMCExpertModeParsers: validate_doublon_holon_4site_terms

const DH_TEST_SEP = "============================================="

function _dh2_content(;
    count::Int = 1,
    complex_type::Int = 1,
    main_rows::Vector{String} = ["0 1 2 0", "1 2 0 0", "2 0 1 0"],
    opt_flags::Vector{Int} = [1, 0, 1, 0, 1, 0],
)
    rows = String[
        DH_TEST_SEP,
        "NDoublonHolon2siteIdx $count",
        "ComplexType $complex_type",
        DH_TEST_SEP,
        DH_TEST_SEP,
    ]
    append!(rows, main_rows)
    for (i, flag) in enumerate(opt_flags)
        push!(rows, "$(100 - i) $flag")
    end
    return join(rows, "\n") * "\n"
end

function _dh4_content(;
    count::Int = 1,
    complex_type::Int = 0,
    main_rows::Vector{String} = ["0 1 0 1 0 0", "1 0 1 0 1 0"],
    opt_flags::Vector{Int} = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
)
    rows = String[
        DH_TEST_SEP,
        "NDoublonHolon4siteIdx $count",
        "ComplexType $complex_type",
        DH_TEST_SEP,
        DH_TEST_SEP,
    ]
    append!(rows, main_rows)
    for (i, flag) in enumerate(opt_flags)
        push!(rows, "$(200 - i) $flag")
    end
    return join(rows, "\n") * "\n"
end

function _orbital_content(;
    complex_type::Int = 1,
    idx_rows::Vector{String} = ["0 0 0", "0 1 1"],
    opt_flags::Vector{Int} = [0, 1],
)
    rows = String[
        DH_TEST_SEP,
        "NOrbitalIdx $(length(opt_flags))",
        "ComplexType $complex_type",
        DH_TEST_SEP,
        DH_TEST_SEP,
    ]
    append!(rows, idx_rows)
    for (idx, flag) in enumerate(opt_flags)
        push!(rows, "$(idx - 1) $flag")
    end
    return join(rows, "\n") * "\n"
end

@testset "DH2 parser uses C index-table layout" begin
    result = parse_doublon_holon_2site_content(_dh2_content(), 3)

    @test result.success
    @test length(result.data.indices) == 1
    @test result.data.indices[1].neighbors == [1 2; 2 0; 0 1]
    @test result.data.opt_flags == Bool[1, 0, 1, 0, 1, 0]
    @test result.data.is_complex

    # The first opt-table column is ignored; row order is the local DH parameter order.
    scrambled = _dh2_content(opt_flags = [0, 1, 1, 0, 0, 1])
    scrambled_result = parse_doublon_holon_2site_content(scrambled, 3)
    @test scrambled_result.success
    @test scrambled_result.data.opt_flags == Bool[0, 1, 1, 0, 0, 1]
end

@testset "DH4 parser uses C index-table layout" begin
    result = parse_doublon_holon_4site_content(_dh4_content(), 2)

    @test result.success
    @test length(result.data.indices) == 1
    @test result.data.indices[1].neighbors == [1 0 1 0; 0 1 0 1]
    @test result.data.opt_flags == Bool[0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
    @test !result.data.is_complex
end

@testset "DH parsers reject malformed C tables" begin
    duplicate_site = _dh2_content(main_rows = ["0 1 2 0", "0 2 1 0", "2 0 1 0"])
    @test !parse_doublon_holon_2site_content(duplicate_site, 3).success

    bad_neighbor = _dh2_content(main_rows = ["0 1 3 0", "1 2 0 0", "2 0 1 0"])
    @test !parse_doublon_holon_2site_content(bad_neighbor, 3).success

    bad_opt_flag = _dh4_content(opt_flags = [0, 1, 2, 1, 0, 1, 0, 1, 0, 1])
    @test !parse_doublon_holon_4site_content(bad_opt_flag, 2).success
end

@testset "Deprecated DH value-term shims remain outside runtime data" begin
    dh2_result = parse_doublon_holon_2site_content("0 1 0.5 0 1\n")
    @test dh2_result.success
    @test dh2_result.data[1] isa DoublonHolon2SiteTerm
    @test dh2_result.data[1].value == 0.5 + 0.0im
    @test dh2_result.data[1].is_complex
    @test validate_doublon_holon_2site_terms(dh2_result.data, 2).is_valid

    dh4_result = parse_doublon_holon_4site_content("0 1 2 3 0.25 1\n")
    @test dh4_result.success
    @test dh4_result.data[1] isa DoublonHolon4SiteTerm
    @test dh4_result.data[1].value == 0.25 + 0.0im
    @test dh4_result.data[1].is_complex
    @test validate_doublon_holon_4site_terms(dh4_result.data, 4).is_valid

    @test !hasfield(ExpertModeData, :doublon_holon_2site_terms)
    @test !hasfield(ExpertModeData, :doublon_holon_4site_terms)
end

@testset "DH projection layout and opt flags are C-compatible" begin
    data = ExpertModeData()
    data.gutzwiller_terms = [GutzwillerTerm(0, 10.0 + 0.0im, false)]
    data.jastrow_terms = [JastrowTerm(0, 1, 20.0 + 0.0im, false)]
    data.doublon_holon_2site_indices = [DoublonHolon2SiteIndex([1 2; 2 0; 0 1])]
    data.doublon_holon_4site_indices = [DoublonHolon4SiteIndex([1 0 1 0; 0 1 0 1])]
    data.doublon_holon_2site_params = [ComplexF64(i, 0) for i = 1:6]
    data.doublon_holon_4site_params = [ComplexF64(10 + i, 0) for i = 1:10]
    data.doublon_holon_2site_opt_flags = Bool[1, 0, 1, 0, 1, 0]
    data.doublon_holon_4site_opt_flags = Bool[0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
    data.doublon_holon_2site_complex = true
    data.doublon_holon_4site_complex = false

    layout = projection_layout(data)
    @test layout.n_gutzwiller == 1
    @test layout.n_jastrow == 1
    @test layout.n_dh2 == 1
    @test layout.n_dh4 == 1
    @test layout.dh2_offset == 2
    @test layout.dh4_offset == 8
    @test layout.n_proj == 18

    params = projection_parameters(data, layout)
    @test params[1:2] == ComplexF64[10.0 + 0.0im, 20.0 + 0.0im]
    @test params[(layout.dh2_offset + 1):(layout.dh2_offset + 6)] ==
          data.doublon_holon_2site_params
    @test params[(layout.dh4_offset + 1):(layout.dh4_offset + 10)] ==
          data.doublon_holon_4site_params

    set_dh_opt_flags!(data)
    dh2_flag_slice = data.optimization_flags[
        (2 * layout.dh2_offset + 1):(2 * (layout.dh2_offset + 6))
    ]
    dh4_flag_slice = data.optimization_flags[
        (2 * layout.dh4_offset + 1):(2 * (layout.dh4_offset + 10))
    ]
    @test dh2_flag_slice == Bool[1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0]
    @test dh4_flag_slice == Bool[0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
                                 1, 0, 0, 0, 1, 0, 0, 0, 1, 0]
end

@testset "DH dispatch and public parse path are fatal-if-present" begin
    mktempdir() do dir
        dh2 = joinpath(dir, "dh2.def")
        write(dh2, _dh2_content())

        data = ExpertModeData()
        data.modpara.nsite = 3
        parse_file_by_type!(data, "DH2", dh2)
        @test length(data.doublon_holon_2site_indices) == 1
        @test length(data.doublon_holon_2site_params) == 6
    end

    mktempdir() do dir
        write(joinpath(dir, "modpara.def"), "Nsite 3\nNCond -1\n")
        write(joinpath(dir, "dh2.def"), _dh2_content())
        write(joinpath(dir, "namelist.def"), "ModPara modpara.def\nDH2 dh2.def\n")

        data = parse_expert_mode_files(joinpath(dir, "namelist.def"))
        @test length(data.doublon_holon_2site_indices) == 1
        @test projection_layout(data).n_proj == 6
    end

    mktempdir() do dir
        write(joinpath(dir, "modpara.def"), "Nsite 3\nNCond -1\n")
        write(joinpath(dir, "namelist.def"), "ModPara modpara.def\nDH2 missing.def\n")
        @test_throws ErrorException parse_expert_mode_files(joinpath(dir, "namelist.def"))
    end

    mktempdir() do dir
        write(joinpath(dir, "modpara.def"), "Nsite 3\nNCond -1\n")
        write(
            joinpath(dir, "dh2.def"),
            _dh2_content(main_rows = ["0 1 2 0", "0 2 1 0", "2 0 1 0"]),
        )
        write(joinpath(dir, "namelist.def"), "ModPara modpara.def\nDH2 dh2.def\n")
        @test_throws ErrorException parse_expert_mode_files(joinpath(dir, "namelist.def"))
    end
end

@testset "DH final layout is used for orbital opt flags regardless of namelist order" begin
    mktempdir() do dir
        write(joinpath(dir, "modpara.def"), "Nsite 3\nNCond -1\n")
        write(joinpath(dir, "orbital.def"), _orbital_content())
        write(joinpath(dir, "dh2.def"), _dh2_content())
        write(
            joinpath(dir, "namelist.def"),
            "ModPara modpara.def\nOrbital orbital.def\nDH2 dh2.def\n",
        )

        data = parse_expert_mode_files(joinpath(dir, "namelist.def"))
        layout = projection_layout(data)

        @test layout.n_proj == 6
        @test layout.dh2_offset == 0
        @test data.optimization_flags[1:12] ==
              Bool[1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0]
        @test data.optimization_flags[13:16] == Bool[0, 0, 1, 1]
    end
end

@testset "SpinJastrow inputs hard-fail before projection layout is used" begin
    mktempdir() do dir
        write(joinpath(dir, "namelist.def"), "SpinJastrow spinjastrow.def\n")
        @test_throws ErrorException parse_expert_mode_files(joinpath(dir, "namelist.def"))
    end
end
