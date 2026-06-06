using Test
using MVMCOptimizers
using MVMCExpertModeParsers
using MVMCExpertModeParsers:
    DoublonHolon2SiteIndex,
    DoublonHolon4SiteIndex,
    GutzwillerTerm,
    JastrowTerm,
    projection_layout

if !isdefined(@__MODULE__, :make_ele_num)
    include(joinpath(@__DIR__, "helpers", "mock_state.jl"))
end
if !isdefined(@__MODULE__, :make_mock_data_for_proj_tests)
    include(joinpath(@__DIR__, "helpers", "mock_data.jl"))
end

@testset "unit/vmc_sampling: make_proj_cnt! / update_proj_cnt! consistency" begin
    data = make_mock_data_for_proj_tests(nsite = 4)
    nsite = data.modpara.nsite
    nproj = data.n_gutzwiller_idx + data.n_jastrow_idx

    # Electron occupation (length = 2*nsite): up[1:nsite], down[nsite+1:2*nsite]
    # up: sites 0,2 ; down: sites 1,3
    ele_num0 = make_ele_num(nsite; up_sites = [0, 2], down_sites = [1, 3])

    proj0 = zeros(Int, nproj)
    MVMCOptimizers.make_proj_cnt!(proj0, ele_num0, data)

    # One hop: up electron 0 -> 1 (s=0); destination up must be empty.
    ri, rj, s = 0, 1, 0
    @test ele_num0[ri + s * nsite + 1] == 1
    @test ele_num0[rj + s * nsite + 1] == 0
    ele_num1 = apply_hop(ele_num0, ri, rj, s, nsite)

    proj1_full = zeros(Int, nproj)
    MVMCOptimizers.make_proj_cnt!(proj1_full, ele_num1, data)

    proj1_inc = similar(proj0)
    MVMCOptimizers.update_proj_cnt!(ri, rj, s, proj1_inc, proj0, ele_num1, data)
    @test proj1_inc == proj1_full

    # fsz variant should match recompute as well (for a single hop).
    proj1_inc_fsz = similar(proj0)
    MVMCOptimizers.update_proj_cnt_fsz!(ri, rj, s, 1 - s, proj1_inc_fsz, proj0, ele_num1, data)
    @test proj1_inc_fsz == proj1_full

    # No-op move (ri==rj) should keep counts unchanged
    proj_noop = similar(proj0)
    MVMCOptimizers.update_proj_cnt!(2, 2, 0, proj_noop, proj0, ele_num0, data)
    @test proj_noop == proj0

    # Second hop: down electron 3 -> 0 (s=1); destination down must be empty.
    ri2, rj2, s2 = 3, 0, 1
    @test ele_num1[ri2 + s2 * nsite + 1] == 1
    @test ele_num1[rj2 + s2 * nsite + 1] == 0
    ele_num2 = apply_hop(ele_num1, ri2, rj2, s2, nsite)

    proj2_full = zeros(Int, nproj)
    MVMCOptimizers.make_proj_cnt!(proj2_full, ele_num2, data)

    proj2_inc = similar(proj0)
    MVMCOptimizers.update_proj_cnt!(ri, rj, s, proj2_inc, proj0, ele_num1, data)
    MVMCOptimizers.update_proj_cnt!(ri2, rj2, s2, proj2_inc, proj2_inc, ele_num2, data)
    @test proj2_inc == proj2_full

    proj2_inc_fsz = similar(proj0)
    MVMCOptimizers.update_proj_cnt_fsz!(ri, rj, s, 1 - s, proj2_inc_fsz, proj0, ele_num1, data)
    MVMCOptimizers.update_proj_cnt_fsz!(ri2, rj2, s2, 1 - s2, proj2_inc_fsz, proj2_inc_fsz, ele_num2, data)
    @test proj2_inc_fsz == proj2_full
end

function _dh_only_data(; nsite::Int, dh2 = nothing, dh4 = nothing)
    data = MVMCExpertModeParsers.ExpertModeData()
    data.modpara = MVMCExpertModeParsers.ModParaParameters(nsite = nsite)
    if dh2 !== nothing
        data.doublon_holon_2site_indices = [DoublonHolon2SiteIndex(dh2)]
        data.doublon_holon_2site_params = [ComplexF64(i, 100 + i) for i = 1:6]
        data.doublon_holon_2site_opt_flags = fill(true, 6)
    end
    if dh4 !== nothing
        data.doublon_holon_4site_indices = [DoublonHolon4SiteIndex(dh4)]
        data.doublon_holon_4site_params = [ComplexF64(10 + i, -10 - i) for i = 1:10]
        data.doublon_holon_4site_opt_flags = fill(true, 10)
    end
    return data
end

@testset "unit/vmc_sampling: DH2/DH4 projection counts" begin
    dh2_neighbors = [
        1 2
        0 2
        0 1
        0 1
    ]
    data2 = _dh_only_data(nsite = 4, dh2 = dh2_neighbors)
    ele_num2 = make_ele_num(4; up_sites = [0, 2], down_sites = [0, 3])
    proj2 = zeros(Int, projection_layout(data2).n_proj)

    MVMCOptimizers.make_proj_cnt!(proj2, ele_num2, data2)

    # site0 is a doublon with one holon neighbor; site1 is a holon with one
    # doublon neighbor. Singly occupied sites are skipped.
    @test proj2 == [0, 0, 1, 1, 0, 0]

    sr_opt_o = zeros(ComplexF64, 2 * (length(proj2) + 1))
    MVMCOptimizers.set_projection_diff!(sr_opt_o, proj2, length(proj2))
    @test [real(sr_opt_o[(i + 1) * 2 + 1]) for i = 0:(length(proj2) - 1)] == proj2

    dh4_neighbors = [
        1 2 3 4
        0 2 3 4
        0 1 3 4
        0 1 2 4
        0 1 2 3
    ]
    data4 = _dh_only_data(nsite = 5, dh4 = dh4_neighbors)
    ele_num4 = make_ele_num(5; up_sites = [0, 2, 3], down_sites = [0, 3])
    proj4 = zeros(Int, projection_layout(data4).n_proj)

    MVMCOptimizers.make_proj_cnt!(proj4, ele_num4, data4)

    # Holon sites 1 and 4 each see two doublons; doublon sites 0 and 3 each
    # see two holons. Site2 is singly occupied and skipped.
    @test proj4 == [0, 0, 0, 0, 2, 2, 0, 0, 0, 0]
end

@testset "unit/vmc_sampling: multiple DH2 indices keep independent strides" begin
    data = MVMCExpertModeParsers.ExpertModeData()
    data.modpara = MVMCExpertModeParsers.ModParaParameters(nsite = 4)
    data.doublon_holon_2site_indices = [
        DoublonHolon2SiteIndex([
            1 2
            0 2
            0 1
            0 1
        ]),
        DoublonHolon2SiteIndex([
            2 3
            2 3
            0 1
            0 1
        ]),
    ]
    data.doublon_holon_2site_params = [ComplexF64(i, 100 + i) for i = 1:12]
    data.doublon_holon_2site_opt_flags = fill(true, 12)

    ele_num = make_ele_num(4; up_sites = [0, 2], down_sites = [0, 3])
    proj = zeros(Int, projection_layout(data).n_proj)

    MVMCOptimizers.make_proj_cnt!(proj, ele_num, data)

    # DH index 1 contributes to xn0=0 slots; DH index 2 contributes to xn0=1 slots.
    @test proj == [0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0]
end

@testset "unit/vmc_sampling: DH update paths match fresh recompute" begin
    data = make_mock_data_for_proj_tests(nsite = 4)
    data.gutzwiller_terms = [
        GutzwillerTerm(0, 0.11 + 0.0im, false),
        GutzwillerTerm(1, 0.12 + 0.0im, false),
    ]
    data.jastrow_terms = [JastrowTerm(0, 1, ComplexF64(0.2 + i / 100), false) for i = 1:data.n_jastrow_idx]
    data.doublon_holon_2site_indices = [
        DoublonHolon2SiteIndex([
            1 2
            0 2
            0 1
            0 1
        ]),
    ]
    data.doublon_holon_2site_params = [ComplexF64(0.3 + i / 10, -i) for i = 1:6]
    data.doublon_holon_2site_opt_flags = fill(true, 6)

    nsite = data.modpara.nsite
    nproj = projection_layout(data).n_proj
    ele_num0 = make_ele_num(nsite; up_sites = [0, 2], down_sites = [1, 3])
    proj0 = zeros(Int, nproj)
    MVMCOptimizers.make_proj_cnt!(proj0, ele_num0, data)

    ri, rj, s = 0, 1, 0
    ele_num1 = apply_hop(ele_num0, ri, rj, s, nsite)
    proj1_full = zeros(Int, nproj)
    MVMCOptimizers.make_proj_cnt!(proj1_full, ele_num1, data)

    proj1_inc = similar(proj0)
    MVMCOptimizers.update_proj_cnt!(ri, rj, s, proj1_inc, proj0, ele_num1, data)
    @test proj1_inc == proj1_full

    proj1_alias = copy(proj0)
    MVMCOptimizers.update_proj_cnt!(ri, rj, s, proj1_alias, proj1_alias, ele_num1, data)
    @test proj1_alias == proj1_full

    proj1_fsz = similar(proj0)
    MVMCOptimizers.update_proj_cnt_fsz!(ri, rj, s, 1 - s, proj1_fsz, proj0, ele_num1, data)
    @test proj1_fsz == proj1_full

    proj1_fsz_alias = copy(proj0)
    MVMCOptimizers.update_proj_cnt_fsz!(ri, rj, s, 1 - s, proj1_fsz_alias, proj1_fsz_alias, ele_num1, data)
    @test proj1_fsz_alias == proj1_full

    # On-site spin flip leaves occupancy class unchanged, so FSZ keeps the old DH tail.
    ele_num_flip = make_ele_num(nsite; up_sites = Int[], down_sites = [0, 3])
    proj_noop = similar(proj0)
    MVMCOptimizers.update_proj_cnt_fsz!(0, 0, 0, 1, proj_noop, proj0, ele_num_flip, data)
    @test proj_noop == proj0
end

@testset "unit/vmc_sampling: DH4 update paths match fresh recompute" begin
    data = make_mock_data_for_proj_tests(nsite = 4)
    data.gutzwiller_terms = [
        GutzwillerTerm(0, 0.11 + 0.0im, false),
        GutzwillerTerm(1, 0.12 + 0.0im, false),
    ]
    data.jastrow_terms = [JastrowTerm(0, 1, ComplexF64(0.2 + i / 100), false) for i = 1:data.n_jastrow_idx]
    data.doublon_holon_4site_indices = [
        DoublonHolon4SiteIndex([
            1 2 3 1
            0 2 3 0
            0 1 3 0
            0 1 2 1
        ]),
    ]
    data.doublon_holon_4site_params = [ComplexF64(0.5 + i / 10, -i) for i = 1:10]
    data.doublon_holon_4site_opt_flags = fill(true, 10)

    nsite = data.modpara.nsite
    nproj = projection_layout(data).n_proj
    ele_num0 = make_ele_num(nsite; up_sites = [0, 2], down_sites = [0, 3])
    proj0 = zeros(Int, nproj)
    MVMCOptimizers.make_proj_cnt!(proj0, ele_num0, data)

    ri, rj, s = 2, 1, 0
    ele_num1 = apply_hop(ele_num0, ri, rj, s, nsite)
    proj1_full = zeros(Int, nproj)
    MVMCOptimizers.make_proj_cnt!(proj1_full, ele_num1, data)

    proj1_inc = similar(proj0)
    MVMCOptimizers.update_proj_cnt!(ri, rj, s, proj1_inc, proj0, ele_num1, data)
    @test proj1_inc == proj1_full

    proj1_alias = copy(proj0)
    MVMCOptimizers.update_proj_cnt!(ri, rj, s, proj1_alias, proj1_alias, ele_num1, data)
    @test proj1_alias == proj1_full

    proj1_fsz = similar(proj0)
    MVMCOptimizers.update_proj_cnt_fsz!(ri, rj, s, 1 - s, proj1_fsz, proj0, ele_num1, data)
    @test proj1_fsz == proj1_full

    proj1_fsz_alias = copy(proj0)
    MVMCOptimizers.update_proj_cnt_fsz!(ri, rj, s, 1 - s, proj1_fsz_alias, proj1_fsz_alias, ele_num1, data)
    @test proj1_fsz_alias == proj1_full
end

@testset "unit/vmc_sampling: DH log projection uses real parameter parts" begin
    data = _dh_only_data(
        nsite = 4,
        dh2 = [
            1 2
            0 2
            0 1
            0 1
        ],
    )
    data.doublon_holon_2site_params .= ComplexF64[
        1 + 10im,
        2 + 20im,
        3 + 30im,
        4 + 40im,
        5 + 50im,
        6 + 60im,
    ]
    cnt_old = [0, 0, 1, 1, 0, 0]
    cnt_new = [0, 1, 0, 2, 0, 1]

    @test MVMCOptimizers.log_proj_val(cnt_old, data) ≈ 7.0
    @test MVMCOptimizers.log_proj_ratio(cnt_new, cnt_old, data) ≈ 9.0
end
