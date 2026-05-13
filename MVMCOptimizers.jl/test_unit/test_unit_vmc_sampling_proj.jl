using Test
using MVMCOptimizers
using MVMCExpertModeParsers

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
