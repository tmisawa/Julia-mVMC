using Test
using MVMCOptimizers

@testset "unit/types: EnergyData defaults" begin
    e = MVMCOptimizers.EnergyData()
    @test e.wc == 0.0 + 0.0im
    @test e.etot == 0.0 + 0.0im
    @test e.etot2 == 0.0 + 0.0im
    @test e.sztot == 0.0 + 0.0im
    @test e.sztot2 == 0.0 + 0.0im
end

@testset "unit/types: SROptData allocation sizes" begin
    sr_opt_size = 3
    n_vmc_sample = 7

    d_c = MVMCOptimizers.SROptData(sr_opt_size, n_vmc_sample, true)
    @test d_c.sr_opt_size == sr_opt_size
    @test length(d_c.sr_opt_oo) == 2 * sr_opt_size * (2 * sr_opt_size + 2)
    @test length(d_c.sr_opt_ho) == 2 * sr_opt_size
    @test length(d_c.sr_opt_o) == 2 * sr_opt_size
    @test length(d_c.sr_opt_o_store) == 2 * sr_opt_size * n_vmc_sample
    @test isempty(d_c.sr_opt_oo_real)
    @test isempty(d_c.sr_opt_ho_real)
    @test isempty(d_c.sr_opt_o_real)
    @test isempty(d_c.sr_opt_o_store_real)

    d_r = MVMCOptimizers.SROptData(sr_opt_size, n_vmc_sample, false)
    @test d_r.sr_opt_size == sr_opt_size
    @test length(d_r.sr_opt_oo) == 2 * sr_opt_size * (2 * sr_opt_size + 2)
    @test length(d_r.sr_opt_ho) == 2 * sr_opt_size
    @test length(d_r.sr_opt_o) == 2 * sr_opt_size
    @test length(d_r.sr_opt_o_store) == 2 * sr_opt_size * n_vmc_sample
    @test length(d_r.sr_opt_oo_real) == sr_opt_size * (sr_opt_size + 2)
    @test length(d_r.sr_opt_ho_real) == sr_opt_size
    @test length(d_r.sr_opt_o_real) == sr_opt_size
    @test length(d_r.sr_opt_o_store_real) == sr_opt_size * n_vmc_sample
end

@testset "unit/types: ElectronConfiguration sizes (fsz vs non-fsz)" begin
    n_sample = 5
    n_site = 4
    n_elec = 3
    n_proj = 11

    n_size = 2 * n_elec
    n_site2 = 2 * n_site

    ec = MVMCOptimizers.ElectronConfiguration(n_sample, n_site, n_elec, n_proj, false)
    @test length(ec.ele_idx) == n_sample * n_size
    @test length(ec.ele_cfg) == n_sample * n_site2
    @test length(ec.ele_num) == n_sample * n_site2
    @test length(ec.ele_proj_cnt) == n_sample * n_proj
    @test isempty(ec.ele_spn)
    @test length(ec.tmp_ele_idx) == n_size
    @test length(ec.tmp_ele_cfg) == n_site2
    @test length(ec.tmp_ele_num) == n_site2
    @test length(ec.tmp_ele_proj_cnt) == n_proj
    @test isempty(ec.tmp_ele_spn)
    @test length(ec.burn_ele_idx) == (n_size + n_site2 + n_site2 + n_proj)
    @test length(ec.burn_ele_cfg) == n_site2
    @test length(ec.burn_ele_num) == n_site2
    @test length(ec.burn_ele_proj_cnt) == n_proj
    @test isempty(ec.burn_ele_spn)
    @test length(ec.counter) == 10

    ec_fsz = MVMCOptimizers.ElectronConfiguration(n_sample, n_site, n_elec, n_proj, true)
    @test length(ec_fsz.ele_idx) == n_sample * n_size
    @test length(ec_fsz.ele_cfg) == n_sample * n_site2
    @test length(ec_fsz.ele_num) == n_sample * n_site2
    @test length(ec_fsz.ele_proj_cnt) == n_sample * n_proj
    @test length(ec_fsz.ele_spn) == n_sample * n_size
    @test length(ec_fsz.tmp_ele_spn) == n_size
    @test length(ec_fsz.burn_ele_idx) == (n_size + n_site2 + n_site2 + n_proj + n_size)
    @test length(ec_fsz.burn_ele_spn) == n_size
    @test length(ec_fsz.counter) == 10
end

@testset "unit/types: SlaterMatrixData size normalization" begin
    # Inputs are clamped to >= 1 in constructor
    sm_c = MVMCOptimizers.SlaterMatrixData(0, 0, 0, true)
    @test length(sm_c.pf_m) == 1
    @test isempty(sm_c.slater_elm_real)
    @test isempty(sm_c.inv_m_real)
    @test isempty(sm_c.pf_m_real)

    n_qp_full = 2
    n_site = 3
    n_elec = 4
    n_site2 = 2 * n_site
    n_size = 2 * n_elec

    sm_r = MVMCOptimizers.SlaterMatrixData(n_qp_full, n_site, n_elec, false)
    @test length(sm_r.slater_elm) == n_qp_full * n_site2 * n_site2
    @test length(sm_r.inv_m) == n_qp_full * (n_size * n_size + 1)
    @test length(sm_r.pf_m) == n_qp_full
    @test length(sm_r.slater_elm_real) == n_qp_full * n_site2 * n_site2
    @test length(sm_r.inv_m_real) == n_qp_full * (n_size * n_size + 1)
    @test length(sm_r.pf_m_real) == n_qp_full
end

