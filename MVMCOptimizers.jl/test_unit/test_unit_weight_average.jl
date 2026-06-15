using Test
using MVMCOptimizers

@testset "unit/weight_average: real SR active range and MPI size1 skip" begin
    state = MVMCOptimizers.VMCOptimizationState(2, 1, 0, 3, 1, 1, false, false)
    sr = state.sr_opt
    sr.sr_opt_oo_real .= collect(1.0:length(sr.sr_opt_oo_real))
    sr.sr_opt_ho_real .= collect(101.0:(100.0 + length(sr.sr_opt_ho_real)))
    state.energy.wc = 2.0 + 0.0im

    ctx_size1 = MVMCOptimizers.ParallelContext(
        true,
        nothing,
        nothing,
        nothing,
        0,
        1,
        0,
        1,
        0,
        1,
        0,
    )

    oo_active_len = sr.sr_opt_size * sr.sr_opt_size
    oo_before = copy(sr.sr_opt_oo_real)
    ho_before = copy(sr.sr_opt_ho_real)

    MVMCOptimizers.weight_average_sr_opt_real!(ctx_size1, state; nsrcg = false)

    @test sr.sr_opt_oo_real[1:oo_active_len] == oo_before[1:oo_active_len] ./ 2
    @test sr.sr_opt_oo_real[(oo_active_len+1):end] == oo_before[(oo_active_len+1):end]
    @test sr.sr_opt_ho_real == ho_before ./ 2
end

@testset "unit/weight_average: real SR-CG active OO length" begin
    state = MVMCOptimizers.VMCOptimizationState(2, 1, 0, 3, 1, 1, false, false)
    sr = state.sr_opt
    sr.sr_opt_oo_real .= collect(1.0:length(sr.sr_opt_oo_real))
    sr.sr_opt_ho_real .= collect(201.0:(200.0 + length(sr.sr_opt_ho_real)))
    state.energy.wc = 4.0 + 0.0im

    oo_active_len = 2 * sr.sr_opt_size
    oo_before = copy(sr.sr_opt_oo_real)
    ho_before = copy(sr.sr_opt_ho_real)

    MVMCOptimizers.weight_average_sr_opt_real!(state; nsrcg = true)

    @test sr.sr_opt_oo_real[1:oo_active_len] == oo_before[1:oo_active_len] ./ 4
    @test sr.sr_opt_oo_real[(oo_active_len+1):end] == oo_before[(oo_active_len+1):end]
    @test sr.sr_opt_ho_real == ho_before ./ 4
end
