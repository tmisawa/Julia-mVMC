"""
Weighted Average Functions

Calculate weighted averages for energy and SR optimization data.
"""

"""
    weight_average_we!(state::VMCOptimizationState)

Calculate weighted average of energy.
Equivalent to C's `WeightAverageWE()`.

Normalizes Wc, Etot, Etot2, Sztot, Sztot2 by Wc.
"""
function weight_average_we!(state::VMCOptimizationState)
    if abs(state.energy.wc) < 1e-15
        @warn "Weight Wc is too small: $(state.energy.wc), etot=$(state.energy.etot)"
        return
    end

    # Debug: Check if etot is zero before normalization
    if abs(state.energy.etot) < 1e-15
        @warn "Etot is zero before normalization: wc=$(state.energy.wc), etot=$(state.energy.etot)" maxlog=3
    end

    inv_w = 1.0 / state.energy.wc
    state.energy.etot *= inv_w
    state.energy.etot2 *= inv_w
    state.energy.sztot *= inv_w
    state.energy.sztot2 *= inv_w
end

"""
    weight_average_we!(ctx, state)

C `WeightAverageWE(comm_parent)`: rank-local `Wc`, `Etot`, `Etot2`, `Sztot`,
`Sztot2` を comm0 で allreduce してから、全 rank が global `Wc` で正規化する。
serial context では既存の `weight_average_we!(state)` と同一経路。
"""
function weight_average_we!(
    ctx::ParallelContext,
    state::VMCOptimizationState,
    diag_timer::CTimer = CTIMER_DISABLED,
)
    ctx.is_mpi || return weight_average_we!(state)
    buf = ComplexF64[
        state.energy.wc,
        state.energy.etot,
        state.energy.etot2,
        state.energy.sztot,
        state.energy.sztot2,
    ]
    ctimer_start!(diag_timer, 961)
    allreduce_sum!(ctx, buf; which = :comm0)
    ctimer_stop!(diag_timer, 961)
    ctimer_start!(diag_timer, 962)
    state.energy.wc = buf[1]
    if abs(state.energy.wc) < 1e-15
        is_output_rank(ctx) &&
            @warn "Weight Wc is too small after MPI allreduce: $(state.energy.wc), etot=$(buf[2])"
        ctimer_stop!(diag_timer, 962)
        return
    end
    inv_w = 1.0 / state.energy.wc
    state.energy.etot = buf[2] * inv_w
    state.energy.etot2 = buf[3] * inv_w
    state.energy.sztot = buf[4] * inv_w
    state.energy.sztot2 = buf[5] * inv_w
    ctimer_stop!(diag_timer, 962)
    return
end

"""
    weight_average_sr_opt!(state::VMCOptimizationState)

Calculate weighted average of SR optimization data (complex version).
Equivalent to C's `WeightAverageSROpt()`.

Normalizes SROptOO and SROptHO by Wc.
Note: C implementation explicitly excludes SROptO from normalization.
"""
function weight_average_sr_opt!(state::VMCOptimizationState)
    if abs(state.energy.wc) < 1e-15
        @warn "Weight Wc is too small: $(state.energy.wc)"
        return
    end

    inv_w = 1.0 / state.energy.wc

    # Normalize SROptOO and SROptHO (not SROptO - per C implementation)
    state.sr_opt.sr_opt_oo .*= inv_w
    state.sr_opt.sr_opt_ho .*= inv_w
end

"""
    weight_average_sr_opt!(ctx, state)

C `WeightAverageSROpt(comm_parent)`: `SROptOO` / `SROptHO` の rank-local sum を
comm0 allreduce し、`WeightAverageWE` 後の global `Wc` で全 rank が正規化する。
"""
function weight_average_sr_opt!(
    ctx::ParallelContext,
    state::VMCOptimizationState,
    diag_timer::CTimer = CTIMER_DISABLED,
)
    ctx.is_mpi || return weight_average_sr_opt!(state)
    if abs(state.energy.wc) < 1e-15
        is_output_rank(ctx) && @warn "Weight Wc is too small: $(state.energy.wc)"
        return
    end
    ctimer_start!(diag_timer, 963)
    allreduce_sum!(ctx, state.sr_opt.sr_opt_oo; which = :comm0)
    ctimer_stop!(diag_timer, 963)
    ctimer_start!(diag_timer, 964)
    allreduce_sum!(ctx, state.sr_opt.sr_opt_ho; which = :comm0)
    ctimer_stop!(diag_timer, 964)
    ctimer_start!(diag_timer, 965)
    inv_w = 1.0 / state.energy.wc
    state.sr_opt.sr_opt_oo .*= inv_w
    state.sr_opt.sr_opt_ho .*= inv_w
    ctimer_stop!(diag_timer, 965)
    return
end

"""
    weight_average_sr_opt_real!(state::VMCOptimizationState)

Calculate weighted average of SR optimization data (real version).
Equivalent to C's `WeightAverageSROpt_real()`.

Normalizes SROptOO_real and SROptHO_real by Wc.
Note: C implementation explicitly excludes SROptO_real from normalization.
"""
@inline function active_sr_opt_oo_real_length(sr::SROptData; nsrcg::Bool = false)
    return nsrcg ? 2 * sr.sr_opt_size : sr.sr_opt_size * sr.sr_opt_size
end

function weight_average_sr_opt_real!(state::VMCOptimizationState; nsrcg::Bool = false)
    if abs(state.energy.wc) < 1e-15
        @warn "Weight Wc is too small: $(state.energy.wc)"
        return
    end

    inv_w = 1.0 / state.energy.wc

    # Normalize SROptOO_real and SROptHO_real (not SROptO_real - per C implementation)
    oo_len = min(active_sr_opt_oo_real_length(state.sr_opt; nsrcg = nsrcg),
                 length(state.sr_opt.sr_opt_oo_real))
    @views state.sr_opt.sr_opt_oo_real[1:oo_len] .*= inv_w
    state.sr_opt.sr_opt_ho_real .*= inv_w
end

"""
    weight_average_sr_opt_real!(ctx, state)

Real-valued counterpart of `WeightAverageSROpt_real(comm_parent)`.
"""
function weight_average_sr_opt_real!(
    ctx::ParallelContext,
    state::VMCOptimizationState;
    nsrcg::Bool = false,
)
    return weight_average_sr_opt_real!(ctx, state, CTIMER_DISABLED; nsrcg = nsrcg)
end

function weight_average_sr_opt_real!(
    ctx::ParallelContext,
    state::VMCOptimizationState,
    diag_timer::CTimer;
    nsrcg::Bool = false,
)
    ctx.is_mpi || return weight_average_sr_opt_real!(state; nsrcg = nsrcg)
    if abs(state.energy.wc) < 1e-15
        is_output_rank(ctx) && @warn "Weight Wc is too small: $(state.energy.wc)"
        return
    end
    ctx.size0 <= 1 && return weight_average_sr_opt_real!(state; nsrcg = nsrcg)

    oo_len = min(active_sr_opt_oo_real_length(state.sr_opt; nsrcg = nsrcg),
                 length(state.sr_opt.sr_opt_oo_real))
    ctimer_start!(diag_timer, 963)
    allreduce_sum!(ctx, @view(state.sr_opt.sr_opt_oo_real[1:oo_len]); which = :comm0)
    ctimer_stop!(diag_timer, 963)
    ctimer_start!(diag_timer, 964)
    allreduce_sum!(ctx, state.sr_opt.sr_opt_ho_real; which = :comm0)
    ctimer_stop!(diag_timer, 964)
    ctimer_start!(diag_timer, 965)
    inv_w = 1.0 / real(state.energy.wc)
    @views state.sr_opt.sr_opt_oo_real[1:oo_len] .*= inv_w
    state.sr_opt.sr_opt_ho_real .*= inv_w
    ctimer_stop!(diag_timer, 965)
    return
end

"""
    weight_average_green_func!(state::VMCOptimizationState)

Calculate weighted average of Green's functions.
Equivalent to C's `WeightAverageGreenFunc()` in average.c.

Normalizes all Green's function arrays by Wc.
"""
function weight_average_green_func!(state::VMCOptimizationState)
    if state.phys_quantities === nothing
        return  # No physical quantities to average
    end

    if abs(state.energy.wc) < 1e-15
        @warn "Weight Wc is too small: $(state.energy.wc), cannot normalize Green's functions"
        return
    end

    phys = state.phys_quantities
    inv_w = 1.0 / state.energy.wc

    # Normalize all Green's function arrays
    phys.phys_cis_ajs .*= inv_w
    phys.phys_cis_ajs_ckt_alt .*= inv_w
    phys.phys_cis_ajs_ckt_alt_dc .*= inv_w
    phys.phys_lanczos_qqqq .*= inv_w
    phys.phys_lanczos_qcisajsq .*= inv_w
    phys.phys_lanczos_qcisajscktaltq .*= inv_w
    phys.phys_lanczos_qcisajscktaltq_dc .*= inv_w
end

"""
    weight_average_green_func!(ctx, state)

C `WeightAverageGreenFunc(comm_parent)`: Green-function accumulators are reduced
to rank0 only and normalized by global `Wc` there. Non-root ranks keep
unspecified local values, matching the C contract that only rank0 writes output.
"""
function weight_average_green_func!(ctx::ParallelContext, state::VMCOptimizationState)
    ctx.is_mpi || return weight_average_green_func!(state)
    state.phys_quantities === nothing && return
    if abs(state.energy.wc) < 1e-15
        is_output_rank(ctx) &&
            @warn "Weight Wc is too small: $(state.energy.wc), cannot normalize Green's functions"
        return
    end

    phys = state.phys_quantities
    reduce_sum_to_root!(ctx, phys.phys_cis_ajs; root = 0, which = :comm0)
    reduce_sum_to_root!(ctx, phys.phys_cis_ajs_ckt_alt; root = 0, which = :comm0)
    reduce_sum_to_root!(ctx, phys.phys_cis_ajs_ckt_alt_dc; root = 0, which = :comm0)
    reduce_sum_to_root!(ctx, phys.phys_lanczos_qqqq; root = 0, which = :comm0)
    reduce_sum_to_root!(ctx, phys.phys_lanczos_qcisajsq; root = 0, which = :comm0)
    reduce_sum_to_root!(ctx, phys.phys_lanczos_qcisajscktaltq; root = 0, which = :comm0)
    reduce_sum_to_root!(ctx, phys.phys_lanczos_qcisajscktaltq_dc; root = 0, which = :comm0)

    if is_output_rank(ctx)
        inv_w = 1.0 / state.energy.wc
        phys.phys_cis_ajs .*= inv_w
        phys.phys_cis_ajs_ckt_alt .*= inv_w
        phys.phys_cis_ajs_ckt_alt_dc .*= inv_w
        phys.phys_lanczos_qqqq .*= inv_w
        phys.phys_lanczos_qcisajsq .*= inv_w
        phys.phys_lanczos_qcisajscktaltq .*= inv_w
        phys.phys_lanczos_qcisajscktaltq_dc .*= inv_w
    end
    return
end
