"""
VMC Parameter Optimization

Main optimization loop using Stochastic Reconfiguration method.
"""

"""
    vmc_para_opt!(data::ExpertModeData;
                  callback::Union{Nothing, Function}=nothing,
                  rng::Union{AbstractRNG, Nothing}=nothing) -> Int

VMC Parameter Optimization using Stochastic Reconfiguration method.

This function performs the same optimization loop as C's `VMCParaOpt()`.

# Arguments
- `data::ExpertModeData`: Expert Mode data with initialized parameters
- `callback::Union{Nothing, Function}`: Optional callback function called at each step
  with signature `callback(step, data, energy, info)`
- `rng::Union{AbstractRNG, Nothing}`: Random number generator (default: `nothing`).
  When `nothing`, a fresh `SFMT19937RNG()` is constructed and seeded via
  `resolve_rnd_seed(ctx, data.modpara.rnd_seed, nothing)` — the same C-parity
  rule as `run_para_opt_from_namelist` (`< 0` → rank0 time seed + bcast,
  `== 0` → `0`, `> 0` → the value, plus the per-group `+ group1` offset
  under MPI; see `resolve_rnd_seed`).
  When a non-`nothing` RNG is passed in, it is used **as-is**; the caller
  is responsible for seeding it. This applies to both `SFMT19937RNG` and
  any other `AbstractRNG` (for example a stable `Random.MersenneTwister`
  for tests).
- `output_dir::Union{String, Nothing}`: If set, output files (zvo_out.dat, zvo_var.dat, zqp_opt.dat, etc.) are written under this directory. The directory is created if it does not exist. Default: `nothing` (current directory).

# Returns
- `info::Int`: Return code (0 = success, non-zero = error)

# Example

For most workflows prefer `run_para_opt_from_namelist`, which wires the
seed / parameter / QP-weight initialization in the C-compatible order
for you. If you do it by hand, seed an `SFMT19937RNG` first because
`initialize_parameters!` needs an initialized RNG (calling
`SFMT19937RNG()` without a `Random.seed!` raises an error on first
draw):

```julia
using Random
using MVMCExpertModeParsers
using MVMCOptimizers
using SFMT

data = parse_expert_mode_files("namelist.def")

rng = SFMT19937RNG()
actual_seed = MVMCOptimizers.resolve_rnd_seed(
    MVMCOptimizers.serial_context(), data.modpara.rnd_seed, nothing)
Random.seed!(rng, actual_seed)

MVMCExpertModeParsers.initialize_parameters!(data; rng=rng)
MVMCExpertModeParsers.read_input_parameters!(data, "namelist.def")
MVMCExpertModeParsers.init_qp_weight!(data)

# `callback` is a keyword argument, so the `do` block syntax does not apply;
# pass an explicit anonymous function instead.
# Pass the same RNG into vmc_para_opt! so it is used as-is (the function
# will not re-seed it).
info = vmc_para_opt!(data;
    rng = rng,
    callback = (step, data, e, info) -> println("Step \$step: Energy = \$e"),
)
```

Letting `vmc_para_opt!` build its own RNG (`rng = nothing`) is also
supported and follows the same seed convention as above; just remember
that any work done before this call (e.g. `initialize_parameters!`)
needs its own seeded RNG.
"""

function vmc_para_opt!(
    data::ExpertModeData;
    callback::Union{Nothing,Function} = nothing,
    rng::Union{AbstractRNG,Nothing} = nothing,
    output_dir::Union{String,Nothing} = nothing,
    skip_sr::Bool = false,
    c_timer::Union{CTimer,Nothing} = nothing,
    ctx::ParallelContext = serial_context(),
)::Int
    # Reject unsupported ModPara inputs (e.g. NSplitSize > 1) before any work.
    validate_supported_modpara(data.modpara)
    validate_supported_para_opt_parallel_modpara(ctx, data.modpara)

    # C-compatible section timer. `nothing` -> disabled singleton (no-op
    # start/stop). When run_para_opt_from_namelist enables it, a concretely
    # typed CTimer{Val{true}} flows in; this assignment is the function-barrier
    # entry, after which the loop-level ctimer_* calls specialise on its type.
    timer = c_timer === nothing ? CTIMER_DISABLED : c_timer

    # If no RNG was passed in, fall back to a SFMT19937RNG seeded with the
    # same C-parity rule as run_para_opt_from_namelist (resolve_rnd_seed:
    # < 0 → rank0 time seed + bcast, == 0 → 0, > 0 → value, + ctx.group1).
    # Review 2026-06-11 F6-(1): the old `> 0 ? : 11272` rule here diverged
    # from the runner and ignored the MPI group1 offset.
    if rng === nothing
        rng = SFMT19937RNG()
        Random.seed!(rng, resolve_rnd_seed(ctx, data.modpara.rnd_seed, nothing))
    end

    # Get optimization parameters
    n_steps = data.modpara.nsr_opt_itr_step
    n_smp = data.modpara.nsr_opt_itr_smp
    # Use NSRCG flag from modpara.def (C implementation uses NSRCG global variable)
    n_sr_cg = data.modpara.nsrcg  # 0 = direct solver (LAPACK), !=0 = CG solver
    all_complex = get_all_complex_flag(data)
    weightavg_diag_timer = ctimer_if_env(timer, "MVMC_WEIGHTAVG_DIAG")
    # i_flg_orbital_general: 0 = sz conserved, non-zero = general (fsz)
    i_flg_orbital_general = data.i_flg_orbital_general
    n_proj_bf = 0  # BackFlow only; DH is a normal projection family in NProj.

    # Calculate NElec from NLocSpin and NCond if not set (C code: readdef.c:593)
    n_elec = data.modpara.nelec
    if n_elec == 0 && data.modpara.ncond != -1
        if data.modpara.ncond % 2 != 0
            @error "NCond must be even, got $(data.modpara.ncond)"
            return 1
        end
        n_elec = (data.modpara.nlocspin + data.modpara.ncond) ÷ 2
        data.modpara.nelec = n_elec
        is_output_rank(ctx) &&
            @info "Calculated NElec = $n_elec from NLocSpin = $(data.modpara.nlocspin) and NCond = $(data.modpara.ncond)"
    end

    # Validate n_elec
    if n_elec <= 0
        @error "NElec must be positive, got $n_elec. Please set NElec in modpara.def or ensure NCond and NLocSpin are set correctly."
        return 1
    end

    # Handle NMPTrans: negative values indicate anti-periodic boundary condition
    # C implementation: if NMPTrans < 0, set APFlag = 1 and NMPTrans *= -1
    if data.modpara.nmp_trans < 0
        # APFlag is handled implicitly by keeping track of the original sign
        # For quantum projection calculations, we need the absolute value
        data.modpara.nmp_trans = abs(data.modpara.nmp_trans)
        @debug "NMPTrans was negative (anti-periodic BC), converted to $(data.modpara.nmp_trans)"
    elseif data.modpara.nmp_trans == 0
        data.modpara.nmp_trans = 1
        is_output_rank(ctx) && @warn "NMPTrans was 0, setting to 1"
    end

    # Initialize optimization state
    n_site = data.modpara.nsite
    counts = _parameter_count_breakdown(data)
    n_proj = counts.n_proj
    n_rbm = counts.n_rbm
    n_orbital_idx = counts.n_orbital_idx
    n_opt_trans = counts.n_opt_trans
    n_para = counts.n_para
    is_output_rank(ctx) &&
        @info "NPara=$n_para (NProj=$n_proj + NRBM=$n_rbm + NOrbitalIdx=$n_orbital_idx + NOptTrans=$n_opt_trans)"

    # Initialize optimization flags if empty
    # C: OptFlag[2*i] = 1 for parameters to optimize
    # Default: all parameters are optimized (OptFlag = 1)
    if isempty(data.optimization_flags)
        # 2 * NPara entries (real and imaginary parts)
        data.optimization_flags = fill(true, 2 * n_para)
        is_output_rank(ctx) &&
            @info "Initialized optimization_flags: all $n_para parameters will be optimized"
    end
    n_qp_full = get_n_qp_full(data)
    n_vmc_sample = data.modpara.nvmc_sample

    state = VMCOptimizationState(
        n_site,
        n_elec,
        n_proj,
        n_para,
        n_qp_full,
        n_vmc_sample,
        all_complex,
        i_flg_orbital_general != 0,
    )

    # Debug: Verify state was created correctly
    @debug "VMCOptimizationState created: ele_idx size=$(length(state.electron_config.ele_idx)), expected=$(n_vmc_sample * 2 * n_elec)"

    # Initialize burn_flag: Step 0 uses new initial sample, Step 1+ uses burn sample
    # Store burn_flag in counter[11] (0 = false, 1 = true)
    if length(state.electron_config.counter) < 11
        resize!(state.electron_config.counter, 11)
    end
    state.electron_config.counter[11] = 0  # Step 0: start with new initial sample

    info = 0

    # [2] VMCParaOpt: wraps the whole optimisation loop and final opt output
    # (matches C's StartTimer(2)/StopTimer(2) around VMCParaOpt()).
    ctimer_start!(timer, 2)

    # Optimization loop
    for step = 0:(n_steps-1)
        # Progress output (rank0 のみ; C の CLI 出力は rank0 限定。plan review F5)
        if is_output_rank(ctx) &&
           (step == 0 || (n_steps < 20) || (step % max(1, n_steps ÷ 20) == 0))
            progress = floor(Int, 100.0 * step / n_steps)
            println("Progress of Optimization: $progress %")
        end

        # [20] UpdateSlaterElm: UpdateSlaterElm_* + UpdateQPWeight (C vmcmain.c:353-362)
        ctimer_start!(timer, 20)
        # 1. Update Slater matrix elements
        if i_flg_orbital_general == 0
            update_slater_elm_fcmp!(data, state)
        else
            update_slater_elm_fsz!(data, state)
        end

        # 2. Update quantum projection weights
        update_qp_weight!(data)
        ctimer_stop!(timer, 20)

        # [3] VMCMakeSample (C vmcmain.c:363-404)
        ctimer_start!(timer, 3)
        # 3. VMC Sampling
        if !all_complex  # real
            # Convert to real arrays if needed ([69] MAll: real/complex copy)
            if !isempty(state.slater_matrix.slater_elm_real)
                ctimer_start!(timer, 69)
                convert_to_real_arrays!(state; threaded = true)
                ctimer_stop!(timer, 69)
            end

            if i_flg_orbital_general == 0
                if n_proj_bf == 0
                    vmc_make_sample_real!(data, state, rng, timer)
                else
                    vmc_bf_make_sample_real!(data, state, rng)
                end
            else
                vmc_make_sample_fsz_real!(data, state, rng, timer)
            end

            # Convert back to complex if needed ([69] MAll: real/complex copy)
            if !isempty(state.slater_matrix.inv_m_real)
                ctimer_start!(timer, 69)
                convert_from_real_arrays!(state; threaded = true)
                ctimer_stop!(timer, 69)
            end
        else  # complex
            if n_proj_bf == 0
                if i_flg_orbital_general == 0
                    vmc_make_sample!(data, state, rng, timer)
                else
                    vmc_make_sample_fsz!(data, state, rng, timer)
                end
            else
                vmc_bf_make_sample!(data, state, rng)
            end
        end
        ctimer_stop!(timer, 3)

        # [4] VMCMainCal (C vmcmain.c:405-418)
        ctimer_start!(timer, 4)
        # 4. Main calculation (energy and SR quantities)
        if n_proj_bf == 0
            if i_flg_orbital_general == 0
                vmc_main_cal!(data, state, timer)
            else
                vmc_main_cal_fsz!(data, state, timer)
            end
        else
            vmc_bf_main_cal!(data, state)
        end
        ctimer_stop!(timer, 4)

        # [21] WeightAverage: WeightAverageWE + [25 SR] + ReduceCounter, with
        # [25] nested inside [21] exactly as in C (vmcmain.c:419-436).
        ctimer_start!(timer, 21)
        ctimer_start!(weightavg_diag_timer, 960)
        # 5. Weighted averages
        weight_average_we!(ctx, state, weightavg_diag_timer)

        ctimer_start!(timer, 25)
        if !all_complex
            weight_average_sr_opt_real!(
                ctx,
                state,
                weightavg_diag_timer;
                nsrcg = data.modpara.nsrcg != 0,
            )
        else
            weight_average_sr_opt!(ctx, state, weightavg_diag_timer)
        end
        ctimer_stop!(timer, 25)

        # 6. Reduce counters (C ReduceCounter(comm_child2); serial is no-op)
        ctimer_start!(weightavg_diag_timer, 966)
        reduce_counter!(ctx, state)
        ctimer_stop!(weightavg_diag_timer, 966)
        ctimer_stop!(weightavg_diag_timer, 960)
        ctimer_stop!(timer, 21)

        # [22] outputData (C vmcmain.c:437-440)
        ctimer_start!(timer, 22)
        # 7. Output data (before optimization, matching C implementation order)
        # C vmcmain.c:441 `if(rank==0) outputData()` — 出力は rank0 のみ（spec §5-9）
        is_output_rank(ctx) && output_data!(data, state, step; output_dir=output_dir)
        ctimer_stop!(timer, 22)

        # If skip_sr is set, only run step 0 and return after output
        if skip_sr
            ctimer_stop!(timer, 2)
            if callback !== nothing
                energy = state.energy.etot
                callback(step, data, energy, 0)
            end
            return 0
        end

        # [5] StochasticOpt: SR solver (C vmcmain.c StochasticOptDiag/CG path)
        ctimer_start!(timer, 5)
        # 8. Stochastic optimization
        if n_sr_cg != 0
            info = stochastic_opt_cg!(data, state, timer)
        else
            info = stochastic_opt!(data, state, timer)
        end
        ctimer_stop!(timer, 5)

        # C stcopt.c:171 は SR の info を rank0 から MPI_Bcast してから return 判定
        # （vmcmain.c:504-506）するため、全 rank が同じ判断で loop を抜ける。
        # rank-local な info のまま early return すると、失敗 rank だけが loop を
        # 抜けて comm0 の collective（sync の Bcast! / readback の barrier）が
        # 不整合になり hang する（review 2026-06-11 F1）。
        info = Int(bcast_scalar(ctx, info))
        if info != 0
            is_output_rank(ctx) &&
                @error "Stochastic optimization error: info=$info at step=$step"
            ctimer_stop!(timer, 2)
            return info
        end

        # [23] SyncModifiedParameter (C vmcmain.c:509-510, comm_parent で bcast)
        ctimer_start!(timer, 23)
        # 9. Sync modified parameters
        sync_modified_parameter!(ctx, data)
        ctimer_stop!(timer, 23)

        # 10. Store optimization data (for averaging)
        if step >= n_steps - n_smp
            store_opt_data!(data, state, step - (n_steps - n_smp))
        end

        # Callback
        if callback !== nothing
            energy = state.energy.etot
            callback(step, data, energy, info)
        end
    end

    # Final output (rank0 のみ; C `if(rank==0) outputData()` 相当)
    is_output_rank(ctx) && println("Start: Output opt params.")
    is_output_rank(ctx) && output_opt_data!(data; output_dir=output_dir)
    is_output_rank(ctx) && println("End: Output opt params.")

    ctimer_stop!(timer, 2)

    return info
end

# Helper functions
function get_all_complex_flag(data::ExpertModeData)::Bool
    # Check if any complex flags are set
    if !isempty(data.complex_flags)
        return any(x -> x != 0, data.complex_flags)
    end
    # Default: check if any terms are complex
    any_complex = false
    for term in data.orbital_terms
        if term.is_complex || imag(term.value) != 0.0
            any_complex = true
            break
        end
    end
    # Also check Gutzwiller and Jastrow terms
    if !any_complex
        for term in data.gutzwiller_terms
            if term.is_complex || imag(term.value) != 0.0
                any_complex = true
                break
            end
        end
    end
    if !any_complex
        for term in data.jastrow_terms
            if term.is_complex || imag(term.value) != 0.0
                any_complex = true
                break
            end
        end
    end
    if !any_complex
        any_complex = (
            data.doublon_holon_2site_complex ||
            data.doublon_holon_4site_complex ||
            any(x -> imag(x) != 0.0, data.doublon_holon_2site_params) ||
            any(x -> imag(x) != 0.0, data.doublon_holon_4site_params)
        )
    end
    return any_complex
end

function get_n_qp_full(data::ExpertModeData)::Int
    # Calculate NQPFull = NSPGaussLeg * NMPTrans * NQPTransOpt
    n_sp_gauss_leg = max(1, data.modpara.nsp_gauss_leg)  # Ensure at least 1
    n_mp_trans = max(1, data.modpara.nmp_trans)  # Ensure at least 1 (handle -1 case)
    n_qp_opt_trans = max(1, data.n_qp_opt_trans)  # Ensure at least 1
    return n_sp_gauss_leg * n_mp_trans * n_qp_opt_trans
end

function convert_to_real_arrays!(state::VMCOptimizationState; threaded::Bool = false)
    # Convert InvM to real arrays
    # Note: slater_elm_real is copied in vmc_make_sample_real! so we skip it here
    # Optimized: only convert if arrays are not empty and have matching sizes
    if !isempty(state.slater_matrix.inv_m) && !isempty(state.slater_matrix.inv_m_real)
        n_copy =
            min(length(state.slater_matrix.inv_m), length(state.slater_matrix.inv_m_real))
        copy_complex_realpart!(
            state.slater_matrix.inv_m_real,
            state.slater_matrix.inv_m,
            n_copy;
            threaded = threaded,
        )
    end
end

function convert_from_real_arrays!(state::VMCOptimizationState; threaded::Bool = false)
    # Convert InvM back from real arrays
    # Optimized: only convert if arrays are not empty and have matching sizes
    if !isempty(state.slater_matrix.inv_m_real) && !isempty(state.slater_matrix.inv_m)
        n_copy =
            min(length(state.slater_matrix.inv_m_real), length(state.slater_matrix.inv_m))
        copy_real_to_complex!(
            state.slater_matrix.inv_m,
            state.slater_matrix.inv_m_real,
            n_copy;
            threaded = threaded,
        )
    end
end
