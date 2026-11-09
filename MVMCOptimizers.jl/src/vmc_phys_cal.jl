"""
VMC Physical Quantity Calculation

Main function for VMCPhysCal mode (NVMCCalMode=1).
Calculates physical quantities like Green's functions with fixed parameters.
"""

"""
    vmc_phys_cal!(data::ExpertModeData;
                  callback::Union{Nothing, Function}=nothing,
                  rng::Union{AbstractRNG, Nothing}=nothing) -> Int

VMC Physical Quantity Calculation mode.
Calculates Green's functions and other physical quantities with fixed parameters.

C implementation: vmcmain.c:VMCPhysCal()

# Arguments
- `data::ExpertModeData`: Initialized Expert Mode data
- `callback`: Optional callback function called at each sampling step
  with signature `callback(ismp, data, energy, info)`
- `rng`: Random number generator (default: `nothing`).
  When `nothing`, a fresh `SFMT19937RNG()` is constructed and seeded
  with the C-parity `resolve_rnd_seed` rule (`0` → `0`, `< 0` → rank0
  time seed + bcast, positive → that value, plus MPI `group1` offset),
  matching `run_phys_cal_from_namelist` as of v0.4.
  When a non-`nothing` RNG is passed in, it is used **as-is**; the
  caller is responsible for seeding it.
- `c_timer::Union{CTimer,Nothing}`: optional C-compatible section timer. When
  `nothing` (default) all `ctimer_*` calls are no-ops and the numerical path is
  unchanged. `run_phys_cal_from_namelist` passes a `CTimer` built from
  `MVMC_C_TIMER`; it relies on the caller to write the report file.

# Returns
- `info::Int`: Return code (0 = success, non-zero = error)

# Output Files
- `zvo_out_XXX.dat`: Energy for sample `XXX = ismp + NDataIdxStart` (one line
  per file, per-sample indexed like C `InitFilePhysCal`)
- `zvo_var_XXX.dat`: Parameters (same per-sample indexing as `zvo_out_XXX.dat`)
- `zvo_cisajs_XXX.dat`: 1-body Green's function (`XXX = ismp + NDataIdxStart`)
- `zvo_cisajscktaltex_XXX.dat`: factored two-body Green (product / `TwoBodyGEx`)
- `zvo_cisajscktalt_XXX.dat`: direct two-body Green (`TwoBodyG`)
"""
function vmc_phys_cal!(
    data::ExpertModeData;
    callback::Union{Nothing,Function} = nothing,
    rng::Union{AbstractRNG,Nothing} = nothing,
    output_dir::Union{String,Nothing} = nothing,
    ctx::ParallelContext = serial_context(),
    c_timer::Union{CTimer,Nothing} = nothing,
)::Int
    # C-compatible section timer. `nothing` -> disabled singleton (no-op
    # dispatch); a concretely typed CTimer flows in from
    # `run_phys_cal_from_namelist` through this function barrier, after which the
    # per-sample ctimer_* calls specialise on its type.
    timer = c_timer === nothing ? CTIMER_DISABLED : c_timer

    # Reject unsupported global / PhysCal combinations before any work.
    validate_supported_modpara(data.modpara)
    validate_supported_phys_cal_modpara(data.modpara)
    validate_supported_phys_cal_data(data)
    # Factored Green support hook. Currently all supported PhysCal modes pass;
    # keep the call here so future unsupported combinations fail before RNG use.
    validate_factored_green_supported(data)

    # Initialize RNG if not provided. Match the C-compatible seed convention
    # used by vmc_para_opt! and run_para_opt_from_namelist.
    # Note: Do NOT re-seed an existing RNG - the caller should manage seeds.
    if rng === nothing
        rng = SFMT19937RNG()
        actual_seed = resolve_rnd_seed(ctx, data.modpara.rnd_seed, nothing)
        Random.seed!(rng, actual_seed)
    end
    # If rng is provided, use it as-is (caller manages the seed)

    # CRITICAL: Consume RNG to match C implementation pattern
    # C's main() calls: init_gen_rand(RndSeed) -> InitParameter() -> ReadInputParameters() -> VMCPhysCal()
    # InitParameter() consumes RNG for Slater initialization, but values are overwritten by ReadInputParameters()
    # To match C's RNG state at VMCPhysCal() entry, we must consume the same amount of RNG
    #
    # Save all declared Para slots, including unmapped storage and RBM/OptTrans.
    # init_parameter! resets retained storage as well as the mapped terms.
    saved_parameters = pack_parameters(data)

    # Call init_parameter! to consume RNG (matches C's InitParameter())
    init_parameter!(data; rng = rng)

    # Restore the fixed parameters without consuming any additional RNG draws.
    unpack_parameters!(data, saved_parameters)

    # Get parameters
    n_data_qty_smp = data.modpara.n_data_qty_smp  # Number of sampling runs
    all_complex = get_all_complex_flag(data)
    i_flg_orbital_general = data.i_flg_orbital_general
    n_proj_bf = 0  # BackFlow only; DH is a normal projection family in NProj.

    # Calculate NElec from NLocSpin and NCond if not set
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
        @error "NElec must be positive, got $n_elec."
        return 1
    end

    # Handle NMPTrans
    if data.modpara.nmp_trans < 0
        data.modpara.nmp_trans = abs(data.modpara.nmp_trans)
        @debug "NMPTrans was negative (anti-periodic BC), converted to $(data.modpara.nmp_trans)"
    elseif data.modpara.nmp_trans == 0
        data.modpara.nmp_trans = 1
        is_output_rank(ctx) && @warn "NMPTrans was 0, setting to 1"
    end

    # Initialize state
    n_site = data.modpara.nsite
    counts = _parameter_count_breakdown(data)
    n_proj = counts.n_proj
    n_para = counts.n_para
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

    # Initialize physical quantities
    initialize_phys_quantities!(state, data)

    # Set calculation mode to measurement
    data.modpara.vmc_calc_mode = 1

    info = 0

    # Initialize quantum projection weights FIRST
    # C: InitQPWeight() is called in main() before VMCPhysCal()
    init_qp_weight!(data)

    # Update Slater matrix elements
    # C: UpdateSlaterElm_fcmp() is called at the start of VMCPhysCal()
    ctimer_start!(timer, 20)
    if i_flg_orbital_general == 0
        update_slater_elm_fcmp!(data, state)
    else
        update_slater_elm_fsz!(data, state)
    end
    ctimer_stop!(timer, 20)

    # Note: C's VMCPhysCal does NOT call UpdateQPWeight() (unlike VMCParaOpt)
    # The QP weights are already initialized above

    # Sampling loop
    for ismp = 0:(n_data_qty_smp-1)
        is_output_rank(ctx) && println(
            "Start: Calculate VMC physical quantities. Sampling: $ismp / $n_data_qty_smp",
        )

        # Reset physical quantities for this sampling run
        reset_phys_quantities!(state)
        clear_phys_quantity!(state; clear_sr = false)

        # Initialize output files
        init_file_phys_cal!(data, ismp)

        # VMC Sampling. C wraps the real Slater copies and the sampling call in
        # StartTimer(3)/StopTimer(3) (vmcmain.c:551-602).
        ctimer_start!(timer, 3)
        if !all_complex  # real
            # CRITICAL: Copy SlaterElm to SlaterElm_real BEFORE VMCMakeSample_real
            # C: for(tmp_i=0;tmp_i<NQPFull*(2*Nsite)*(2*Nsite);tmp_i++) SlaterElm_real[tmp_i]= creal(SlaterElm[tmp_i]);
            # This is done in VMCPhysCal BEFORE calling VMCMakeSample_real
            n_copy_slater = min(
                length(state.slater_matrix.slater_elm),
                length(state.slater_matrix.slater_elm_real),
            )
            copy_complex_realpart!(
                state.slater_matrix.slater_elm_real,
                state.slater_matrix.slater_elm,
                n_copy_slater;
                threaded = true,
            )

            # Also copy InvM to InvM_real
            # C: for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++) InvM_real[tmp_i]= creal(InvM[tmp_i]);
            n_copy_inv = min(
                length(state.slater_matrix.inv_m),
                length(state.slater_matrix.inv_m_real),
            )
            copy_complex_realpart!(
                state.slater_matrix.inv_m_real,
                state.slater_matrix.inv_m,
                n_copy_inv;
                threaded = true,
            )

            # DEBUG: Verify SlaterElm_real has been copied
            @debug "vmc_phys_cal!: Before vmc_make_sample_real!, SlaterElm_real sum = $(sum(abs, state.slater_matrix.slater_elm_real))"

            if i_flg_orbital_general == 0
                if n_proj_bf == 0
                    with_vmc_real_sampler_team(data,state,ctx) do
                        vmc_make_sample_real!(data, state, rng, timer; ctx = ctx)
                    end
                else
                    vmc_bf_make_sample_real!(data, state, rng)
                end
            else
                vmc_make_sample_fsz_real!(data, state, rng, timer; ctx = ctx)
            end
        else  # complex
            if n_proj_bf == 0
                if i_flg_orbital_general == 0
                    vmc_make_sample!(data, state, rng, timer; ctx = ctx)
                else
                    vmc_make_sample_fsz!(data, state, rng, timer; ctx = ctx)
                end
            else
                vmc_bf_make_sample!(data, state, rng)
            end
        end
        ctimer_stop!(timer, 3)

        # Main calculation (energy + Green's functions)
        ctimer_start!(timer, 4)
        if n_proj_bf == 0
            if i_flg_orbital_general == 0
                vmc_main_cal!(data, state, timer, ctx)  # Will calculate Green's functions if mode=1
            else
                vmc_main_cal_fsz!(data, state, timer, ctx)
            end
        else
            vmc_bf_main_cal!(data, state)
        end
        ctimer_stop!(timer, 4)

        # Weighted averages (C StartTimer(21)/StopTimer(21) spans WE, the Green
        # average and ReduceCounter; vmcmain.c:618-624).
        ctimer_start!(timer, 21)
        weight_average_we!(ctx, state)
        weight_average_green_func!(ctx, state)

        # Reduce counters (C ReduceCounter(comm_child2); serial is no-op)
        reduce_counter!(ctx, state)
        ctimer_stop!(timer, 21)

        # Output data. Pass the 0-based sample counter; output_data_phys! drives
        # the energy/param write mode from it (first sample truncates) and numbers
        # the Green files with ismp + NDataIdxStart internally.
        ctimer_start!(timer, 22)
        is_output_rank(ctx) && output_data_phys!(data, state, ismp; output_dir = output_dir)

        # Close files
        is_output_rank(ctx) && close_file_phys_cal!(data, ismp)
        ctimer_stop!(timer, 22)

        # Callback
        if callback !== nothing
            energy = state.energy.etot
            callback(ismp, data, energy, info)
        end

        is_output_rank(ctx) && println(
            "End  : Calculate VMC physical quantities. Sampling: $ismp / $n_data_qty_smp",
        )
    end

    return info
end

"""
    init_file_phys_cal!(data::ExpertModeData, ismp::Int)

Initialize output files for physical quantity calculation.
Equivalent to C's `InitFilePhysCal()`.
"""
function init_file_phys_cal!(data::ExpertModeData, ismp::Int)
    # Files are opened in output_data_phys! when needed
    # This function is a placeholder for future file initialization
end

"""
    output_data_phys!(data::ExpertModeData, state::VMCOptimizationState, ismp::Int)

Output physical quantity data to files. `ismp` is the 0-based sample index.
Equivalent to C's `outputData()` in VMCPhysCal mode.

All per-sampling files use the C-visible file index `XXX = ismp + NDataIdxStart`:
the energy/parameter files (`zvo_out_XXX.dat` / `zvo_var_XXX.dat`) and the Green
files. C `InitFilePhysCal` opens each sample's `zvo_out_XXX.dat` /
`zvo_var_XXX.dat` with `"w"`, so every file holds exactly one sample.
"""
function output_data_phys!(
    data::ExpertModeData,
    state::VMCOptimizationState,
    ismp::Int;
    output_dir::Union{String,Nothing} = nothing,
)
    file_idx = physcal_output_file_index(data, ismp)

    # Per-sample energy/parameter files, C-indexed and truncated (see data_io.jl).
    output_data!(data, state, ismp; output_dir = output_dir, file_index = file_idx)

    # Output Green's functions, numbered with the C-visible per-sampling index.
    if state.phys_quantities !== nothing
        output_green_func!(data, state, file_idx; output_dir = output_dir)
    end
end

"""
    physcal_output_file_index(data::ExpertModeData, ismp::Int) -> Int

C-visible PhysCal per-sampling file index: `ismp + NDataIdxStart`
(so the usual first files are `*_001.dat`).
"""
physcal_output_file_index(data::ExpertModeData, ismp::Int)::Int =
    ismp + data.modpara.n_data_idx_start

"""
    close_file_phys_cal!(data::ExpertModeData, ismp::Int)

Close output files for physical quantity calculation.
Equivalent to C's `CloseFilePhysCal()`.
"""
function close_file_phys_cal!(data::ExpertModeData, ismp::Int)
    # Files are closed automatically when using `open()` with do block
    # This function is a placeholder for future file closing logic
end
