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
  with `data.modpara.rnd_seed`, using the C-compatible fallback of
  `11272` when the recorded seed is `<= 0` (matches `vmc_para_opt!`
  and `run_para_opt_from_namelist`).
  When a non-`nothing` RNG is passed in, it is used **as-is**; the
  caller is responsible for seeding it.

# Returns
- `info::Int`: Return code (0 = success, non-zero = error)

# Output Files
- `zvo_out.dat`: Energy (one line per sample, truncated on the first sample;
  not per-sample indexed, unlike C — see `output_data_phys!`)
- `zvo_var.dat`: Parameters (same convention as `zvo_out.dat`)
- `zvo_cisajs_XXX.dat`: 1-body Green's function (`XXX = ismp + NDataIdxStart`)
- `zvo_cisajscktaltex_XXX.dat`: factored two-body Green (product / `TwoBodyGEx`)
- `zvo_cisajscktalt_XXX.dat`: direct two-body Green (`TwoBodyG`)
"""
function vmc_phys_cal!(
    data::ExpertModeData;
    callback::Union{Nothing,Function} = nothing,
    rng::Union{AbstractRNG,Nothing} = nothing,
    output_dir::Union{String,Nothing} = nothing,
)::Int
    # Reject unsupported ModPara inputs (e.g. NSplitSize > 1) before any work.
    validate_supported_modpara(data.modpara)
    if MVMCExpertModeParsers.has_doublon_holon(data)
        @error "DoublonHolon (DH2/DH4) inputs are parsed but not executable until DH-2 connects projection counts/logs/loaders"
        return 1
    end
    # Reject TwoBodyGEx in FSZ / general-orbital mode before any sampling or RNG
    # side effects (its Green measurement path is not yet wired).
    validate_factored_green_supported(data)

    # Initialize RNG if not provided. Match the C-compatible seed convention
    # used by vmc_para_opt! and run_para_opt_from_namelist: when
    # data.modpara.rnd_seed <= 0, fall back to 11272.
    # Note: Do NOT re-seed an existing RNG - the caller should manage seeds.
    if rng === nothing
        rng = SFMT19937RNG()
        actual_seed = data.modpara.rnd_seed > 0 ? data.modpara.rnd_seed : 11272
        Random.seed!(rng, actual_seed)
    end
    # If rng is provided, use it as-is (caller manages the seed)

    # CRITICAL: Consume RNG to match C implementation pattern
    # C's main() calls: init_gen_rand(RndSeed) -> InitParameter() -> ReadInputParameters() -> VMCPhysCal()
    # InitParameter() consumes RNG for Slater initialization, but values are overwritten by ReadInputParameters()
    # To match C's RNG state at VMCPhysCal() entry, we must consume the same amount of RNG
    #
    # Save current parameter values (they were already loaded from InOrbital etc.)
    saved_orbital_values = [term.value for term in data.orbital_terms]
    saved_gutzwiller_values = [term.value for term in data.gutzwiller_terms]
    saved_jastrow_values = [term.value for term in data.jastrow_terms]

    # Call init_parameter! to consume RNG (matches C's InitParameter())
    init_parameter!(data; rng = rng)

    # Restore the original parameter values (simulates C's ReadInputParameters() overwriting)
    for (i, term) in enumerate(data.orbital_terms)
        term.value = saved_orbital_values[i]
    end
    for (i, term) in enumerate(data.gutzwiller_terms)
        term.value = saved_gutzwiller_values[i]
    end
    for (i, term) in enumerate(data.jastrow_terms)
        term.value = saved_jastrow_values[i]
    end

    # Get parameters
    n_data_qty_smp = data.modpara.n_data_qty_smp  # Number of sampling runs
    all_complex = get_all_complex_flag(data)
    i_flg_orbital_general = data.i_flg_orbital_general
    n_proj_bf = 0  # BackFlow only; DH is a normal projection family and is guarded above in DH-1.

    # Calculate NElec from NLocSpin and NCond if not set
    n_elec = data.modpara.nelec
    if n_elec == 0 && data.modpara.ncond != -1
        if data.modpara.ncond % 2 != 0
            @error "NCond must be even, got $(data.modpara.ncond)"
            return 1
        end
        n_elec = (data.modpara.nlocspin + data.modpara.ncond) ÷ 2
        data.modpara.nelec = n_elec
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
        @warn "NMPTrans was 0, setting to 1"
    end

    # Initialize state
    n_site = data.modpara.nsite
    n_proj = MVMCExpertModeParsers.projection_layout(data).n_proj

    n_orbital_idx = if data.modpara.n_orbital_idx > 0
        data.modpara.n_orbital_idx
    elseif !isempty(data.orbital_terms)
        maximum(t.idx for t in data.orbital_terms) + 1
    else
        0
    end

    n_rbm = has_rbm_terms(data) ? MVMCExpertModeParsers.count_rbm_parameters(data) : 0
    n_para = n_proj + n_rbm + n_orbital_idx
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
    if i_flg_orbital_general == 0
        update_slater_elm_fcmp!(data, state)
    else
        update_slater_elm_fsz!(data, state)
    end

    # Note: C's VMCPhysCal does NOT call UpdateQPWeight() (unlike VMCParaOpt)
    # The QP weights are already initialized above

    # Sampling loop
    for ismp = 0:(n_data_qty_smp-1)
        println(
            "Start: Calculate VMC physical quantities. Sampling: $ismp / $n_data_qty_smp",
        )

        # Reset physical quantities for this sampling run
        reset_phys_quantities!(state)
        clear_phys_quantity!(state)

        # Initialize output files
        init_file_phys_cal!(data, ismp)

        # VMC Sampling
        if !all_complex  # real
            # CRITICAL: Copy SlaterElm to SlaterElm_real BEFORE VMCMakeSample_real
            # C: for(tmp_i=0;tmp_i<NQPFull*(2*Nsite)*(2*Nsite);tmp_i++) SlaterElm_real[tmp_i]= creal(SlaterElm[tmp_i]);
            # This is done in VMCPhysCal BEFORE calling VMCMakeSample_real
            n_copy_slater = min(
                length(state.slater_matrix.slater_elm),
                length(state.slater_matrix.slater_elm_real),
            )
            @inbounds for i = 1:n_copy_slater
                state.slater_matrix.slater_elm_real[i] =
                    real(state.slater_matrix.slater_elm[i])
            end

            # Also copy InvM to InvM_real
            # C: for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++) InvM_real[tmp_i]= creal(InvM[tmp_i]);
            n_copy_inv = min(
                length(state.slater_matrix.inv_m),
                length(state.slater_matrix.inv_m_real),
            )
            @inbounds for i = 1:n_copy_inv
                state.slater_matrix.inv_m_real[i] = real(state.slater_matrix.inv_m[i])
            end

            # DEBUG: Verify SlaterElm_real has been copied
            @debug "vmc_phys_cal!: Before vmc_make_sample_real!, SlaterElm_real sum = $(sum(abs, state.slater_matrix.slater_elm_real))"

            if i_flg_orbital_general == 0
                if n_proj_bf == 0
                    vmc_make_sample_real!(data, state, rng)
                else
                    vmc_bf_make_sample_real!(data, state, rng)
                end
            else
                vmc_make_sample_fsz_real!(data, state, rng)
            end
        else  # complex
            if n_proj_bf == 0
                if i_flg_orbital_general == 0
                    vmc_make_sample!(data, state, rng)
                else
                    vmc_make_sample_fsz!(data, state, rng)
                end
            else
                vmc_bf_make_sample!(data, state, rng)
            end
        end

        # Main calculation (energy + Green's functions)
        if n_proj_bf == 0
            if i_flg_orbital_general == 0
                vmc_main_cal!(data, state)  # Will calculate Green's functions if mode=1
            else
                vmc_main_cal_fsz!(data, state)
            end
        else
            vmc_bf_main_cal!(data, state)
        end

        # Weighted averages
        weight_average_we!(state)
        weight_average_green_func!(state)

        # Reduce counters (no-op in single process)
        reduce_counter!(state)

        # Output data. Pass the 0-based sample counter; output_data_phys! drives
        # the energy/param write mode from it (first sample truncates) and numbers
        # the Green files with ismp + NDataIdxStart internally.
        output_data_phys!(data, state, ismp; output_dir = output_dir)

        # Close files
        close_file_phys_cal!(data, ismp)

        # Callback
        if callback !== nothing
            energy = state.energy.etot
            callback(ismp, data, energy, info)
        end

        println(
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

The energy/parameter files (`zvo_out.dat` / `zvo_var.dat`) use `ismp` directly so
the first sample (`ismp == 0`) truncates and later samples append, matching
optimization-mode semantics (and so a re-run does not accumulate stale lines).
Unlike C, these two files are not per-sample indexed; that parity is deferred to
the fixture/e2e work. The Green files are numbered `ismp + NDataIdxStart`
(`physcal_output_file_index`) to match C's per-sampling file index.
"""
function output_data_phys!(
    data::ExpertModeData,
    state::VMCOptimizationState,
    ismp::Int;
    output_dir::Union{String,Nothing} = nothing,
)
    # Output energy and parameters (same as optimization mode): the 0-based sample
    # counter selects the write mode (ismp == 0 -> truncate, else append).
    output_data!(data, state, ismp; output_dir = output_dir)

    # Output Green's functions, numbered with the C-visible per-sampling index.
    if state.phys_quantities !== nothing
        file_idx = physcal_output_file_index(data, ismp)
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
