"""
Parameter Initialization Utilities

Functions for initializing variational parameters similar to C implementation's InitParameter().
Based on mVMC/src/mVMC/parameter.c
"""

# Constants from C implementation
const D_AMP_MAX = 4.0  # Maximum amplitude for Slater parameters (D_AmpMax)

"""
    init_parameter!(data::ExpertModeData; rng::AbstractRNG=SFMT19937RNG())

Initialize variational parameters in ExpertModeData similar to C's InitParameter().

This function:
- Initializes Proj (Gutzwiller + Jastrow) parameters to 0.0
- Initializes RBM parameters randomly if FlagRBM > 0
- Initializes Slater (Orbital) parameters randomly based on OptFlag
- Initializes OptTrans parameters from ParaQPOptTrans

# Arguments
- `data::ExpertModeData`: Expert Mode data structure to initialize
- `rng::AbstractRNG`: Random number generator (default: SFMT19937RNG())

# Notes
- Uses the same random initialization logic as C implementation
- AllComplexFlag determines whether complex or real parameters are used
- Only parameters with OptFlag > 0 are initialized randomly
- Parameters with OptFlag == 0 are set to 0.0
"""
function init_parameter!(data::ExpertModeData; rng::AbstractRNG = SFMT19937RNG())
    modpara = data.modpara

    # Note: Do NOT re-seed the RNG here!
    # The RNG should be seeded by the caller (e.g., vmc_para_opt!) before calling this function.
    # C implementation seeds once in main() and then calls InitParameter() without re-seeding.
    # Re-seeding here would reset the RNG state and cause mismatch with C implementation.

    # Initialize Proj parameters (Gutzwiller + Jastrow) to 0.0
    # In C: for(i=0;i<NProj;i++) Proj[i] = 0.0;
    for term in data.gutzwiller_terms
        term.value = ComplexF64(0.0)
    end
    for term in data.jastrow_terms
        term.value = ComplexF64(0.0)
    end
    fill!(data.doublon_holon_2site_params, 0.0 + 0.0im)
    fill!(data.doublon_holon_4site_params, 0.0 + 0.0im)

    # C uses AllComplexFlag across variational-factor groups (not only ModPara.ComplexType).
    gutzwiller_complex_flag = any(term -> term.is_complex, data.gutzwiller_terms)
    jastrow_complex_flag = any(term -> term.is_complex, data.jastrow_terms)
    dh_complex_flag = data.doublon_holon_2site_complex || data.doublon_holon_4site_complex
    orbital_complex_flag = any(term -> term.is_complex, data.orbital_terms)
    all_complex_flag =
        modpara.complex_flag != 0 ||
        gutzwiller_complex_flag ||
        jastrow_complex_flag ||
        dh_complex_flag ||
        orbital_complex_flag

    # Initialize RBM parameters if FlagRBM > 0
    # Note: FlagRBM is determined by whether RBM terms exist
    flag_rbm = (
        length(data.charge_rbm_phys_layer_terms) > 0 ||
        length(data.spin_rbm_phys_layer_terms) > 0 ||
        length(data.general_rbm_phys_layer_terms) > 0 ||
        length(data.charge_rbm_hidden_layer_terms) > 0 ||
        length(data.spin_rbm_hidden_layer_terms) > 0 ||
        length(data.general_rbm_hidden_layer_terms) > 0 ||
        length(data.charge_rbm_phys_hidden_terms) > 0 ||
        length(data.spin_rbm_phys_hidden_terms) > 0 ||
        length(data.general_rbm_phys_hidden_terms) > 0
    )

    if flag_rbm
        n_proj = projection_layout(data).n_proj
        rbm_sections = (
            data.charge_rbm_phys_layer_terms,
            data.spin_rbm_phys_layer_terms,
            data.general_rbm_phys_layer_terms,
            data.charge_rbm_hidden_layer_terms,
            data.spin_rbm_hidden_layer_terms,
            data.general_rbm_hidden_layer_terms,
            data.charge_rbm_phys_hidden_terms,
            data.spin_rbm_phys_hidden_terms,
            data.general_rbm_phys_hidden_terms,
        )
        rbm_section_sizes = [_rbm_section_nparam(terms) for terms in rbm_sections]
        n_rbm = sum(rbm_section_sizes)
        rbm_values = fill(ComplexF64(0.0), n_rbm)

        if !all_complex_flag
            # Real RBM: RBM[i] = 0.01*(genrand_real2() - 0.5)/(double)Nneuron
            # In C: OptFlag[2*i+2*NProj] > 0
            # C implementation: Nneuron = Nneuron + NneuronCharge + NneuronSpin + NneuronGeneral
            # (see readdef.c:720-724)
            nneuron_total =
                modpara.nneuron +
                modpara.nneuron_charge +
                modpara.nneuron_spin +
                modpara.nneuron_general
            # If Nneuron is 0, use 1 to avoid division by zero (C implementation would also handle this)
            nneuron_divisor = nneuron_total > 0 ? Float64(nneuron_total) : 1.0

            for i in 0:(n_rbm - 1)
                # Check OptFlag: opt_flag_idx = 2*i + 2*NProj
                opt_flag_idx = 2 * i + 2 * n_proj + 1
                should_optimize = (
                    opt_flag_idx <= length(data.optimization_flags) &&
                    data.optimization_flags[opt_flag_idx]
                )
                if should_optimize
                    # C implementation: RBM[i] = 0.01*(genrand_real2() - 0.5)/(double)Nneuron
                    rbm_values[i + 1] = ComplexF64(0.01 * (rand(rng) - 0.5) / nneuron_divisor)
                else
                    rbm_values[i + 1] = ComplexF64(0.0)
                end
            end
        else
            # Complex RBM: RBM[i] = 1e-2*genrand_real2()*cexp(2.0*I*M_PI*genrand_real2())
            for i = 0:(n_rbm-1)
                opt_flag_idx = 2 * i + 2 * n_proj + 1
                should_optimize = (
                    opt_flag_idx <= length(data.optimization_flags) &&
                    data.optimization_flags[opt_flag_idx]
                )
                if should_optimize
                    r1 = rand(rng)
                    r2 = rand(rng)
                    rbm_values[i + 1] = ComplexF64(1e-2 * r1 * exp(2.0im * π * r2))
                else
                    rbm_values[i + 1] = ComplexF64(0.0)
                end
            end
        end

        # Scatter RBM values by (section_offset + local idx), matching C RBM layout.
        section_offset = 0
        for (terms, n_section) in zip(rbm_sections, rbm_section_sizes)
            for term in terms
                local_idx = term.idx
                if 0 <= local_idx < n_section
                    term.value = rbm_values[section_offset + local_idx + 1]
                else
                    term.value = ComplexF64(0.0)
                end
            end
            section_offset += n_section
        end
    end

    # Initialize Slater (Orbital) parameters
    # In C: OptFlag[2*i+2*NProj + 2*FlagRBM*NRBM] determines if Slater[i] should be optimized
    # C: AllComplexFlag = iComplexFlgGutzwiller + iComplexFlgJastrow + iComplexFlgDH2 + iComplexFlgDH4 + iComplexFlgOrbital
    # iComplexFlgOrbital comes from ComplexType header in orbital files (orbitalidx.def, orbitalidxpara.def)
    # If any term-group is complex, AllComplexFlag != 0.
    n_proj = projection_layout(data).n_proj
    n_rbm = count_rbm_parameters(data)

    # C implementation: for(i=0;i<NSlater;i++) { Slater[i] = 2*(genrand_real2()-0.5); }
    # NSlater = NOrbitalIdx = n_orbital_idx
    # Important: C implementation loops over NSlater (unique indices), not over all orbital_terms!
    # Julia's orbital_terms may have multiple entries with the same idx (for different i,j pairs)

    # Get NSlater (number of unique orbital parameters)
    # Since MVMCExpertModeParsers.jl now pre-offsets parallel orbital indices,
    # we can simply find the maximum idx across all orbital_terms.
    # For Kitaev FSZ: indices range from 0 to 275, so n_slater = 276.
    n_slater = 0
    if !isempty(data.orbital_terms)
        n_slater = maximum(t.idx for t in data.orbital_terms) + 1
    elseif modpara.n_orbital_idx > 0
        n_slater = modpara.n_orbital_idx
    end

    # First, generate Slater values for each unique idx (matching C's loop over NSlater)
    slater_values = Vector{ComplexF64}(undef, n_slater)

    if !all_complex_flag
        # Real Slater: Slater[i] = 2*(genrand_real2()-0.5) for OptFlag > 0
        # Uniform distribution [-1, 1)
        # In C: OptFlag[2*i+2*NProj + 2*FlagRBM*NRBM] > 0
        for i = 0:(n_slater-1)
            # Check OptFlag: opt_flag_idx = 2*i + 2*NProj + 2*FlagRBM*NRBM
            # Julia uses 1-based indexing, so add 1
            opt_flag_idx = 2 * i + 2 * n_proj + 2 * flag_rbm * n_rbm + 1
            # If optimization_flags is empty OR index out of range, default to initializing
            # (FSZ mode may have more Slater params than the flag array covers)
            should_optimize = (
                isempty(data.optimization_flags) ||
                opt_flag_idx > length(data.optimization_flags) ||
                data.optimization_flags[opt_flag_idx]
            )
            if should_optimize
                # OptFlag > 0 or out of range: initialize randomly
                # C: genrand_real2() returns [0, 1), so rand(rng) is equivalent
                slater_values[i+1] = ComplexF64(2.0 * (rand(rng) - 0.5))
            else
                # OptFlag == 0: set to 0.0
                slater_values[i+1] = ComplexF64(0.0)
            end
        end
    else
        # Complex Slater: Slater[i] = 2*(genrand_real2()-0.5) + 2*I*(genrand_real2()-0.5)
        # Then divide by sqrt(2.0)
        for i = 0:(n_slater-1)
            opt_flag_idx = 2 * i + 2 * n_proj + 2 * flag_rbm * n_rbm + 1
            # If optimization_flags is empty OR index out of range, default to initializing
            # (FSZ mode may have more Slater params than the flag array covers)
            should_optimize = (
                isempty(data.optimization_flags) ||
                opt_flag_idx > length(data.optimization_flags) ||
                data.optimization_flags[opt_flag_idx]
            )
            if should_optimize
                # OptFlag > 0 or out of range: initialize randomly
                real_part = 2.0 * (rand(rng) - 0.5)
                imag_part = 2.0 * (rand(rng) - 0.5)
                slater_values[i+1] = ComplexF64(real_part + imag_part * im) / sqrt(2.0)
            else
                # OptFlag == 0: set to 0.0
                slater_values[i+1] = ComplexF64(0.0)
            end
        end
    end

    # Now, apply Slater values to all orbital_terms based on their idx
    # Since indices are pre-offset in MVMCExpertModeParsers.jl, use term.idx directly
    for term in data.orbital_terms
        idx = term.idx
        if idx >= 0 && idx < n_slater
            term.value = slater_values[idx+1]
        end
    end

    # Initialize OptTrans parameters from ParaQPOptTrans.
    # In C: for(i=0;i<NOptTrans;i++) OptTrans[i] = ParaQPOptTrans[i];
    if !isempty(data.para_qp_opt_trans)
        data.opt_trans = copy(data.para_qp_opt_trans)
    end

    return data
end

"""
    sync_modified_parameter!(data::ExpertModeData)

Sync and modify variational parameters similar to C's SyncModifiedParameter().

This function:
- Rescales Slater parameters to ensure maximum amplitude <= D_AMP_MAX (4.0)
- Similar to C implementation's SyncModifiedParameter() after InitParameter()

# Arguments
- `data::ExpertModeData`: Expert Mode data structure to modify

# Notes
- Only rescales Slater (Orbital) parameters
- Finds maximum absolute value and scales all parameters proportionally
- Ensures max(|Slater[i]|) <= D_AMP_MAX
"""
function sync_modified_parameter!(data::ExpertModeData)
    # Find maximum absolute value of Slater parameters
    xmax = 0.0
    for term in data.orbital_terms
        abs_val = abs(term.value)
        if abs_val > xmax
            xmax = abs_val
        end
    end

    # C実装では、常に ratio = D_AmpMax/xmax を計算し、常にスケーリングする
    # (parameter.c:159-161行目を参照)
    # これにより、xmax < D_AMP_MAX の場合でもスケーリングされる
    if xmax > 0.0
        ratio = D_AMP_MAX / xmax
        for term in data.orbital_terms
            term.value *= ratio
        end
    end

    return data
end

"""
    initialize_parameters!(data::ExpertModeData; rng::AbstractRNG=SFMT19937RNG())

Complete parameter initialization workflow: InitParameter() followed by SyncModifiedParameter().

This is the equivalent of calling both init_parameter!() and sync_modified_parameter!() in sequence,
matching the C implementation's workflow.

# Arguments
- `data::ExpertModeData`: Expert Mode data structure to initialize.
- `rng::AbstractRNG`: Random number generator. The default constructs a
  fresh `SFMT19937RNG()` but **does not seed it**; using that default
  directly raises `ArgumentError: SFMT19937 not initialized` on the
  first random draw. Either pass in an already-seeded RNG (preferred,
  see the example below) or use the higher-level
  `MVMCOptimizers.run_para_opt_from_namelist`, which wires seeding for
  you.

# Example
```julia
using Random
using SFMT  # provides SFMT19937RNG

data = parse_expert_mode_files("namelist.def")

rng = SFMT19937RNG()
actual_seed = data.modpara.rnd_seed > 0 ? data.modpara.rnd_seed : 11272
Random.seed!(rng, actual_seed)

initialize_parameters!(data; rng=rng)
```
"""
function initialize_parameters!(data::ExpertModeData; rng::AbstractRNG = SFMT19937RNG())
    init_parameter!(data; rng = rng)
    sync_modified_parameter!(data)
    return data
end
