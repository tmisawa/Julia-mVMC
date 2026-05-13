"""
ModPara.def Parser

Parser for modpara.def files containing main simulation parameters.
"""

"""
    parse_modpara_def(filepath::String) -> ParseResult{ModParaParameters}

Parse modpara.def file from file path.
"""
function parse_modpara_def(filepath::String)::ParseResult{ModParaParameters}
    try
        content = read_def_file(filepath)
        return parse_modpara_content(content)
    catch e
        return ParseResult{ModParaParameters}(false, nothing, "Error reading file: $e", 0)
    end
end

"""
    parse_modpara_content(content::String) -> ParseResult{ModParaParameters}

Parse modpara.def content from string.
"""
function parse_modpara_content(content::String)::ParseResult{ModParaParameters}
    context = ParsingContext("modpara.def")
    params = ModParaParameters()

    lines = split(content, '\n')

    for (line_num, line) in enumerate(lines)
        context.line_number = line_num
        clean_line_str = clean_line(line)

        if isempty(clean_line_str)
            continue
        end

        tokens = split_def_line(clean_line_str)
        if length(tokens) < 2
            continue
        end

        # Handle "param = value" format
        if length(tokens) >= 3 && tokens[2] == "="
            param_name = tokens[1]
            param_value = tokens[3]
        else
            param_name = tokens[1]
            param_value = tokens[2]
        end

        try
            parse_modpara_parameter!(params, param_name, param_value, context)
        catch e
            push!(
                context.errors,
                "Line $line_num: Error parsing parameter '$param_name': $e",
            )
        end
    end

    success = length(context.errors) == 0
    return ParseResult{ModParaParameters}(
        success,
        success ? params : nothing,
        join(context.errors, "; "),
        context.line_number,
    )
end

"""
    parse_modpara_parameter!(params::ModParaParameters, name::String, value::String, context::ParsingContext)

Parse a single parameter and update the ModParaParameters struct.
"""
function parse_modpara_parameter!(
    params::ModParaParameters,
    name::String,
    value::String,
    context::ParsingContext,
)
    # Basic system parameters
    # Support both "NSite" and "Nsite" (C implementation accepts both)
    if name == "NSite" || name == "Nsite"
        params.nsite = safe_parse_int(value, DEFAULT_NSITE)
    elseif name == "NElec" || name == "Nelec"
        params.nelec = safe_parse_int(value, DEFAULT_NELEC)
    elseif name == "NLocSpin" || name == "NlocalSpin"
        params.nlocspin = safe_parse_int(value, DEFAULT_NLOCSPIN)
    elseif name == "NCond" || name == "Ncond"
        params.ncond = safe_parse_int(value, -1)

        # Calculation modes
        # Support both "NVMCCalMode" and "VMCCalMode" (C implementation uses "NVMCCalMode")
    elseif name == "VMCCalMode" || name == "NVMCCalMode"
        params.vmc_calc_mode = safe_parse_int(value, DEFAULT_VMC_CALC_MODE)
    elseif name == "LanczosMode" || name == "NLanczosMode"
        params.lanczos_mode = safe_parse_int(value, DEFAULT_LANCZOS_MODE)

        # VMC parameters
    elseif name == "NSROptItrStep"
        params.nsr_opt_itr_step = safe_parse_int(value, DEFAULT_NSR_OPT_ITR_STEP)
    elseif name == "NSROptItrSmp"
        params.nsr_opt_itr_smp = safe_parse_int(value, DEFAULT_NSR_OPT_ITR_SMP)
    elseif name == "NSROptFixSmp"
        params.nsr_opt_fix_smp = safe_parse_int(value, 0)
    elseif name == "NVMCWarmUp"
        params.nvmc_warmup = safe_parse_int(value, DEFAULT_NVMC_WARMUP)
    elseif name == "NVMCInterval"
        params.nvmc_interval = safe_parse_int(value, DEFAULT_NVMC_INTERVAL)
    elseif name == "NVMCSample"
        params.nvmc_sample = safe_parse_int(value, DEFAULT_NVMC_SAMPLE)

        # SR parameters
    elseif name == "DSROptRedCut"
        params.dsr_opt_red_cut = safe_parse_float(value, DEFAULT_DSR_OPT_RED_CUT)
    elseif name == "DSROptStaDel"
        params.dsr_opt_sta_del = safe_parse_float(value, DEFAULT_DSR_OPT_STA_DEL)
    elseif name == "DSROptStepDt"
        params.dsr_opt_step_dt = safe_parse_float(value, DEFAULT_DSR_OPT_STEP_DT)
    elseif name == "DSROptCGTol"
        params.dsr_opt_cg_tol = safe_parse_float(value, DEFAULT_DSR_OPT_CG_TOL)
    elseif name == "NSROptCGMaxIter"
        params.nsr_opt_cg_max_iter = safe_parse_int(value, DEFAULT_NSR_OPT_CG_MAX_ITER)

        # Random number generation
    elseif name == "RndSeed"
        params.rnd_seed = safe_parse_int(value, DEFAULT_RND_SEED)
    elseif name == "NSplitSize"
        params.nsplit_size = safe_parse_int(value, DEFAULT_NSPLIT_SIZE)

        # Quantum projection
    elseif name == "NSPGaussLeg"
        params.nsp_gauss_leg = safe_parse_int(value, DEFAULT_NSP_GAUSS_LEG)
    elseif name == "NSPStot"
        params.nsp_stot = safe_parse_int(value, DEFAULT_NSP_STOT)
    elseif name == "NMPTrans"
        params.nmp_trans = safe_parse_int(value, DEFAULT_NMP_TRANS)
    elseif name == "2Sz"
        params.two_sz = safe_parse_int(value, DEFAULT_TWO_SZ)

        # Data output
    elseif name == "NDataIdxStart"
        params.n_data_idx_start = safe_parse_int(value, DEFAULT_N_DATA_IDX_START)
    elseif name == "NDataQtySmp"
        params.n_data_qty_smp = safe_parse_int(value, DEFAULT_N_DATA_QTY_SMP)
    elseif name == "CDataFileHead"
        params.c_data_file_head = value
    elseif name == "CParaFileHead"
        params.c_para_file_head = value

        # File control
    elseif name == "NFileFlushInterval"
        params.n_file_flush_interval = safe_parse_int(value, DEFAULT_N_FILE_FLUSH_INTERVAL)
    elseif name == "NStore"
        # NStoreO: 0 = normal, !=0 = store O samples for BLAS calculation
        params.nstore_o = safe_parse_int(value, 1)  # Default is 1 (from C SetDefaultValuesModPara)
    elseif name == "NSRCG"
        # NSRCG: 0 = direct solver (LAPACK), !=0 = CG solver
        params.nsrcg = safe_parse_int(value, 0)  # Default is 0 (from C SetDefaultValuesModPara)

    # Complex flag
    elseif name == "ComplexType"
        params.complex_flag = safe_parse_int(value, DEFAULT_COMPLEX_FLAG)

        # RBM parameters
    elseif name == "Nneuron"
        params.nneuron = safe_parse_int(value, DEFAULT_NNEURON)
    elseif name == "NneuronGeneral"
        params.nneuron_general = safe_parse_int(value, DEFAULT_NNEURON)
    elseif name == "NneuronCharge"
        params.nneuron_charge = safe_parse_int(value, DEFAULT_NNEURON)
    elseif name == "NneuronSpin"
        params.nneuron_spin = safe_parse_int(value, DEFAULT_NNEURON)
    elseif name == "NBlockSize_RBMRatio"
        params.nblock_size_rbm_ratio = safe_parse_int(value, DEFAULT_NBLOCK_SIZE_RBM_RATIO)

        # Lanczos parameters
    elseif name == "NOneBodyG"
        params.n_one_body_g = safe_parse_int(value, DEFAULT_N_ONE_BODY_G)
    elseif name == "NTwoBodyG"
        params.n_two_body_g = safe_parse_int(value, DEFAULT_N_TWO_BODY_G)
    elseif name == "NTwoBodyGEx"
        params.n_two_body_g_ex = safe_parse_int(value, DEFAULT_N_TWO_BODY_G_EX)

        # Exchange update
    elseif name == "NExUpdatePath"
        params.nex_update_path = safe_parse_int(value, DEFAULT_NEX_UPDATE_PATH)

    else
        push!(context.warnings, "Unknown parameter: $name")
    end
end
