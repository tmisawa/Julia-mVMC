"""
Runtime compatibility contract for unsupported ModPara inputs.

The parser (`MVMCExpertModeParsers.jl`) still parses every ModPara field,
including ones the Julia runtime does not yet implement. The rejection of
unsupported inputs lives here, in the optimizer/runtime layer, so that the
parser stays a faithful reader of the input format and the runtime owns the
"what is actually supported" contract.
"""

"""
    validate_supported_modpara(modpara::ModParaParameters)

Reject ModPara inputs that the Julia runtime does not yet support.

Non-mutating. Called at the independent runtime entry points
(`vmc_para_opt!`, `vmc_phys_cal!`) before any optimization or physical
calculation proceeds, so unsupported runs fail fast with a clear message.

Currently rejected:

- `NSplitSize < 1` (i.e. `0` or negative): invalid value. NSplitSize is a
  process-split count and must be at least 1. Thrown as `ArgumentError`,
  matching the codebase convention for invalid parameter *values*.
- `NLanczosMode < 0` or `NLanczosMode > 2`: invalid / unknown values. Entry
  point validators below define the currently supported subset.
- `NSRCG >= 2`: Julia currently supports the direct SR solver (`NSRCG = 0`)
  and the standard SR-CG solver (`NSRCG = 1`) only.
- `useDiagScale != 0`: C's preconditioned CG mode is not ported yet.
- `RescaleSmat != 0`: C's S-matrix rescaling mode is not ported yet.

`NSplitSize >= 1`, `NLanczosMode = 0/1/2`, and `NSRCG <= 1` are the globally
valid settings. Entry-point-specific validators below reject combinations that
are still unsupported for parameter optimization or physical measurement.
"""
function validate_supported_modpara(modpara::ModParaParameters)
    if modpara.nsplit_size < 1
        throw(
            ArgumentError(
                "NSplitSize must be >= 1; got NSplitSize = $(modpara.nsplit_size).",
            ),
        )
    end
    if modpara.lanczos_mode < 0 || modpara.lanczos_mode > 2
        error(
            "NLanczosMode must be 0, 1, or 2; got NLanczosMode = " *
            "$(modpara.lanczos_mode).",
        )
    end
    if modpara.nsrcg >= 2
        error(
            "NSRCG >= 2 is not supported: Julia-mVMC currently implements " *
            "the direct SR solver (NSRCG = 0) and standard SR-CG solver " *
            "(NSRCG = 1) only (got NSRCG = $(modpara.nsrcg)). " *
            "Set NSRCG = 0 or 1, or fall back to the C reference at " *
            "https://github.com/issp-center-dev/mVMC.",
        )
    end
    if modpara.use_diag_scale != 0
        error(
            "useDiagScale != 0 is not supported: C's preconditioned CG " *
            "mode is not implemented in Julia-mVMC " *
            "(got useDiagScale = $(modpara.use_diag_scale)). " *
            "Set useDiagScale = 0, or fall back to the C reference at " *
            "https://github.com/issp-center-dev/mVMC.",
        )
    end
    if modpara.rescale_smat != 0
        error(
            "RescaleSmat != 0 is not supported: C's S-matrix rescaling " *
            "mode is not implemented in Julia-mVMC " *
            "(got RescaleSmat = $(modpara.rescale_smat)). " *
            "Set RescaleSmat = 0, or fall back to the C reference at " *
            "https://github.com/issp-center-dev/mVMC.",
        )
    end
    return nothing
end

"""
    validate_supported_para_opt_modpara(modpara)

Validate parameter-optimization-only ModPara combinations.

`NSplitSize > 1` is supported for the direct SR solver (`NSRCG = 0`) by
splitting VMC samples inside each comm1 group. SR-CG remains restricted to
`NSplitSize = 1` until its sampled matrix-vector product is made comm1-aware.
"""
function validate_supported_para_opt_modpara(modpara::ModParaParameters)
    if modpara.lanczos_mode > 0
        error(
            "NLanczosMode > 0 is not supported for parameter optimization in " *
            "Julia-mVMC (got NLanczosMode = $(modpara.lanczos_mode)). " *
            "Use NLanczosMode = 0 for ParaOpt, or run physical measurement " *
            "with NLanczosMode = 1 or 2.",
        )
    end
    if modpara.nsplit_size > 1 && modpara.nsrcg != 0
        error(
            "NSplitSize > 1 with SR-CG is not supported: sample splitting by " *
            "NSplitSize is currently implemented only for the direct SR solver " *
            "(NSRCG = 0), got NSplitSize = $(modpara.nsplit_size), " *
            "NSRCG = $(modpara.nsrcg). Set NSRCG = 0 or NSplitSize = 1.",
        )
    end
    return nothing
end

function _has_opt_trans_derived_qp_sectors(data::ExpertModeData)
    n_qp_opt_trans = max(1, data.n_qp_opt_trans)
    return n_qp_opt_trans > 1 || length(data.opt_trans) > 1 ||
           length(data.qp_opt_trans) > 1
end

"""
    validate_supported_para_opt_data(data)

Validate parameter-optimization settings that require parsed data, not just
ModPara. `NSplitSize > 1` supports standard-projection `NQPFull > 1` when
`NQPOptTrans == 1` on sz-conserved paths; OptTrans-derived QP sectors and FSZ
standard-projection `NQPFull > 1` stay rejected with split.
"""
function validate_supported_para_opt_data(data::ExpertModeData)
    if data.modpara.nsplit_size > 1
        if data.i_flg_orbital_general != 0 &&
           (data.modpara.nsp_gauss_leg > 1 || abs(data.modpara.nmp_trans) > 1)
            error(
                "NSplitSize > 1 with FSZ standard-projection NQPFull > 1 " *
                "is not supported: " *
                "split standard-projection sampling currently supports " *
                "sz-conserved NQPFull > 1 paths only, got " *
                "NSplitSize = $(data.modpara.nsplit_size), " *
                "NSPGaussLeg = $(data.modpara.nsp_gauss_leg), " *
                "NMPTrans = $(data.modpara.nmp_trans). Use NSplitSize = 1 " *
                "or NSPGaussLeg = 1 and |NMPTrans| = 1 for FSZ inputs.",
            )
        end
        if _has_opt_trans_derived_qp_sectors(data)
            error(
                "NSplitSize > 1 with NQPOptTrans > 1 / OptTrans is not " *
                "supported: grouped QP-split sampling currently supports " *
                "standard-projection NQPFull only (NQPOptTrans = 1), got " *
                "NSplitSize = $(data.modpara.nsplit_size), " *
                "NQPOptTrans = $(data.n_qp_opt_trans). Use NSplitSize = 1 " *
                "for OptTrans-derived QP sectors.",
            )
        end
    end
    return nothing
end

"""
    validate_supported_phys_cal_modpara(modpara)

Validate physical-measurement-only ModPara combinations.

Physical-measurement support for `NSplitSize > 1` is scoped by parsed data:
sz-conserved normal Green paths are supported, while FSZ/general-orbital,
Full Lanczos, and OptTrans-derived split combinations are rejected by
`validate_supported_phys_cal_data`.
"""
function validate_supported_phys_cal_modpara(modpara::ModParaParameters)
    return nothing
end

"""
    validate_supported_phys_cal_data(data)

Validate physical-measurement settings that require parsed data.
`NSplitSize > 1` is currently limited to sz-conserved normal Green paths.
Full Lanczos is limited to `NSplitSize = 1` on the sz-conserved path;
FSZ/general-orbital Lanczos remains separate work.
"""
function validate_supported_phys_cal_data(data::ExpertModeData)
    if data.modpara.nsplit_size > 1
        if data.i_flg_orbital_general != 0
            error(
                "NSplitSize > 1 is not supported for FSZ / general-orbital " *
                "PhysCal in Julia-mVMC (got NSplitSize = " *
                "$(data.modpara.nsplit_size), i_flg_orbital_general = " *
                "$(data.i_flg_orbital_general)). Current PhysCal split " *
                "support is limited to sz-conserved normal Green paths.",
            )
        end
        if data.modpara.lanczos_mode > 0
            error(
                "NSplitSize > 1 with NLanczosMode > 0 is not supported for " *
                "PhysCal in Julia-mVMC (got NSplitSize = " *
                "$(data.modpara.nsplit_size), NLanczosMode = " *
                "$(data.modpara.lanczos_mode)). Current PhysCal split " *
                "support is limited to sz-conserved normal Green paths.",
            )
        end
        if _has_opt_trans_derived_qp_sectors(data)
            error(
                "NSplitSize > 1 with NQPOptTrans > 1 / OptTrans is not " *
                "supported for PhysCal in Julia-mVMC: grouped PhysCal split " *
                "currently supports standard-projection QP sectors only, got " *
                "NSplitSize = $(data.modpara.nsplit_size), " *
                "NQPOptTrans = $(data.n_qp_opt_trans). Use NSplitSize = 1 " *
                "for OptTrans-derived QP sectors.",
            )
        end
    end
    if data.modpara.lanczos_mode > 0 && data.i_flg_orbital_general != 0
        error(
            "NLanczosMode > 0 is not supported in FSZ / general-orbital mode " *
            "(i_flg_orbital_general = $(data.i_flg_orbital_general)). " *
            "Use sz-conserved inputs for Full Lanczos, or set NLanczosMode = 0.",
        )
    end
    if data.modpara.lanczos_mode > 0
        for term in data.transfer_terms
            if term.spin1 != term.spin2
                error(
                    "NLanczosMode > 0 does not support spin-flip Transfer terms " *
                    "in Julia-mVMC (site1=$(term.site1), spin1=$(term.spin1), " *
                    "site2=$(term.site2), spin2=$(term.spin2)).",
                )
            end
        end
        for term in data.inter_all_terms
            if term.spin0 != term.spin1 || term.spin2 != term.spin3
                error(
                    "NLanczosMode > 0 does not support spin-changing InterAll " *
                    "terms in Julia-mVMC.",
                )
            end
        end
    end
    return nothing
end

"""
    validate_supported_para_opt_parallel_modpara(ctx, modpara)

Validate parameter-optimization settings whose MPI path has extra constraints.

`NSRCG = 1` is supported under MPI for `NSplitSize = 1`: Julia's CG
`operate_by_s!` follows C's comm0 search-vector broadcast and sampled
matrix-vector allreduce. Unsupported global ModPara combinations are rejected
by `validate_supported_modpara`.
"""
function validate_supported_para_opt_parallel_modpara(
    ctx::ParallelContext,
    modpara::ModParaParameters,
)
    validate_supported_para_opt_modpara(modpara)
    return nothing
end
