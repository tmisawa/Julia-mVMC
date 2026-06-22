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
- `NSplitSize > 1`: C's grouped MPI/QP split by `NSplitSize` is not yet
  implemented in Julia-mVMC. Multi-rank MPI sample parallel runs are supported
  only with `NSplitSize = 1`. Thrown as `error(...)` / `ErrorException`,
  matching the codebase's unsupported-*feature* convention (the BackFlow stubs).
- `NLanczosMode > 0`: full Lanczos is not supported (only the step-0
  comparison matches C; the post-Lanczos eigenvector / overlap pipeline is
  not ported — see `docs/manual/05_compatibility.md`). Rejected early because
  C also builds the indirect one-body Green list when `NLanczosMode > 1`
  (`readdef.c`), which the Julia canonical-list path does not reproduce.
  Thrown as `error(...)` / `ErrorException` (unsupported-feature convention).
- `NSRCG >= 2`: Julia currently supports the direct SR solver (`NSRCG = 0`)
  and the standard SR-CG solver (`NSRCG = 1`) only.
- `useDiagScale != 0`: C's preconditioned CG mode is not ported yet.
- `RescaleSmat != 0`: C's S-matrix rescaling mode is not ported yet.

`NSplitSize = 1`, `NLanczosMode = 0`, and `NSRCG <= 1` are the supported
settings. The fields are parsed for input-format fidelity; `NSplitSize = 1`
supports both serial execution and MPI sample-parallel execution.
"""
function validate_supported_modpara(modpara::ModParaParameters)
    if modpara.nsplit_size < 1
        throw(
            ArgumentError(
                "NSplitSize must be >= 1; got NSplitSize = $(modpara.nsplit_size).",
            ),
        )
    elseif modpara.nsplit_size > 1
        error(
            "NSplitSize > 1 is not supported: grouped MPI/QP splitting by " *
            "NSplitSize is not implemented in Julia-mVMC " *
            "(got NSplitSize = $(modpara.nsplit_size)). " *
            "Set NSplitSize = 1 for serial or supported multi-rank MPI " *
            "sample-parallel runs, or fall back to the C reference at " *
            "https://github.com/issp-center-dev/mVMC.",
        )
    end
    if modpara.lanczos_mode > 0
        error(
            "NLanczosMode > 0 is not supported: full Lanczos is not implemented " *
            "in Julia-mVMC (got NLanczosMode = $(modpara.lanczos_mode)). " *
            "Only the step-0 (variational) quantities match the C reference. " *
            "Set NLanczosMode = 0, or fall back to the C reference at " *
            "https://github.com/issp-center-dev/mVMC.",
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
    return nothing
end
