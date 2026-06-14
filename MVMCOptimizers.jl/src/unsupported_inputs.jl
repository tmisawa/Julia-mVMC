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

`NSplitSize = 1` and `NLanczosMode = 0` are the only supported settings. The
fields are parsed for input-format fidelity; `NSplitSize = 1` supports both
serial execution and v0.4 MPI sample-parallel execution.
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
    return nothing
end

"""
    validate_supported_para_opt_parallel_modpara(ctx, modpara)

Reject parameter-optimization settings whose MPI path is not C-compatible yet.

`NSRCG != 0` selects the CG SR solver. C's CG `operate_by_S` broadcasts the
search vector and allreduces the sampled matrix-vector product inside each CG
iteration. Julia's current CG solver is still rank-local, so under MPI it would
silently solve a non-C-equivalent system. Serial CG remains available.
"""
function validate_supported_para_opt_parallel_modpara(
    ctx::ParallelContext,
    modpara::ModParaParameters,
)
    if ctx.is_mpi && modpara.nsrcg != 0
        error(
            "NSRCG != 0 is not supported under MPI in Julia-mVMC v0.4: " *
            "the CG SR solver does not yet implement C-compatible " *
            "operate_by_S broadcast/allreduce. Set NSRCG = 0 for MPI " *
            "runs, run this input without MPI, or fall back to the C " *
            "reference at https://github.com/issp-center-dev/mVMC.",
        )
    end
    return nothing
end
