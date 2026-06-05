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
- `NSplitSize > 1`: MPI parallelization is not supported in Julia-mVMC.
  Thrown as `error(...)` / `ErrorException`, matching the codebase's
  unsupported-*feature* convention (the BackFlow stubs).
- `NLanczosMode > 0`: full Lanczos is not supported (only the step-0
  comparison matches C; the post-Lanczos eigenvector / overlap pipeline is
  not ported — see `docs/manual/05_compatibility.md`). Rejected early because
  C also builds the indirect one-body Green list when `NLanczosMode > 1`
  (`readdef.c`), which the Julia canonical-list path does not reproduce.
  Thrown as `error(...)` / `ErrorException` (unsupported-feature convention).

`NSplitSize = 1` (single process) and `NLanczosMode = 0` are the only
supported settings. The fields are parsed for input-format fidelity but have
no runtime effect at those values.
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
            "NSplitSize > 1 is not supported: MPI parallelization is not supported " *
            "in Julia-mVMC (got NSplitSize = $(modpara.nsplit_size)). " *
            "Set NSplitSize = 1, or fall back to the C reference at " *
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
