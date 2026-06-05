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

`NSplitSize = 1` (single process) is the only supported setting. The field
is parsed for input-format fidelity but has no runtime effect when `== 1`.
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
    return nothing
end
