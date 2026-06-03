# 4. Physical-quantity measurement (VMCPhysCal) — experimental

> ⚠️ `VMCPhysCal` (`NVMCCalMode = 1`) is **experimental** in this release. Major
> components are implemented and exercised against simple reference
> data, but parts of the parsing and product-side accumulation are still
> in flux. Treat numbers as a sanity check, not as production output,
> until a later release (v0.3 or later).

## What is implemented

| Function | Status | Notes |
|----------|--------|-------|
| `vmc_phys_cal!` (entry point, `src/vmc_phys_cal.jl`) | ✅ wired | Drives sampling and Green-function accumulation. |
| 1-body Green function `<c†_i c_j>` | ✅ | Matches C reference for HeisenbergChain. |
| 2-body Green function (direct), `<c†_i c_j c†_k c_l>` | ✅ | Direct path only; the product path is not yet implemented. |
| Weighted average over QP weights (`weight_average_green_func!`) | ✅ | Same convention as C (Σ w_i G_i / Σ w_i). |
| Output to `zvo_cisajs.dat` and `zvo_cisajscktaltdc.dat` | ✅ | C-compatible per-row format. |

## Known limitations

- **Product-side two-body Green functions** (`cisajscktalt.def`)
  — the factored `<c†c>×<c†c>` accumulator (C path
  `CalculateGreenFunc_BF` and friends) is not yet ported.
- **Backflow correlation factor** — the `vmc_bf_*` entry points raise
  an error in this release. Inputs that activate Back Flow (`n_proj_bf > 0`, i.e.
  any `BackFlow*` keyword in `namelist.def`) are not supported; remove
  those keywords or fall back to the C reference at
  <https://github.com/issp-center-dev/mVMC>.
- **MPI parallelisation** — `reduce_counter!` is a no-op in this release, and
  `NSplitSize > 1` is rejected with an unsupported-MPI error until MPI
  support is implemented. `NSplitSize = 1` is the only supported setting.
- **`InterAllTerm` spin metadata** — when the input does not provide
  spin information, `vmc_main_cal.jl` substitutes default values (see
  the TODO at `src/vmc_main_cal.jl` near the InterAll loop).

## When to fall back to C-mVMC

For published physics results, the safest path in this release is:

1. Use Julia-mVMC for `VMCParaOpt` (parameter optimisation) — verified
   bit-level for the modes listed in
   [`03_optimization.md`](03_optimization.md).
2. Hand off the optimised `zqp_opt.dat` to C-mVMC for `VMCPhysCal`
   (physical-quantity measurement) until the product-side path lands.

The output formats of `zqp_opt.dat` are byte-compatible (same column
layout), so the hand-off requires no conversion script.
