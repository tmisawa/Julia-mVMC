# 4. Physical-quantity measurement (VMCPhysCal) — experimental

> ⚠️ `VMCPhysCal` (`NVMCCalMode = 1`) is **experimental** in this release. Major
> components are implemented and exercised against C reference data — the
> one-body, direct two-body, and factored/product two-body Green functions now
> all match C-mVMC to the per-quantity gate tolerance (see
> [`../../test/integration/phys_cal_equivalent.jl`](../../test/integration/phys_cal_equivalent.jl)).
> Treat numbers as a sanity check, not as production output, until a later
> release, and note the FSZ caveat on the factored path below.

## What is implemented

| Function | Status | Notes |
|----------|--------|-------|
| `vmc_phys_cal!` (entry point, `src/vmc_phys_cal.jl`) | ✅ wired | Drives sampling and Green-function accumulation. |
| 1-body Green function `<c†_i c_j>` | ✅ | Matches C reference for HeisenbergChain. |
| 2-body Green function (direct), `<c†_i c_j c†_k c_l>` (`TwoBodyG`/`greentwo.def`) | ✅ | Output `zvo_cisajscktalt_*.dat`. |
| 2-body Green function (factored/product), `<c†_i c_j>·<c†_k c_l>` (`TwoBodyGEx`/`greentwoex.def`) | ✅ (non-FSZ) | Output `zvo_cisajscktaltex_*.dat`. Matches C to the gate tolerance for real, complex (cmp) and Kondo systems. FSZ is rejected at runtime (see below). |
| Doublon-holon projection (`DH2`/`DH4`) | ✅ | DH-present PhysCal fixture is gated against C reference data. |
| Weighted average over QP weights (`weight_average_green_func!`) | ✅ | Same convention as C (Σ w_i G_i / Σ w_i). |
| Output to `zvo_cisajs_*.dat`, `zvo_cisajscktalt_*.dat`, `zvo_cisajscktaltex_*.dat` | ✅ | C-compatible per-row / value-only format. |

## Known limitations

- **Factored two-body Green under FSZ** — the product-side `TwoBodyGEx`
  (`greentwoex.def` → `zvo_cisajscktaltex_*.dat`) path is supported for the
  spin-conserving (`mode = :real` / `:cmp`) sector and is rejected at runtime
  (`validate_factored_green_supported`) when combined with the FSZ generalised
  orbital. Use the C reference for factored Green under FSZ.
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
2. For `VMCPhysCal` (physical-quantity measurement), the one-body, direct and
   factored two-body Green functions are gated against C references
   ([`phys_cal_equivalent.jl`](../../test/integration/phys_cal_equivalent.jl)) and
   can be run via
   [`run_phys_cal_from_namelist`](../../MVMCOptimizers.jl/src/run_phys_cal_from_namelist.jl);
   still fall back to C-mVMC for FSZ factored Green, Backflow, Lanczos
   (`NLanczosMode > 0`), or MPI (`NSplitSize > 1`).

The output formats of `zqp_opt.dat` are byte-compatible (same column
layout), so the hand-off requires no conversion script.
