# 4. Physical-quantity measurement (VMCPhysCal) — verified for supported paths

> `VMCPhysCal` (`NVMCCalMode = 1`) is **C-verified for supported paths** in this
> release. The
> one-body, direct two-body, and factored/product two-body Green functions now
> all match C-mVMC to the per-quantity gate tolerance (see
> [`../../test/integration/phys_cal_equivalent.jl`](../../test/integration/phys_cal_equivalent.jl)).
> `NLanczosMode = 1/2` output is also C-referenced on the sz-conserved
> `NSplitSize = 1` path (see
> [`../../test/integration/lanczos_equivalent.jl`](../../test/integration/lanczos_equivalent.jl)).
> Unsupported combinations still raise errors; see the limitations below.

## What is implemented

| Function | Status | Notes |
|----------|--------|-------|
| `vmc_phys_cal!` (entry point, `src/vmc_phys_cal.jl`) | ✅ wired | Drives sampling and Green-function accumulation. |
| 1-body Green function `<c†_i c_j>` | ✅ | Matches C reference for HeisenbergChain. |
| 2-body Green function (direct), `<c†_i c_j c†_k c_l>` (`TwoBodyG`/`greentwo.def`) | ✅ | Output `zvo_cisajscktalt_*.dat`. |
| 2-body Green function (factored/product), `<c†_i c_j>·<c†_k c_l>` (`TwoBodyGEx`/`greentwoex.def`) | ✅ | Output `zvo_cisajscktaltex_*.dat`. Matches C to the gate tolerance for real, complex (cmp), Kondo, and FSZ/general-orbital fixtures. |
| Doublon-holon projection (`DH2`/`DH4`) | ✅ | DH-present PhysCal fixture is gated against C reference data. |
| Weighted average over QP weights (`weight_average_green_func!`) | ✅ | Same convention as C (Σ w_i G_i / Σ w_i). |
| Output to `zvo_cisajs_*.dat`, `zvo_cisajscktalt_*.dat`, `zvo_cisajscktaltex_*.dat` | ✅ | C-compatible per-row / value-only format. |

## Known limitations

- **Backflow correlation factor** — the `vmc_bf_*` entry points raise
  an error in this release. Inputs that activate Back Flow (`n_proj_bf > 0`, i.e.
  any `BackFlow*` keyword in `namelist.def`) are not supported; remove
  those keywords or fall back to the C reference at
  <https://github.com/issp-center-dev/mVMC>.
- **MPI parallelisation** — physical measurement supports multi-rank execution
  through MPI.jl-compatible launchers. `NSplitSize > 1` is supported for
  sz-conserved normal-Green PhysCal runs (`NLanczosMode = 0`), including
  standard-projection `NQPFull > 1` when `NQPOptTrans = 1`. FSZ/general-orbital
  PhysCal split, PhysCal Lanczos split, and OptTrans-derived QP sectors with
  split remain rejected at runtime.
  The CI MPI smoke gate covers rank0 Green output, reduce-to-root paths, and
  the same-chain-count `NSplitSize = 1` vs `NSplitSize = 2` self-consistency
  path, but it is not a site performance benchmark.
- **`InterAllTerm` spin metadata** — when the input does not provide
  spin information, `vmc_main_cal.jl` substitutes default values (see
  the TODO at `src/vmc_main_cal.jl` near the InterAll loop).

## When to fall back to C-mVMC

For published physics results, the recommended path in this release is:

1. Use Julia-mVMC for `VMCParaOpt` (parameter optimisation) — verified
   bit-level for the modes listed in
   [`03_optimization.md`](03_optimization.md).
2. For `VMCPhysCal` (physical-quantity measurement), use the supported paths:
   the one-body, direct and factored two-body Green functions are gated against C references
   ([`phys_cal_equivalent.jl`](../../test/integration/phys_cal_equivalent.jl)) and
   can be run via
   [`run_phys_cal_from_namelist`](../../MVMCOptimizers.jl/src/run_phys_cal_from_namelist.jl).
   `NLanczosMode = 1/2` also writes Full-Lanczos files on the sz-conserved
   `NSplitSize = 1` path with C-reference gates. Still fall back to C-mVMC for
   Backflow, FSZ/general-orbital Lanczos, or Lanczos PhysCal runs that need
   `NSplitSize > 1`.

The output formats of `zqp_opt.dat` are byte-compatible (same column
layout), so the hand-off requires no conversion script.
