# Physics Calculation

[日本語](../ja/physics_calc.md)

`VMCPhysCal` (`NVMCCalMode = 1`) is C-verified for supported paths in v0.5.0.
Use `run_phys_cal_from_namelist` for file-based workflows.

## Green functions

| Input | Output | Status |
|-------|--------|--------|
| `greenone.def` (`OneBodyG`) | `zvo_cisajs_*.dat` | C-referenced |
| `greentwo.def` (`TwoBodyG`) | `zvo_cisajscktalt_*.dat` | C-referenced |
| `greentwoex.def` (`TwoBodyGEx`) | `zvo_cisajscktaltex_*.dat` | C-referenced for supported real, complex, Kondo, DH, and FSZ fixtures |

The factored/product two-body quantity follows the C weighted-average
convention:

```math
G_{\mathrm{ex}} =
\frac{\sum_q w_q\,G^{(1)}_q\,\overline{G^{(2)}_q}}
     {\sum_q w_q}.
```

## Full Lanczos PhysCal

`NLanczosMode = 1/2` output is supported for sz-conserved `NSplitSize = 1`
PhysCal runs. The gate compares:

- `zvo_ls_out_001.dat`;
- `zvo_ls_qqqq_001.dat`;
- `zvo_ls_cisajs_001.dat`;
- `zvo_ls_cisajscktalt_001.dat`;
- `zvo_ls_cisajscktaltex_001.dat` when `TwoBodyGEx` is present.

## MPI PhysCal

`NSplitSize > 1` is supported for sz-conserved normal-Green runs
(`NLanczosMode = 0`). The MPI smoke gate checks rank0 Green output,
reduce-to-root paths, and `NSplitSize = 1` versus `NSplitSize = 2`
self-consistency for the same chain count.

## Limitations

Fall back to C-mVMC for:

- BackFlow;
- FSZ/general-orbital Lanczos;
- Lanczos PhysCal runs requiring `NSplitSize > 1`;
- FSZ/general-orbital PhysCal split;
- OptTrans-derived QP sectors with split.

Unsupported combinations raise errors rather than silently producing output.
