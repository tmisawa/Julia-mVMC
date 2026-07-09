# 物理量計算

[English](../en/physics_calc.md)

`VMCPhysCal` (`NVMCCalMode = 1`) は v0.5.0 の対応 path で C reference と検証済みです。
file-based workflow では `run_phys_cal_from_namelist` を使います。

## Green functions

| Input | Output | Status |
|-------|--------|--------|
| `greenone.def` (`OneBodyG`) | `zvo_cisajs_*.dat` | C-referenced |
| `greentwo.def` (`TwoBodyG`) | `zvo_cisajscktalt_*.dat` | C-referenced |
| `greentwoex.def` (`TwoBodyGEx`) | `zvo_cisajscktaltex_*.dat` | 対応 real / complex / Kondo / DH / FSZ fixture で C-referenced |

factored/product two-body quantity は C と同じ weighted-average convention に従います。

```math
G_{\mathrm{ex}} =
\frac{\sum_q w_q\,G^{(1)}_q\,\overline{G^{(2)}_q}}
     {\sum_q w_q}.
```

## Full Lanczos PhysCal

`NLanczosMode = 1/2` output は sz-conserved `NSplitSize = 1` PhysCal run で
support されています。gate は次を比較します。

- `zvo_ls_out_001.dat`;
- `zvo_ls_qqqq_001.dat`;
- `zvo_ls_cisajs_001.dat`;
- `zvo_ls_cisajscktalt_001.dat`;
- `TwoBodyGEx` がある場合の `zvo_ls_cisajscktaltex_001.dat`。

## MPI PhysCal

`NSplitSize > 1` は sz-conserved normal-Green run (`NLanczosMode = 0`) で
support されています。MPI smoke gate は rank0 Green output、reduce-to-root path、
同じ chain count での `NSplitSize = 1` と `NSplitSize = 2` の self-consistency を確認します。

## 制限

以下は C-mVMC に fallback してください。

- BackFlow;
- FSZ/general-orbital Lanczos;
- `NSplitSize > 1` が必要な Lanczos PhysCal run;
- FSZ/general-orbital PhysCal split;
- split を伴う OptTrans-derived QP sector。

未対応の組合せは、静かに output を出すのではなく error になります。
