# 出力ファイル

[English](../en/output_files.md)

Julia-mVMC は、対応済み範囲で C-mVMC の output file naming に従います。

## Optimization output

| File | Contents |
|------|----------|
| `zvo_out.dat` | step ごとの energy, variance, spin diagnostic |
| `zvo_var.dat` | step ごとの parameter snapshot |
| `zqp_opt.dat` | final optimized parameters |
| `zqp_gutzwiller_opt.dat` | optimized Gutzwiller block |
| `zqp_jastrow_opt.dat` | optimized Jastrow block |
| `zqp_orbital_opt.dat` and related files | optimized orbital blocks |
| `zvo_CalcTimer.dat` | optional C-compatible timer output |

main energy row は次の列を持ちます。

```text
real(<H>) imag(<H>) <H^2> variance <Sz> <Sz^2>
```

variance は cancellation を含む式で計算されます。

```math
\mathrm{variance} =
\frac{\langle H^2\rangle - \langle H\rangle^2}
     {\langle H\rangle^2}.
```

そのため ``\langle H\rangle`` の小さな BLAS summation-order difference が
variance column では増幅されることがあります。

## PhysCal Green output

| File | Input |
|------|-------|
| `zvo_cisajs_001.dat` | `greenone.def` |
| `zvo_cisajscktalt_001.dat` | `greentwo.def` |
| `zvo_cisajscktaltex_001.dat` | `greentwoex.def` |

one-sample PhysCal run では、通常 suffix は `_001.dat` です。

## Full Lanczos output

| File | Contents |
|------|----------|
| `zvo_ls_out_001.dat` | Lanczos energy, norm, alpha |
| `zvo_ls_qqqq_001.dat` | QQQQ moments |
| `zvo_ls_cisajs_001.dat` | one-body Lanczos Green output |
| `zvo_ls_cisajscktalt_001.dat` | direct two-body Lanczos Green output |
| `zvo_ls_cisajscktaltex_001.dat` | requested 時の factored two-body Lanczos Green output |

## MPI での rank0 output

MPI run では user-facing file は rank0 だけが書きます。non-root rank は sampling と
reduction に参加しますが、output file を重複して書きません。
