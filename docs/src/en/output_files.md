# Output Files

[日本語](../ja/output_files.md)

Julia-mVMC follows C-mVMC output naming where supported.

## Optimization output

| File | Contents |
|------|----------|
| `zvo_out.dat` | per-step energy, variance, and spin diagnostics |
| `zvo_var.dat` | per-step parameter snapshot |
| `zqp_opt.dat` | final optimized parameters |
| `zqp_gutzwiller_opt.dat` | optimized Gutzwiller block |
| `zqp_jastrow_opt.dat` | optimized Jastrow block |
| `zqp_orbital_opt.dat` and related files | optimized orbital blocks |
| `zvo_CalcTimer.dat` | optional C-compatible timer output |

The main energy row contains:

```text
real(<H>) imag(<H>) <H^2> variance <Sz> <Sz^2>
```

The variance is computed by cancellation:

```math
\mathrm{variance} =
\frac{\langle H^2\rangle - \langle H\rangle^2}
     {\langle H\rangle^2}.
```

Small BLAS summation-order differences in ``\langle H\rangle`` can therefore be
amplified in the variance column.

## PhysCal Green output

| File | Input |
|------|-------|
| `zvo_cisajs_001.dat` | `greenone.def` |
| `zvo_cisajscktalt_001.dat` | `greentwo.def` |
| `zvo_cisajscktaltex_001.dat` | `greentwoex.def` |

For a one-sample PhysCal run, the suffix is usually `_001.dat`.

## Full Lanczos output

| File | Contents |
|------|----------|
| `zvo_ls_out_001.dat` | Lanczos energy, norm, and alpha |
| `zvo_ls_qqqq_001.dat` | QQQQ moments |
| `zvo_ls_cisajs_001.dat` | one-body Lanczos Green output |
| `zvo_ls_cisajscktalt_001.dat` | direct two-body Lanczos Green output |
| `zvo_ls_cisajscktaltex_001.dat` | factored two-body Lanczos Green output when requested |

## Rank0 output under MPI

MPI runs write user-facing files from rank0 only. Non-root ranks participate in
sampling and reductions but do not duplicate output files.
