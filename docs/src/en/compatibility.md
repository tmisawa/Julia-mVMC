# Compatibility and Verification

[日本語](../ja/compatibility.md)

Julia-mVMC treats C-mVMC compatibility as the default correctness target.

## Function correspondence

| C function | Julia function | Status |
|------------|----------------|--------|
| `VMCParaOpt` | `vmc_para_opt!` | verified |
| `VMCPhysCal` | `vmc_phys_cal!` | verified for supported paths |
| `VMCMainCal` | `vmc_main_cal!` | verified |
| `VMCMainCal_fsz` | `vmc_main_cal_fsz!` | verified |
| `VMCMakeSample` | `vmc_make_sample!` | verified |
| `InitParameter` | `init_parameter!` | verified |
| `ReadInitParameter` | `read_initial_def!` | verified |
| `ReadInputParameters` | `read_input_parameters!` | parser verified; block consumption depends on feature |
| `CalculateGreenFunc` | `calculate_green_func!` | verified for supported direct/factored paths |

## Public test gates

Run from the repository root:

```bash
julia --project=@. test/integration/runtests.jl
julia --project=@. test/integration/ctest_equivalent.jl
julia --project=@. test/integration/phys_cal_equivalent.jl
julia --project=@. test/integration/lanczos_equivalent.jl
```

MPI smoke tests:

```bash
OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --project=@. test/mpi/run_mpi_smoke.jl
```

## Reference provenance

Committed reference fixtures live under `test/integration/reference/`. They
record the C commit, build flags, and regeneration commands in metadata files
where applicable. Normal CI does not require a local C-mVMC build.

## Floating-point tolerances

Energy and spin accumulators generally match the C reference within BLAS
summation-order noise. Derived quantities such as variance can amplify small
differences through cancellation.

`NSRCG = 1` parameter updates use a tolerance gate rather than bit parity
because truncated CG is sensitive to FMA and reduction-order differences.

## Threading compatibility

Julia-mVMC uses `JULIA_NUM_THREADS` for shared-memory execution. The
C-compatible policy is conservative:

- sample-level Markov-chain splitting is disabled;
- `VMCMakeSample` keeps sequential rank-local chain semantics;
- selected inner-loop threading is opt-in through `JULIA_MVMC_INNER_THREADS=1`;
- PfaPack QP-level threading remains a debug/benchmark triage mode.

## Not supported in v0.5.0

- BackFlow;
- spin Jastrow;
- `NSplitSize > 1` with SR-CG;
- `NSRCG >= 2`;
- `useDiagScale != 0`;
- `RescaleSmat != 0`;
- FSZ/general-orbital PhysCal split;
- PhysCal Lanczos split;
- FSZ/general-orbital Lanczos;
- OptTrans-derived QP sectors with split.
