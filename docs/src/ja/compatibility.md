# 互換性と検証

[English](../en/compatibility.md)

Julia-mVMC は C-mVMC compatibility を標準の correctness target として扱います。

## Function correspondence

| C function | Julia function | Status |
|------------|----------------|--------|
| `VMCParaOpt` | `vmc_para_opt!` | verified |
| `VMCPhysCal` | `vmc_phys_cal!` | supported paths verified |
| `VMCMainCal` | `vmc_main_cal!` | verified |
| `VMCMainCal_fsz` | `vmc_main_cal_fsz!` | verified |
| `VMCMakeSample` | `vmc_make_sample!` | verified |
| `InitParameter` | `init_parameter!` | verified |
| `ReadInitParameter` | `read_initial_def!` | verified |
| `ReadInputParameters` | `read_input_parameters!` | parser verified; block consumption depends on feature |
| `CalculateGreenFunc` | `calculate_green_func!` | supported direct/factored paths verified |

## Public test gates

repository root から実行します。

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

committed reference fixture は `test/integration/reference/` 配下にあります。必要な
fixture では C commit、build flags、regeneration command を metadata に記録しています。
通常の CI では local C-mVMC build は不要です。

## Floating-point tolerances

energy / spin accumulator は通常、BLAS summation-order noise の範囲で C reference と
一致します。variance のような derived quantity は cancellation により小さな差を増幅します。

`NSRCG = 1` の parameter update は bit parity ではなく tolerance gate です。truncated CG は
FMA と reduction-order difference に敏感なためです。

## Threading compatibility

Julia-mVMC は shared-memory execution に `JULIA_NUM_THREADS` を使います。C-compatible
policy は conservative です。

- sample-level Markov-chain splitting は無効。
- `VMCMakeSample` は sequential rank-local chain semantics を保つ。
- selected inner-loop threading は `JULIA_MVMC_INNER_THREADS=1` で opt-in。
- PfaPack QP-level threading は debug / benchmark triage mode のまま。

## v0.5.0 で未対応のもの

- BackFlow;
- spin Jastrow;
- `NSplitSize > 1` with SR-CG;
- `NSRCG >= 2`;
- `useDiagScale != 0`;
- `RescaleSmat != 0`;
- FSZ/general-orbital PhysCal split;
- PhysCal Lanczos split;
- FSZ/general-orbital Lanczos;
- OptTrans-derived QP sector with split.
