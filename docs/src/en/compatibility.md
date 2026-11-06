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

## C reference version

The C authority is mVMC 1.4.0 (`v1.4.0`, `e83aa6cb`; nested StdFace `6fa4ef1f`,
BLIS `a32257ee`, PFAPACK `ebaa11b7`, the same pins as 1.3.0). C-native fixtures
were re-verified or regenerated against it; see
`test/integration/reference/README.md` and the fixture READMEs under
`MVMCOptimizers.jl/test_unit/fixtures/`. Adopted 1.4.0 contracts:

- SR-CG returns the C status: the iteration count on convergence, `-(iter+1)` on
  the breakdown guards (non-finite `delta`, `d.q`, `alpha`, `beta`; `|d.q|` or the
  old `|delta|` below `1e-20`) and `-(max_iter+1)` at the iteration limit. A
  failure aborts by default (`NSRCGAbortOnFail 1`, a Julia error with C's
  message); `NSRCGAbortOnFail 0` continues with the approximate solution, and
  `NSRCGFallback 1` first retries with the other CG variant from `x = 0`.
- `NSRCG 2` (or `useDiagScale 1` with `NSRCG 1`) selects DiagScale-preconditioned
  CG; `RescaleSmat 1` rescales the Slater block before SR-CG and requires
  `NSRCG 1`. SR-CG with a negative `DSROptStepDt` is rejected.
- `zvo_SRinfo.dat` CG rows end with `"%5d %d"` (the C status) and the header has
  the ` info` title.
- The Jastrow reader requires every directed pair exactly once with the same
  class for `(i,j)` and `(j,i)`.
- The real two-electron Pfaffian update reads the old b-site (`rbOld`); the
  PhysCal timer report has the `addPhysTwist [54]` line.
- `reweight 1` is rejected (reweighting is not ported).

With the default CG tolerance some ordinary inputs reach the iteration limit or
the absolute `1e-20` guards and therefore abort by default under 1.4.0; set
`NSRCGAbortOnFail 0` for the previous continue-anyway behavior.

Kept deviations: a real-parameter RBM keeps the correct RBM factor (C 1.4.0
rejects RBM with real parameters, Julia-mVMC#59); OptTrans flags and real-mode
OptTrans derivative slots (#94, #55); an SR-CG step with no active parameter is a
no-op (C 1.4.0 aborts and reads `r[0]` out of bounds). New 1.4.0 features
(reweighting, twist, spin Jastrow, t-J, `NExUpdatePath = 6`, backflow) are not
ported. The GPU path is unchanged: it supports direct SR only.

## Floating-point tolerances

Energy and spin accumulators generally match the C reference within BLAS
summation-order noise. Derived quantities such as variance can amplify small
differences through cancellation.

`NSRCG = 1` parameter updates use a tolerance gate rather than bit parity
because truncated CG is sensitive to FMA and reduction-order differences.
The discrete C 1.4.0 SR-CG status, fallback and abort decisions and the integer
`zvo_SRinfo.dat` fields are compared exactly (`test/integration/srcg_option_equivalent.jl`).

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
- `reweight = 1` (C 1.4.0 reweighting);
- FSZ/general-orbital PhysCal split;
- PhysCal Lanczos split;
- FSZ/general-orbital Lanczos;
- OptTrans-derived QP sectors with split.
