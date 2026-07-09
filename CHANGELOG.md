# Changelog

## v0.5.0 - 2026-07-08

This release extends the experimental PhysCal and MPI split coverage after
v0.4.2. It adds C-referenced Full-Lanczos PhysCal output, split normal-Green
PhysCal runs, sz-conserved standard-projection sample splitting, pair-hop
fixtures, and FSZ/general-orbital `TwoBodyGEx` measurements.

### Added

- Added `VMCPhysCal` Full-Lanczos output support for `NLanczosMode = 1/2` on
  the sz-conserved `NSplitSize = 1` path, with C-reference fixtures for
  Lanczos energy, QQQQ, one-body, direct two-body, and factored two-body Green
  output where supported.
- Added C-reference integration coverage for pair-hop Hubbard fixtures in real
  and FSZ/general-orbital modes.
- Added direct-SR `VMCParaOpt` `NSplitSize > 1` support for sz-conserved
  standard-projection `NQPFull > 1` inputs when `NQPOptTrans = 1`
  (`NSPGaussLeg > 1` and/or `NMPTrans > 1`), including MPI self-consistency
  smoke coverage.
- Added `VMCPhysCal` `NSplitSize > 1` support for sz-conserved normal-Green
  runs (`NLanczosMode = 0`), including MPI reduce-to-root and rank0 output
  smoke coverage.
- Added FSZ/general-orbital `TwoBodyGEx` / `greentwoex.def` support for
  factored two-body Green output and a C-referenced FSZ PhysCal fixture.
- Added RBM-bearing `initial.def` / `zqp_opt.dat` parameter loading and promoted
  `GeneralRBM_cmp` to the C ctest-equivalent integration gate.

### Changed

- Updated the public version metadata for the in-repo packages
  `MVMCOptimizers` and `MVMCExpertModeParsers` to `0.5.0`.
- Tightened unsupported-input validation for split PhysCal, split Lanczos,
  FSZ/general-orbital standard-projection splitting, OptTrans-derived QP
  sectors, and additional SR-CG modes.
- Updated README and manual compatibility notes for the v0.5.0 PhysCal,
  Lanczos, MPI split, and unsupported-input scope.

### Notes

- `VMCPhysCal` remains experimental. Split PhysCal is limited to
  sz-conserved normal-Green runs (`NLanczosMode = 0`); FSZ/general-orbital
  split, Lanczos split, and OptTrans-derived QP-sector split runs remain
  rejected.
- `VMCPhysCal` Lanczos support is limited to the sz-conserved
  `NSplitSize = 1` path. FSZ/general-orbital Lanczos, split Lanczos, and
  ParaOpt Lanczos remain unsupported.
- `VMCParaOpt` split standard projection is limited to sz-conserved inputs.
  FSZ standard-projection `NQPFull > 1`, OptTrans-derived QP sectors with
  `NSplitSize > 1`, SR-CG with `NSplitSize > 1`, `NSRCG >= 2`,
  `useDiagScale != 0`, and `RescaleSmat != 0` remain rejected.
- GitHub-generated source ZIP/TAR archives do not include submodule contents.
  Use `git clone --recurse-submodules https://github.com/tmisawa/Julia-mVMC`
  for a functional checkout.
- `PfaPack.jl` and `SFMT.jl` are submodules and remain at their own package
  version `0.1.0` in this release.

## v0.4.2 - 2026-06-23

This patch release adds the next MPI compatibility slice after v0.4.1. It adds
MPI-compatible standard SR-CG (`NSRCG = 1`) for `NSplitSize = 1` and direct-SR
grouped sample splitting for `NSplitSize > 1` with `NStore = 0/1` when
`NQPFull = 1`.

### Added

- Added MPI-compatible standard SR-CG support for `VMCParaOpt` with
  `NSRCG = 1` and `NSplitSize = 1`, including the C-compatible
  broadcast/allreduce ordering in `operate_by_s`, rank0 parameter update, and
  parameter synchronization across ranks.
- Added direct-SR `NSplitSize > 1` support for `VMCParaOpt` with
  `NStore = 0/1`, non-FSZ and FSZ paths, and the `NQPFull = 1` scope.
- Added MPI smoke coverage for serial-vs-MPI `NSRCG = 1` first-step
  C-reference checks and four-way direct-SR `NSplitSize`/`NStore`
  self-consistency gates under `mpiexec -n 2/-n 4`.

### Changed

- Updated the public version metadata for the in-repo packages
  `MVMCOptimizers` and `MVMCExpertModeParsers` to `0.4.2`.
- Tightened unsupported-input validation for SR-CG split runs, additional CG
  modes, grouped QP splitting, and `VMCPhysCal` `NSplitSize > 1` before
  MPI context construction or long-running computation.
- Updated README and manual compatibility notes for the v0.4.2 MPI scope.

### Notes

- `NSplitSize > 1` support is currently the R1 direct-SR scope:
  `VMCParaOpt`, `NSRCG = 0`, `NStore = 0/1`, non-FSZ/FSZ, and `NQPFull = 1`.
  Grouped QP splitting (`NSplitSize > 1` with `NQPFull > 1`) is not ported.
- Standard SR-CG (`NSRCG = 1`) remains limited to `NSplitSize = 1`.
  `NSplitSize > 1` with SR-CG, `NSRCG >= 2`, `useDiagScale != 0`, and
  `RescaleSmat != 0` remain rejected.
- `VMCPhysCal` still supports only `NSplitSize = 1`.
- GitHub-generated source ZIP/TAR archives do not include submodule contents.
  Use `git clone --recurse-submodules https://github.com/tmisawa/Julia-mVMC`
  for a functional checkout.
- `PfaPack.jl` and `SFMT.jl` are submodules and remain at their own package
  version `0.1.0` in this release.

## v0.4.1 - 2026-06-19

This patch release hardens the experimental MPI sample-parallel path introduced
in v0.4.0 and improves the real Hubbard optimization kernels. The supported MPI
scope is unchanged: `VMCParaOpt` with `NSplitSize = 1`, `NSRCG = 0`, direct SR,
and one software thread per MPI rank.

### Added

- Added CI/site smoke coverage for the Hubbard MPI path alongside the existing
  Heisenberg and `VMCPhysCal` MPI checks.
- Added focused MPI failure-mode checks for unsupported `NSRCG != 0` and
  `NSplitSize > 1` inputs before long computation or rank-divergent collectives.
- Added gated C-style diagnostic timers for Julia `VMCParaOpt` runs to support
  kernel-level performance triage without changing the default execution path.

### Changed

- Updated the public version metadata for the in-repo packages
  `MVMCOptimizers` and `MVMCExpertModeParsers` to `0.4.1`.
- Reused main-calculation workspaces and optimized real Hubbard kernels,
  including Slater/Green-function update paths, `CalHamiltonian1`, and
  WeightAverage active-range handling.
- Tightened the release benchmark scaffold for Genkai/Ohtaka MPI timing,
  correctness, and launcher checks.
- Documented the current MPI performance status in
  <https://github.com/tmisawa/Julia-mVMC/issues/37>.

### Notes

- MPI support remains experimental and limited to `NSplitSize = 1`.
  C's grouped MPI/QP split (`NSplitSize > 1`) still raises an unsupported-input
  error.
- `VMCParaOpt` under MPI still supports the direct SR solver only
  (`NSRCG = 0`). `NSRCG != 0` under MPI remains rejected pending a separate
  MPI SR-CG design and C-parity gate.
- GitHub-generated source ZIP/TAR archives do not include submodule contents.
  Use `git clone --recurse-submodules https://github.com/tmisawa/Julia-mVMC`
  for a functional checkout.
- `PfaPack.jl` and `SFMT.jl` are submodules and remain at their own package
  version `0.1.0` in this release.

## v0.4.0 - 2026-06-14

This release adds experimental MPI sample-parallel execution, hardens the
parameter-layout and input-overlay paths, and brings the serial `NSRCG = 1`
SR-CG solver onto a documented C-reference gate. The default serial path
remains C-parity oriented; grouped MPI/QP splitting and MPI CG solver runs are
still rejected.

### Added

- Added MPI.jl-based sample-parallel infrastructure for `VMCParaOpt` and
  `VMCPhysCal` with `NSplitSize = 1`, including rank0 output/readback guards,
  comm0 reductions, seed offset handling, and MPI smoke gates for
  `mpiexec -n 2/-n 4`.
- Added fail-fast MPI detection policy for common launcher environments,
  including explicit `JULIA_MVMC_MPI` handling.
- Added serial `NSRCG = 1` SR-CG parity fixes and a first-step C-reference
  integration fixture (`heisenberg_chain_real_nsrcg`).
- Added MPI-aware parameter synchronization helpers and shared flat-layout
  walkers for variational parameter blocks.

### Changed

- Updated the public version metadata for the in-repo packages
  `MVMCOptimizers` and `MVMCExpertModeParsers` to `0.4.0`.
- Applied `InOrbital` overlays by variational-parameter index instead of array
  position, matching duplicate-index semantics.
- Tightened rank0-only logging/warning behavior under MPI so multi-rank runs do
  not duplicate user-visible output.
- Documented the `NSRCG = 1` serial tolerance gate and the remaining MPI
  limitations in the README and manual.

### Notes

- MPI support is experimental and currently limited to `NSplitSize = 1`.
  C's grouped MPI/QP split (`NSplitSize > 1`) still raises an unsupported-input
  error.
- `VMCParaOpt` under MPI supports the direct SR solver only (`NSRCG = 0`).
  `NSRCG != 0` under MPI raises an unsupported-input error because C's CG
  `operate_by_S` broadcast/allreduce path is not ported yet.
- Serial `NSRCG = 1` is supported, but the post-CG parameter update is a
  tolerance gate (`NSRCG_PARAM_TOL = 1e-2`) rather than bit parity because
  truncated SR-CG amplifies FMA and reduction-order differences.
- GitHub-generated source ZIP/TAR archives do not include submodule contents.
  Use `git clone --recurse-submodules https://github.com/tmisawa/Julia-mVMC`
  for a functional checkout.
- `PfaPack.jl` and `SFMT.jl` are submodules and remain at their own package
  version `0.1.0` in this release.

## v0.3.0 - 2026-06-10

This release expands the C-reference coverage for physical-quantity
calculation, tightens unsupported-input handling, and adds conservative
single-process threading infrastructure. The default execution path remains
C-parity oriented; MPI and independent-chain sampling are still unsupported.

### Added

- Added `VMCPhysCal` runner support from `namelist.def` and C-referenced
  end-to-end fixtures for one-body, direct two-body, and factored/product
  two-body Green-function output.
- Added `greentwoex.def` / `TwoBodyGEx` parsing and factored two-body Green
  output (`zvo_cisajscktaltex`) for the supported non-FSZ path.
- Added DH2/DH4 projection-layout parsing and runtime support, including
  C-reference regression fixtures.
- Added `InOrbitalParallel` offset handling, warn-only `OptTrans` input
  handling, and additional `ReadInputParameters` / parameter-sync regression
  tests.
- Added conservative threading helpers for selected inner loops, local
  accumulator construction/merge, and CI guards for `JULIA_NUM_THREADS > 1`.
- Added explicit unsupported-input gates for `NSplitSize > 1` and
  `NLanczosMode > 0`.

### Changed

- Updated the public version metadata for the in-repo packages
  `MVMCOptimizers` and `MVMCExpertModeParsers` to `0.3.0`.
- Expanded the C ctest-equivalent and PhysCal reference gates used by CI.
- Kept `VMCMainCal` sample processing sequential by default and removed the
  sample-level threading opt-in to preserve C-compatible stochastic semantics.
- Documented the threading policy, PhysCal coverage, and unsupported-input
  contract in the manual and package README files.

### Notes

- `JULIA_MVMC_INNER_THREADS=1` enables only selected C OpenMP-equivalent
  inner-loop threading. Leave it unset for the default C-parity path.
- `JULIA_MVMC_PFAPACK_THREADS=1` remains a debug/benchmark triage mode and is
  not yet treated as C-compatible.
- `NSplitSize > 1` still raises an unsupported-MPI error. MPI support is
  planned as a separate design/review track.
- GitHub-generated source ZIP/TAR archives do not include submodule contents.
  Use `git clone --recurse-submodules https://github.com/tmisawa/Julia-mVMC`
  for a functional checkout.
- `PfaPack.jl` and `SFMT.jl` are submodules and remain at their own package
  version `0.1.0` in this release.

## v0.2.0 - 2026-06-03

This release is an infrastructure, verification, and diagnostics update.
The VMC solver path is unchanged in intent, and the committed C-reference
checks continue to verify numerical compatibility.

### Added

- Added a C ctest-equivalent integration runner covering 12 standard C-mVMC
  fixtures with the same `ref_mean.dat` / `ref_std.dat` acceptance rule as
  C-mVMC's `test/python/runtest.py`.
- Added committed ctest-equivalent reference data and regeneration tooling
  for the supported standard fixtures.
- Added a C-compatible opt-in section timer for Julia VMCParaOpt runs.
  Enable it with `MVMC_C_TIMER=1` to write `zvo_CalcTimer.dat`.
- Added CI coverage for Ubuntu and macOS on Julia 1.11 and 1.12.
- Added Julia 1.11 manifest snapshot `Manifest-v1.11.toml` as a
  reproducibility aid.

### Changed

- Updated the public version metadata for the in-repo packages
  `MVMCOptimizers` and `MVMCExpertModeParsers` to `0.2.0`.
- Updated README, citation metadata, and manual pages for the v0.2 release.
- Documented timer behavior, ctest-equivalent integration tests, and
  reference-data provenance.
- Removed the old `TimerOutputs` dependency from `MVMCOptimizers`.

### Notes

- `MVMC_C_TIMER` is disabled when unset or exactly `0`; any other value,
  including `false`, `off`, `no`, or an empty string, enables timing.
- Timer-enabled runs include measurement overhead. Use them for bottleneck
  analysis, not as primary elapsed-time benchmarks.
- GitHub-generated source ZIP/TAR archives do not include submodule contents.
  Use `git clone --recurse-submodules https://github.com/tmisawa/Julia-mVMC`
  for a functional checkout.
- `PfaPack.jl` and `SFMT.jl` are submodules and remain at their own package
  version `0.1.0` in this release.

## v0.1.0 - 2026-05-14

- Initial public release of Julia-mVMC.
