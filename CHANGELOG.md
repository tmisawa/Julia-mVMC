# Changelog

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
