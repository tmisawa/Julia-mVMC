# Lanczos PhysCal performance

Issue [#89](https://github.com/tmisawa/Julia-mVMC/issues/89) reports a slow
single-thread mode-1 PhysCal measurement. A profile of the reported six-site
input with `NVMCSample=2000`, `NDataQtySmp=1` attributed 799 of 821 sampling
snapshots to `calculate_lanczos_h2!`; most were in intermediate-state Slater
factorization and inversion. A warm, in-process API call with profiling took
0.902 s and allocated 277 MiB on the host below. The updated path took 0.137 s
and allocated 166 MiB with the same profiler and shortened input. These times
include profiling overhead and are not the full issue benchmark.

Measurements use Linux x86-64 on an AMD Ryzen 9 PRO 8945HS, Julia 1.13.1,
ILP64 OpenBLAS, one Julia and one BLAS thread, and the fixed input and seed from
[the issue](https://github.com/tmisawa/Julia-mVMC/issues/89). The Julia API
is measured inside one process, excluding the first call and explicit GC;
Rust commit `5f5a490e3f4a4b77149b147f684ed2bbae252d77` is a release CLI
process measured from launch to exit. Both use one computational thread.
For the full issue workload (`NVMCSample=20000`, `NDataQtySmp=20`), the
updated Julia calls took 22.505, 22.596, and 22.548 s (median 22.548 s)
after one warmup; Rust took 13.09, 13.01, and 13.12 s (median 13.09 s)
after one warmup. The timing scopes differ as described above. For the first
measured run, all 320 QQQQ values in the 20 output files were finite, with
maximum Julia/Rust absolute difference `1.14e-13`.

For the separate short `native_c_lanczos_576` InterAll fixtures in the Rust
repository, the nodal mode-1/2 C outputs agree, but the ordinary real
mode-1/2 fixtures have a pre-existing Julia/C discrepancy: their maximum QQQQ
difference is 1.164 on both the baseline and updated Julia implementations.
This comparison is distinct from the issue workload and is not a regression
introduced by this change.

The real overlap path reuses the original inverse for up to two changed
electrons. It falls back to direct factorization near a node or for more
changes, and never forms an unused intermediate inverse. The complex path
retains the original calculation. `test_unit_lanczos_pfaffians.jl` compares
updates with direct factorization in multiple projection sectors, including a
singular intermediate state; `test/integration/lanczos_equivalent.jl` checks
C-referenced mode-2 outputs.
