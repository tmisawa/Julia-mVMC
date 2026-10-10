# Calibrated partial-pool matrix-work gate

The initial uncapped partial-pool candidate (`e86730a`) improved Opt but
regressed L32 PhysCal, so it was not adopted:

| PhysCal sites | Before (s) | Ungated candidate (s) | Change |
|---:|---:|---:|---:|
| 32 | 4.654 | 5.834 | +25.4% |
| 64 | 22.394 | 16.978 | -24.2% |

These are warmed 100-group, 320-total-sample, 1 rank × 16 thread runs,
three alternating pairs. All numerical output differences were observed as
zero. The timing and output checks are retained in the accompanying JSON.

An independent kernel probe alternated serial/parallel calls for eight QPs
in a sixteen-thread pool, validating normalized inverse residuals. Each
repetition averaged 100 calls after warming both paths. Busy serial gaps
were outside the reported kernel durations. Median microseconds:

| Matrix order | Serial, no gap | Parallel, no gap | Serial, 1ms gap | Parallel, 1ms gap |
|---:|---:|---:|---:|---:|
| 24 | 34.40 | 20.81 | 33.67 | 129.41 |
| 32 | 56.07 | 27.97 | 54.64 | 64.71 |
| 40 | 87.55 | 34.08 | 86.25 | 167.83 |
| 48 | 126.87 | 49.87 | 129.01 | 100.59 |
| 64 | 241.63 | 78.46 | 241.46 | 153.17 |
| 96 | 675.04 | 173.23 | 676.07 | 176.87 |

With serial gaps, small matrices no longer amortize static dispatch. Require
`qp_num * n_size^3 >= 32768 * nthreads()` when QPs are fewer than pool
threads. At 8 QPs / 16 threads this retains the serial path at n=32 and
n=40, while allowing n=48 and n=64. Full-pool paths keep their prior behavior
and the caller retains its existing work gate. This is an empirical dispatch
estimate, not a numerical tolerance or a universal hardware crossover.

Validation of the calibrated source at this commit: the 18 explicit worker
execution assertions and the complete sixteen-thread unit suite passed
(45,350 optimizer assertions plus base and Slater suites). Matrix arithmetic,
RNG operations and numerical tolerances are unchanged. The first worker-test
fixture incorrectly expected success for a block-diagonal matrix whose zero
intermediate column gives INFO=K-1 in C DSKTF2 normal mode (dsktf2.f:216–228).
The corrected execution test uses the existing dense residual-test family.

End-to-end timing and RNG audits of the calibrated source remain required
before merging. The PR body records subsequent verification; the final mvmc-rs
reference-update report will retain the complete evidence. Earlier Opt and
PhysCal timings are for the ungated source, not the calibrated candidate.

`qp-crossover.jl` reproduces the historical calibration on commit `e86730a`
with the same Julia 1.13.1 project, one BLAS thread and JULIA_NUM_THREADS=16,0.
Do not run it on the calibrated source and label those timings as ungated.
Standalone kernel times are not a whole-VMC speedup claim.
