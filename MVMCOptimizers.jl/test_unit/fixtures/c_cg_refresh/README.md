# Independent C sampled CG refresh records

## mVMC 1.4.0 records (current)

The six `*.txt` files are C mVMC 1.4.0 (`v1.4.0`, `e83aa6cb`) outputs of the
optional fixed-input kernel adapter `ctest_cg_refresh.c` (imported unchanged from
the Rust port's `c_toolbox/ctest_cg_refresh.c`, https://github.com/AtelierArith/mvmc-rs,
adapter SHA256 `9df41611459585db474b016c5905d753c158194e5516dbf9f1d82baf11d3ea11`).
`ctest_cg_main_upstream.inc` is the byte-for-byte concatenation of 1.4.0
`src/mVMC/stcopt_cg_impl.c` lines 1–21 (license) and 350–668
(`StochasticOptPreconCG_DiagScale_Main`, `StochasticOptCG_Main`, sampled operator);
upstream SHA256 `84a1bb1867a957f25d8a822f67213ff044bf52448d96263c2efdf085bb815dd5`,
extraction SHA256 `63d35a521f1b3ddcea93b7009881594fb20f091a2968dc4aecc4709621482cb6`.
`ctest_cg_dot_upstream.inc` (`stcopt_cg.c:26–34`) is unchanged from 1.3.0
(extraction SHA256 `f36c3f7d0aa0fc22e7e3cbffde3c95b79808112979e1e0ebc8fdd853d995d214`).

The status column is the C 1.4.0 return value: tolerance zero never converges,
so every limit records `-(limit+1)`, or `-(iter+1)` at an earlier breakdown guard.
The `*_diag_scale.txt` files run the DiagScale-preconditioned Main (`diag` argument).
Inputs are the operand rows of the previous (1.3.0) records, unchanged.

Regenerated 2026-11-06, Linux x86_64, GCC 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1),
`-O0 -ffp-contract=off`, no OpenMP/SIMD, LP64 OpenBLAS pthread 0.3.26+ds-1ubuntu0.1
(library SHA256 `bfc7492adbf84a8f567720a9e1fae2afc18f3d817da233e7f4d453683485308e`),
one thread. The output hashes equal the independent mvmc-rs #682 regeneration:

| File | C 1.4.0 output SHA256 |
| --- | --- |
| real.txt | f5b6e8340e886943423f3e9d85d8d3ca298c0c61d4bb6fc3d9fd649c934f3fa2 |
| complex.txt | 7009acebae2fc64bdfbd58d2b8fd471248c0475ba627f844bf4e80f6f0fc3088 |
| sampled_complex.txt | eec31746b75c09e58d44a973ee8d1f0dca9e48670ad4e7662c5ccec4f5990091 |
| real_diag_scale.txt | 3a3a7007db4e4c0830dbbff76ba2e0b06f10398e8a2ee23ba1a649784b2a4e45 |
| complex_diag_scale.txt | d9f0f158b16d23971ebde5efdb044b6070709a26943f56e47cb4016736f0dd6d |
| sampled_complex_diag_scale.txt | 5551ab36a2db88ea1a364f5f3b63552949505d261913a02e4d90386ca2755d96 |

Against the 1.3.0 records, the real and complex plain trajectories are
bit-identical and only the status column changed (1.3.0 returned the iteration
count). The sampled complex trajectory stops at iteration 5 for every limit above
5 (status `-6`): its old `delta` falls below the new absolute `1e-20` guard,
where 1.3.0 continued (108 data rows differ for that reason).

```sh
cc -O0 -ffp-contract=off -DMVMC_SRCG_REAL ctest_cg_refresh.c -lopenblas -lm -o /tmp/real-cg
cc -O0 -ffp-contract=off ctest_cg_refresh.c -lopenblas -lm -o /tmp/complex-cg
OPENBLAS_NUM_THREADS=1 /tmp/real-cg real.txt /tmp/real.txt
OPENBLAS_NUM_THREADS=1 /tmp/complex-cg complex.txt /tmp/complex.txt
OPENBLAS_NUM_THREADS=1 /tmp/complex-cg sampled_complex.txt /tmp/sampled_complex.txt
# DiagScale: append `diag`, e.g. /tmp/real-cg real.txt /tmp/real_diag_scale.txt diag
```

## Historical mVMC 1.3.0 generation

The text below describes the original 1.3.0 records, which these files replace.

These three files are imported without modification from the Rust port's
optional C kernel oracle, documented in `c_toolbox/ctest_cg_refresh.md` of
https://github.com/AtelierArith/mvmc-rs . They are NOT Rust-generated results
or re-labelled historical Julia solver expectations. Input operand records
originated in archived Julia1.13.1 records; only means, diagonal, stored
samples and gradients are read by the C adapter, never archived solver output.

Authority: mVMC1.3.0 `src/mVMC/stcopt_cg_impl.c`, unchanged Main/operator
definitions (lines255–426) plus original license (lines1–21), SHA256
41452de5fe766409431c6e1cf73cd2b12485faaeb4bfbcdec1e53444e9af1cf9.
Extracted include SHA256
9f941cc5f9208eed0efc80937f292837c48b7bc7d078f467715ca001d9e9f43f.
Unchanged dot definition `stcopt_cg.c:26–34`, upstream SHA256
66a36cdff6f22e23367e12c08fee505daff6194447a667f5b89a65eb30130700;
include SHA256 f36c3f7d0aa0fc22e7e3cbffde3c95b79808112979e1e0ebc8fdd853d995d214.
Upstream https://github.com/issp-center-dev/mVMC; GPL-3.0-or-later,
copyright2016 University of Tokyo. Adapter SHA256
5e41bb737a018fda5d1ae014672954d5544f4f78fe0882d76be38e81ba77548e.
The post-allocation-guard adapter regenerated all three files with identical
hashes. No numerical C function body is patched.

Serial fixed-input kernel, NOT C initialization/sampling or MPI execution.
Generation: Linux x86_64, GCC13.3.0,
-O0 -ffp-contract=off, LP64 OpenBLAS pthread0.3.26, one thread;
library SHA256 bfc7492adbf84a8f567720a9e1fae2afc18f3d817da233e7f4d453683485308e.
Real executable c72f45946fab01d6de2b2afa66357b5149439428672a8a29c2add4762102783a;
complex executable969413856298b78745fcf74c3882fa5d9863bbc844a36bee91e0a31df77f1631.

| File | Input SHA256 | C output SHA256 |
| --- | --- | --- |
| real.txt | 95765162416342c86bace224de36dd315e5a0efa9ee1dd501556b8074df7c0a4 | e9c1372e2aa7b1ff12e15c44c74f5875030e622b80c985eacdee35c853122844 |
| complex.txt | 8258518a952cfe34db3a409ee9034a5c25905b7cf87ca0799b71ed4360c64eb5 | bd48dd44525c79082f6db667ff6dc6103f7c363328198044c8c839df230b80ea |
| sampled_complex.txt | 87829b292761dada69121c4d1fd175c0b4a71e029b5c6ae6dded2dcd2b80894f | 6036d436653df6c02846827194cfe82f43597e93f392f5dd52b932f2e166b590 |

Each file has shape n/sample-count/complex, six hexadecimal operand/product
rows, then solution/residual/direction for restarted limits1–41. Hex records
retain acquisition precision; computed assertions are numerical, not bitwise.
Tolerance0, shift1e-5, weight=sample count, both20/40 residual refreshes.

The exact adapter and both licensed upstream extracts are bundled alongside
these records, with the hashes above. This makes regeneration independent of
an unpublished Rust toolbox commit. Each output contains the original input
operands as its first six rows; the adapter consumes those operands only.
Optional regeneration from this directory in the recorded LP64 environment:

```sh
cc -O0 -ffp-contract=off -DMVMC_SRCG_REAL ctest_cg_refresh.c -lopenblas -lm -o /tmp/real-cg
cc -O0 -ffp-contract=off ctest_cg_refresh.c -lopenblas -lm -o /tmp/complex-cg
OPENBLAS_NUM_THREADS=1 /tmp/real-cg real.txt /tmp/real.txt
OPENBLAS_NUM_THREADS=1 /tmp/complex-cg complex.txt /tmp/complex.txt
OPENBLAS_NUM_THREADS=1 /tmp/complex-cg sampled_complex.txt /tmp/sampled_complex.txt
```

Normal Julia tests read these checked-in fixtures only, without a C compiler,
oracle invocation or Rust checkout. The original validation used Julia1.13.1
and OpenBLAS0.3.30 ILP64/Haswell with one thread.

The operator and first CG step use the original component budget
`8*(n+samples)*eps`. Later truncated iterates are not compared componentwise
across BLAS providers: the real, complex and sampled-complex covariance
condition numbers are approximately `3.6e9`, `5.0e8` and `9.6e6`.
An iteration-count-times-epsilon budget cannot bound their forward sensitivity.
On macOS aarch64 the archived trajectories first exceed that budget at
iterations 7, 7 and 6, even though the independent residual checks pass.

All limits 1--41 still exercise the production solver and validate the archived
C records. Both native and Julia residuals are checked against the explicit
covariance with the original budget
`16*(n+samples)*limit*eps*(abs(g)+sum(abs(A*x)))`. Limits 20 and 40 also
verify that the solver replaces its residual by `g-S*x` using the same local
backend. This checks the refresh boundary without requiring cross-provider
bitwise equality. Residual consistency alone is not a convergence test.

Separate analytic SPD tests provide the convergence gate: 64 paired sample
directions give a known diagonal covariance, including nonzero means, a
diagonal shift and real/imaginary sample contributions. Its condition number
is exactly 4096. The solver must cross both refresh boundaries, converge,
and agree with the independent analytic solution and residual within
dimension/condition-scaled roundoff bounds. The tiny-SPD native regression
still detects the former absolute denominator cutoff. Existing fixture bytes
and the production CG recurrence are unchanged.
