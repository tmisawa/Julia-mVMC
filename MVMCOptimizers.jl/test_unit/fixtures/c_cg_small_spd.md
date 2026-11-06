# C CG small-SPD kernel fixture

Upstream mVMC 1.4.0 (`v1.4.0`, `e83aa6cb`) `src/mVMC/stcopt_cg_impl.c:484–591`,
the complete `fn_StochasticOptCG_Main` definition copied unchanged into the
companion probe. Source SHA256
84a1bb1867a957f25d8a822f67213ff044bf52448d96263c2efdf085bb815dd5.
Probe SHA256 4d50e858c1b894db9077ce29ff725f53d19b7154c1b4fb689cb3d58959ff643f.
Origin https://github.com/issp-center-dev/mVMC; GPL-3.0-or-later, copyright
2016 University of Tokyo (inline notice included).

Only MPI rank/size and a one-dimensional sampled operator are adapters.
O=g=2^-54, mean O=0, W=1, stabilization=0, tolerance1e-30, limit1.
S=O*O=2^-108 is strictly SPD, condition number1. delta=2^-108 is above the
convergence threshold1e-60, but d*S*d=2^-216 is below the absolute
`fabs(dq) < 1.0e-20` breakdown guard that 1.4.0 added (C:540). C 1.4.0
therefore returns status `-(0+1)` before the first update: x stays 0 and
r=d=g=2^-54. Output columns: status, solution, residual, search direction.

History: under mVMC 1.3.0 (`d73d06bd`, lines 258–351, source SHA256
41452de5fe766409431c6e1cf73cd2b12485faaeb4bfbcdec1e53444e9af1cf9, probe SHA256
d7a0c22250081c9130a358380493afabfa865b05b682bde48c7c0aea481a4be1) the same
input had no denominator guard and returned `1 18014398509481984 0 0`
(x=2^54, the analytical solution). The fixture now records the 1.4.0 contract.

Optional reproduction from repository root:

```sh
cc -std=c11 -O0 -ffp-contract=off \
  MVMCOptimizers.jl/test_unit/fixtures/c_cg_small_spd_probe.c -lm -o /tmp/c_cg_small_spd
/tmp/c_cg_small_spd
```

Regenerated 2026-11-06 on Linux x86_64, GCC 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1),
`-O0 -ffp-contract=off`, no OpenMP/SIMD. The probe invokes no BLAS, no sampling,
no RNG and no MPI runtime; this is not a complete C sampled SR/MPI execution
claim. The status is compared exactly; x is the untouched initial zero and r/d
are copies of g, compared within eps(Float64)*2^-54.
