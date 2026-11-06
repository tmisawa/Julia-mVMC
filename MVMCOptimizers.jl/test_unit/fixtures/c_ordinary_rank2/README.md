# Ordinary real C-order rank2 regression (#176/#184)

Static native.txt is byte-for-byte copied from the independently acquired
Rust repository fixture tests/fixtures/issue176_initializer/rank2-boundary/native.txt,
SHA256 8492f90b077e0559907e5c75a82d4ef8864d48c90675ed95861c7287b084e411.
It is native C output, not Rust-generated or Julia-generated expectations.
original-README.md, environment.txt and source-sha256.txt preserve original
acquisition provenance, including historical absolute paths (not runtime dependencies).
LapackLicence preserves the authoritative PFAPACK license and copyright notices.

Authoritative unchanged source: mVMC1.3.0 src/pfapack/fortran/dsktf2.f
SHA256 ea1591aa942377add5cd42597c282d1f282700251fa4b825f026cad1c6c2f65b,
upper/normal lines208-260; dskr2.f lines162-178 SHA256
57bf85cfabc22055af6c95d5765200ef5c27c9b5dff074da264fc76bbf3a56f2.
The new MVMC private Julia port retains pivot first-max, swaps, INFO, temporary
products, zero-vector skip, scale and left-associated update. Only the ordinary
real child is routed to it; PfaPack gitlink/API, complex, FSZ and inverse solver
remain unchanged. No method piracy, new FFI or runtime oracle.

Optional independent reproduction from the Rust repository:
`bash c_toolbox/issue176_rank2/reproduce.sh` (README source/extraction review).
Reviewed actual handle93753 exit0, container73c57e563c61,
/tmp/mvmc-issue176-rank2.RORK4x, GNU13.3.0 -O0 -ffp-contract=off,
Linux x86_64 OpenBLAS0.3.26 Haswell pthreads LP64, actual threads1.
Probe SHA256 1c383e056f02f7fc82911a8d971da647c5673b9a8150b9cca804d587a1703ab6;
script SHA256 75b6fdeabae28ac9278fd72cacd53ad014e5d53a2a03fac516448460d3c0fef5.
Native factor INFO0/PF-1 proves this dyadic cancellation boundary. The fixture
inverse is diagnostic only; no forward accuracy claim for this ill-conditioned
matrix. PF uses the reviewed relative4EPS scalar policy for two dyadic factors
and sign, not a whole-matrix budget. The native child is FSZ, but factor input
is explicitly identical to the ordinary full-occupancy child's four spin-sites.
This does not establish native ordinary runner/first-walker or SR parity.

Focused developer test (Julia1.13.1, existing Manifest-v1.13.toml):
```
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia +1.13.1 --project=. \
  -e 'using LinearAlgebra; BLAS.set_num_threads(1); include("MVMCOptimizers.jl/test_unit/test_unit_ordinary_real_rank2.jl")'
```
Parent reviewed the private port, caller, tests, license and original native
provenance. New candidate actual Julia1.13.1/ILP64 OpenBLAS threads1 runs all
terminated0: focused13, optimizer subpackage standard24270 (15+25+24230),
filtered root PhysCal6 and original PairHop prefix1 controls16. All883 checkout
files matched before/after. Proof: /tmp/mvmc-julia-ordinary-rank2-proof.YagFur;
commands.md and results.md record exact commands, runtime and source/log hashes.
These are bounded checks, not whole Julia workspace or native runner parity.
Parent/C owner confirmed independent native N6 first-walker handle26130
terminal0, stage bHW3fT/input40e, source/library integrity checks0. Native
PF119.26565040127034 matches the rank2-route first saved QP0. This comparison is
bounded to that checkpoint, not full trajectory, solver, SR or complex/FSZ parity.
Final packaging/publication is controlled by the parent. No fixture regeneration
or tolerance change.

After the frozen runs, only the caller's obsolete performance comment and this
README were corrected. Executable algorithm and test code are unchanged; the
post-comment package manifest has new hashes, not retroactively assigned to the
historical frozen-source manifest. No scalar port performance claim is made.

## mVMC 1.4.0 re-verification (Julia-mVMC#95, 2026-11-06)

The C reference moved to mVMC 1.4.0 (`v1.4.0`, `e83aa6cb`). The authoritative
PFAPACK `dsktf2.f`/`dskr2.f` (nested pin `ebaa11b7`) and `matrix.c` are
byte-identical in 1.4.0. `vmcmake_fsz_real.c` and `projection.c`, which the
original acquisition harness also compiled (`source-sha256.txt`, historical 1.3.0
paths), changed only by guarded t-J/doublon-only/spin-Jastrow branches that are
inactive for this input (`NExUpdatePath` 0, no spin Jastrow). The fixture is
unchanged and is a valid 1.4.0 reference.
