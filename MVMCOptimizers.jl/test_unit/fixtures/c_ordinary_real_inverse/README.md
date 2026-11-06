# Ordinary real direct inverse, private MVMC route (#184/#176)

Base PR54 bd7f3da33d5a785ffdb646ed32eed9369a679a41. No dependency gitlink,
complex/FSZ factor/inverse, new FFI, method piracy or RNG change.
The private helper retains existing TRTRI, upper copy, negated vT, column
permutation, BLAS TRMM and row permutation; only real scalar solver uses
the four numerator/denominator branches from unchanged invert.tcc.
Original source SHA256 d477c7d2bb6d2b0c29872d06ec75fd7e804e1e30ace00888b77d31d2c74b821c;
MPL-2.0 source-adjacent license and NOTICE preserve the C/Julia port origin.

Real4 input/output/operator are static copies of independently acquired
native C fixtures, NOT Rust/Julia-generated replacements. Original acquisition
README and real-direct-probe README are retained verbatim, with historical
absolute paths/scopes intact. Normal tests do not invoke C or read the Rust
repository/toolbox. Native C real4 max difference0 vs historical Julia; condition
estimate14.2, residual2.7755575615628914e-16, backward error1.826024711554534e-17.
Existing256eps abs+rel plus scaled residual is unchanged; it is not a generic
ill-conditioned matrix or SR/runner parameter tolerance.

Static payload SHA256:
real4_identity.input.txt 4ff4aaa135d7653d1ba42bcd63c7b88aa92e883f61ef46b1a84e6a03070a126d
real4_identity.c.txt 84c8f2f28e34d492156c80e063c26ac404e97d5b6a4d3c012ad164d5f1a66f45
real4_identity.operator.txt 89eb1b48cf0eb01ca3ec45d8b3e299a0a5385bf728346d40c1770801aa53c046

Analytic dyadic4 tridiagonal system exercises all four solve branches. Public
2x2 finite divisors2^1023 and2^-1023 have off-diagonal +/- reciprocal and zero
diagonals. Internal smallest-subnormal/Inf/zero classifications are diagnostics,
not accepted finite public models or provider-independent NaN TRMM assertions.
The original optional native real_direct_solver_probe.cc SHA256
ad44d79244708a1a067381d76c683bd0b52ac0a4cea728466787617f5f63446d,
GNU13.3 -O0 -ffp-contract=off, native LP64 OpenBLAS acquisition exited0;
compiler/library closure limitations remain explicit in original README.

BEFORE candidate, actual Julia diagnostic75234 terminal0: analytic4, native4,
finite2^+/-1023 all finite/maxerror0. This is not a fabricated Julia RED.
Internal smallest-subnormal zero numerator gave NaN, equal numerator gave Inf.
Old script/log/status preserved /tmp/mvmc-julia-real-inverse-proof.gjTnCL.
Rust-only old public NaN evidence is not claimed as Julia public failure.
Julia1.13.1 ILP64 OpenBLAS threads1, existing Manifest09ebd06d.
No native first-walker/full SR/trajectory adoption or tolerance waiver is inferred.

This helper deliberately retains the existing FULL TRMM pipeline. C invert.tcc
also has an ILAENV-selected !full panel/trmmt branch with skew reconstruction;
that branch is not newly ported or universally covered by this local patch.
Nonidentity pivot regression is a separate independent analytic4 control:
tridiagonal upper entries2,1,4, axes permutation[3,2,1,4], observed factor
pivots[1,2,1,4] (genuine net swap). The literal dyadic inverse is checked by
S*X=I and the ordinary child's row-major storage. All arithmetic is exactly
representable; no N6 forward budget or new tolerance follows. Existing native
N6 first-walker diagnostics are separate evidence, not imported golden values.

## mVMC 1.4.0 re-verification (Julia-mVMC#95, 2026-11-06)

The C reference moved to mVMC 1.4.0 (`v1.4.0`, `e83aa6cb`). `ltl2inv/invert.tcc`
and the rest of `ltl2inv` are byte-identical in 1.4.0; the fixture is unchanged
and is a valid 1.4.0 reference. `extern/mVMC-1.3.0` paths in the original
acquisition notes are historical.
