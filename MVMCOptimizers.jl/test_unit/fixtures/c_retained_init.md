# Independent retained-slot initialization fixture

`c_retained_init.txt`: first eight lines are C InitParameter Para (real/imag),
then exactly 624 unsigned SFMT words immediately after initialization.
Seed 11272; real mode; NProj=2, NRBM=3, NSlater=3, Nneuron=1.
Active real flags are zero-based Para 2,4,5,7: exactly four draws.
Mapped RBM idx=2 and Slater idx=1 deliberately leave active unmapped slots.
This is a standalone initialization-kernel probe, NOT a full C/MPI trajectory.

The probe copies ONLY the complete `void InitParameter()` definition from
mVMC 1.3.0 src/mVMC/parameter.c lines 35–94, without algorithm modifications.
Origin: https://github.com/issp-center-dev/mVMC (vendored mVMC-1.3.0).
Original source SHA256: 46ad04622f4475337028cee633bd76ce318a6d5058d03f202b55204500399fb0.
Probe SHA256: c96d020a36195e21a67be019c0c92eccae974e119b20d6463991fc7018a1601c.
Fixture SHA256: 36f6af91291c9d902a423b18eee554dcf0a1933b6e5dbcf92ed35ac875ce43eb.
Original function is GPL-3.0-or-later, copyright 2016 University of Tokyo.

Linux x86_64, GCC 13.3.0, C11, -O0 -ffp-contract=off; no BLAS in probe.
SFMT submodule 1526553009f318ae78338151460fda78beadddc2:
SFMT.c SHA256 4eed94fb587cefa3022d259379fdd492ddd76a4709d23aae9c7dd6ebc52db321;
SFMT-real.c SHA256 1e283c215c48ec94a1083c963314bd4a1c6d4498fe4ffbbbb0b80240178de36e.

Optional reproduction from repository root (tests never invoke this probe):

```sh
make -C SFMT.jl/deps/sfmt
cc -D_DEFAULT_SOURCE -std=c11 -O0 -ffp-contract=off \
  MVMCOptimizers.jl/test_unit/fixtures/c_retained_init_probe.c \
  -L SFMT.jl/deps/sfmt -lsfmt -lm -o /tmp/c_retained_init_probe
LD_LIBRARY_PATH="$PWD/SFMT.jl/deps/sfmt" /tmp/c_retained_init_probe
```

Julia regression: Julia 1.13.1, local Manifest-v1.13.toml matching the pinned
reference workspace; OpenBLAS ILP64, one thread. Parameters use explicit
eps(Float64)*0.01 absolute / 2*eps(Float64) relative tolerance; next624 and flags are exact.
Measured maximum absolute C–Julia parameter difference was 0.0 on both
same-width initializations: no arithmetic first divergence occurred. The bound
covers rounding of the real RBM scaled draw (magnitude <= .005, subtract/multiply
then divide by Nneuron=1); real Slater draws use binary-exact subtract/scale.
This small local operation bound is not a tolerance for sampling drift.
The Manifest is not part of this patch. No Rust-generated expected values.
## mVMC 1.4.0 re-verification (Julia-mVMC#95, 2026-11-06)

The C reference moved to mVMC 1.4.0 (`v1.4.0`, `e83aa6cb`). `InitParameter` is
byte-identical in 1.4.0 `src/mVMC/parameter.c` (still lines 35–92; whole-file
SHA256 `801273e47a6f7b3390244a6bdcb74fa09825791e8b8283d2d6d7e36d25470cad`, changed
only in the DH offsets for spin Jastrow, `NSpinJastrowIdx = 0` here), so the probe
and fixture are unchanged and are valid 1.4.0 references.
