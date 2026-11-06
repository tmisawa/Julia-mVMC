# Private ordinary complex quotient/inverse (#184/#176)

Base PR54 79485f4ba64092e0cdd0105fd0973293e369fe41. New private GNU/LLVM
scalar quotient and four direct solver branches are used only by the ordinary
complex child. Existing complex turbo factor, FSZ routes, real helpers,
TRTRI/copy/vT/column-permutation/FULL TRMM/row-permutation remain unchanged.
C ILAENV !full panel/trmmt branch is not newly ported. No dependency gitlink,
method piracy, new FFI, global mutable state or RNG change.

Linux+glibc HostPlatform libc selects GNU13.3 L_divdc3; Linux-musl and other
platforms select LLVM17. macOS aarch64 selects the existing observed FMA product
sums, others unfused product/sum. No generic Sys.islinux-only GNU assumption.
Both explicit kernels are tested against their own native datasets on Linux;
testing the archived ARM arithmetic is NOT native macOS execution proof.

Origin source/attribution: origin-gcc13.inc preserves the GCC13.3 whole-source
SHA97d364b0d74d6bf073eb3829ab12884ce243c68af8f7601d84baade12ee15599,
double-mode branch and complete divide plus original GPL/runtime-exception
notice. origin-llvm17.c is unchanged LLVM17 compiler-rt divdc3.c with original
Apache/LLVM notice. Source-adjacent licenses and NOTICES preserve this lineage.
The reviewed Rust scalar operational port additionally supplied binary scaling
and platform FMA selection; it is not mislabelled as upstream C authored code.
original-rust-validation-README.md records that separate bounded validation.

quotients-*.txt are byte-for-byte copies of existing independent native C
datasets tests/fixtures/interall/c_complex_division{,_linux_gnu,_macos_arm}.txt.
Headers retain compiler/runtime/source/library hashes. GNU/LLVM have373 rows;
Mac has375 including two archived equal-operand cases. Every input and
expected component is represented by hexadecimal IEEE payload. Expectations
are native, not Rust/Julia-generated. Directory name proves no InterAll scope.
Normal tests only read these static fixtures; never compile C or read toolbox.
Retained component policy16eps relative+4minimum-subnormal quanta; exact signed
infinity and NaN classification. Signed-zero inputs are preserved; zero outputs
use the retained numerical policy, not a newly introduced computed-bit gate.

complex6 input/output/operator copied from independently acquired native
issue184_inverse/complex6_pair_pivots.* fixtures. Original acquisition README
retains GNU13.3/nativeLP64OpenBLAS0.3.26 and Julia1.13.1ILP64threads1 metadata,
source/compiler flags, condition estimate11.7445, backward error/residual.
All A/M/vT planes use existing256eps abs+rel plus scaled backward residual,
not a new range/model/SR bound. Genuine net nonidentity control is independently
analytic dyadic4 with literal inverse; no nativeN6 forward tolerance is added.

FIRST old Julia diagnostic59767 terminal0, Julia1.13.1ILP64threads1,glibc:
373 cases each. inv-multiply firstdifference row2: z=1,w=-.1-1.7i,
delta6.94e-18-1.11e-16i, inside existing component budget. Julia / firstdifference
row5 zero divisor: NaN+NaNi vs nativeInf+NaNi (classification difference).
Budget/classification failures inv-multiply47, Julia /23. This is not a claim
that all supported finite models fail. Script/log/status retained externally
/tmp/mvmc-julia-complex-inverse-proof.T1V4st. No new C capture or golden regen.

Focused test requires Julia1.13.1 and existing Manifest09ebd06d, actual BLAS1.
No SR fixture adoption, full sampling/MPI/wholeworkspace or nativeMac claim.
Candidate publication and full regression require parent source review.

## mVMC 1.4.0 re-verification (Julia-mVMC#95, 2026-11-06)

The C reference moved to mVMC 1.4.0 (`v1.4.0`, `e83aa6cb`). `ltl2inv/invert.tcc`
(SHA256 `d477c7d2bb6d2b0c29872d06ec75fd7e804e1e30ace00888b77d31d2c74b821c`),
`ilaenv_lauum.cc`, `ilaenv_wrap.f90` and `ltl2inv.cc` are byte-identical in 1.4.0;
the fixture is unchanged and is a valid 1.4.0 reference. `extern/mVMC-1.3.0` paths
in the original acquisition notes are historical.
