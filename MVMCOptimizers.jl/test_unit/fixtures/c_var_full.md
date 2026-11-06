# Fixed-value C var formatting reference

Generated separately on Linux x86_64 with GCC13.3.0,
`cc -std=c11 -O0`, from mVMC1.3.0 `vmcmain.c:655–657`'s literal
fprintf templates and full `for(i=0;i<NPara;i++)` loop. No sampling, BLAS,
computed floating-point normalization or Julia/Rust oracle is involved.
The test only reads the checked-in bytes; it does not compile/invoke C.

Input: Etot=-3+0im, Etot2=9+0im, NPara33. Slots are Gutz2, Jast2,
DH2 six, DH4 ten, nine RBM slots, Slater2, OptTrans2. Zero-based slots
1 and 3 are zero; all others hold `(i+1)/8 - (i+1)/16*im`.
These are binary-exact constants. Each complex value has a literal third
field `0.0`, followed by a space; the whole row ends with one newline.
This is a writer fixture, not a full accepted-input/trajectory reference.

Authoritative `vmcmain.c` SHA-256:
`fdcd661c4eb028786fc58ae5e1f5232426c943b547f95a0d74e888e602256d63`.
Original standalone formatting driver SHA-256:
`60fca5b77c96d7e33119a39efbd604a9a4af7e9bf5f28ab5fd2f17d7cbb2eee2`.
Driver: `AtelierArith/mvmc-rs`'s optional `c_toolbox/physcal_181_format.c`.
## mVMC 1.4.0 re-verification (Julia-mVMC#95, 2026-11-06)

The C reference moved to mVMC 1.4.0 (`v1.4.0`, `e83aa6cb`). The `zvo_var.dat`
templates and loop are byte-identical at 1.4.0 `vmcmain.c:657–659` (whole-file
SHA256 `da14885da14fd5b98feb27f7c53efd6fa229a4a33d27c30fe44848e6514394f0`); the
fixture is unchanged.
