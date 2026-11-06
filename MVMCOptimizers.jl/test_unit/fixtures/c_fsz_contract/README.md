# Independent native FSZ contracts

These are offline native C reference results, not Julia- or Rust-generated
expectations. Normal Julia tests read only these fixtures; no oracle is run.

Original authority: mVMC-1.3.0 `src/mVMC/matrix.c`, SHA256
`849488176375102fbc3bab7001e0dc27dfa7c8e049fd90ee1d6248ebfba10764`.
The GPL header (lines 1–21) and complex/real children (125–179, 223–278)
were extracted verbatim. The driver compiled the original native factorization
and utu2pfa/utu2inv implementations, not replacement arithmetic.

`rows.txt` SHA256
`d9eb1fd0d556d24accade5794c146f276548c699db35f76ce3d5bb3d8741512b`
contains fourteen real/complex cases. Each row has 41 fields: mode, case,
zero-based QP start, local QP, child status, factor INFO, inverse calls,
PF real/imag, then sixteen column-major inverse real/imag pairs.
The regular upper triangle is (.5,.125,.25,.375,.625,1.25), electron sites
(0,1,0,1), spins (0,0,1,1), Nsite=2, Nsize=4. Its analytic PF is .640625.
`complex_regular` replaces entry 01 with .5+.25i only for the complex case.
Zero/Inf/NaN and offset cases retain the native status rather than inventing
a Julia global-index status. Failure rows contain native dirty scratch:
Julia deliberately retains its stronger all-QP staged-publication boundary.

`sum-overflow.txt` SHA256
`049bb082853ea86b877c82364cf2f1985790285632787d95e66cfd7ff2e0716c`
sets entry 01 to 1e308 (real) or 1e308+1e308i (complex). Native complex PF
components are both 1.25e308, but their sum is nonfinite: status=1,
INFO=0, inverse calls=0. The real case succeeds.

Acquisition used Linux x86_64, GNU C/C++/Fortran 13.3, -O0
-ffp-contract=off, original Fortran SKTRF/SKTF2/LASKTRF/SKR2/SKR2K,
original ltl2inv with BLAS_EXTERNAL, LP64 OpenBLAS 0.3.26 Haswell, one thread.
This is a standalone child check, not a full native executable/MPI trajectory.
The verbatim child extraction and independent acquisition driver are bundled
as `children.inc` and `rows_probe.c`. Optional reproduction:
`bash generate_rows.sh /absolute/mVMC-1.3.0/src`. Native source identity is
checked against the SHA256 above; no unpublished Rust repository is required.
Compilation uses the original native sources supplied explicitly by developer.

`real-initializer-101.txt` SHA256
`b2a98556d7fde8964a8eb63f68408190917c6e32b45652888d8707eeef88482c`
is a mechanical extraction from independent acquisition
`/tmp/mvmc-issue176-initializer.rtEpR2/final-checkpoints.json`, profile
`success-attempt101`, final checkpoint. Lines: attempt/status/abort/cursor/count,
tmp indices/config/occupations/spins, raw624, next624.
The original real initializer (`vmcmake_fsz_real.c` SHA256
`ef84d79cc8a59e8da611bf1118c76d67954b865dd3d0ee8a1b685384d3457f9b`,
lines 421–502 plus GPL header) ran the native child with a finite sparse
Nsite=128, Ne=1 table, entries (83,220)=1 and (220,83)=-1, SFMT seed1.
No scripted kernel status: first success is attempt101. C then aborts after
incrementing/checking the limit; Julia instead returns its existing error Int.
The native probe SHA256 is
`f2208728396e76e53879434ccda5f2508e899c1a05105545f7469ceb84244353`.
Julia's existing SFMT public wrapper exposes next words, not raw state/cursor
or a draw counter; tests do not claim those unobserved fields were measured.
Independent complex-initializer acquisition completed separately at
`/tmp/mvmc-julia-fsz-complex101.4NALRe` (exit0). The optional
`generate_complex.sh /absolute/mVMC-1.3.0/src` reproduces it using the original
complex initializer lines435–517, complex matrix lines77–179, projection
lines58–152 and their GPL headers. Complex source SHA256:
`2b43a47a85ec13833d453a010962334147024a6d6aae01b99bb8b2a0b80c15db`.
Its original SFMT seed1, finite sparse matrix and natural 101st success were
not scripted. Actual GNU13.3/OpenBLAS0.3.26 Haswell, threads1; source checks0.
Binary SHA256 `158ac5e7b45ac10e4965635e7ef0fae7a2b91d008f007cd77620739c871836f9`;
result SHA256 `80904237b85afb4ca9f49ce7bc807bab55ba35dc0ff9b5da2f345d3f4f43eee7`.
Driver SHA256 `894592fb1400ad95f803089fdb5581de9932dac7a91cde1a0a6c163fa863363e`;
generator SHA256 `556e45cd48772e96d55a56c01aef1716a5fd6147356994d3f0cae74242510cbc`.
Exact acquisition command: `bash generate_complex.sh /workspaces/mvmc-rs/extern/mVMC-1.3.0/src`
in the Linux container; handle68078 exited0. The source/environment/linked-path
records are bundled in `complex-acquisition/`. `libraries-post.sha256` hashes
the six resolved file-backed dependencies after acquisition (including loader);
these are explicitly post-acquisition closure, not hashes captured before run.
OpenBLAS resolved to `openblas-pthread/libopenblasp-r0.3.26.so`, SHA256
`bfc7492adbf84a8f567720a9e1fae2afc18f3d817da233e7f4d453683485308e`.
The virtual linux-vdso has no file hash. This validates only serial native
kernel/initializer behavior: MPI is a one-rank adapter, not native MPI evidence.
`complex-initializer-101.txt` has the same SHA256 as the real fixture: this is
independently acquired identity of discrete checkpoints, not copied expectations.
The earlier stage `.aGXi7h` failed compilation on an extraction boundary;
its failure is preserved, not counted as successful acquisition. Collective
per-attempt retry agreement is separate and is not repaired by this serial patch.

Validation command (Julia 1.13.1, local Manifest-v1.13.toml):

```sh
OPENBLAS_NUM_THREADS=1 julia +1.13.1 --compiled-modules=existing --project=. \
  -e 'using LinearAlgebra; BLAS.set_num_threads(1); include("MVMCOptimizers.jl/test_unit/test_unit_fsz_native_contract.jl")'
```

The local Manifest and compiled dependency libraries are not patch files.
Computed values use explicit small-dimension numerical comparisons, not bits;
statuses, range/shadow/padding preservation and next624 remain exact.

The component policy matches the reviewed native child test: divide both
components by their maximum magnitude before subtraction, then allow 64eps
plus four minimum-subnormal units divided by that scale. There is no unit
absolute floor: zero cannot replace a roughly 1e-308 inverse component.
Regular PF .640625 uses exact equality (binary-exact inputs/result).
Each successful inverse also satisfies the independently reconstructed
`Slater * storedInv = -I` residual. Four products and three additions per entry
are bounded against the sum of absolute products, evaluated at 256-bit precision
to avoid overflow in the huge input. This is a small-matrix regression budget,
not a general condition-independent solver policy.

Real-kernel complex shadows deliberately remain unchanged at the low-level
call, matching Julia's separate real/complex state architecture rather than
Rust's refresh convenience. The real sampler's proposals, updates and log-IP
read real arrays. `vmc_para_opt!` explicitly converts from real arrays after
sampling; PhysCal's `vmc_main_cal_fsz!` recomputes the complex matrix/PF for each
saved sample before using complex measurements. Neither path requires the
low-level real kernel to write complex shadows. The existing real sampler
Slater-table synchronization is not changed by this patch.

### Public real-FSZ analytic control and retained RED checkpoints

The focused test additionally creates two separate output directories and runs
the existing public `vmc_phys_cal!` API after the namelist wrapper's
parse -> fixed load -> input overlays -> Slater sync preparation sequence.
The convenience `run_phys_cal_from_namelist` wrapper itself is not executed by
this control: it exposes neither a callback nor its internally created RNG.
The lower public API exposes both without instrumentation or production edits.

Input is two itinerant sites, four electrons (`NCond=4`, `NElec=2` pairs),
General orbitals, identity translation, `NMPTrans=1`, `NSPGaussLeg=1`, and
`ComplexType=0` in both declared parameter families (Gutzwiller and Slater).
Their active real optimization flags are zero; inactive imaginary storage is
not required to be zero, but the complete stored flag vector must remain
unchanged during the run. Fixed input values are respectively 0 and 1; C
`parameter.c:SyncModifiedParameter` scales Slater by `D_AmpMax/xmax`, so the
post-sync fixed vector is `[0,4]`. This is not an optimization update.
The one declared Gutzwiller coefficient is mapped to both sites (`[0,0]`).
The namelist uses `TransSym` for the identity translation mapping; Julia's
`QPTrans` keyword instead denotes its momentum-style parser. The parsed
translation indices, signs and weight are asserted explicitly.
DH/Jastrow/RBM/OptTrans are absent, not implicitly claimed exercised. There
are six upper-triangle General mappings, all using the same Slater coefficient.
Julia's parser accepts flattened records `I J idx sign`. Native C
`readdef.c:GetInfoOrbitalGeneral` instead accepts
`i spin_i j spin_j idx sign`, with `I=i+2*spin_i`, `J=j+2*spin_j`.
The six mappings are exactly the six C-supported `I<J` pairs; the generated
Julia file is **not** a native C reader fixture. Parsed General mode, term
complex flags, and `get_all_complex_flag=false` determine the real-FSZ branch;
no `:real`/`:fsz` label selects execution.

Each run has seed 1, warmup 1, interval 1, one frame, one saved walker,
`NSROptItrStep=1`, and Lanczos 0. PhysCal does not execute optimization steps.
Only onsite Coulomb terms are present, strength 2 at each site. C
`calham_fsz_real.c:79-83` accumulates `U*n_up*n_down`; full occupancy therefore
gives E=4, H2=16, variance=0, Sz=Sz2=0. All four diagonal one-body
occupations are 1, and every off-diagonal entry is 0 by Pauli exclusion.
These immutable analytic expectations are not sampled or Rust-generated
goldens. The dimension-four normalized observable comparison uses 64eps times
its physical scale (16 for H2; unit scale for zero variance/occupation), the
same small-operation coefficient as the existing focused comparisons; the
subnormal inverse comparator and all existing reference budgets are unchanged.

Callback count/index/status/energy, unchanged fixed parameters and flags, and
ordered 16 Green records are checked. Next624 is exact between the two
same-input/same-seed executions only. No raw SFMT state, cursor, draw count,
saved state, hidden PF/inverse, C sampling trajectory, or MPI parity is claimed
through this public API. Native child fixtures above cover kernel numerics
separately. At the historical prelaunch proposal checkpoint, actual terminal
results and backend/thread metadata were pending parent review and explicit
launch authorization. Completed results are recorded in the final section below.

First public-control execution, owner handle `98415`, exited 1 with 8 passes,
4 failures and 2 errors (14 assertions). Its preserved log is
`/tmp/mvmc-julia-realfsz-public-proof.vOBRlV/focused.log`, SHA256
`ee3847d893ac5c58a18997bd1656f6bbb71e3e931eed7494d35e47fda9aa5a6b`.
The 868-file pre/post source manifests were identical, SHA256
`58610d3f8edd227cf00e3d8cae0fafe532baf2bdd62e8f231d6db4ebee228090`.
This was a harness failure, not evidence of a Julia numerical defect: it used
the momentum `QPTrans` keyword for translation records, expected two mapped
Gutzwiller terms rather than one declared coefficient, omitted C Slater
normalization in its expectation, and incorrectly required inactive imaginary
flag storage to be zero. The missing translation mapping aborted the public
call before sampling/callback/output; subsequent native testsets did not run.
Only those four harness issues are corrected. Production and every numerical
budget remain unchanged. The new full focused rerun is separately authorized;
its actual terminal results are not yet available.

The corrected harness execution `9869` also exited 1: 14 passes, 1 failure,
2 errors (17 assertions), 58.9 seconds. Unlike `98415`, all setup assertions
passed and the public call reached `vmc_make_sample_fsz_real!`. At saving its
first sample, the real caller passed `log_ip_old::Float64` to
`save_ele_config_fsz!`, whose sole method required `ComplexF64`. This is an
actual Julia public real-FSZ dispatch defect, not an input/roundoff defect.
Preserved log `/tmp/mvmc-julia-realfsz-corrected-proof.JYL9H0/focused.log`
has SHA256 `1cfcafa1e0a89721cd3fc61227c1d147362d0589281a233d9649c3b8979df9b5`;
the identical 868-file pre/post manifests have SHA256
`a9a41d765ee879a5b5488208d7a2a1df9547181f380ee1a9936e4e6937eabf0e`.

The dedicated PR54 append repairs only the save-helper argument signature to
`Union{Float64,ComplexF64}`. The helper does not read that argument: its body
copies five saved electron/configuration/projection/spin planes. Both actual
callers keep their real/complex arithmetic and existing arguments unchanged.
A literal, sampler-independent regression calls the helper with each log-IP
type and checks all five saved planes and unchanged input arrays exactly.
The unchanged public analytic control remains the end-to-end real-path
regression. No mode forcing, complex conversion of real arithmetic, tolerance
change, or separate upstream PR is introduced. At the historical prelaunch
union-signature checkpoint, the fresh full focused result and subsequent
standard/bounded-root revalidation were pending review. Completed results are
recorded below; earlier GREEN results are not relabelled as proof of this new
production hash.

### Final union-signature validation

The pending statements above describe historical acquisition/review stages;
the following final executions supersede them without relabelling either RED
run. Production SHA256 for all final executions is
`6cf602be483f967de71b253fead910d92d3583bea7f9b02caccabad2c7c4deb0`,
focused test SHA256 is
`26cb92a5b023d67ab0fda9187bb1dd5c120b723f04de3571ecfe324b3174b3a8`.

| Execution | Actual owner handle / exit | Passed assertions |
| --- | --- | --- |
| Full focused file | 54915 / 0 | 385: helper 12, public real 143, native child 15, initializers 20, native rows 195 |
| Standard optimizer subpackage | 13634 / 0 | 24,257: main 15, Slater 25, unit 24,217 |
| Original filtered PhysCal FSZ root | 39307 / 0 | 6, 59.4 seconds |
| Original two PairHop prefix-1 root cases | 71208 / 0 | 16, 50.7 seconds |

Focused artifacts are `/tmp/mvmc-julia-realfsz-dispatch-proof.t9G4p5`;
its log SHA256 is
`e1c4300b6fcb437a90314e9ed3f11a20254c0939a62027e5ebae8b72b34ac09c`.
The identical 868-file source manifests have SHA256
`04af935d80cb0c23b2b5dff000f0b53e56c6743afcf6704784fdb6feb93dfe39`.
Both real public calls reached the end of measurement and passed the analytic
observables, fixed-parameter/full-flag retention, existing callback, and exact
next624 repeatability checks. Raw SFMT state/cursor/count remain unobserved.

Final standard/root artifacts are
`/tmp/mvmc-julia-fsz54-final-union-proof.dpNaPE`, including exact commands,
logs, and identical 873-file source manifests (including ignored native
dependency files), SHA256
`b5685fe1b742a1f28ccfbb23068359ae668322c74a3cdc461b9c527b97d90ccb`.
Standard log SHA256 is
`d1bd2086a6e710b88fd8a13a9e6d5b4ae5924b61c8f259e0b8977ec64f7ca5ed`;
PhysCal log SHA256 is
`355be13c01ae03acbbc74fd92126c12dd23bb977b22c8a501e96196cf57e85ca`;
PairHop log SHA256 is
`59f690426fd07ec6fe4dce0c809e2a0b500f3a85a6e3c21fc856b1dc4c7bd2c9`.
All logs print Julia 1.13.1, Julia threads 1, actual ILP64
`libopenblas64_.so` BLAS threads 1, and the loaded dedicated-fork package
paths/runtime libraries. This is not a prehashed system/compiled-cache closure.
Manifest and both original root tests retained the hashes recorded above.
This result summary was added after terminal verification; it does not change
production/tests or pretend the final README bytes were compiled test inputs.

Standard means the optimizer subpackage, not the whole Julia workspace.
Bounded root inputs/counters/seeds and every existing numerical budget were
unchanged. These original FSZ root inputs use ComplexType 1; the new lower
public API control separately verifies ComplexType 0. No full C sampling
trajectory, MPI, collective-per-attempt retry, or whole #179/#185 closure is
claimed. No reference golden was regenerated by these Julia executions.

Before this addition, the untouched root tests had SHA256
`d08a09d36b8673aabbd95bfe95165c2da5e1d3434d5292366ebf887235068e44`
(`test/integration/phys_cal_equivalent.jl`) and
`18c5d366437752b0dadce4a762e315c0637426165476dcf8cc5fc531d08a70f8`
(`test/integration/pairhop_equivalent.jl`). The uncommitted local
`Manifest-v1.13.toml` SHA256 is
`09ebd06dab244510094b99fe7c6efa2fe7a3d22221d1336b951123a5a6e8befc`.
The authorized execution must record source closure before/after, Julia
1.13.1, actual loaded package paths, actual BLAS provider and thread count.
Manifest and generated runtime/cache artifacts must not be committed.
## mVMC 1.4.0 re-verification (Julia-mVMC#95, 2026-11-06)

The C reference moved to mVMC 1.4.0 (`v1.4.0`, `e83aa6cb`). Re-run on Linux
x86_64, GNU 13.3.0, `-O0 -ffp-contract=off`, LP64 OpenBLAS 0.3.26 (pthread,
library SHA256 `bfc7492a…5308e`), one thread:

- `rows.txt` / `sum-overflow.txt`: `bash generate_rows.sh <mVMC-1.4.0>/src` exits 0;
  `matrix.c`, PFAPACK and `ltl2inv` are byte-identical in 1.4.0 and both outputs
  compare equal (`cmp`).
- `complex-initializer-101.txt`: 1.4.0 changed `vmcmake_fsz.c`
  (SHA256 `f6322ed48eb186c097074ea9287906fcf58129a44d91b932d0a6750ea9daba66`) and
  `projection.c` (`9459b9819d08bd67dc40bc9a2103a43851667af3775e4a32752f1f7f559e5794`):
  `makeInitialSample_fsz` moved to lines 461–556 and gained a t-J branch for
  `NExUpdatePath` 4/5, and `MakeProjCnt` (58–166) a spin-Jastrow block. The 1.4.0
  extraction (`sed -n '1,21p;461,556p'`, `'1,21p;58,166p'`, `matrix.c`
  `'1,21p;77,179p'`) with `complex_probe.c` plus the 1.4.0 globals
  `int NExUpdatePath=0,NSpinJastrowIdx=0,**SpinJastrowIdx;` reproduces
  `complex-initializer-101.txt` byte for byte (SHA256 `b2a98556…8482c`).
- `real-initializer-101.txt`: its acquisition input is not bundled, so it was not
  re-run; in 1.4.0 `makeInitialSample_fsz_real` (`vmcmake_fsz_real.c` 447–541)
  keeps the 1.3.0 loop verbatim in the branch taken for `NExUpdatePath` other than
  4/5, so the recorded checkpoint is unchanged. It equals the re-run complex
  record as before.

The `generate_*.sh` scripts still check the 1.3.0 hashes they were written for;
they document the historical acquisition.

