# Third-Party Licenses

Julia-mVMC depends on two external Julia packages which are vendored as
git submodules (clone with `git clone --recurse-submodules`, or run
`git submodule update --init --recursive` after cloning). The licenses
for those submodules — and for the upstream sources they themselves
bundle — are documented inside each submodule.

Julia-mVMC's own source is licensed under GPL-3.0-or-later; see
[LICENSE](LICENSE).

---

## PfaPack.jl (git submodule)

- **Repository**: <https://github.com/tmisawa/PfaPack.jl>
- **License**: mixed BSD-3-Clause (primary) + MPL-2.0 (utu2 / C++
  wrapper sub-component) + LapackLicence (Wimmer Fortran kernels)
  + BSD-3-Clause (BLIS headers)
- **Per-file license map and full upstream license texts**:
  [PfaPack.jl/THIRD_PARTY_LICENSES.md](PfaPack.jl/THIRD_PARTY_LICENSES.md)
- **Acknowledgments**: PfaPack.jl was authored primarily by Satoshi
  Terasaki (AtelierArith, <s.terasaki@atelier-arith.jp>) and Takahiro
  Misawa. It bundles upstream sources from Michael Wimmer's PfaPack
  distribution (Fortran kernels), the xrq-phys/Pfaffine C++ utu2
  routines (RuQing Xu), and a header subset from the BLIS project.

## SFMT.jl (git submodule)

- **Repository**: <https://github.com/tmisawa/SFMT.jl>
- **License**: BSD-3-Clause
- **Full upstream license text**:
  [SFMT.jl/THIRD_PARTY_LICENSES.md](SFMT.jl/THIRD_PARTY_LICENSES.md)
- **Acknowledgments**: SFMT.jl was authored primarily by Satoshi
  Terasaki (AtelierArith) and Takahiro Misawa. The bundled SFMT C
  source is from Mutsuo Saito and Makoto Matsumoto (Hiroshima
  University).

**License compatibility**: BSD-3-Clause and MPL-2.0 are both
compatible with GPL-3.0-or-later (the latter via MPL §3.3 Secondary
License provisions). The combined work (Julia-mVMC + submodules) is
distributable under GPL-3.0-or-later, with the upstream BSD-3 / MPL
file headers preserved inside each submodule.

---

## C mVMC reference (test fixtures, NOT source)

Julia-mVMC ships the following fixtures in `test/integration/reference/`
to verify bit-level agreement with the C mVMC reference implementation:

- `**/zvo_out_first10.dat` — the first 10 SR steps of the C reference
  run's `zvo_out_001.dat` (output snapshots).
- `**/inputs/*.def` — the full set of expert-mode input definition
  files (e.g. `modpara.def`, `locspn.def`, `trans.def`, `coulombintra.def`,
  `orbitalidx.def`, `namelist.def`, ...) used to drive the C reference
  run. These are model / Hamiltonian definitions in mVMC's plain-text
  expert-mode format, not C source code.
- `**/inputs/initial.def` (HeisenbergChain_fsz and HubbardChain only) —
  initial variational parameters used by the C run, mirrored so
  Julia-mVMC starts from the same state.
- `**/physcal_ref/` (PhysCal e2e gate, Plan 3b) — a separate PhysCal input
  set (`inputs/*.def` incl. a hand-authored `greentwoex.def`, `NVMCCalMode=1`),
  the fixed `physcal_ref/zqp_opt.dat`, and the C-mVMC Green-function outputs
  `physcal_ref/expected/zvo_cisajs_001.dat`, `zvo_cisajscktalt_001.dat`,
  `zvo_cisajscktaltex_001.dat`. Output snapshots of a C reference run.

The full C mVMC source itself is **not** bundled.

- **Upstream**: <https://github.com/issp-center-dev/mVMC>
- **C reference commit used to regenerate these fixtures**:
  - opt / ctest-equivalent fixtures (`zvo_out_first10.dat`, `ctest_ref/`,
    `inputs/`): `5e7ea400ae35b566cfa2de6e342efe962f179a41` (master), gcc-15.
  - PhysCal e2e fixtures (`physcal_ref/`):
    `66f17422968009f8cc70f1dec94b2f52e562d344` (`develop`), Apple Clang 15 +
    gfortran (`USE_GEMMT=OFF`); the DH fixture was generated from local branch
    `feature/omp-simd-pfupdate` @ `622166afe33c6be3402d7c926db7e9c0003a47c4`,
    based on that develop commit plus benchmark/test-data commits. See each
    `physcal_ref/metadata.txt`.

  See [test/integration/reference/README.md](test/integration/reference/README.md)
  for the regeneration procedure.
- **License**: GPL-3.0-or-later. The bundled fixtures inherit the
  upstream license; combination with the rest of Julia-mVMC
  (also GPL-3.0-or-later) is straightforward.

---

## Credit consent

Satoshi Terasaki (AtelierArith) has confirmed his consent to be
credited as primary author of the PfaPack.jl and SFMT.jl submodules
(see each submodule's `Project.toml` `authors` field and `CITATION` /
README acknowledgments). The original SFMT19937.jl wrapper and the
ltl2inv / blalink C++ wrappers (originally under
`MVMCPfaPack.jl/deps/` in the mono-repo) were authored by Terasaki and
have been moved to the respective external submodules under the
licenses he chose (BSD-3-Clause for SFMT.jl; BSD-3-Clause + MPL-2.0
for PfaPack.jl).
