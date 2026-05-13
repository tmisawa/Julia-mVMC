# 2. Input files (`.def` format)

Julia-mVMC reads the same C-mVMC **expert-mode** `.def` files. The
top-level entry point is `namelist.def`, which lists every other
`.def` file by keyword. All paths in `namelist.def` are resolved
relative to its own directory.

## Supported `.def` files

| File / namelist keyword | Purpose | Status |
|------|------|------|
| `namelist.def` | top-level dispatch, lists all other files | ✅ |
| `modpara.def` (`ModPara`) | system size, sampling, SR, RNG seed | ✅ |
| `locspn.def` (`LocSpin`) | local-spin site flags | ✅ |
| `trans.def` (`Trans`) | one-body hopping / transfer | ✅ |
| `coulombintra.def` (`CoulombIntra`) | on-site U | ✅ |
| `coulombinter.def` (`CoulombInter`) | nearest-neighbour V | ✅ |
| `hund.def` (`Hund`) | Hund coupling | ✅ |
| `exchange.def` (`Exchange`) | spin-spin exchange | ✅ |
| `interall.def` (`InterAll`) | general 4-fermion interaction | ✅ |
| `gutzwilleridx.def` (`Gutzwiller`) | Gutzwiller correlator indexing | ✅ |
| `jastrowidx.def` (`Jastrow`) | Jastrow indexing | ✅ |
| `orbitalidx.def` (`Orbital`) | Slater orbital indexing | ✅ |
| `orbitalidxgen.def` (`OrbitalGeneral`) | fsz mode generalised orbital | ✅ |
| `orbitalidxpara.def` (`OrbitalParallel`) | parallel-orbital block | ✅ |
| `qptransidx.def` (`TransSym`) | quantum projection symmetry | ✅ |
| `greenone.def` (`OneBodyG`) | one-body Green-function targets | ✅ |
| `greentwo.def` (`TwoBodyG`) | two-body Green-function targets (direct mode) | ✅ |
| `initial.def` (CLI 2nd argument in C) | starting variational parameters | ✅ auto-detected from namelist dir |

## In* overlay files (consumed by `read_input_parameters!`)

These optional files overwrite parameter values *after* `initial.def`
has been applied (matching C's `vmcmain.c` order
`InitParameter → ReadInitParameter → ReadInputParameters`). When both
`initial.def` and an `In*.def` carry the same parameter, the `In*.def`
value wins.

| File / namelist keyword | Target field | Status |
|------|------|------|
| `InGutzwiller.def` (`InGutzwiller`) | `gutzwiller_terms[i].value` | ✅ |
| `InJastrow.def` (`InJastrow`) | `jastrow_terms[i].value` | ✅ |
| `InOrbital.def` / `InOrbitalAntiParallel.def` | `orbital_terms[i].value` | ✅ |
| `InOrbitalGeneral.def` | `orbital_terms[i].value` (fsz layout) | ✅ |
| `InChargeRBM_PhysLayer.def`, `_HiddenLayer.def`, `_PhysHidden.def` | corresponding RBM term arrays | ✅ |
| `InSpinRBM_*` and `InGeneralRBM_*` (3 layers each) | corresponding RBM term arrays | ✅ |

## Not yet supported

These are recognised by the parser but **not** consumed by the
optimisation or sampling pipeline. Reference data generated with them
active is not reproducible bit-for-bit.

| File | Status / reason |
|------|--------|
| `pairhop.def` (`PairHop`) | Pair-hopping interaction not wired into `CalculateHamiltonian`. |
| `dh2.def` / `dh4.def` (`DH2`, `DH4`) | Doublon-holon correlator stubs only. |
| `InDH2.def`, `InDH4.def` | `read_input_parameters!` emits a warning and skips (not yet wired into the `Proj` array layout). |
| `InOrbitalParallel.def` | `read_input_parameters!` emits a warning and skips (the C-side `iNOrbitalAntiParallel` offset path is not yet implemented). |
| `InOptTrans.def` | `read_input_parameters!` emits a warning and skips. The `FlagOptTrans` gate and `OptTrans[]` storage do not exist on the Julia side yet. |
| `cisajscktalt.def` (TwoBodyG product side) | Direct two-body Green is supported; the `<c†c>×<c†c>` factorisation path is not. |
| `OptTrans` block in `initial.def` | Refused with a warning by `read_initial_def!` because `FlagOptTrans` / `OptTrans[]` are not implemented yet. |

## Compatibility notes

- Section headers (`NGutzwillerIdx`, `NTransfer`, etc.) follow the C
  convention: 5 ignored lines (header + counter + 3 separators) before
  the data block.
- `qptransidx.def` accepts the 3-column form (`itmpsgn = 1` defaulted)
  in addition to the canonical 4-column form.
- For `interall.def`, fsz-style spin-flip terms are honoured.
- The RNG layout matches C: `RndSeed` from `modpara.def` seeds an
  SFMT19937 stream; values <= 0 fall back to `11272`. See
  `MVMCOptimizers.run_para_opt_from_namelist` for the precise phase
  ordering (`init_parameter!` → `read_initial_def!` →
  `read_input_parameters!` → `sync_modified_parameter!` →
  `init_qp_weight!`), which mirrors C's `vmcmain.c:264-281`.

## Generating inputs

The samples in [`examples/inputs/`](../../examples/inputs/) and
[`test/integration/reference/`](../../test/integration/reference/) were
all produced by C-mVMC's StdFace. To create new inputs, use C-mVMC's
StdFace generator and copy the resulting `.def` files into a directory;
Julia-mVMC will read them as-is.
