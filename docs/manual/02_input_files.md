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
| `dh2.def` / `dh4.def` (`DH2`, `DH4`) | doublon-holon index tables and projection slices | ✅ |
| `orbitalidx.def` (`Orbital`) | Slater orbital indexing | ✅ |
| `orbitalidxgen.def` (`OrbitalGeneral`) | fsz mode generalised orbital | ✅ |
| `orbitalidxpara.def` (`OrbitalParallel`) | parallel-orbital block | ✅ |
| `qptransidx.def` (`TransSym`) | quantum projection symmetry | ✅ |
| `opttrans.def` (`OptTrans`) | optional QPOptTrans mapping and OptTrans weights | ✅ |
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
| `InDH2.def`, `InDH4.def` | doublon-holon projection parameter vectors | ✅ |
| `InOrbital.def` / `InOrbitalAntiParallel.def` | `orbital_terms[i].value` | ✅ |
| `InOrbitalParallel.def` | `orbital_terms[i].value` after the anti-parallel offset | ✅ |
| `InOrbitalGeneral.def` | `orbital_terms[i].value` (fsz layout) | ✅ |
| `InOptTrans.def` | `opt_trans` | ✅ |
| `InChargeRBM_PhysLayer.def`, `_HiddenLayer.def`, `_PhysHidden.def` | corresponding RBM term arrays | ✅ |
| `InSpinRBM_*` and `InGeneralRBM_*` (3 layers each) | corresponding RBM term arrays | ✅ |

## Not yet supported

These are recognised by the parser but **not** consumed by the
optimisation or sampling pipeline. Reference data generated with them
active is not reproducible bit-for-bit.

| File | Status / reason |
|------|--------|
| `spinjastrow.def` (`SpinJastrow`) | Not implemented; parser hard-fails if the keyword is present because it would change projection offsets. |

> Note: the factored/product two-body Green (`TwoBodyGEx` / `greentwoex.def` →
> `zvo_cisajscktaltex`) **is** supported for sz-conserved and FSZ/general-orbital
> `NSplitSize = 1` PhysCal runs and gated against C; see
> [`04_physics_calc.md`](04_physics_calc.md).

## Compatibility notes

- Section headers (`NGutzwillerIdx`, `NTransfer`, etc.) follow the C
  convention: 5 ignored lines (header + counter + 3 separators) before
  the data block.
- `qptransidx.def` accepts the 3-column form (`itmpsgn = 1` defaulted)
  in addition to the canonical 4-column form.
- `pairhop.def` / `pairhopp.def` (`PairHop`) is consumed by the Hamiltonian in
  both non-FSZ and FSZ paths. Each input row is expanded to `(i,j)` and `(j,i)`,
  matching C-mVMC's internal `PairHopping` list.
- For `interall.def`, fsz-style spin-flip terms are honoured.
- The RNG layout matches C: `RndSeed` from `modpara.def` seeds an
  SFMT19937 stream with the C-parity rule (v0.4): missing line → `11272`
  (C `readdef.c` default), `0` → `0`, negative → a time-derived seed
  (rank 0, broadcast under MPI), positive → the value; under MPI the
  per-group offset `+ group1` is added (C `vmcmain.c:257`). See
  `MVMCOptimizers.run_para_opt_from_namelist` for the precise phase
  ordering (`init_parameter!` → `read_initial_def!` →
  `read_input_parameters!` → `sync_modified_parameter!` →
  `init_qp_weight!`), which mirrors C's `vmcmain.c:264-281`.

## Generating inputs

Most samples in [`examples/inputs/`](../../examples/inputs/) and
[`test/integration/reference/`](../../test/integration/reference/) were
produced by C-mVMC's StdFace. DH2/DH4 fixtures add hand-authored DH index
tables on top because StdFace does not emit DH files. To create new non-DH
inputs, use C-mVMC's StdFace generator and copy the resulting `.def` files into
a directory; Julia-mVMC will read them as-is.
