# Input Files

[日本語](../ja/input_files.md)

Julia-mVMC reads C-mVMC expert-mode `.def` files. The entry point is
`namelist.def`, and paths listed there are resolved relative to the directory
containing `namelist.def`.

This page summarizes the Julia-mVMC v0.5.0 support status. For detailed input
syntax and keyword definitions, see the C-mVMC manual:
[Input files for Expert mode](https://issp-center-dev.github.io/mVMC/docs/expert.html)
and [Input files for Standard mode](https://issp-center-dev.github.io/mVMC/docs/standard.html).
Julia-mVMC consumes expert-mode `.def` files directly; Standard-mode files can
be converted with C-mVMC's `vmcdry.out` workflow when needed.

## Core files

| File / keyword | Purpose | v0.5.0 status |
|----------------|---------|---------------|
| `namelist.def` | Lists all other input files | Supported |
| `modpara.def` (`ModPara`) | system size, sampling, SR settings, RNG seed | Supported |
| `locspn.def` (`LocSpin`) | local-spin site flags | Supported |
| `trans.def` (`Trans`) | hopping / transfer terms | Supported |
| `coulombintra.def` (`CoulombIntra`) | on-site Coulomb interaction | Supported |
| `coulombinter.def` (`CoulombInter`) | inter-site Coulomb interaction | Supported |
| `hund.def` (`Hund`) | Hund coupling | Supported |
| `exchange.def` (`Exchange`) | spin exchange | Supported |
| `interall.def` (`InterAll`) | general four-fermion interaction | Supported with documented spin-metadata fallback |
| `pairhop.def` / `pairhopp.def` (`PairHop`) | pair hopping | Supported for non-FSZ and FSZ paths |

## Variational-parameter and projection files

| File / keyword | Purpose | v0.5.0 status |
|----------------|---------|---------------|
| `gutzwilleridx.def` (`Gutzwiller`) | Gutzwiller correlator indexing | Supported |
| `jastrowidx.def` (`Jastrow`) | Jastrow correlator indexing | Supported |
| `dh2.def`, `dh4.def` (`DH2`, `DH4`) | doublon-holon projection slices | Supported |
| `orbitalidx.def` (`Orbital`) | Slater orbital indexing | Supported |
| `orbitalidxgen.def` (`OrbitalGeneral`) | FSZ/generalized orbital indexing | Supported |
| `orbitalidxpara.def` (`OrbitalParallel`) | parallel-orbital block | Supported |
| `qptransidx.def` (`TransSym`) | quantum-number projection symmetry | Supported |
| `opttrans.def` (`OptTrans`) | OptTrans mapping and weights | Supported for `NSplitSize = 1`; split restrictions apply |
| `initial.def` | initial variational parameters, including RBM triples | Auto-detected when present |

## Measurement files

| File / keyword | Output | v0.5.0 status |
|----------------|--------|---------------|
| `greenone.def` (`OneBodyG`) | `zvo_cisajs_*.dat` | Supported |
| `greentwo.def` (`TwoBodyG`) | `zvo_cisajscktalt_*.dat` | Supported |
| `greentwoex.def` (`TwoBodyGEx`) | `zvo_cisajscktaltex_*.dat` | Supported for verified PhysCal paths |

## `In*.def` overlays

`In*.def` files overwrite parameter values after `initial.def` is applied,
matching the C ordering `InitParameter -> ReadInitParameter ->
ReadInputParameters`. If both `initial.def` and an overlay provide the same
parameter, the overlay wins.

Supported overlays include `InGutzwiller.def`, `InJastrow.def`, `InDH2.def`,
`InDH4.def`, `InOrbital*.def`, `InOptTrans.def`, and RBM layer overlays.

## Recognized but not implemented

`spinjastrow.def` is recognized by the parser but not implemented in the
sampling and optimization pipeline. Inputs that activate it are rejected.

## RNG seed rule

`RndSeed` in `modpara.def` follows the C-parity rule:

| Value | Meaning |
|-------|---------|
| missing | C default `11272` |
| `0` | seed `0` |
| negative | time-derived seed |
| positive | that seed value |

Under MPI, Julia-mVMC also applies the C-compatible per-group offset.
