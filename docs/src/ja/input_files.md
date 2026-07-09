# 入力ファイル

[English](../en/input_files.md)

Julia-mVMC は C-mVMC expert-mode の `.def` file を読み込みます。entry point は
`namelist.def` で、そこに書かれた path は `namelist.def` のある directory からの
相対 path として解決します。

この page では Julia-mVMC v0.5.0 の対応状況をまとめます。入力構文や keyword の
詳細は C-mVMC manual の
[Input files for Expert mode](https://issp-center-dev.github.io/mVMC/docs/expert.html)
および [Input files for Standard mode](https://issp-center-dev.github.io/mVMC/docs/standard.html)
を参照してください。Julia-mVMC は expert-mode `.def` file を直接読み込みます。
Standard-mode file が必要な場合は、C-mVMC の `vmcdry.out` workflow で expert-mode
file に変換して利用します。

## Core files

| File / keyword | Purpose | v0.5.0 status |
|----------------|---------|---------------|
| `namelist.def` | 他の入力 file の一覧 | Supported |
| `modpara.def` (`ModPara`) | system size, sampling, SR setting, RNG seed | Supported |
| `locspn.def` (`LocSpin`) | local-spin site flags | Supported |
| `trans.def` (`Trans`) | hopping / transfer terms | Supported |
| `coulombintra.def` (`CoulombIntra`) | on-site Coulomb interaction | Supported |
| `coulombinter.def` (`CoulombInter`) | inter-site Coulomb interaction | Supported |
| `hund.def` (`Hund`) | Hund coupling | Supported |
| `exchange.def` (`Exchange`) | spin exchange | Supported |
| `interall.def` (`InterAll`) | general four-fermion interaction | spin metadata fallback 付きで supported |
| `pairhop.def` / `pairhopp.def` (`PairHop`) | pair hopping | non-FSZ / FSZ path で supported |

## Variational-parameter / projection files

| File / keyword | Purpose | v0.5.0 status |
|----------------|---------|---------------|
| `gutzwilleridx.def` (`Gutzwiller`) | Gutzwiller correlator indexing | Supported |
| `jastrowidx.def` (`Jastrow`) | Jastrow correlator indexing | Supported |
| `dh2.def`, `dh4.def` (`DH2`, `DH4`) | doublon-holon projection slices | Supported |
| `orbitalidx.def` (`Orbital`) | Slater orbital indexing | Supported |
| `orbitalidxgen.def` (`OrbitalGeneral`) | FSZ/generalized orbital indexing | Supported |
| `orbitalidxpara.def` (`OrbitalParallel`) | parallel-orbital block | Supported |
| `qptransidx.def` (`TransSym`) | quantum-number projection symmetry | Supported |
| `opttrans.def` (`OptTrans`) | OptTrans mapping and weights | `NSplitSize = 1` で supported。split restriction あり。 |
| `initial.def` | initial variational parameters, RBM triples | 存在すれば auto-detect |

## Measurement files

| File / keyword | Output | v0.5.0 status |
|----------------|--------|---------------|
| `greenone.def` (`OneBodyG`) | `zvo_cisajs_*.dat` | Supported |
| `greentwo.def` (`TwoBodyG`) | `zvo_cisajscktalt_*.dat` | Supported |
| `greentwoex.def` (`TwoBodyGEx`) | `zvo_cisajscktaltex_*.dat` | verified PhysCal path で supported |

## `In*.def` overlay

`In*.def` file は `initial.def` 適用後に parameter value を上書きします。これは C の
`InitParameter -> ReadInitParameter -> ReadInputParameters` の順序に対応します。
`initial.def` と overlay が同じ parameter を持つ場合、overlay が優先されます。

対応 overlay には `InGutzwiller.def`, `InJastrow.def`, `InDH2.def`,
`InDH4.def`, `InOrbital*.def`, `InOptTrans.def`, RBM layer overlay が含まれます。

## 認識するが未実装の file

`spinjastrow.def` は parser が認識しますが、sampling / optimization pipeline では
未実装です。有効化された入力は reject されます。

## RNG seed rule

`modpara.def` の `RndSeed` は C-parity rule に従います。

| Value | Meaning |
|-------|---------|
| missing | C default `11272` |
| `0` | seed `0` |
| negative | time-derived seed |
| positive | その seed value |

MPI では C-compatible な per-group offset も適用します。
