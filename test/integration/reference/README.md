# Integration Test Reference Data

This directory contains C-mVMC reference outputs and the corresponding `.def`
input fixtures used by Julia-mVMC integration tests.

- `test/integration/runtests.jl` verifies the first 10 SR steps against
  selected C outputs with tight tolerances.
- `test/integration/ctest_equivalent.jl` mirrors C-mVMC's standard ctest
  criterion: compare the first two final optimisation summary values against
  `ref_mean.dat` / `ref_std.dat` with the same `3 * std` and `1e-8` gates.

## Provenance

- **C source**: [`issp-center-dev/mVMC`](https://github.com/issp-center-dev/mVMC)
  at commit `5e7ea400ae35b566cfa2de6e342efe962f179a41`
  (`update cmakelists`, the head of `v1.3.0-3-g5e7ea40`).
- **Build flags**: `cmake -DCMAKE_BUILD_TYPE=Release ..`
- **Toolchain**: `gcc-15` (Homebrew), macOS arm64.
- **OpenMP threads**: `OMP_NUM_THREADS=1` for C-mVMC `runtest.py` runs.
- **Run command** (per model):

  ```bash
  cd mVMC/build/test/python
  OMP_NUM_THREADS=1 python3 runtest.py <ModelName>
  ```

  which writes `work/<ModelName>/output/zvo_out_001.dat` and
  `work/<ModelName>/output/zqp_opt.dat`.
- **Strict first-10 reference**: first 10 lines of `zvo_out_001.dat`,
  saved as `<model>/zvo_out_first10.dat`.
- **C ctest reference**: C's `data/<ModelName>/ref/ref_mean.dat` and
  `ref_std.dat`, saved under `<model>/ctest_ref/`.

## Models

| Subdirectory | Mode | C runner argument | Strict first-10 | C ctest-equivalent |
|--------------|------|-------------------|-----------------|-------------------|
| `heisenberg_chain_real` | real | `HeisenbergChain` | yes | yes |
| `heisenberg_chain_real_nsrcg` | real, NSRCG=1 | hand-authored from `HeisenbergChain` | 1-step only | no |
| `hubbard_chain_real` | real | `HubbardChain` | yes | yes |
| `hubbard_tetragonal_real` | real | `HubbardTetragonal` | no | yes |
| `hubbard_tetragonal_momentum_projection_real` | real | `HubbardTetragonal_MomentumProjection` | no | yes |
| `kondo_chain_real` | real | `KondoChain` | no | yes |
| `heisenberg_chain_cmp` | complex | `HeisenbergChain_cmp` | yes | yes |
| `hubbard_chain_cmp` | complex | `HubbardChain_cmp` | no | yes |
| `kondo_chain_cmp` | complex | `KondoChain_cmp` | no | yes |
| `kondo_chain_stot1_cmp` | complex | `KondoChain_Stot1_cmp` | no | yes |
| `general_rbm_cmp` | complex | `GeneralRBM_cmp` | no | yes |
| `heisenberg_chain_fsz` | fsz (generalized orbital) | `HeisenbergChain_fsz` | yes | yes |
| `hubbard_chain_fsz` | fsz (generalized orbital) | `HubbardChain_fsz` | no | yes |
| `kondo_chain_fsz` | fsz (generalized orbital) | `KondoChain_fsz` | no | yes |
| `hubbard_chain_pairhop_real` | real | hand-authored from `HubbardChain` | 1-step PairHop only | no |
| `hubbard_chain_pairhop_fsz` | fsz (generalized orbital) | hand-authored from `HubbardChain_fsz` | 1-step PairHop only | no |

Each subdirectory contains:

- `zvo_out_first10.dat` — first 10 SR steps of the C reference run, for the
  strict first-10 fixtures only.
- `zvo_out_first1.dat` — first optimisation row for the `NSRCG=1` first-step
  fixture.
- `zqp_opt_1step.dat` — C's raw `NSROptItrSmp=1` one-step output for the
  `NSRCG=1` fixture, used to gate the post-CG parameter update.
- `zvo_out_mpi2_first1.dat` / `zvo_SRinfo_mpi2_1step.dat` /
  `zqp_opt_mpi2_1step.dat` — the corresponding `NSRCG=1` `mpiexec -n 2`
  C reference files used by the MPI smoke gate.
- `ctest_ref/ref_mean.dat` and `ctest_ref/ref_std.dat` — reference vectors
  used by the C ctest-equivalent runner.
- `inputs/*.def` — the full set of expert-mode input files used to drive the run.
- `inputs/zqp_opt.dat` / `inputs/initial.def` (when applicable) — initial
  variational parameters used by the C run, mirrored here so that
  Julia-mVMC starts from the same state.

### PairHop first-step fixtures

`hubbard_chain_pairhop_real` and `hubbard_chain_pairhop_fsz` are one-step
fixtures derived from the corresponding Hubbard references with:

- `PairHop pairhopp.def` added to `namelist.def`
- one `pairhopp.def` data row, which C expands internally to both `(i,j)` and
  `(j,i)`
- `NSROptItrStep = 1`
- `NSROptItrSmp = 1`

The committed `zvo_out_first1.dat` files were generated with C-mVMC
`622166afe33c6be3402d7c926db7e9c0003a47c4`, using:

```bash
OMP_NUM_THREADS=1 vmc.out -e namelist.def initial.def
```

## Regenerating

To regenerate this reference data from a fresh clone of the C mVMC tree:

1. Clone `issp-center-dev/mVMC` and check out the commit recorded above.
2. Build:

   ```bash
   mkdir -p build && cd build
   cmake -DCMAKE_BUILD_TYPE=Release ..
   make -j
   ```

3. Run the C tests to populate `work/`:

   ```bash
   cd test/python
   export OMP_NUM_THREADS=1
   python3 runtest.py HeisenbergChain
   python3 runtest.py HubbardChain
   python3 runtest.py HubbardTetragonal
   python3 runtest.py HubbardTetragonal_MomentumProjection
   python3 runtest.py KondoChain
   python3 runtest.py HeisenbergChain_cmp
   python3 runtest.py HubbardChain_cmp
   python3 runtest.py KondoChain_cmp
   python3 runtest.py KondoChain_Stot1_cmp
   python3 runtest.py GeneralRBM_cmp
   python3 runtest.py HeisenbergChain_fsz
   python3 runtest.py HubbardChain_fsz
   python3 runtest.py KondoChain_fsz
   ```

4. Regenerate the C ctest-equivalent fixture data:

   ```bash
   julia --project=@. test/integration/tools/generate_ctest_fixtures.jl \
     --c-test-dir <mVMC-build>/test/python
   ```

5. For strict first-10 fixtures, also refresh `zvo_out_first10.dat`:

   - Copy `work/<Model>/*.def` and optional
     `data/<Model>/zqp_opt.dat` / `data/<Model>/initial.def` into
     `<model>/inputs/`.
   - Take the first 10 lines of `work/<Model>/output/zvo_out_001.dat` and
     save them as `<model>/zvo_out_first10.dat`.

### NSRCG=1 first-step fixture

`heisenberg_chain_real_nsrcg` is a hand-authored first-step fixture derived
from `heisenberg_chain_real` with:

- `NSROptItrStep = 1`
- `NSROptItrSmp = 1`
- `NSRCG = 1`
- `DSROptCGTol = 1.0e-10`
- `NSROptCGMaxIter = 0`
- `NSplitSize = 1`

The committed `zvo_out_first1.dat` and `zqp_opt_1step.dat` were generated with
C-mVMC `mVMC/build-openblas-debug/src/mVMC/vmc.out`, linked against Homebrew
OpenBLAS (`/opt/homebrew/opt/openblas/lib/libopenblas.dylib`), using:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 vmc.out -e namelist.def
```

The committed `zvo_out_mpi2_first1.dat`, `zvo_SRinfo_mpi2_1step.dat`, and
`zqp_opt_mpi2_1step.dat` were generated from the same input and binary with:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 2 vmc.out -e namelist.def
```

This fixture deliberately uses a C OpenBLAS reference rather than the older
Accelerate-backed mac reference. Truncated SR-CG is sensitive to BLAS and
reduction order; C(Accelerate) and C(OpenBLAS) differ by `~3e-3` over 10 steps
for this input, so this fixture is intentionally gated with a coarse parameter
tolerance while keeping the first `zvo_out` row tight.

## PhysCal e2e fixtures (`<model>/physcal_ref/`)

`test/integration/phys_cal_equivalent.jl` runs `run_phys_cal_from_namelist`
(`NVMCCalMode = 1`) for the systems below and compares the three produced Green
files against committed C references with the per-quantity tolerances of
`tools/green_compare.jl` (one-body `rtol = 1e-10`; direct + factored
`rtol = 1e-9`; `atol = 1e-12`).

| Subdirectory | Mode | Fixed params (`zqp_opt.dat`) |
|--------------|------|------------------------------|
| `heisenberg_chain_real/physcal_ref` | real | freshly optimised (NVMCCalMode=0) |
| `heisenberg_chain_cmp/physcal_ref`  | cmp  | freshly optimised (NVMCCalMode=0) |
| `heisenberg_chain_fsz/physcal_ref`  | fsz  | reused `../inputs/zqp_opt.dat` |
| `hubbard_chain_real/physcal_ref`    | real | reused `../inputs/zqp_opt.dat` |
| `hubbard_chain_dh_real/physcal_ref` | real | reused Hubbard params + hand-authored DH2/DH4 slice |
| `kondo_chain_real/physcal_ref`      | real | reused `../inputs/zqp_opt.dat` |

Each `physcal_ref/` contains:

- `inputs/` — the PhysCal input set: the optimisation `.def` files plus a
  hand-authored `greentwoex.def` (factored `TwoBodyGEx` terms over constituents
  already present in `greenone.def`), `modpara.def` with `NVMCCalMode = 1`, and a
  `namelist.def` that adds the `TwoBodyGEx greentwoex.def` line. The DH fixture
  also adds hand-authored `dh2.def` / `dh4.def` ring-neighbour index tables.
  The FSZ fixture keeps the same `Orbital` + `OrbitalParallel` AP+P layout as
  its optimisation fixture so the reused `zqp_opt.dat` parameter order matches C.
  Isolated from the committed optimisation `inputs/` (`NVMCCalMode = 0`), so a
  PhysCal fixture cannot regress the optimisation suite.
- `zqp_opt.dat` — the fixed variational parameters fed to the runner (`opt_para`).
- `expected/zvo_cisajs_001.dat`, `zvo_cisajscktalt_001.dat`,
  `zvo_cisajscktaltex_001.dat` — the committed C Green references.
- `metadata.txt` — per-system provenance (commands, params, Julia-vs-C result).

### PhysCal reference provenance

- **C source**: `issp-center-dev/mVMC`, branch `develop` @
  `66f17422968009f8cc70f1dec94b2f52e562d344` (the canonical integration head).
  The DH and FSZ fixtures were generated with a reference build at
  `622166afe33c6be3402d7c926db7e9c0003a47c4`, which is based on that commit
  plus the OpenMP SIMD benchmark branch and test-data-only changes.
- **Build**: `cmake -DCONFIG=apple -DCMAKE_BUILD_TYPE=Release -DUSE_GEMMT=OFF
  -DUSE_SCALAPACK=OFF -DTesting=OFF` — no BLIS (reference `dskr2k`/`zskr2k`).
- **Toolchain**: Apple Clang 15.0.0 (C/C++) + gfortran 15.2.0 (Homebrew) + libomp,
  Apple Accelerate BLAS, macOS arm64. (The opt-side `ctest_ref` above used gcc-15
  @ master `5e7ea40`; this differs only at ~1e-12, far below the gate tolerance.
  Apple Clang was used because gcc-15.2's fixincludes headers had drifted ahead of
  the installed CommandLineTools SDK on the build host.)
- **OpenMP**: `OMP_NUM_THREADS=1`, single MPI rank.

### Regenerating PhysCal references (per model)

```bash
export OMP_NUM_THREADS=1
VMC=<mVMC-develop-build>/src/mVMC/vmc.out
# (A) fixed params: reuse <model>/inputs/zqp_opt.dat, or optimise (NVMCCalMode=0):
( cd opt_stage && "$VMC" -e namelist.def )          # -> output/zqp_opt.dat
# (B) PhysCal: NVMCCalMode=1 input set incl. greentwoex.def, params as 2nd arg:
( cd phys_stage && "$VMC" -e namelist.def zqp_opt.dat )
# -> output/{zvo_cisajs,zvo_cisajscktalt,zvo_cisajscktaltex}_001.dat
```

## Full Lanczos R1 PhysCal fixtures (`<model>/physcal_ref/`)

`test/integration/lanczos_equivalent.jl` runs `run_phys_cal_from_namelist`
(`NVMCCalMode = 1`, `NLanczosMode = 1`) for the systems below and compares the
produced Lanczos files against committed C references.

| Subdirectory | Mode | Source C fixture | Fixed params (`zqp_opt.dat`) |
|--------------|------|------------------|------------------------------|
| `hubbard_chain_lanczos/physcal_ref` | real | `HubbardChainLanczos` | copied from the C fixture |
| `spin_chain_lanczos/physcal_ref`    | real | `SpinChainLanczos`    | copied from the C fixture |

Each Lanczos `physcal_ref/` contains:

- `inputs/` — expert-mode input files generated from the C standard-mode fixture.
- `zqp_opt.dat` — fixed variational parameters fed to the runner (`opt_para`).
- `expected/zvo_ls_out_001.dat` — C R1 energy, norm, and alpha output.
- `expected/zvo_ls_qqqq_001.dat` — C 16-value flattened QQQQ output.
- `metadata.txt` — per-system provenance and generation-time Julia-vs-C result.

### Regenerating Full Lanczos R1 references (per model)

```bash
export OMP_NUM_THREADS=1
VMC=<mVMC-build>/src/mVMC/vmc.out
# Standard-mode fixture expansion + PhysCal run:
( cd stage && "$VMC" -s StdFace.def zqp_opt.dat )
# -> expert .def files and output/{zvo_ls_out,zvo_ls_qqqq}_001.dat
#
# The committed expert input set can also be replayed directly:
( cd phys_stage && "$VMC" -e namelist.def zqp_opt.dat )
```
