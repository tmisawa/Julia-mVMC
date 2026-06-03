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
| `hubbard_chain_real` | real | `HubbardChain` | yes | yes |
| `hubbard_tetragonal_real` | real | `HubbardTetragonal` | no | yes |
| `hubbard_tetragonal_momentum_projection_real` | real | `HubbardTetragonal_MomentumProjection` | no | yes |
| `kondo_chain_real` | real | `KondoChain` | no | yes |
| `heisenberg_chain_cmp` | complex | `HeisenbergChain_cmp` | yes | yes |
| `hubbard_chain_cmp` | complex | `HubbardChain_cmp` | no | yes |
| `kondo_chain_cmp` | complex | `KondoChain_cmp` | no | yes |
| `kondo_chain_stot1_cmp` | complex | `KondoChain_Stot1_cmp` | no | yes |
| `heisenberg_chain_fsz` | fsz (generalized orbital) | `HeisenbergChain_fsz` | yes | yes |
| `hubbard_chain_fsz` | fsz (generalized orbital) | `HubbardChain_fsz` | no | yes |
| `kondo_chain_fsz` | fsz (generalized orbital) | `KondoChain_fsz` | no | yes |

Each subdirectory contains:

- `zvo_out_first10.dat` — first 10 SR steps of the C reference run, for the
  strict first-10 fixtures only.
- `ctest_ref/ref_mean.dat` and `ctest_ref/ref_std.dat` — reference vectors
  used by the C ctest-equivalent runner.
- `inputs/*.def` — the full set of expert-mode input files used to drive the run.
- `inputs/zqp_opt.dat` / `inputs/initial.def` (when applicable) — initial
  variational parameters used by the C run, mirrored here so that
  Julia-mVMC starts from the same state.

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
   python3 runtest.py HeisenbergChain_fsz
   python3 runtest.py HubbardChain_fsz
   python3 runtest.py KondoChain_fsz
   ```

4. Regenerate the C ctest-equivalent fixture data:

   ```bash
   julia --project=@. test/integration/tools/generate_ctest_fixtures.jl \
     --c-test-dir ../private-mVMC/mVMC/build/test/python
   ```

5. For strict first-10 fixtures, also refresh `zvo_out_first10.dat`:

   - Copy `work/<Model>/*.def` and optional
     `data/<Model>/zqp_opt.dat` / `data/<Model>/initial.def` into
     `<model>/inputs/`.
   - Take the first 10 lines of `work/<Model>/output/zvo_out_001.dat` and
     save them as `<model>/zvo_out_first10.dat`.
