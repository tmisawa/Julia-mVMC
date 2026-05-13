# Integration Test Reference Data

This directory contains C-mVMC reference outputs and the corresponding `.def`
input fixtures used by `test/integration/runtests.jl` to verify that
Julia-mVMC reproduces the C implementation at the bit level.

## Provenance

- **C source**: [`issp-center-dev/mVMC`](https://github.com/issp-center-dev/mVMC)
  at commit `5e7ea400ae35b566cfa2de6e342efe962f179a41`
  (`update cmakelists`, the head of `v1.3.0-3-g5e7ea40`).
- **Build flags**: `cmake -DCMAKE_BUILD_TYPE=Release ..`
- **Toolchain**: `gcc-15` (Homebrew), macOS arm64.
- **Run command** (per model):

  ```bash
  cd mVMC/build/test/python
  python3 runtest.py <ModelName>
  ```

  which writes `work/<ModelName>/output/zvo_out_001.dat`.
- **Reference value bundled here**: first 10 lines of `zvo_out_001.dat`,
  saved as `<model>/zvo_out_first10.dat`.

## Models

| Subdirectory | Mode | C runner argument |
|--------------|------|-------------------|
| `heisenberg_chain_real` | real | `HeisenbergChain` |
| `heisenberg_chain_cmp` | complex | `HeisenbergChain_cmp` |
| `heisenberg_chain_fsz` | fsz (generalized orbital) | `HeisenbergChain_fsz` |
| `hubbard_chain_real` | real | `HubbardChain` |

Each subdirectory contains:

- `zvo_out_first10.dat` — first 10 SR steps of the C reference run.
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
   python3 runtest.py HeisenbergChain
   python3 runtest.py HeisenbergChain_cmp
   python3 runtest.py HeisenbergChain_fsz
   python3 runtest.py HubbardChain
   ```

4. For each model `<Model>` corresponding to `<model>` in the table above:

   - Copy `work/<Model>/*.def` and `work/<Model>/zqp_opt.dat`/`initial.def`
     (when present) into `<model>/inputs/`.
   - Take the first 10 lines of `work/<Model>/output/zvo_out_001.dat` and
     save them as `<model>/zvo_out_first10.dat`.
