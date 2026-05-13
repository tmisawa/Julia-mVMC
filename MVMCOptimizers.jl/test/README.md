# MVMCOptimizers.jl tests

Subpackage-local test suite. Runs in seconds and has **no external file
dependencies** — fixtures live in `test_unit/helpers/` as in-memory builders,
and integration tests against the C reference live separately at the
workspace root (`test/integration/runtests.jl`).

## Run

From the `MVMCOptimizers.jl/` directory:

```bash
julia --project=@. -e 'using Pkg; Pkg.test()'
```

`runtests.jl` exercises:

1. **Legendre polynomial / Gauss-Legendre quadrature** — sanity checks on
   the quadrature primitives used by `init_qp_weight!`.
2. **Slater Update Tests** (`test_slater_update.jl`) — single-update / inverse
   refresh consistency.
3. **Unit tests** under `../test_unit/` (covered by `test_unit/INDEX.md`):
   - `test_unit_stochastic_opt.jl` — SR matrix / gradient / parameter update
     map (RBM 9-section layout included).
   - `test_unit_vmc_sampling_*.jl` — incremental ↔ full re-evaluation parity
     for RBM / projection counters / log-ratio computations.
   - `test_unit_slater_update.jl` — orbital index ↔ qp_trans matrices.
   - `test_unit_vmc_main_cal_sr.jl` — SR diff,
     `calculate_oo!` / `calculate_oo_real!`.
   - `test_unit_parameter_sync.jl` — Gutzwiller/Jastrow shift, Slater
     rescale, and an explicit assertion that `data.para_qp_trans` is
     **not** rescaled by `sync_modified_parameter!` (C normalises the
     dedicated `OptTrans[]` array, not the QPTrans phase factors).
   - `test_unit_types.jl` — size invariants.

## Integration tests against the C reference

Bit-level comparison against zvo_out outputs from C-mVMC runs lives at the
**workspace root**, not here:

```bash
# from the workspace root
julia --project=@. test/integration/runtests.jl
```

Inputs and reference outputs are bundled under
`test/integration/reference/<model>/`, so no external C build is needed.
See `test/integration/reference/README.md` for provenance details.

## Layout

```
MVMCOptimizers.jl/test/
├── README.md             (this file)
├── runtests.jl           (entry point — Pkg.test() runs this)
├── test_slater_update.jl (Slater single-update tests)
└── samples/              (legacy in-tree sample inputs; see workspace
                          test/integration/reference/ for the actively
                          maintained fixture set)
```
