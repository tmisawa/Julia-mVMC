# MVMCOptimizers.jl

> Part of [**Julia-mVMC**](https://github.com/tmisawa/Julia-mVMC) — the Julia port of [issp-center-dev/mVMC](https://github.com/issp-center-dev/mVMC). See the top-level repository for installation guidance and the full user manual.

The main VMC optimization package: Stochastic-Reconfiguration parameter optimization (equivalent to C mVMC's `VMCParaOpt()`) and physical-quantity calculation (equivalent to `VMCPhysCal()`). Works against `ExpertModeData` produced by [`MVMCExpertModeParsers.jl`](../MVMCExpertModeParsers.jl), with Pfaffian / inverse routines from [`PfaPack.jl`](../PfaPack.jl) and the C-compatible RNG from [`SFMT.jl`](../SFMT.jl).

## Scope (v0.4)

- **`VMCParaOpt` (parameter optimization)** — verified against the C reference by strict first-10-step integration checks, C ctest-equivalent gates for supported standard fixtures, and a serial `NSRCG = 1` first-step tolerance gate; see the integration tests at `../test/integration/`.
- **`VMCPhysCal` (physical-quantity calculation)** — experimental. One-body (`zvo_cisajs`), direct two-body (`TwoBodyG` → `zvo_cisajscktalt`), and factored/product two-body (`TwoBodyGEx`/`greentwoex.def` → `zvo_cisajscktaltex`) Green functions are supported (factored is non-FSZ), including DH2/DH4-present fixtures, and gated against C references via [`../test/integration/phys_cal_equivalent.jl`](../test/integration/phys_cal_equivalent.jl); run through [`run_phys_cal_from_namelist`](src/run_phys_cal_from_namelist.jl). See [`../docs/manual/04_physics_calc.md`](../docs/manual/04_physics_calc.md).
- **Threading** — conservative `JULIA_MVMC_INNER_THREADS=1` opt-ins cover selected inner loops. Sample-level `VMCMainCal` threading is intentionally disabled for C-parity; `JULIA_MVMC_PFAPACK_THREADS=1` remains a debug/benchmark triage mode only.
- **MPI** — v0.4 supports multi-rank sample-parallel execution with `NSplitSize = 1` through MPI.jl-compatible launchers. `VMCParaOpt` MPI currently supports the direct SR solver only (`NSRCG = 0`); C's grouped MPI/QP split (`NSplitSize > 1`) and MPI CG solver runs (`NSRCG != 0`) are still rejected.
- **Not supported in this release**: BackFlow (`vmc_bf_*` entry points raise an error), full Lanczos (only step-0 is comparable), grouped MPI/QP splitting (`NSplitSize > 1`), and MPI CG solver runs (`NSRCG != 0` under MPI).

## Installation

This subpackage is **not** published as a standalone registered package and `Pkg.add("MVMCOptimizers")` will not work in this release. It is intended to be used as part of the Julia-mVMC workspace via the root [`Project.toml`](../Project.toml) `[sources]` block. From a clone of [tmisawa/Julia-mVMC](https://github.com/tmisawa/Julia-mVMC):

```bash
git clone --recurse-submodules https://github.com/tmisawa/Julia-mVMC
cd Julia-mVMC
julia --project=@. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
```

After that, `using MVMCOptimizers` works inside any code activated against the root workspace.

## Quick start

Most workflows go through `run_para_opt_from_namelist`, which mirrors the C `vmc.out -s namelist.def` driver. The four runnable examples in [`../examples/`](../examples/) all use it:

```bash
julia --project=@. ../examples/heisenberg_chain_real.jl
```

Inline form:

```julia
using MVMCOptimizers

result = run_para_opt_from_namelist("path/to/namelist.def";
    nsteps = 50,
    nsmp = 50,                   # NSROptItrSmp override; must satisfy nsteps >= nsmp
    mode = :real,                # :real, :cmp, or :fsz
    output_dir = "out",
    seed = 11272,                # SFMT19937 seed (overrides modpara.def's RndSeed)
    initial_def = :auto,         # :auto picks up initial.def if present
)

println("final energy / site = ", result.final_energy_per_site)
```

Lower-level entry points are also exported for callers that want to manage the
state themselves:

```julia
data = MVMCExpertModeParsers.parse_expert_mode_files("namelist.def")
# ... seed the RNG, run initialize_parameters! / read_input_parameters! /
#     read_initial_def! / init_qp_weight! in the same order as
#     run_para_opt_from_namelist (see C `vmcmain.c:256-281`) ...
info = vmc_para_opt!(data;
    callback = (step, data, energy, info) -> step % 10 == 0 &&
        println("step $step  E = ", real(energy)),
)
```

## Public API

| Function | Purpose |
|----------|---------|
| `run_para_opt_from_namelist(path; nsteps, mode, nsmp=nothing, ...)` | High-level driver. Parses `.def` files, seeds the RNG, runs `nsteps` SR steps, writes `zvo_out.dat` / `zqp_opt.dat`. `nsmp` overrides `NSROptItrSmp` for final-sample averaging and must satisfy `nsteps >= nsmp`. |
| `vmc_para_opt!(data; callback=nothing)` | Lower-level SR optimization loop (no I/O, no parser). |
| `vmc_phys_cal!(data; rnd_seed=0)` | Physical-quantity calculation (experimental in this release). |
| `read_initial_def!(data, path)` | Overlay variational parameters from an `initial.def` file. |

## C mVMC compatibility

| C function | Julia entry point |
|------------|-------------------|
| `VMCParaOpt` (in `vmcmain.c`) | `vmc_para_opt!` |
| `VMCPhysCal` | `vmc_phys_cal!` |
| `VMC_MakeSample` / `_real` / `_fsz` / `_fsz_real` | `vmc_make_sample!` / `vmc_make_sample_real!` / `vmc_make_sample_fsz!` / `vmc_make_sample_fsz_real!` |
| `VMC_MainCal` / `_fsz` | `vmc_main_cal!` / `vmc_main_cal_fsz!` |
| `StochasticOpt` / `StochasticOptCG` | `stochastic_opt!` / `stochastic_opt_cg!` |
| `WeightAverageWE` / `WeightAverageSROpt` | `weight_average_we!` / `weight_average_sr_opt!` |
| `SyncModifiedParameter` | `sync_modified_parameter!` |
| `UpdateQPWeight` | `update_qp_weight!` |
| `VMC_BF_MakeSample` / `_real` / `VMC_BF_MainCal` | stubs — **raise an error in this release** |

The SR step state is held in `VMCOptimizationState` (defined in `src/types.jl`),
which bundles `EnergyData`, `SlaterMatrixData`, `ElectronConfiguration`, and
`SROptData`. These mirror the C globals on a one-to-one basis.

## Testing

```bash
cd MVMCOptimizers.jl
julia --project=@. -e 'using Pkg; Pkg.test()'
```

The integration tests in the root workspace (`test/integration/runtests.jl` from the repository root) additionally exercise this package end-to-end against four fixture models with bit-level comparison to the C reference.

## License

This subpackage is part of [Julia-mVMC](https://github.com/tmisawa/Julia-mVMC) and is licensed under **GPL-3.0-or-later**. See the subpackage [`LICENSE`](LICENSE) (a verbatim copy of GPLv3) and the workspace-level [`THIRD_PARTY_LICENSES.md`](../THIRD_PARTY_LICENSES.md) for bundled-source notices.

## References

- [mVMC](https://github.com/issp-center-dev/mVMC) — original C reference implementation.
- [`../docs/manual/`](../docs/manual/) — installation, input files, optimization, physics calculation, C/Julia compatibility notes.
