# Julia-mVMC

Julia port of the [mVMC](https://github.com/issp-center-dev/mVMC) (many-variable Variational Monte Carlo) solver for quantum lattice models.

## Status (v0.4)

| Component | Status | Notes |
|-----------|--------|-------|
| VMCParaOpt (parameter optimization) | ✅ Verified | Strict first-10-step C-reference checks plus C ctest-equivalent gates for supported standard fixtures; see `test/integration/`. |
| VMCPhysCal (physical quantities) | 🚧 Experimental | C-referenced one-body, direct two-body, and factored/product two-body Green-function output for supported fixtures. |
| Shared-memory threading | 🚧 Experimental | Conservative inner-loop opt-ins only; sample-level `VMCMainCal` threading is intentionally disabled for C-parity. |
| Lanczos | ⚠️ Step-0 only | Single-step comparison verified; full Lanczos not yet ported. |
| BackFlow | ❌ Not supported | Planned for a future release. |
| MPI parallelization | 🚧 Experimental (v0.4 R1) | `VMCParaOpt` (`NSRCG = 0` direct SR solver only) and `VMCPhysCal` run under MPI.jl-compatible launchers with `NSplitSize = 1`; rank0 output/readback and comm0 reductions are smoke-tested for `mpiexec -n 2/-n 4`. C's grouped MPI/QP split (`NSplitSize > 1`) and MPI CG solver (`NSRCG != 0`) are still rejected. |

## Installation

Requires **Julia 1.11+**, `gfortran`, `g++`, `make`, BLAS/LAPACK.

The supported install path in this release is to clone the repo **with submodules** and activate the workspace project:

```bash
git clone --recurse-submodules https://github.com/tmisawa/Julia-mVMC
cd Julia-mVMC
julia --project=@. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
```

If you cloned without `--recurse-submodules`, run `git submodule update --init --recursive` before `Pkg.instantiate()`.

`Pkg.add(url=..., subdir=...)` is **not** a supported install path in this release because Julia-mVMC and its dependencies are coupled via the workspace `[sources]` block in the root `Project.toml` (relative paths into the submodules, not committed URLs). See [docs/manual/01_install.md](docs/manual/01_install.md) for native prerequisites and detailed setup.

## Quickstart

After completing the [Installation](#installation) step above (`Pkg.instantiate()` and `Pkg.build()` must have run at least once on a clean clone), run:

```bash
julia --project=@. examples/heisenberg_chain_real.jl
```

Expected output ends with a line like `Final energy / site = -0.44...`. Each example runs a 50-step VMCParaOpt by default (override with `JULIA_MVMC_EXAMPLE_STEPS`).

## Documentation

[`docs/manual/`](docs/manual/) — installation, input files, optimization, physics calculations, and C/Julia compatibility notes.

[`CHANGELOG.md`](CHANGELOG.md) — release notes.

## Repository structure

```
Julia-mVMC/
├── MVMCOptimizers.jl/        # main package (VMCParaOpt, VMCPhysCal)
├── MVMCExpertModeParsers.jl/ # .def file parsers
├── PfaPack.jl/               # submodule: Pfaffian (Fortran/C++ wrapper)
├── SFMT.jl/                  # submodule: SFMT random number generator (C wrapper)
├── examples/                 # 4 runnable examples
├── docs/manual/              # user manual
└── test/integration/         # C reference comparison tests
```

`PfaPack.jl/` and `SFMT.jl/` are external Julia packages developed in
their own repositories (<https://github.com/tmisawa/PfaPack.jl>,
<https://github.com/tmisawa/SFMT.jl>) and vendored here as git
submodules under non-GPL open-source licenses (BSD-3-Clause + MPL-2.0
for PfaPack, BSD-3-Clause for SFMT). MPL-2.0 is a file-scoped weak
copyleft, not a permissive license — see
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) for the per-license
notices.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE) and [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).

## Citation

See [CITATION.cff](CITATION.cff). Plain-text:

> Misawa, T. (2026). Julia-mVMC v0.3.0. https://github.com/tmisawa/Julia-mVMC

## Acknowledgments

The [PfaPack.jl](https://github.com/tmisawa/PfaPack.jl) and [SFMT.jl](https://github.com/tmisawa/SFMT.jl) Julia wrappers used here as submodules were primarily authored by **Satoshi Terasaki** ([AtelierArith](https://atelier-arith.jp/)). Their original placement inside Julia-mVMC as `MVMCPfaPack.jl/` and `SFMT19937.jl/` predates the initial v0.1 release; they were extracted to standalone repositories under non-GPL open-source licenses (BSD-3-Clause for SFMT; BSD-3-Clause + MPL-2.0 for PfaPack) with Terasaki's consent.

The mVMC C reference implementation is developed at ISSP, University of Tokyo: <https://github.com/issp-center-dev/mVMC>.

## Contributing

Issues and pull requests are welcome.
