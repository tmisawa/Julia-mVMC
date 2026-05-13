# 1. Installation

## Verified environment

- **Julia**: 1.11+ (CI tests 1.11 and 1.12)
- **OS**: macOS, Linux. Windows is **not verified**; reports welcome.

## Native prerequisites

`PfaPack.jl` and `SFMT.jl` (vendored as git submodules) build native
libraries on first `Pkg.add` / `Pkg.build`. The system needs a Fortran
compiler, a C++ compiler, `make`, and BLAS/LAPACK.

Ubuntu / Debian:

```bash
sudo apt-get install -y gfortran g++ make libblas-dev liblapack-dev
```

macOS (with Homebrew):

```bash
brew install gcc gfortran openblas
# Apple's libblas / liblapack already cover LAPACK at runtime; gfortran
# from Homebrew is required because Apple ships only Clang.
```

## Install (clone-based, with submodules)

The supported install path for v0.1 is to clone the repo **with
submodules** and activate the workspace project:

```bash
git clone --recurse-submodules https://github.com/tmisawa/Julia-mVMC
cd Julia-mVMC
julia --project=@. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
```

If you already cloned without `--recurse-submodules`, run
`git submodule update --init --recursive` from the repo root before
`Pkg.instantiate()`.

The Julia-mVMC workspace consists of two in-repo subpackages
(`MVMCOptimizers.jl/`, `MVMCExpertModeParsers.jl/`) plus two git
submodules (`PfaPack.jl/`, `SFMT.jl/`). They are tied together by the
`[sources]` block in the workspace root `Project.toml` using **relative
paths**. Running `julia --project=@.` from the workspace root puts all
four on the load path simultaneously.

## Why no `Pkg.add(url=..., subdir=...)`?

`Pkg.add(url=..., subdir="MVMCOptimizers.jl")` does **not** work for v0.1
and is intentionally not documented as an install path. The reason is
that `MVMCOptimizers.jl/Project.toml` declares its siblings via
relative paths in `[sources]`:

```toml
[sources]
MVMCExpertModeParsers = {path = "../MVMCExpertModeParsers.jl"}
PfaPack = {path = "../PfaPack.jl"}
SFMT = {path = "../SFMT.jl"}
```

When `Pkg` clones the repo and looks only at the `MVMCOptimizers.jl/`
subdirectory, it cannot resolve `../MVMCExpertModeParsers.jl` against
the larger checkout, and the submodules are not pulled in either, so
the dependency resolver fails (the subpackages are not on a registry).

The clone-based workflow above sidesteps this entirely. Adding a
URL-based install path is a v0.2+ concern and would require committing
the GitHub URL into each subpackage's `[sources]` block.

## Smoke test

```bash
JULIA_MVMC_EXAMPLE_STEPS=5 julia --project=@. examples/heisenberg_chain_real.jl
```

Expected: a `Final energy / site = ...` line and exit 0 within seconds.
Run all four [examples](../../examples/) the same way.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `Pkg.build` fails on `PfaPack` with `gfortran: not found` | Fortran compiler missing | Install `gfortran` (see above). |
| Link error on `dgemm_` / `zgemm_` | BLAS/LAPACK not on the linker path | Install `libblas-dev liblapack-dev` (or equivalent), then `Pkg.build("PfaPack")`. |
| `UndefVarError: SFMT19937RNG` | Build artifacts stale after Julia upgrade | `using Pkg; Pkg.build("SFMT")` and restart Julia. |
| `PfaPack.jl/` or `SFMT.jl/` directories are empty after clone | Submodules not initialised | `git submodule update --init --recursive`. |
| `julia: command not found` | Julia not installed or not on PATH | Install from <https://julialang.org/downloads/>. v1.11+ required. |
| Tests pass but examples fail with `LoadError` referencing `[sources]` | Running with the wrong project | Use `--project=@.` from the repo root, **not** `--project=MVMCOptimizers.jl`. |
