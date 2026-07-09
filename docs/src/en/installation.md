# Installation

[日本語](../ja/installation.md)

## Verified environment

| Item | Requirement |
|------|-------------|
| Julia | 1.11 or newer. CI tests Julia 1.11 and 1.12. |
| OS | Linux and macOS. Windows is not verified in this release. |
| Native tools | `gfortran`, `g++`, `make`, BLAS, and LAPACK. |

`PfaPack.jl` and `SFMT.jl` are git submodules and build native libraries during
the first package build.

## Native prerequisites

Ubuntu / Debian:

```bash
sudo apt-get update
sudo apt-get install -y gfortran g++ make libblas-dev liblapack-dev
```

macOS with Homebrew:

```bash
brew install gcc
```

Apple's system BLAS/LAPACK is enough at runtime. Homebrew `gfortran` is required
because Apple does not ship a Fortran compiler.

## Clone and instantiate

The supported installation path is a workspace checkout with submodules:

```bash
git clone --recurse-submodules https://github.com/tmisawa/Julia-mVMC
cd Julia-mVMC
julia --project=@. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
```

If the repository was cloned without submodules, initialize them first:

```bash
git submodule update --init --recursive
```

Run Julia from the workspace root with `--project=@.`. The root project ties
together `MVMCOptimizers.jl`, `MVMCExpertModeParsers.jl`, `PfaPack.jl`, and
`SFMT.jl` through relative `[sources]` entries.

## Why clone-based installation?

`Pkg.add(url=..., subdir="MVMCOptimizers.jl")` is not supported in v0.5.0. The
subpackage projects refer to sibling packages and submodules with relative
paths, which are only valid in the full workspace checkout.

## Smoke test

```bash
JULIA_MVMC_EXAMPLE_STEPS=5 julia --project=@. examples/heisenberg_chain_real.jl
```

The run should finish with a `Final energy / site = ...` line.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `gfortran: not found` | Fortran compiler missing | Install `gfortran`. |
| macOS linker reports missing `System` | Homebrew `gfortran` cannot find the SDK | Set `SDKROOT="$(xcrun --sdk macosx --show-sdk-path)"` and rebuild. |
| BLAS symbols such as `dgemm_` are missing | BLAS/LAPACK development libraries missing | Install BLAS/LAPACK packages and rebuild `PfaPack`. |
| `PfaPack.jl/` or `SFMT.jl/` is empty | Submodules were not initialized | Run `git submodule update --init --recursive`. |
| Examples fail to resolve workspace packages | Wrong Julia project | Use `--project=@.` from the repository root. |
