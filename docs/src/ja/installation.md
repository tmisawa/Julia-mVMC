# インストール

[English](../en/installation.md)

## 検証済み環境

| Item | Requirement |
|------|-------------|
| Julia | 1.11 以降。CI では Julia 1.11 / 1.12 を検証。 |
| OS | Linux と macOS。Windows はこの release では未検証。 |
| Native tools | `gfortran`, `g++`, `make`, BLAS, LAPACK。 |

`PfaPack.jl` と `SFMT.jl` は git submodule で、初回 package build 時に native
library を build します。

## Native prerequisites

Ubuntu / Debian:

```bash
sudo apt-get update
sudo apt-get install -y gfortran g++ make libblas-dev liblapack-dev
```

Homebrew を使う macOS:

```bash
brew install gcc
```

runtime の BLAS/LAPACK は Apple system library で十分です。Apple は Fortran
compiler を同梱しないため、Homebrew `gfortran` が必要です。

## Clone と instantiate

対応している install path は、submodule を含む workspace checkout です。

```bash
git clone --recurse-submodules https://github.com/tmisawa/Julia-mVMC
cd Julia-mVMC
julia --project=@. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
```

submodule なしで clone した場合は、先に次を実行します。

```bash
git submodule update --init --recursive
```

Julia は repository root から `--project=@.` で起動してください。root project は
relative `[sources]` entry によって `MVMCOptimizers.jl`,
`MVMCExpertModeParsers.jl`, `PfaPack.jl`, `SFMT.jl` をまとめます。

## なぜ clone-based installation か

v0.5.0 では `Pkg.add(url=..., subdir="MVMCOptimizers.jl")` は未対応です。
subpackage project が sibling package と submodule を relative path で参照しており、
これは full workspace checkout でのみ有効だからです。

## Smoke test

```bash
JULIA_MVMC_EXAMPLE_STEPS=5 julia --project=@. examples/heisenberg_chain_real.jl
```

最後に `Final energy / site = ...` が表示されれば正常です。

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `gfortran: not found` | Fortran compiler がない | `gfortran` を install する。 |
| macOS linker が `System` を見つけられない | Homebrew `gfortran` が SDK を見つけられない | `SDKROOT="$(xcrun --sdk macosx --show-sdk-path)"` を設定して rebuild。 |
| `dgemm_` など BLAS symbol が missing | BLAS/LAPACK development library がない | BLAS/LAPACK package を install し、`PfaPack` を rebuild。 |
| `PfaPack.jl/` または `SFMT.jl/` が空 | submodule が未初期化 | `git submodule update --init --recursive` を実行。 |
| example が workspace package を解決できない | Julia project が違う | repository root から `--project=@.` を使う。 |
