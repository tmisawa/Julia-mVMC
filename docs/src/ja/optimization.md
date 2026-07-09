# パラメータ最適化

[English](../en/optimization.md)

`VMCParaOpt` (`NVMCCalMode = 0`) は Julia-mVMC で最も成熟している path です。
通常は high-level wrapper `run_para_opt_from_namelist` を使います。

## API

```julia
using MVMCOptimizers

result = run_para_opt_from_namelist(
    namelist_path;
    nsteps,
    mode,
    nsmp = nothing,
    output_dir = tempname(),
    seed = nothing,
    initial_def = :auto,
)
```

| Argument | Meaning |
|----------|---------|
| `namelist_path` | `namelist.def` への path |
| `nsteps` | SR optimization step 数 |
| `mode` | sanity label: `:real`, `:cmp`, `:fsz` |
| `nsmp` | ctest-style averaging 用 final-sample count |
| `output_dir` | `zvo_*` / `zqp_*` file の出力先 |
| `seed` | SFMT19937 seed override |
| `initial_def` | `:auto`, `:none`, または明示的な `initial.def` path |

## 検証済み model

strict public integration gate では、最初の 10 SR step を committed C-mVMC
reference と比較します。

| Fixture | Mode | Coverage |
|---------|------|----------|
| `heisenberg_chain_real` | real | spin-only Hamiltonian, real Slater |
| `heisenberg_chain_cmp` | complex | complex orbital path |
| `heisenberg_chain_fsz` | FSZ/general orbital | spin を個別追跡する path |
| `hubbard_chain_real` | real | Gutzwiller/Jastrow を含む charge fluctuation |

より広い C ctest-equivalent gate では、対応する standard fixture と
`GeneralRBM_cmp` を C の `ref_mean.dat` / `ref_std.dat` acceptance rule で確認します。

## MPI scope

v0.5.0 では次の sample-parallel MPI を support します。

- `NQPFull = 1` の direct SR (`NSRCG = 0`) with `NSplitSize >= 1`;
- `NQPOptTrans = 1` の sz-conserved standard projection `NQPFull > 1`;
- `NSplitSize = 1` の standard SR-CG (`NSRCG = 1`)。

Smoke gate では rank0-only output、comm0 reduction、direct-SR
`NSplitSize/NStore` self-consistency、standard-projection self-consistency、
`mpiexec` 下の one-step `NSRCG = 1` C-reference behavior を確認します。

## reject される組合せ

この release では以下を意図的に reject します。

- `NSplitSize > 1` with SR-CG;
- `NSRCG >= 2`;
- `useDiagScale != 0`;
- `RescaleSmat != 0`;
- OptTrans-derived QP sector with `NSplitSize > 1`;
- FSZ standard projection with `NQPFull > 1`.

未対応の stochastic / MPI semantics が C-compatible path に混ざることを避けるためです。

## Timer output

C-compatible section timing data を書くには `MVMC_C_TIMER=1` を設定します。

```bash
MVMC_C_TIMER=1 julia --project=@. examples/heisenberg_chain_real.jl
```

timer は output directory に `zvo_CalcTimer.dat` を書きます。
