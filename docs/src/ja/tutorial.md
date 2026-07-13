# チュートリアル

[English](../en/tutorial.md)

この tutorial では、同梱された expert-mode 入力を使って小さな parameter
optimization を実行します。[インストール](installation.md)が完了している前提です。

## Heisenberg-chain example を実行する

repository root で次を実行します。

```bash
JULIA_MVMC_EXAMPLE_STEPS=5 julia --project=@. examples/heisenberg_chain_real.jl
```

この script は `examples/inputs/heisenberg_chain_real/namelist.def` を読み、
`VMCParaOpt` を実行して final energy per site を表示します。

## Example の中身

example は次の wrapper を呼び出します。

```julia
using MVMCOptimizers

result = MVMCOptimizers.run_para_opt_from_namelist(
    "examples/inputs/heisenberg_chain_real/namelist.def";
    nsteps = 5,
    nsmp = 5,
    mode = :real,
)
```

high-level wrapper は C-mVMC の `vmcmain.c` と同じ phase order を辿ります。

1. `namelist.def` から expert-mode file を parse する。
2. variational parameter を初期化する。
3. 必要なら `initial.def` を読む。
4. `In*.def` overlay を適用する。
5. derived parameter と quantum-projection weight を同期する。
6. SR optimization loop を実行する。

## 同梱 mode を試す

```bash
JULIA_MVMC_EXAMPLE_STEPS=5 julia --project=@. examples/heisenberg_chain_cmp.jl
JULIA_MVMC_EXAMPLE_STEPS=5 julia --project=@. examples/heisenberg_chain_fsz.jl
JULIA_MVMC_EXAMPLE_STEPS=5 julia --project=@. examples/hubbard_chain.jl
```

| Example | Mode | Purpose |
|---------|------|---------|
| `heisenberg_chain_real.jl` | real | real Slater spin-chain path |
| `heisenberg_chain_cmp.jl` | complex | complex orbital path |
| `heisenberg_chain_fsz.jl` | FSZ/general orbital | spin を個別追跡する path |
| `hubbard_chain.jl` | real | charge fluctuation を含む Hubbard path |

## Output directory

example script は wrapper の default output directory を使います。永続的な directory
へ出力したい場合は wrapper を直接呼びます。

```julia
result = MVMCOptimizers.run_para_opt_from_namelist(
    "examples/inputs/heisenberg_chain_real/namelist.def";
    nsteps = 5,
    nsmp = 5,
    mode = :real,
    output_dir = "run_heisenberg_real",
)
```

主な file は[出力ファイル](output_files.md)で説明します。
