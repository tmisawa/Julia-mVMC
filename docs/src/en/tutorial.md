# Tutorial

[日本語](../ja/tutorial.md)

This tutorial runs a small parameter optimization from bundled expert-mode input
files. It assumes that [installation](installation.md) has completed.

## Run the Heisenberg-chain example

From the repository root:

```bash
JULIA_MVMC_EXAMPLE_STEPS=5 julia --project=@. examples/heisenberg_chain_real.jl
```

The script loads `examples/inputs/heisenberg_chain_real/namelist.def`, runs
`VMCParaOpt`, and prints the final energy per site.

## What the example does

The example calls:

```julia
using MVMCOptimizers

result = MVMCOptimizers.run_para_opt_from_namelist(
    "examples/inputs/heisenberg_chain_real/namelist.def";
    nsteps = 5,
    nsmp = 5,
    mode = :real,
)
```

The high-level wrapper mirrors C-mVMC's `vmcmain.c` phase ordering:

1. parse expert-mode files through `namelist.def`;
2. initialize variational parameters;
3. optionally read `initial.def`;
4. apply `In*.def` overlays;
5. synchronize derived parameters and quantum-projection weights;
6. run the SR optimization loop.

## Try the bundled modes

```bash
JULIA_MVMC_EXAMPLE_STEPS=5 julia --project=@. examples/heisenberg_chain_cmp.jl
JULIA_MVMC_EXAMPLE_STEPS=5 julia --project=@. examples/heisenberg_chain_fsz.jl
JULIA_MVMC_EXAMPLE_STEPS=5 julia --project=@. examples/hubbard_chain.jl
```

| Example | Mode | Purpose |
|---------|------|---------|
| `heisenberg_chain_real.jl` | real | real Slater spin-chain path |
| `heisenberg_chain_cmp.jl` | complex | complex orbital path |
| `heisenberg_chain_fsz.jl` | FSZ/general orbital | individually tracked spin path |
| `hubbard_chain.jl` | real | charge-fluctuating Hubbard path |

## Output location

The example scripts use the wrapper default output directory. For a persistent
directory, call the wrapper directly:

```julia
result = MVMCOptimizers.run_para_opt_from_namelist(
    "examples/inputs/heisenberg_chain_real/namelist.def";
    nsteps = 5,
    nsmp = 5,
    mode = :real,
    output_dir = "run_heisenberg_real",
)
```

Important files are described in [Output files](output_files.md).
