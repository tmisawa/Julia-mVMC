# Examples

Each example runs `MVMCOptimizers.run_para_opt_from_namelist` and prints the final energy per site. Inputs are bundled under `examples/inputs/<model>/`.

| File | Model | Mode |
|------|-------|------|
| `heisenberg_chain_real.jl` | 16-site Heisenberg chain | real |
| `heisenberg_chain_cmp.jl` | 16-site Heisenberg chain | complex |
| `heisenberg_chain_fsz.jl` | 16-site Heisenberg chain | fsz (generalized orbital) |
| `hubbard_chain.jl` | Hubbard chain | real |

## Run

From the repo root with the workspace project activated:

```bash
julia --project=@. examples/heisenberg_chain_real.jl
```

Default is 50 SR steps. Override with the `JULIA_MVMC_EXAMPLE_STEPS` env var:

```bash
JULIA_MVMC_EXAMPLE_STEPS=10 julia --project=@. examples/heisenberg_chain_real.jl
```

Each script takes seconds (5–10 steps) to a few minutes (50 steps) on a single core.
