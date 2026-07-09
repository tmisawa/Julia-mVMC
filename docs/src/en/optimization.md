# Parameter Optimization

[日本語](../ja/optimization.md)

`VMCParaOpt` (`NVMCCalMode = 0`) is the most mature Julia-mVMC path. Use the
high-level wrapper `run_para_opt_from_namelist`.

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
| `namelist_path` | path to `namelist.def` |
| `nsteps` | number of SR optimization steps |
| `mode` | sanity label: `:real`, `:cmp`, or `:fsz` |
| `nsmp` | final-sample count for ctest-style averaging |
| `output_dir` | directory for `zvo_*` and `zqp_*` files |
| `seed` | optional SFMT19937 seed override |
| `initial_def` | `:auto`, `:none`, or an explicit `initial.def` path |

## Verified models

The strict public integration gate compares the first 10 SR steps against
committed C-mVMC references.

| Fixture | Mode | Coverage |
|---------|------|----------|
| `heisenberg_chain_real` | real | spin-only Hamiltonian, real Slater |
| `heisenberg_chain_cmp` | complex | complex orbital path |
| `heisenberg_chain_fsz` | FSZ/general orbital | individually tracked spin path |
| `hubbard_chain_real` | real | charge fluctuations with Gutzwiller/Jastrow |

A broader C ctest-equivalent gate covers supported standard fixtures and
`GeneralRBM_cmp` with C's `ref_mean.dat` / `ref_std.dat` acceptance rule.

## MPI scope

v0.5.0 supports sample-parallel MPI for:

- direct SR (`NSRCG = 0`) with `NSplitSize >= 1` and `NQPFull = 1`;
- direct SR with sz-conserved standard projection `NQPFull > 1` when
  `NQPOptTrans = 1`;
- standard SR-CG (`NSRCG = 1`) with `NSplitSize = 1`.

Smoke gates cover rank0-only output, comm0 reductions, direct-SR
`NSplitSize/NStore` self-consistency, standard-projection self-consistency, and
one-step `NSRCG = 1` C-reference behavior under `mpiexec`.

## Rejected combinations

The following combinations are intentionally rejected in this release:

- `NSplitSize > 1` with SR-CG;
- `NSRCG >= 2`;
- `useDiagScale != 0`;
- `RescaleSmat != 0`;
- `NSplitSize > 1` with OptTrans-derived QP sectors;
- FSZ standard projection with `NQPFull > 1`.

These rejections avoid silently mixing unsupported stochastic or MPI semantics
with C-compatible paths.

## Timer output

Set `MVMC_C_TIMER=1` to write C-compatible section timing data:

```bash
MVMC_C_TIMER=1 julia --project=@. examples/heisenberg_chain_real.jl
```

The timer writes `zvo_CalcTimer.dat` in the output directory.
