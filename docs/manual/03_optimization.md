# 3. Parameter optimization (VMCParaOpt)

The Julia port of C-mVMC's `VMCParaOpt` (`NVMCCalMode = 0`) is the most
mature component. Use it via the high-level wrapper
`run_para_opt_from_namelist`, which mirrors `vmcmain.c:264–281` exactly.

## API

```julia
using MVMCOptimizers

result = run_para_opt_from_namelist(
    namelist_path::AbstractString;
    nsteps::Integer,
    mode::Symbol,                         # :real | :cmp | :fsz
    output_dir::AbstractString = tempname(),
    seed::Union{Integer,Nothing} = nothing,
    initial_def::Union{AbstractString,Symbol,Nothing} = :auto,
)
```

| Argument | Behaviour |
|----------|-----------|
| `namelist_path` | Path to `namelist.def`. Other `.def` files are resolved relative to it. |
| `nsteps` | Number of SR steps. Overrides `NSROptItrStep` from `modpara.def`. |
| `mode` | Sanity label (one of `:real`, `:cmp`, `:fsz`); the actual run mode is determined by `complex_flag` and orbital files in the parsed input. |
| `output_dir` | Where `zvo_out.dat`, `zqp_opt.dat`, etc. are written. Defaults to a fresh tempdir. |
| `seed` | SFMT19937 seed. `nothing` → use `RndSeed` from `modpara.def` (fallback `11272` when ≤ 0). |
| `initial_def` | `:auto` (default) loads `inputs/initial.def` if present and aborts on a present-but-broken file; `:none` / `nothing` skips entirely; an explicit path errors if loading fails. |

Returns a `NamedTuple` with `status`, `output_dir`, `zvo_first_n` (first
`nsteps` raw lines of `zvo_out.dat`), and `final_energy_per_site`.

## Verified models (bit-level vs C reference)

The public CI continuously verifies four models — one per major
category — against committed C-mVMC reference output, comparing the
first 10 SR steps within BLAS summation-order noise (≤ 1e-10 for the
linear accumulators, ≤ 1e-9 for the squared/derived columns; see
[`test/integration/runtests.jl`](../../test/integration/runtests.jl)
and [`test/integration/reference/README.md`](../../test/integration/reference/README.md)
for provenance):

| Fixture | Mode | Coverage |
|---------|------|----------|
| `heisenberg_chain_real` | real | spin-only Hamiltonian, real Slater |
| `heisenberg_chain_cmp` | complex | spin-only, complex orbital path |
| `heisenberg_chain_fsz` | fsz (generalised orbital) | individually-tracked spin per electron |
| `hubbard_chain_real` | real | charge fluctuations + Gutzwiller/Jastrow |

Internal smoke runs (not part of the public CI, not bundled as
fixtures) cover additional models in the same categories — Kondo
chains, HubbardTetragonal momentum-projection, GeneralRBM_cmp, and
the Lanczos chains at step 0 — and have historically reproduced C
output within the same tolerances. Treat them as "expected to work
but not continuously verified in this repo": if you hit a regression
on an unbundled model, please open an issue with a reduced fixture so
it can be added to the public CI.

`SpinChainLanczos` and `HubbardChainLanczos` reproduce C only at
step 0 (the SR section is bypassed there); the full Lanczos section
is **not** ported — see [`05_compatibility.md`](05_compatibility.md).

## Output files

Written under `output_dir`:

| File | Contents |
|------|----------|
| `zvo_out.dat` | Per-step row: `real(<H>) imag(<H>) <H²> variance <Sz> <Sz²>`. Step 0 truncates the file; subsequent steps append. |
| `zvo_var.dat` | Per-step row: parameter snapshot (Gutzwiller, Jastrow, Slater real/imag). |
| `zqp_opt.dat` | Final optimised parameters (one column per parameter). |
| `zqp_gutzwiller_opt.dat`, `zqp_jastrow_opt.dat`, `zqp_orbital*_opt.dat` | Per-block optimised parameter dumps (mirrors C). |

## Reading variance with care

Column 4 (`variance = (<H²> - <H>²) / <H>²`) is computed by a cancelling
subtraction. With the integration fixtures this picks up ~1e-10 of
amplified BLAS noise on top of the ~1e-13 noise in `<H>` itself. Treat
`variance` as a derived diagnostic, not a primary observable.

## End-to-end example

See [`examples/heisenberg_chain_real.jl`](../../examples/heisenberg_chain_real.jl)
and the other three example scripts. The 50-step default takes seconds
on a single core; CI runs them at `JULIA_MVMC_EXAMPLE_STEPS=5`.
