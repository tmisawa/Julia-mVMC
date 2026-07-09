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
    nsmp::Union{Integer,Nothing} = nothing,
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
| `nsmp` | Number of final optimisation samples used for C ctest-style averaging. `nothing` preserves `NSROptItrSmp`; an integer overrides it and must satisfy `nsteps >= nsmp`. |
| `output_dir` | Where `zvo_out.dat`, `zqp_opt.dat`, etc. are written. Defaults to a fresh tempdir. |
| `seed` | SFMT19937 seed. `nothing` → resolve `RndSeed` from `modpara.def` with the C-parity rule (v0.4): missing → `11272`, `0` → `0`, negative → time-derived seed, positive → the value; `+ group1` under MPI. |
| `initial_def` | `:auto` (default) loads `inputs/initial.def` if present and aborts on a present-but-broken file; `:none` / `nothing` skips entirely; an explicit path errors if loading fails. |

## MPI limitations

v0.5.0 supports `VMCParaOpt` sample-parallel MPI with the direct SR
solver (`NSRCG = 0`) for `NSplitSize >= 1` when each sample uses either
`NQPFull = 1` or sz-conserved standard-projection `NQPFull > 1` with
`NQPOptTrans = 1` (`NSPGaussLeg > 1` and/or `NMPTrans > 1`). Non-FSZ, FSZ
`NQPFull = 1`, and sz-conserved standard-projection paths are gated by
self-consistency comparisons:
`(NSplitSize=1, NStore=0/1)` under `mpiexec -n 2` and
`(NSplitSize=2, NStore=0/1)` under `mpiexec -n 4` must agree for the relevant
parameter-optimization output files.

The standard SR-CG solver (`NSRCG = 1`) remains limited to `NSplitSize = 1`.
That path broadcasts the rank0 search vector and allreduces the sampled
matrix-vector product inside `operate_by_s`, matching C's collective structure
for `NSplitSize = 1`.

Additional C CG modes are not ported yet: `NSRCG >= 2`, `useDiagScale != 0`,
and `RescaleSmat != 0` are rejected with clear errors. `NSplitSize > 1` with
SR-CG is also rejected because the comm1-aware CG path is not implemented.
OptTrans-derived QP sectors (`NQPOptTrans > 1` or active `OptTrans`) remain
unsupported with `NSplitSize > 1`; use `NSplitSize = 1` for those inputs.
FSZ standard-projection `NQPFull > 1` (`NSPGaussLeg > 1` or
`abs(NMPTrans) > 1`) is not supported.

The CI MPI smoke gate runs under `mpiexec` and checks rank0-only output,
comm0 reductions, `VMCPhysCal` Green reductions, and failure-mode regressions:
SR-CG `operate_by_s` collectives and a one-step `NSRCG = 1` C-reference
comparison are exercised under `mpiexec -n 2`, while `NSRCG >= 2` and
`NSplitSize > 1` with SR-CG or active OptTrans are rejected before MPI context
construction. This is a correctness smoke gate, not a performance or
site-compatibility benchmark.

## NSRCG=1 coverage

`NSRCG = 1` runs are supported and covered by one-step C reference fixtures
(`heisenberg_chain_real_nsrcg`) for serial execution and `mpiexec -n 2`. The
first `zvo_out.dat` row uses the same tight energy/Sz tolerances as the
standard fixtures. The post-CG parameter update is a tolerance gate
(`NSRCG_PARAM_TOL = 1e-2` in
[`test/integration/runtests.jl`](../../test/integration/runtests.jl) and the
matching MPI smoke), not a bit-parity gate, because truncated SR-CG amplifies
FMA and reduction-order differences between C and Julia.

Returns a `NamedTuple` with `status`, `output_dir`, `zvo_first_n` (first
`nsteps` raw lines of `zvo_out.dat`), `ctest_values` (the first two
final-sample averages used by the C ctest-equivalent runner),
`final_energy_per_site`, `effective_nsteps`, and `effective_nsmp`.

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

A broader C ctest-equivalent CI gate covers the 12 supported standard C
fixtures plus `GeneralRBM_cmp` using C's `ref_mean.dat` / `ref_std.dat`
criterion rather than a first-10-step comparison. See
[`05_compatibility.md`](05_compatibility.md) for the runner and model
selection details.

Internal smoke runs (not part of the public CI, not bundled as
fixtures) cover additional models in the same categories — Kondo
chains, HubbardTetragonal momentum-projection, and the Lanczos chains at step
0 — and have historically reproduced C output within the same tolerances.
Treat them as "expected to work but not continuously verified in this repo":
if you hit a regression on an unbundled model, please open an issue with a
reduced fixture so it can be added to the public CI.

`SpinChainLanczos` and `HubbardChainLanczos` reproduce C only at
step 0 in the ParaOpt path (the SR section is bypassed there).
`NLanczosMode > 0` is still rejected for ParaOpt; the R1 Lanczos
support is limited to PhysCal energy/QQQQ output — see
[`05_compatibility.md`](05_compatibility.md).

## Output files

Written under `output_dir`:

| File | Contents |
|------|----------|
| `zvo_out.dat` | Per-step row: `real(<H>) imag(<H>) <H²> variance <Sz> <Sz²>`. Step 0 truncates the file; subsequent steps append. |
| `zvo_var.dat` | Per-step row: parameter snapshot (Gutzwiller, Jastrow, Slater real/imag). |
| `zqp_opt.dat` | Final optimised parameters (one column per parameter). |
| `zqp_gutzwiller_opt.dat`, `zqp_jastrow_opt.dat`, `zqp_orbital*_opt.dat` | Per-block optimised parameter dumps (mirrors C). |

## C-compatible timing output

The C-compatible section timer is disabled by default. Enable it explicitly
with the `MVMC_C_TIMER` environment variable:

```bash
MVMC_C_TIMER=1 julia --project=@. examples/heisenberg_chain_real.jl
```

When enabled, `run_para_opt_from_namelist` writes `zvo_CalcTimer.dat` under
`output_dir`. The file uses the same timer id / label layout as C-mVMC's
`zvo_CalcTimer.dat`, so the same comparison scripts can be used for
section-level bottleneck analysis.

Unset `MVMC_C_TIMER`, or set it to `0`, to keep the timer disabled. The
legacy `MVMC_TIMER` environment variable is accepted as a deprecated alias
for the same C-compatible timer:

```bash
MVMC_C_TIMER=0 MVMC_TIMER=0 julia --project=@. examples/heisenberg_chain_real.jl
```

The implementation treats an unset variable or the exact string `0` as
disabled; **any other value enables timing — including intuitively "off"
looking values such as `false`, `off`, `no`, or an empty string.** Only
an unset variable or the literal `0` keeps the timer off; use
`MVMC_C_TIMER=1` to enable it. Timer-enabled runs include timing overhead
and should be used for bottleneck breakdowns, not as the primary
elapsed-time benchmark.

## Reading variance with care

Column 4 (`variance = (<H²> - <H>²) / <H>²`) is computed by a cancelling
subtraction. With the integration fixtures this picks up ~1e-10 of
amplified BLAS noise on top of the ~1e-13 noise in `<H>` itself. Treat
`variance` as a derived diagnostic, not a primary observable.

## End-to-end example

See [`examples/heisenberg_chain_real.jl`](../../examples/heisenberg_chain_real.jl)
and the other three example scripts. The 50-step default takes seconds
on a single core; CI runs them at `JULIA_MVMC_EXAMPLE_STEPS=5`.
