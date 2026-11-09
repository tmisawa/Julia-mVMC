# L=128 CPU performance, 2026-11-10

Target workload from the 2026-11-05 Rust benchmark: half-filled Hubbard
chain, L=128, real, U=4, t=1, seed 79, 8 QP planes, 320 samples/walker.
Optimization is 20 SR steps; PhysCal is one fixed-parameter window without
Lanczos. Julia 1.13.3, Xeon E5-2699 v3, physical cores 0..N-1, one BLAS
thread, `JULIA_MVMC_INNER_THREADS=1`, `JULIA_MVMC_PFAPACK_THREADS=1`.

## Reproduction

```sh
julia +1.13 --project=. benchmark/cpu_l128/measure.jl \
  /path/to/mvmc-rs/benchmark/cpu_scaling_l128/inputs \
  /tmp/julia-l128 20 3 both
```

Runs use the same input files and the same `measure.jl` entry point as the
Rust benchmark harness (`scripts/benchmark/julia.jl`). One warmup plus three
measured repetitions; the table is the median of the three (1-core/4-core
rows use one measured repetition because the run time is dominated by the
serial kernel and the spread is small).

## Results (seconds)

| Workload | Cores | Before (2026-11-05) | After | Rust | After/Rust |
|---|---:|---:|---:|---:|---:|
| Optimize | 1 | 78.035 | 66.27 | 69.036 | 0.96 |
| Optimize | 4 | 40.257 | 24.37 | 24.377 | 1.00 |
| Optimize | 8 | 34.204 | 16.88 | 16.721 | 1.01 |
| Optimize | 16 | 29.170 | 16.12 | 14.779 | 1.09 |
| PhysCal | 1 | 4.377 | 3.752 | 3.728 | 1.01 |
| PhysCal | 4 | 2.841 | 1.473 | 1.382 | 1.07 |
| PhysCal | 8 | 2.304 | 1.006 | 0.946 | 1.06 |
| PhysCal | 16 | 2.006 | 0.981 | 0.901 | 1.09 |

All eight target cells are now within **0.96–1.09x** of the Rust port; the
8/16-core PhysCal cells are 1.06x and 1.09x. The remaining difference is the
residual per-thread factor/inverse throughput and the ~0.01 s parser floor.
The last two rows carry roughly ±4% run-to-run spread on this shared host.

## Main changes

- Ported the Rust `utu2inv` real path's **64-column panelized upper-triangle
  `Mᵀ A`** (`_ordinary_trmmt_upper_real!`) plus skew-symmetric lower-triangle
  restoration. The product touches half the entries and stays cache resident.
- The LTL rank-2 update processes two adjacent columns per pass, sharing every
  `A[i,k]`/`A[i,kk]` load (the Rust AVX2 two-column pairing), and the
  skew-tridiagonal solve runs in 16-column blocks so its three strided row
  slices stay in L1/L2. Profiling at 16 threads showed the solve was the
  worst-scaling phase; blocking made it 1.6x faster serially and 1.3x at
  8 workers.
- `CalculateMAll` saved-sample batches run on the persistent, CPU-pinned hop
  team instead of a fresh `@threads` region per batch; this removed task
  migration and improved the 16-core stage.
- Accepted-hop rank-1 update uses 16 register-held accumulators (same
  per-element addition order, no FMA).
- Real Slater derivative accumulates in `Float64`; measurement batches swap
  inversion vectors instead of copying; batch capacity is memory-bounded.
- CPU main-calculation loops use dynamic default-pool scheduling.
- PhysCal skips all SR-accumulator allocation/reset/merge work (the 60 MB
  `sr_opt` arrays are not used in measurement mode) and avoids clearing them
  in `clear_phys_quantity!`.
- PhysCal one-body Green functions run in parallel with per-thread reusable
  scratch.
- Specialized byte-level integer scanner for the 16k-line `orbitalidx.def` /
  `jastrowidx.def` tables: Jastrow parse **34 ms → 2.1 ms**, Orbital
  **11 ms → 2.0 ms**, and the separate Jastrow contract pass is folded into
  the same scan. A backward OptFlag-only orbital scan replaces reparsing the
  site-pair table.
- Optimization-mode MainCal keeps one sleeping worker team for the whole
  saved-sample pass instead of creating/pinning a team for every batch
  (`CalculateMAll` 16-core 0.26 s → 0.21 s). PhysCal keeps per-batch teams
  because its Green kernels share the default-pool threads.
- The final row permutation walks contiguous columns (Rust
  `permute_rows_forward` order) instead of striding through rows.

## Validation

```sh
julia +1.13 --project=MVMCOptimizers.jl -e 'using Pkg; Pkg.test()'
# 54681/54681 pass

julia +1.13 --project=. test/integration/runtests.jl
# 120/120 pass vs C 1.4.0 reference (including hubbard_chain_real)

julia +1.13 --project=MVMCExpertModeParsers.jl -e 'using Pkg; Pkg.test()'
# pass

julia +1.13 --threads=4,0 --project=. -e \
  'using MVMCExpertModeParsers; include("MVMCOptimizers.jl/test_unit/test_unit_real_hop_threading.jl"); include("MVMCOptimizers.jl/test_unit/test_unit_real_maincal_batch.jl"); include("MVMCOptimizers.jl/test_unit/test_unit_green_threading.jl")'
# pass (including RNG tape and worker-error teardown)
```

The C-reference integration tests exercise the real-Slater, batch, and
PhysCal paths. The PhysCal one-body Green output changes only at the
2.8e-17 roundoff level relative to the serial kernel; direct two-body Green
and energy output are bit-identical.
