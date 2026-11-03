# CPU L128 threading (#92)

The CPU sampler formerly updated accepted-hop inverse matrices serially even
with inner threading enabled. Launching static Julia tasks for every hop is
slower than the serial update on this input. A sampling-call-scoped team now
reuses at most `min(QPs, Julia threads) - 1` workers while the caller owns all
RNG draws, proposals, decisions and ordered QP sums. Workers use private scratch
and disjoint matrix planes. The team is joined before returning, including on
errors; periodic full factorizations use the same team.

MainCal factors bounded groups of saved samples before consuming their results
in the original sample order. This fills a 16-thread pool with eight QPs without
copying one optimizer state per worker or regrouping sample reductions. Each
matrix retains the existing C-order factor/inverse implementation. Direct Green
operators use cyclic worker assignment with independent output slots; the sample
accumulation and factored-Green dependency barrier remain ordered. Partial QP
pools also parallelize the existing Slater derivative contraction.

`JULIA_MVMC_INNER_THREADS=1` remains the explicit inner-threading opt-in;
`--threads` alone does not change that existing policy. The zero/unset policy,
small-workload paths, nested callers and complex/FSZ sampling are retained.
There is no dependency or lock-file change.

## Reproduce

Use **`julia +1.13`**, the committed `Manifest-v1.13.toml`, normal native
PfaPack/SFMT build prerequisites, and the common opt/phys inputs from
[mvmc-rs](https://github.com/AtelierArith/mvmc-rs/tree/ab10344c499872923ccbb96695d3a99c64d61e7a/benchmark/cpu_scaling_l128/inputs).
The version resolved by this channel is recorded, not assumed. Source identifiers,
input/manifest/driver hashes, actual loaded libraries, affinity and all repetitions
are retained by the driver. Output directories must be fresh.

```bash
julia +1.13 --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.build("PfaPack"); Pkg.build("SFMT")'
bench_inputs=/path/to/mvmc-rs/benchmark/cpu_scaling_l128/inputs
for n in 1 4 8 16; do
  JULIA_MVMC_INNER_THREADS=1 JULIA_MVMC_PFAPACK_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    taskset -c "0-$((n-1))" julia +1.13 --startup-file=no --threads="$n" --project=. \
    benchmark/cpu_l128/measure.jl "$bench_inputs" "/tmp/julia-l128-after-t$n" 20 3 both
done
```

Check physical-core numbering before using this host's affinity. Run the same
script against a separate clean checkout of before revision
`611031147f1bad9704c55084f6817feeba481b13` using that checkout's `--project`;
keep its native builds and inputs pinned. BLAS is one thread. Each public API
includes parsing and file output, with one warmup and three measurements;
process startup, initial JIT, explicit pre-call GC and benchmark serialization
are excluded. Opt is 20 SR steps, 320 samples; PhysCal is one fixed-parameter
320-sample window without Lanczos. CPU/GPU sampler trajectories differ.

```bash
julia +1.13 --project=. benchmark/cpu_l128/verify.jl \
  /tmp/julia-l128-before-t8 /tmp/julia-l128-after-t8
OPENBLAS_NUM_THREADS=1 julia +1.13 --threads=4,0 --project=. -e '
  include("MVMCOptimizers.jl/test_unit/test_unit_real_hop_threading.jl")
  include("MVMCOptimizers.jl/test_unit/test_unit_real_maincal_batch.jl")
  include("MVMCOptimizers.jl/test_unit/test_unit_green_threading.jl")'
```

The output verifier compares every retained warmup/measured output family with
an explicit per-value `1e-11 + 1e-11*max(abs(x),abs(y))` bound, separately checking
Green integer columns exactly. The supplied-word-tape regression checks exact
draw count, the configuration/counter history at every draw, all saved
configuration arrays, and that no worker consumes RNG. Independent inverse
residuals, unchanged C fixtures, mixed kernel statuses, batch tails and repeated
reuse are separate checks. Performance tables are recorded after final validation.
