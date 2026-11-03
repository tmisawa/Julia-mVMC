# L128 CPU Opt and PhysCal (#92)

Same-host full-workload medians improve by **2.18x for 16-core Opt** and
**2.01x for 16-core PhysCal**. All retained numerical output pairs and RNG
primitive streams agree between revisions. C/Rust remain faster at 4/8/16 cores.

## Conditions and sources

Linux x86_64, Xeon E5-2699 v3 (18 physical cores, 36 logical CPUs), Juliaup
`+1.13` resolving to **Julia 1.13.3**, ILP64 OpenBLAS **0.3.30 / Haswell**, one
BLAS thread, GNU Fortran 13.3.0. The native dependency hashes and actual loaded
library paths are retained. The first N physical cores (CPU 0..N-1) are used;
SMT siblings are excluded. No CPU frequency control or exclusive host ownership
is claimed. GPU work ran concurrently with its host on CPU17; development
builds/correctness checks used CPU16/17 and their SMT siblings. CPU timing cases
ran sequentially. These are three-observation point medians, not confidence bounds.

Before: `611031147f1bad9704c55084f6817feeba481b13`.
After: `032de0b34f29d646b457a0c2b8175082e799d8ca` (runtime `8d29a60`, followed by
an admission-test-only commit). Later report files do not change runtime code.
`Manifest-v1.13.toml` remains byte-identical. PfaPack is `9161cda4`, SFMT `15265530`.

Half-filled real Hubbard chain L=128, U=4, t=1, seed 79, eight QP planes,
320 samples, NVMCWarmUp=1 and NVMCInterval=1. Opt runs 20 SR steps; PhysCal runs
one fixed-parameter window without Lanczos. Inputs are copied unchanged from
[mvmc-rs ab10344c](https://github.com/AtelierArith/mvmc-rs/tree/ab10344c499872923ccbb96695d3a99c64d61e7a/benchmark/cpu_scaling_l128/inputs).
Fixed parameter SHA-256:
`9130219b4af78438bee2f30171833e7a8f86a01ff11e23276af439652230fd90`.
All eight timing cohorts have identical input, manifest and driver hashes.

Each case uses one warmup and three measured calls: 64 complete API calls in
total. Timings include public CPU API input parsing, initialization, sampling,
observables/SR, file output and API readback. Process startup, initial JIT,
explicit pre-call GC and benchmark serialization are excluded. Natural GC during
the API is included. Inner threading is explicitly enabled with
`JULIA_MVMC_INNER_THREADS=1`; the existing default/zero policy remains unchanged.

## Results

Seconds; each cell is the median of all three measured repetitions.
`summary.toml` retains every value and range; raw CSV files include warmups.

| Physical cores | Opt before | Opt after | Before/after | PhysCal before | PhysCal after | Before/after |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 77.476 | 77.078 | 1.005x | 4.640 | 4.384 | 1.058x |
| 4 | 73.914 | 39.793 | 1.857x | 4.460 | 2.618 | 1.704x |
| 8 | 61.539 | 33.702 | 1.826x | 3.776 | 2.262 | 1.669x |
| 16 | 62.233 | 28.504 | 2.183x | 3.821 | 1.900 | 2.011x |

Single-core differences are observations across sequential cohorts without
frequency control; the changes target parallel work.

For context, the latest separately recorded
[C/Rust #629 comparison](https://github.com/AtelierArith/mvmc-rs/blob/ab10344c499872923ccbb96695d3a99c64d61e7a/benchmark/cpu_scaling_l128/results/20261103-issue629/report.md)
uses the same definitions and parameters:

| Workload | Cores | C | Rust | This Julia |
|---|---:|---:|---:|---:|
| Opt | 8 | 19.565 | 17.660 | 33.702 |
| Opt | 16 | 21.817 | 16.417 | 28.504 |
| PhysCal | 8 | 1.261 | 1.026 | 2.262 |
| PhysCal | 16 | 1.291 | 0.961 | 1.900 |

C/Rust use process wall including startup/exit (C includes `mpirun`), excluding
Docker startup/exec launch. Julia uses the warm API scope above. The C/Rust
session and BLAS builds differ; this contextual table is not a newly paired
three-language comparison. The remaining 16-core wall gap is about 1.31x/1.74x
versus C/Rust for Opt and 1.47x/1.98x for PhysCal under those respective scopes.

## Threading changes and remaining work

- Accepted-hop inverse updates were serial. Sampling now reuses at most
  `min(QPs, threads)-1` workers with private scratch. The caller retains RNG,
  proposals, acceptance and QP reductions. Periodic full factorizations share
  this team; its lifetime ends before returning to the owner loop.
- MainCal factors bounded groups of saved samples before consuming them in the
  original order. Two samples provide 16 independent planes for this eight-QP
  input. Only matrix buffers are replicated, with no whole optimizer per worker.
- Direct Green operators use cyclic work assignment and separate result slots;
  sample accumulation stays ordered. Partial QP pools also run the existing
  independent Slater contraction at 16 threads.

A separate final-source two-step timer run at 16 cores measured All=2.980s,
sampling=1.481s, MainCal=1.236s, accepted-hop updates=0.667s,
MainCal full Pfaffian/inverse=0.695s, local energy=0.206s and derivatives=0.256s.
These inclusive, nested timers must not be summed or substituted for the
20-step performance table. Sampling and repeated full factorizations remain
substantial work; the measurements do not attribute the residual gap to one
instruction, allocation or synchronization mechanism. Spinning worker CPU time
is not a measurement of useful arithmetic utilization.

## Numerical and RNG validation

- **112 output-file pairs**, including every warmup and measured repetition at
  every thread count: finite values and shapes, exact integer Green indices,
  per-value `1e-11 + 1e-11*max(abs(x),abs(y))` bounds. Observed maximum absolute
  difference and error/bound ratio are both **0**. This is not a bitwise
  comparison of computed floats. The separately instrumented audit outputs add
  seven passing file pairs with maximum numerical difference 0.
- The full Opt/PhysCal RNG audits at eight threads forward the actual native
  SFMT primitives. Kind-tagged word traces, seed calls, draw counts and the
  next 624 words match exactly. Opt consumes **4,110,359** words
  (4,109,847 UInt32 calls  + 512 real2 calls); PhysCal **206,174**
  (205,662 + 512). Both seed 79. These are before/after Julia audits, not a new
  full C-process state oracle.
- Supplied-word-tape tests check every draw's configuration/counter history,
  exact cursor and saved arrays, and that workers never consume RNG.
  Independent inverse residuals, batch tails/reuse, failure publication,
  Green/factored-Green values and worker teardown on exceptions are separate
  regressions. Other proposal families retain their original scheduling.
- Local optimizer unit suite: **46,085** assertions, plus three later admission
  assertions; four-thread regressions: **397** assertions including those
  admission cases. All three optimizer/C-ctest/PhysCal integration suites pass,
  including 65 Lanczos, 16 PairHop, 66 C-ctest and 31 PhysCal e2e assertions.
  The historical sampler-body hashes remain unchanged. Final Linux/macOS and
  Julia 1.11/1.12/1.13 status is available on [PR #93](https://github.com/tmisawa/Julia-mVMC/pull/93).

## Evidence and reproduction

See [the driver documentation](../../README.md) for isolated before/after
checkout preparation and commands. Use `julia +1.13`; the resolved version is
recorded in each cohort. `raw.tar.gz` contains every timing/input/output cohort,
metadata, validation logs and final timer diagnostic. Its `SHA256SUMS` covers
all included files. Large raw RNG binaries are retained locally and omitted
from this archive; their complete stream hashes, primitive counts, seeds and
next 624 words are included in `audit-*/rng.toml`. `audit.jl` regenerates them
without replacing the RNG. Opt trace SHA-256:
`cb7b6f387d449a8cd5079b668c95176c7bdcaa6ef74a374d38e04cce4f041d37`;
PhysCal: `790f05190122864cbaff0c5bb2eb351c2cc587060111f0c5da259b2131838213`.

```bash
mkdir /tmp/julia-l128-evidence-check
tar -xzf benchmark/cpu_l128/results/20261103/raw.tar.gz -C /tmp/julia-l128-evidence-check
(cd /tmp/julia-l128-evidence-check && sha256sum -c SHA256SUMS)
julia +1.13 --project=. benchmark/cpu_l128/summarize.jl /tmp/julia-l128-evidence-check/raw
for n in 1 4 8 16; do
  julia +1.13 --project=. benchmark/cpu_l128/verify.jl \
    "/tmp/julia-l128-evidence-check/raw/before-t$n" "/tmp/julia-l128-evidence-check/raw/after-t$n"
done
```
