# Historical ungated real-QP dispatch investigation

**This initial candidate is superseded.** Its L32 PhysCal regression prevented
adoption. See [WORK_GATE.md](WORK_GATE.md) for the calibrated gate and subsequent
validation. The timings below belong to ungated commit `e86730a`.

Related to [mvmc-rs issue #496](https://github.com/AtelierArith/mvmc-rs/issues/496). The prior SR-store milestone is upstream Julia PR #77,
merged as `d493f113ecc009d44e70327ffcbb500a566a1380`.

## Change and numerical contract

The ordinary real sampling wrapper already selects a threaded Pfaffian workspace
when inner threading is enabled and `qp_num * n_size^3 >= 131072`. However,
the threaded inverse routine then returns to serial execution whenever QP count
is less than the thread pool size. Eight QPs therefore run serially in a
sixteen-thread pool despite exceeding the existing work gate.

The candidate changes only this latter QP-count guard from `qp_num < nthreads()`
to `qp_num < 2`. Existing static scheduling assigns independent QP planes and
private scratch to useful workers. Nested/worker callers still take the serial
path. The work-size gate, complex path, arithmetic within each matrix, serial
reductions, RNG calls, MPI ownership and native locking policy are unchanged.

Baseline commit: `85f9ea12c49e46639b8284334388116d5f16a337`, whose tracked tree
is identical to merged upstream `d493f113ecc009d44e70327ffcbb500a566a1380`.
Candidate is the isolated `julia-issue496-real-qp-budget` worktree based on that
baseline. Each full-run directory records the exact source diff and input hashes.

## Protocol

Linux x86_64 Dev Container on Ryzen 9 PRO 8945HS (8 physical cores, 16 logical
CPUs). Julia 1.13.1 with its pinned Manifest-v1.13.toml, MPI.jl 0.20.27 / MPICH
4.2.0. Julia OpenBLAS 0.3.30 ILP64 and native PfaPack helper OpenBLAS 0.3.26
LP64 both use one thread. No benchmarks overlap with another CPU-intensive job.

Half-filled periodic Hubbard chain, t=1, U=4, Lsub=4, eight QPs. L32 and L64,
320 total samples, NSplitSize=1 (320 samples per rank at one rank, 80 at four).
Opt uses 300 steps; PhysCal uses 100 groups with fixed C-generated optimized
parameters from the existing rank/thread sweep.

Both paths receive a full warmup. Three pairs alternate baseline/candidate,
candidate/baseline, baseline/candidate. The typed Ref selector is installed
before warmups in both variants, and no Profile run occurs during this primary
comparison. Timing covers the maximum rank's production API plus completion
barrier, including input parsing, initialization and output; process startup,
compilation, selector installation and warmups are excluded.

Run `run-paired-opt.sh 1 4` or `run-paired-physcal.sh 1 4` inside the reference
container with exclusive CPU access. Scripts deliberately refuse existing result
directories. The analyzers require six complete records and validate every energy
output row. Results do not replace historical C/Rust measurements or represent a
fresh three-language comparison.

## Verification status

Completed warmed 300-step Opt results, one MPI rank × sixteen threads, median
of three alternating pairs (seconds):

| Sites | Baseline | Candidate | Change |
|---:|---:|---:|---:|
| 32 | 11.628508343 | 11.023358190 | -5.20% |
| 64 | 57.268995731 | 37.152654182 | -35.13% |

All six pairs favored the candidate. All 300 × 6 output fields across the six
runs per size were finite and within the repository's long-run repeatability
bound (absolute and relative `1e-11`); observed maximum absolute difference was
zero. This output comparison does not replace the pending exact RNG/control
audit. Raw repetitions and observations are in `paired-summary.json`.

- Full Julia unit suite at 16 threads passed 45,338 optimizer assertions, plus
  the base and Slater suites; see `julia-full-tests-t16.log`.
- Expanded independent inverse-residual coverage includes zero QPs, uneven
  QP ranges and more QPs than workers. Integer pivot sentinels prove that workers
  actually execute matrix kernels when there are fewer QPs than pool threads.
- Twenty-step L32/L64, 1 rank × 16 thread audits matched all native SFMT words,
  index, draw/reseed counts, saved integer configurations and acceptance counters.
  Computed output differences were observed as zero; floating values are not
  checked bitwise.
- Full 300-step audits and PhysCal results are pending. Do not infer full
  validation from the short pilot. `audit-full.sh` checks actual unmodified
  baseline/candidate sources, without the timing selector.

The short L64 20-step pilot measured 4.5983 → 3.1097 seconds (-32.37%); this is
not the requested full 300-step benchmark. `failed-parse.log` preserves an initial
harness parse failure that occurred before timing; trimming the extracted method
before `Meta.parse` corrected it.

## Portable paired-run command

With the candidate checkout active in a Julia 1.13.1 project and the same
320-sample Expert input, invoke the supplied `paired-julia.jl` as follows:

```sh
OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
JULIA_NUM_THREADS=16,0 JULIA_NUM_GC_THREADS=1 JULIA_MVMC_MPI=1 \
JULIA_MVMC_INNER_THREADS=1 JULIA_MVMC_PFAPACK_THREADS=0 JULIA_MVMC_PFAPACK_LOCK=1 \
UCX_ERROR_SIGNALS=SIGILL,SIGBUS,SIGFPE UCX_MEMTYPE_CACHE=no \
mpiexec -n 1 julia +1.13.1 --project=. --startup-file=no \
  benchmark/results/real_qp_dispatch_20261010/paired-julia.jl \
  /absolute/path/to/opt-inputs/namelist.def 300 3 /absolute/path/to/output 1
```

The full Opt logs and all measured energy rows accompany this report as gzip
files. The analysis script expects the original run-directory structure;
`paired-summary.json` contains each raw timing and the verified output differences.
The earlier short-pilot timing is labelled separately in that summary. Source,
test, manifest and native-helper hashes are recorded in `provenance.json`.
