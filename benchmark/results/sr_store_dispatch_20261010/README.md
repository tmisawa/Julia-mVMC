# SR sample-store dispatch (Linux x86_64, 2026-10-10)

Small NStore vectors previously started a static Julia thread region for every sample.
The candidate retains the elementwise C storage arithmetic and uses a vector-size gate:
65,536 real values or 32,768 complex values. Large stores retain threaded execution.

Periodic half-filled Hubbard chain, t=1, U=4, Lsub=4, eight QPs. Opt300, total320 samples.
Rank4 inputs use 80 samples per rank (NSplitSize=1). Ryzen 9 PRO 8945HS: eight physical
cores / sixteen logical CPUs. Julia 1.13.1, MPI.jl 0.20.27 / MPICH4.2; BLAS1 in both
Julia OpenBLAS0.3.30 ILP64 and the native helper OpenBLAS0.3.26 LP64.

| MPI ranks | Threads/rank | Sites | Baseline s | Candidate s | Change |
|---:|---:|---:|---:|---:|---:|
| 1 | 16 | 32 | 17.803 | 11.627 | -34.7% |
| 1 | 16 | 64 | 60.181 | 56.652 | -5.9% |
| 4 | 4 | 32 | 3.634 | 3.494 | -3.8% |
| 4 | 4 | 64 | 13.509 | 13.311 | -1.5% |

Each process warmed both paths, then measured three alternating pairs. Every candidate
pair was faster (12/12). Times are the maximum MPI-rank production API duration;
startup, JIT, build and warmups are excluded. No profiler runs during this confirmation.
A typed Ref gate selector was installed before both warmups: baseline uses the original
64-item gate; candidate uses the new thresholds. Its dispatch overhead is common.
Separate unmodified-production pilots gave L32 / rank1 / threads16 medians of
18.139 s baseline and 11.653 s candidate. Small four-rank changes should be interpreted
conservatively; they are not a claim of universal scaling or optimal thresholds.

Baseline source: 043ec52b9c9a1542083c233bd1923fe24f36d03c, same reviewed tree as upstream
main 02afdae0a19732c727e0d09132d9761c57d68fcb. Exact source/library hashes and container
identity are in provenance.json and container.txt. Source patches and full optional
oracle evidence are retained in the Rust workspace investigation for issue #496.

Verification: 41,392 optimizer unit assertions plus 25 Slater assertions and 15 base
assertions passed. The new 53-assertion SR-store suite also passed at one and sixteen
threads. Native SFMT 624 words/index/reseed/draw count, all saved and burn/temp electron
configurations and acceptance/rejection counters matched for all eight size/layout/
20-or-300-step cases (20 paired ranks). Floating output maximum absolute difference
was observed as zero; numerical tests use explicit justified tolerances rather than
bitwise floating assertions. Do not regenerate historical fixtures from these runs.

Reproduce the kernel crossover with benchmark/sr_store_dispatch.jl and the command
in its header. Whole-run confirmation uses paired-julia.jl and confirm-julia.sh;
adjust the isolated project/input paths to the local checked-out baseline/candidate.
Use the pinned Manifest-v1.13.toml and the same native provider/thread settings.

Related to https://github.com/AtelierArith/mvmc-rs/issues/496.
