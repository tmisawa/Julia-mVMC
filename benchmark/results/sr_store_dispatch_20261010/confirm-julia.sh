#!/usr/bin/env bash
set -euo pipefail
cd /workspaces/mvmc-rs
root=/workspaces/mvmc-rs/bench-out/issue496-thread-granularity-20261010
cache=/home/vscode/.cache/mvmc
export OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export JULIA_NUM_GC_THREADS=1 JULIA_MVMC_MPI=1 JULIA_MVMC_INNER_THREADS=1 JULIA_MVMC_PFAPACK_THREADS=0
export UCX_ERROR_SIGNALS=SIGILL,SIGBUS,SIGFPE UCX_MEMTYPE_CACHE=no
git -C "$cache/julia-issue496-thread-budget" diff > "$root/paired-source-pre.diff"
sha256sum "$cache/julia-issue496-thread-budget/MVMCOptimizers.jl/src/"*.jl > "$root/paired-source-pre.sha256"
for pair in '1 16' '4 4'; do
    read -r ranks threads <<< "$pair"
    for size in 32 64; do
        inputs="/workspaces/mvmc-rs/bench-out/rank-thread-sweep-20261010/r$ranks-t$threads/L$size/opt-inputs"
        expected_samples=$((320/ranks))
        rg -q "^NVMCSample[[:space:]]+$expected_samples$" "$inputs/modpara.def"
        dest="$root/paired/L$size-r$ranks-t$threads"
        mkdir -p "$dest"
        sha256sum "$inputs/"*.def > "$dest/input-sha256.txt"
        JULIA_NUM_THREADS="$threads,0" /opt/mpich/bin/mpiexec -n "$ranks" \
            "$cache/tools/julia-1.13.1/bin/julia" --project="$cache/julia-issue496-thread-budget-project" \
            --startup-file=no "$root/paired-julia.jl" "$inputs/namelist.def" 300 3 "$dest" "$ranks" \
            > "$dest/run.log" 2>&1
        rg '^PAIRED|^WORLD|^THREADS' "$dest/run.log"
    done
done
sha256sum --check "$root/paired-source-pre.sha256" > "$root/paired-source-post-check.log"
git -C "$cache/julia-issue496-thread-budget" diff > "$root/paired-source-post.diff"
cmp "$root/paired-source-pre.diff" "$root/paired-source-post.diff"
