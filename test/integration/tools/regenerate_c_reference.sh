#!/bin/bash
# Regenerate the C-native integration references under test/integration/reference
# from an UNMODIFIED C mVMC vmc.out (Julia-mVMC#95: mVMC 1.4.0, v1.4.0 / e83aa6cb).
#
# Optional developer command; no Julia test runs it.
#
# Usage (Linux x86_64 with MPICH; the reference environment is recorded in
# reference/C_1.4.0_PROVENANCE.txt):
#   test/integration/tools/regenerate_c_reference.sh <path/to/vmc.out> <output-dir>
# then compare <output-dir> with reference/ and copy the files over.
#
# Inputs are read unchanged from reference/ (copied into a staging directory).
# One MPI rank unless noted, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1.
# - PhysCal (`<model>/physcal_ref/expected/`): `vmc.out -e namelist.def zqp_opt.dat`.
# - First optimization rows: `vmc.out -e namelist.def [initial.def]` on
#   `<model>/inputs` with NSROptItrStep set to the number of rows kept (rows of
#   earlier steps do not depend on the total step count); `initial.def` is passed
#   when present, as the Julia runner loads it by default.
# - NSRCG=1 first step (`heisenberg_chain_real_nsrcg`): serial and `mpiexec -n 2`
#   zvo_out row, raw one-step zqp_opt.dat and (MPI) zvo_SRinfo.dat.
set -euo pipefail
vmc=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
out=${2:?output directory}
ref=$(cd "$(dirname "$0")/../reference" && pwd)
mkdir -p "$out"
out=$(cd "$out" && pwd)
stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run_c() { # nproc dir args...
  local np=$1 dir=$2
  shift 2
  (cd "$dir" && mpiexec -n "$np" "$vmc" "$@" > stdout.txt 2>&1)
}

for model in heisenberg_chain_real heisenberg_chain_cmp heisenberg_chain_fsz kondo_chain_real \
  hubbard_chain_real hubbard_chain_dh_real hubbard_chain_lanczos spin_chain_lanczos; do
  work="$stage/$model-physcal"
  cp -r "$ref/$model/physcal_ref/inputs" "$work"
  cp "$ref/$model/physcal_ref/zqp_opt.dat" "$work/zqp_opt.dat"
  run_c 1 "$work" -e namelist.def zqp_opt.dat
  dest="$out/$model/physcal_ref/expected"
  mkdir -p "$dest"
  for f in "$ref/$model/physcal_ref/expected"/*; do
    cp "$work/output/$(basename "$f")" "$dest/"
  done
done

first_rows() { # model rows file
  local model=$1 rows=$2 file=$3
  work="$stage/$model-opt$rows"
  cp -r "$ref/$model/inputs" "$work"
  sed -i -E "s/^(NSROptItrStep[[:space:]]+).*/\1$rows/" "$work/modpara.def"
  if [ -f "$work/initial.def" ]; then
    run_c 1 "$work" -e namelist.def initial.def
  else
    run_c 1 "$work" -e namelist.def
  fi
  mkdir -p "$out/$model"
  head -n "$rows" "$work/output/zvo_out_001.dat" > "$out/$model/$file"
}
for model in heisenberg_chain_real heisenberg_chain_cmp heisenberg_chain_fsz hubbard_chain_real; do
  first_rows "$model" 10 zvo_out_first10.dat
done
for model in hubbard_chain_pairhop_real hubbard_chain_pairhop_fsz; do
  first_rows "$model" 1 zvo_out_first1.dat
done

model=heisenberg_chain_real_nsrcg
mkdir -p "$out/$model"
for np in 1 2; do
  work="$stage/$model-np$np"
  cp -r "$ref/$model/inputs" "$work"
  run_c "$np" "$work" -e namelist.def
  if [ "$np" = 1 ]; then
    head -n 1 "$work/output/zvo_out_001.dat" > "$out/$model/zvo_out_first1.dat"
    cp "$work/output/zqp_opt.dat" "$out/$model/zqp_opt_1step.dat"
  else
    head -n 1 "$work/output/zvo_out_001.dat" > "$out/$model/zvo_out_mpi2_first1.dat"
    cp "$work/output/zqp_opt.dat" "$out/$model/zqp_opt_mpi2_1step.dat"
    cp "$work/output/zvo_SRinfo.dat" "$out/$model/zvo_SRinfo_mpi2_1step.dat"
  fi
done

# hubbard_chain_real_nsrcg: hand-authored NSRCG=1 first step of hubbard_chain_real
# (its projection block has a nonzero S diagonal, so RescaleSmat acts); its
# references are the variants below (`plain` = the fixture modpara).

# C 1.4.0 SR-CG option variants of the NSRCG=1 first steps (Julia-mVMC#95):
# `variants/<name>/modpara.def` replaces the fixture modpara (later lines win, as
# in C). Serial; the C stderr and exit status are kept because a variant may abort.
for model in heisenberg_chain_real_nsrcg hubbard_chain_real_nsrcg; do
  for vdir in "$ref/$model"/variants/*/; do
    name=$(basename "$vdir")
    work="$stage/$model-variant-$name"
    cp -r "$ref/$model/inputs" "$work"
    cp "$vdir/modpara.def" "$work/modpara.def"
    args=(-e namelist.def)
    [ -f "$work/initial.def" ] && args+=(initial.def)
    set +e
    (cd "$work" && mpiexec -n 1 "$vmc" "${args[@]}" > stdout.txt 2> stderr.txt)
    status=$?
    set -e
    dest="$out/$model/variants/$name"
    mkdir -p "$dest"
    echo "$status" > "$dest/exit_status.txt"
    grep -E '^(warning|error|remark):' "$work/stderr.txt" > "$dest/stderr.txt" || true
    cp "$work/output/zvo_SRinfo.dat" "$dest/zvo_SRinfo.dat"
    if [ "$status" = 0 ]; then
      head -n 1 "$work/output/zvo_out_001.dat" > "$dest/zvo_out_first1.dat"
      cp "$work/output/zqp_opt.dat" "$dest/zqp_opt_1step.dat"
    fi
  done
done

{
  echo "generator=test/integration/tools/regenerate_c_reference.sh"
  echo "generated_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "source=${SOURCE_DESC:-unknown}"
  echo "vmc_out_sha256=$(sha256sum "$vmc" | cut -d' ' -f1)"
  echo "arch=$(uname -m) kernel=$(uname -sr)"
  echo "gcc=$(gcc --version | head -1)"
  echo "gfortran=$(gfortran --version | head -1)"
  echo "mpich=$(mpichversion 2>/dev/null | head -1 | tr -s ' ')"
  echo "openblas=$(dpkg-query -W -f='${Version}' libopenblas0-pthread 2>/dev/null || echo unknown)"
  echo "cmake_flags=${CMAKE_FLAGS:-unknown}"
  echo "threads=OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1, mpiexec -n 1 (nsrcg MPI files: -n 2)"
  echo "image=${IMAGE_DESC:-unknown}"
} > "$out/provenance.txt"
