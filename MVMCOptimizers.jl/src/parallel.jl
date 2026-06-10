"""
Parallel (MPI) infrastructure for Julia-mVMC v0.4.

Mirrors C-mVMC's communicator structure (mVMC/src/mVMC/vmcmain.c:239-257):
comm0 = whole run, comm1 = QP/sample split within a group (size NSplitSize),
comm2 = cross-group counter reduction. Serial runs (no mpiexec) use a no-op
SerialContext and never touch MPI APIs, preserving v0.3 bit parity.

Design doc: docs/specs/2026-06-10-julia-mvmc-v0.4-mpi-design.md (v3.1)
"""

using MPI: MPI

# C safempi.c:29  #define D_MpiSendMax 1048576 (elements per allreduce chunk)
const D_MPI_SEND_MAX = 1_048_576

"""
    split_loop(loop_length, mpi_rank, mpi_size) -> (ist, ien)

C `SplitLoop` (mVMC/src/mVMC/splitloop.c:33-63) と同一の分割。0-based half-open
`[ist, ien)` を返す。`mpi_size >= loop_length` のとき `mpi_rank >= loop_length` の
rank は空 range を受け取る。`loop_length == 0` も C と同じく全 rank 空 range。
"""
function split_loop(loop_length::Int, mpi_rank::Int, mpi_size::Int)
    if mpi_size < loop_length
        imod = loop_length % mpi_size
        if imod == 0
            idiv = loop_length ÷ mpi_size
            ist = idiv * mpi_rank
            ien = ist + idiv
        else
            idiv = (loop_length - imod) ÷ mpi_size
            if mpi_rank < mpi_size - imod
                ist = idiv * mpi_rank
                ien = ist + idiv
            else
                ist = idiv * mpi_rank + mpi_rank - (mpi_size - imod)
                ien = ist + idiv + 1
            end
        end
    else
        if mpi_rank < loop_length
            ist = mpi_rank
            ien = mpi_rank + 1
        else
            ist = loop_length
            ien = loop_length
        end
    end
    return ist, ien
end

"""
    split_range(loop_length, mpi_rank, mpi_size) -> UnitRange{Int}

`split_loop` の 1-based Julia range 版（`(ist+1):ien`）。
"""
function split_range(loop_length::Int, mpi_rank::Int, mpi_size::Int)
    ist, ien = split_loop(loop_length, mpi_rank, mpi_size)
    return (ist + 1):ien
end

"""
    ParallelContext

C-mVMC の comm0/comm1/comm2 構造（mVMC/src/mVMC/vmcmain.c:239-246）の直写し。
serial 実行では `serial_context()`（is_mpi=false、全 size=1）を使い、通信 wrapper は
すべて no-op になる。
"""
struct ParallelContext
    is_mpi::Bool
    comm0::Union{Nothing,MPI.Comm}
    comm1::Union{Nothing,MPI.Comm}   # QP / sample 分担（グループ内、size=NSplitSize）
    comm2::Union{Nothing,MPI.Comm}   # counter 集約（グループ横断）
    rank0::Int
    size0::Int
    rank1::Int
    size1::Int
    rank2::Int
    size2::Int
    group1::Int                      # = rank0 ÷ NSplitSize（seed offset, vmcmain.c:257）
end

serial_context() =
    ParallelContext(false, nothing, nothing, nothing, 0, 1, 0, 1, 0, 1, 0)

"出力（zvo_*.dat / readback / 進捗 print）を担う rank か（C: `if(rank==0)`)。"
is_output_rank(ctx::ParallelContext) = ctx.rank0 == 0

# mpiexec 起動の検出に使う環境変数（spec §4.1）。
const MPI_ENV_KEYS = ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMI_RANK")

mpi_env_detected() = any(k -> haskey(ENV, k), MPI_ENV_KEYS)

"""
    resolve_mpi_mode() -> :serial | :mpi | :mpi_guarded_serial

`JULIA_MVMC_MPI` の 3 値 semantics（spec §4.1、F12+A7）:
auto/unset → 検出時のみ :mpi。`0` → 未検出なら :serial、検出時は
:mpi_guarded_serial（MPI.Init 後に size>1 なら error+Abort、size==1 なら serial）。
`1` → 常に :mpi（Init 失敗はそのまま error）。
"""
function resolve_mpi_mode()
    raw = strip(get(ENV, "JULIA_MVMC_MPI", ""))
    detected = mpi_env_detected()
    if raw == "1"
        return :mpi
    elseif raw == "0"
        return detected ? :mpi_guarded_serial : :serial
    elseif raw == "" || raw == "auto"
        return detected ? :mpi : :serial
    else
        error("JULIA_MVMC_MPI must be \"0\", \"1\", \"auto\", or unset; got \"$raw\"")
    end
end

"""
    build_parallel_context(nsplit_size) -> ParallelContext

C vmcmain.c:239-257 の comm split を再現する。MPI mode では `MPI.Initialized()`
guard 付きで `MPI.Init()` し（二重 Init 回避、F13）、library 側からは
`MPI.Finalize()` を呼ばない（MPI.jl が Julia exit 時に自動 finalize）。
`size0 % nsplit_size != 0` は C と同じく warning のみで継続する（F8）。
"""
function build_parallel_context(nsplit_size::Int)
    mode = resolve_mpi_mode()
    mode === :serial && return serial_context()

    MPI.Initialized() || MPI.Init()
    comm0 = MPI.Comm_dup(MPI.COMM_WORLD)
    rank0 = MPI.Comm_rank(comm0)
    size0 = MPI.Comm_size(comm0)

    if mode === :mpi_guarded_serial
        # F12: mpiexec 配下の JULIA_MVMC_MPI=0。size>1 は出力破壊防止のため abort。
        if size0 > 1
            rank0 == 0 && @error "JULIA_MVMC_MPI=0 under a multi-rank mpiexec run is " *
                                 "not allowed (each rank would write the same output " *
                                 "files). Unset JULIA_MVMC_MPI or run without mpiexec."
            MPI.Abort(comm0, 1)
        end
        return serial_context()
    end

    nsplit = max(nsplit_size, 1)
    if size0 % nsplit != 0 && rank0 == 0
        @warn "load imbalance. MPI size0=$size0 NSplitSize=$nsplit"   # C vmcmain.c:248-250
    end
    group1 = rank0 ÷ nsplit
    comm1 = MPI.Comm_split(comm0, group1, rank0)
    rank1 = MPI.Comm_rank(comm1)
    size1 = MPI.Comm_size(comm1)
    group2 = rank1
    comm2 = MPI.Comm_split(comm0, group2, rank0)
    rank2 = MPI.Comm_rank(comm2)
    size2 = MPI.Comm_size(comm2)
    return ParallelContext(true, comm0, comm1, comm2,
                           rank0, size0, rank1, size1, rank2, size2, group1)
end
