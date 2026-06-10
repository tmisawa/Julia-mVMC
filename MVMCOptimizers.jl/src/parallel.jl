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
