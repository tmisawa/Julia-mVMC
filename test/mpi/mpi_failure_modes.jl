# MPI launcher failure-mode worker.
# Usage:
#   julia --project=<workspace> test/mpi/mpi_failure_modes.jl <mode> <output_dir>
#
# Modes:
#   nsrcg  - NSRCG != 0 must be rejected under a real MPI context.
#   nsplit - NSplitSize > 1 must be rejected before MPI.Init().
using MVMCOptimizers
using MVMCExpertModeParsers
using MPI

length(ARGS) == 2 || error("usage: mpi_failure_modes.jl <nsrcg|nsplit> <output_dir>")

const mode = ARGS[1]
const outdir = ARGS[2]
const fixture = joinpath(@__DIR__, "..", "integration", "reference",
                         "hubbard_chain_real", "inputs", "namelist.def")

function expect_error_contains(f, pieces)
    try
        f()
    catch err
        msg = sprint(showerror, err)
        for piece in pieces
            occursin(piece, msg) || error("expected error message to contain '$piece'; got: $msg")
        end
        return msg
    end
    error("expected an error containing $(join(pieces, ", "))")
end

function run_nsrcg_rejection()
    data = MVMCExpertModeParsers.parse_expert_mode_files(fixture)
    data.modpara.nsrcg = 1
    data.modpara.nsr_opt_itr_step = 1
    data.modpara.nsr_opt_itr_smp = 1
    ctx = MVMCOptimizers.build_parallel_context(data.modpara.nsplit_size)

    expect_error_contains(
        () -> MVMCOptimizers.vmc_para_opt!(data; ctx = ctx, output_dir = outdir),
        ("NSRCG != 0", "operate_by_S broadcast/allreduce"),
    )
    println("failure-mode worker: nsrcg expected rejection ok")
end

function run_nsplit_rejection()
    data = MVMCExpertModeParsers.parse_expert_mode_files(fixture)
    data.modpara.nsplit_size = 2

    expect_error_contains(
        () -> MVMCOptimizers.validate_supported_modpara(data.modpara),
        ("NSplitSize > 1", "grouped MPI/QP splitting by NSplitSize is not implemented"),
    )
    MPI.Initialized() && error("NSplitSize > 1 rejection should happen before MPI.Init()")
    println("failure-mode worker: nsplit expected rejection ok")
end

try
    if mode == "nsrcg"
        run_nsrcg_rejection()
    elseif mode == "nsplit"
        run_nsplit_rejection()
    else
        error("unknown mode '$mode'; expected nsrcg or nsplit")
    end
catch err
    @error "mpi_failure_modes worker failed" exception = (err, catch_backtrace())
    if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD, 1)
    end
    rethrow()
end
