# MPI launcher failure-mode worker.
# Usage:
#   julia --project=<workspace> test/mpi/mpi_failure_modes.jl <mode> <output_dir>
#
# Modes:
#   nsrcg2 - NSRCG >= 2 must be rejected before MPI.Init().
#   nsplit_srcg - NSplitSize > 1 with SR-CG must be rejected before MPI.Init().
using MVMCOptimizers
using MVMCExpertModeParsers
using MPI

length(ARGS) == 2 || error("usage: mpi_failure_modes.jl <nsrcg2|nsplit_srcg> <output_dir>")

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

function run_nsrcg2_rejection()
    data = MVMCExpertModeParsers.parse_expert_mode_files(fixture)
    data.modpara.nsrcg = 2

    expect_error_contains(
        () -> MVMCOptimizers.validate_supported_modpara(data.modpara),
        ("NSRCG >= 2", "standard SR-CG solver"),
    )
    MPI.Initialized() && error("NSRCG >= 2 rejection should happen before MPI.Init()")
    println("failure-mode worker: nsrcg2 expected rejection ok")
end

function run_nsplit_srcg_rejection()
    data = MVMCExpertModeParsers.parse_expert_mode_files(fixture)
    data.modpara.nsplit_size = 2
    data.modpara.nsrcg = 1

    expect_error_contains(
        () -> MVMCOptimizers.validate_supported_para_opt_modpara(data.modpara),
        ("NSplitSize > 1 with SR-CG", "NSRCG = 1"),
    )
    MPI.Initialized() && error("NSplitSize > 1 with SR-CG rejection should happen before MPI.Init()")
    println("failure-mode worker: nsplit_srcg expected rejection ok")
end

try
    if mode == "nsrcg2"
        run_nsrcg2_rejection()
    elseif mode == "nsplit_srcg"
        run_nsplit_srcg_rejection()
    else
        error("unknown mode '$mode'; expected nsrcg2 or nsplit_srcg")
    end
catch err
    @error "mpi_failure_modes worker failed" exception = (err, catch_backtrace())
    if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD, 1)
    end
    rethrow()
end
