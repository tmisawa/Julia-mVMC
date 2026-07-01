# MPI launcher failure-mode worker.
# Usage:
#   julia --project=<workspace> test/mpi/mpi_failure_modes.jl <mode> <output_dir>
#
# Modes:
#   nsrcg2 - NSRCG >= 2 must be rejected before MPI.Init().
#   nsplit_srcg - NSplitSize > 1 with SR-CG must be rejected before MPI.Init().
#   nsplit_opttrans - NSplitSize > 1 with OptTrans must be rejected before MPI.Init().
#   nsplit_fsz_projection - FSZ standard-projection NQPFull > 1 must be rejected before MPI.Init().
using MVMCOptimizers
using MVMCExpertModeParsers
using MPI

length(ARGS) == 2 ||
    error("usage: mpi_failure_modes.jl <nsrcg2|nsplit_srcg|nsplit_opttrans|nsplit_fsz_projection> <output_dir>")

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

function run_nsplit_opttrans_rejection()
    data = MVMCExpertModeParsers.parse_expert_mode_files(fixture)
    data.modpara.nsplit_size = 2
    data.n_qp_opt_trans = 2
    data.opt_trans = ComplexF64[1.0 + 0.0im, 0.5 + 0.0im]
    data.qp_opt_trans = [[0], [0]]

    expect_error_contains(
        () -> MVMCOptimizers.validate_supported_para_opt_data(data),
        ("NSplitSize > 1 with NQPOptTrans > 1", "OptTrans"),
    )
    MPI.Initialized() && error("OptTrans split rejection should happen before MPI.Init()")
    println("failure-mode worker: nsplit_opttrans expected rejection ok")
end

function run_nsplit_fsz_projection_rejection()
    data = MVMCExpertModeParsers.ExpertModeData()
    data.modpara.nsplit_size = 2
    data.modpara.nsp_gauss_leg = 1
    data.modpara.nmp_trans = 2
    data.n_qp_opt_trans = 1
    data.i_flg_orbital_general = 1

    expect_error_contains(
        () -> MVMCOptimizers.validate_supported_para_opt_data(data),
        ("NSplitSize > 1 with FSZ standard-projection NQPFull > 1", "NMPTrans = 2"),
    )
    MPI.Initialized() && error("FSZ projection split rejection should happen before MPI.Init()")
    println("failure-mode worker: nsplit_fsz_projection expected rejection ok")
end

try
    if mode == "nsrcg2"
        run_nsrcg2_rejection()
    elseif mode == "nsplit_srcg"
        run_nsplit_srcg_rejection()
    elseif mode == "nsplit_opttrans"
        run_nsplit_opttrans_rejection()
    elseif mode == "nsplit_fsz_projection"
        run_nsplit_fsz_projection_rejection()
    else
        error("unknown mode '$mode'; expected nsrcg2, nsplit_srcg, nsplit_opttrans, or nsplit_fsz_projection")
    end
catch err
    @error "mpi_failure_modes worker failed" exception = (err, catch_backtrace())
    if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD, 1)
    end
    rethrow()
end
