# mpiexec 配下（または serial）で run_phys_cal_from_namelist を 1 回実行する worker。
# 使い方: julia --project=<workspace> test/mpi/mpi_physcal_smoke.jl <output_dir>
using MVMCOptimizers
using MPI

const refdir = joinpath(@__DIR__, "..", "integration", "reference",
                        "heisenberg_chain_real", "physcal_ref")
const namelist = joinpath(refdir, "inputs", "namelist.def")
const opt_para = joinpath(refdir, "zqp_opt.dat")
const outdir = ARGS[1]
const expected_files = (
    "zvo_out.dat",
    "zvo_var.dat",
    "zvo_cisajs_001.dat",
    "zvo_cisajscktalt_001.dat",
    "zvo_cisajscktaltex_001.dat",
)

# top-level abort guard: 片 rank だけが exception で死ぬと他 rank が collective /
# barrier で hang するため、worker 全体を try/catch で包む。
try
    result = MVMCOptimizers.run_phys_cal_from_namelist(
        namelist; opt_para = opt_para, mode = :real, output_dir = outdir)

    rank = MPI.Initialized() && !MPI.Finalized() ? MPI.Comm_rank(MPI.COMM_WORLD) : 0
    @assert result.status == 0
    @assert result.n_para_consumed == 14

    if rank == 0
        for name in expected_files
            @assert isfile(joinpath(outdir, name)) "missing PhysCal output: $name"
        end
        println("physcal worker: root rank ok")
    else
        println("physcal worker: non-root rank ok")
    end
catch err
    @error "mpi_physcal_smoke worker failed" exception = (err, catch_backtrace())
    if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD, 1)
    end
    rethrow()
end
