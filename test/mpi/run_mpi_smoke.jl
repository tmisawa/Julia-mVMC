# MPI smoke gate:
# - R0: rank0-only output/readback and launch-policy regressions.
# - R1: `mpiexec -n 2` runs independent chains and comm0 allreduce, so MPI output
#   is intentionally not bit-identical to serial output.
# 使い方: JULIA_NUM_THREADS=1 julia --project=. test/mpi/run_mpi_smoke.jl
using MPI: mpiexec
using Test

const worker = joinpath(@__DIR__, "mpi_smoke.jl")
const hubbard_worker = joinpath(@__DIR__, "mpi_hubbard_smoke.jl")
const physcal_worker = joinpath(@__DIR__, "mpi_physcal_smoke.jl")
const weight_average_worker = joinpath(@__DIR__, "mpi_weight_average_smoke.jl")
const srcg_operate_worker = joinpath(@__DIR__, "mpi_srcg_operate_smoke.jl")
const srcg_e2e_worker = joinpath(@__DIR__, "mpi_srcg_e2e_smoke.jl")
const failure_worker = joinpath(@__DIR__, "mpi_failure_modes.jl")
const project = abspath(joinpath(@__DIR__, "..", ".."))
const para_opt_files = (
    "zqp_gutzwiller_opt.dat",
    "zqp_jastrow_opt.dat",
    "zqp_opt.dat",
    "zqp_orbital_opt.dat",
    "zvo_out.dat",
    "zvo_var.dat",
)
const physcal_files = (
    "zvo_out.dat",
    "zvo_var.dat",
    "zvo_cisajs_001.dat",
    "zvo_cisajscktalt_001.dat",
    "zvo_cisajscktaltex_001.dat",
)

function assert_files_present(dir::AbstractString, names)
    for name in names
        @test isfile(joinpath(dir, name))
    end
end

@testset "R1 mpiexec -n 2 smoke (rank0 output + allreduce path)" begin
    serial_dir = mktempdir()
    mpi_dir = mktempdir()

    # 1) serial sanity（MPI env なし → SerialContext）。
    run(`$(Base.julia_cmd()) --project=$project $worker $serial_dir`)
    @test !isempty(readdir(serial_dir))

    # 2) mpiexec -n 2（両 rank が同じ output_dir を受け取る; rank0 のみ書く）。
    # R1 では comm0 allreduce で rank0 output は serial と別値になり得る。
    run(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $worker $mpi_dir`)

    # 3) rank0-only output: ファイル集合が serial と同一（重複 write なし）。
    @test sort(readdir(serial_dir)) == sort(readdir(mpi_dir))
    zvo_lines = readlines(joinpath(mpi_dir, "zvo_out.dat"))
    @test length(zvo_lines) == 4
    @test all(line -> length(split(strip(line))) >= 2, zvo_lines)
    println("R1 smoke: mpiexec -n 2 rank0 output completed through allreduce path")
end

@testset "R1 mpiexec -n 2 WeightAverageWE numeric smoke" begin
    out = read(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $weight_average_worker`,
               String)
    @test count("weight-average worker: root rank ok", out) == 1
    @test count("weight-average worker: non-root rank ok", out) == 1
end

@testset "R1 mpiexec -n 2 Hubbard para-opt smoke" begin
    mpi_dir = mktempdir()
    out = read(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $hubbard_worker $mpi_dir`,
               String)
    assert_files_present(mpi_dir, para_opt_files)
    zvo_lines = readlines(joinpath(mpi_dir, "zvo_out.dat"))
    @test length(zvo_lines) == 4
    @test all(line -> length(split(strip(line))) >= 2, zvo_lines)
    @test count("hubbard worker: root rank ok", out) == 1
    @test count("hubbard worker: non-root rank ok", out) == 1
    println("R1 smoke: mpiexec -n 2 Hubbard para-opt completed")
end

@testset "R2 mpiexec -n 2 WeightAverage warning rank0 gate" begin
    cmd = addenv(
        `$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $weight_average_worker`,
        "JULIA_MVMC_SMOKE_LOG_STDOUT" => "1",
        "JULIA_MVMC_SMOKE_TINY_WC" => "1",
    )
    out = read(cmd, String)
    @test count("Weight Wc is too small after MPI allreduce", out) == 1
    @test count("weight-average worker: root rank ok", out) == 1
    @test count("weight-average worker: non-root rank ok", out) == 1
end

@testset "v0.5 mpiexec -n 2 SR-CG operate_by_s collective smoke" begin
    out = read(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $srcg_operate_worker`,
               String)
    @test count("srcg-operate worker: root rank ok", out) == 1
    @test count("srcg-operate worker: non-root rank ok", out) == 1
end

@testset "v0.5 mpiexec -n 2 SR-CG e2e C-reference smoke" begin
    mpi_dir = mktempdir()
    out = read(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $srcg_e2e_worker $mpi_dir`,
               String)
    assert_files_present(mpi_dir, para_opt_files)
    @test isfile(joinpath(mpi_dir, "zvo_SRinfo.dat"))
    @test count("srcg-e2e worker: root rank ok", out) == 1
    @test count("srcg-e2e worker: non-root rank ok", out) == 1
end

@testset "v0.5 failure mode: NSRCG >= 2 rejects before MPI.Init" begin
    outdir = mktempdir()
    out = read(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $failure_worker nsrcg2 $outdir`,
               String)
    @test count("failure-mode worker: nsrcg2 expected rejection ok", out) == 2
    @test isempty(readdir(outdir))
end

@testset "v0.5 failure mode: NSplitSize > 1 rejects before MPI.Init" begin
    outdir = mktempdir()
    out = read(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $failure_worker nsplit $outdir`,
               String)
    @test count("failure-mode worker: nsplit expected rejection ok", out) == 2
    @test isempty(readdir(outdir))
end

@testset "R1 mpiexec -n 4 para-opt smoke" begin
    mpi_dir = mktempdir()
    run(`$(mpiexec()) -n 4 $(Base.julia_cmd()) --project=$project $worker $mpi_dir`)
    assert_files_present(mpi_dir, para_opt_files)
    zvo_lines = readlines(joinpath(mpi_dir, "zvo_out.dat"))
    @test length(zvo_lines) == 4
    @test all(line -> length(split(strip(line))) >= 2, zvo_lines)
    println("R1 smoke: mpiexec -n 4 para-opt completed")
end

@testset "R1 mpiexec -n 2 PhysCal smoke (rank0 Green output)" begin
    mpi_dir = mktempdir()
    out = read(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $physcal_worker $mpi_dir`,
               String)
    assert_files_present(mpi_dir, physcal_files)
    @test count("Start: Calculate VMC physical quantities.", out) == 1
    @test count("End  : Calculate VMC physical quantities.", out) == 1
    @test count("physcal worker: root rank ok", out) == 1
    @test count("physcal worker: non-root rank ok", out) == 1
    println("R1 smoke: mpiexec -n 2 PhysCal completed through Green reduce path")
end

# 任意追加 smoke（plan Step 3.5）: stdout の rank0 gate と JULIA_MVMC_MPI policy の
# end-to-end 検証。R0 gate の必須条件ではないが、F5/F12/A7 の回帰を防ぐ。
@testset "rank0-only stdout under mpiexec -n 2 (plan review F5)" begin
    outdir = mktempdir()
    cmd = addenv(
        `$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $worker $outdir`,
        "JULIA_MVMC_SMOKE_LOG_STDOUT" => "1",
    )
    out = read(cmd, String)
    @test count("Progress of Optimization: 0 %", out) == 1
    @test count("Start: Output opt params.", out) == 1
    @test count("NPara=14", out) == 1
    @test count("worker: non-root rank ok", out) == 1
    @test count("worker: root rank ok", out) == 1
end

@testset "JULIA_MVMC_MPI=0 under mpiexec -n 2 aborts (F12)" begin
    outdir = mktempdir()
    cmd = addenv(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $worker $outdir`,
                 "JULIA_MVMC_MPI" => "0")
    proc = run(pipeline(ignorestatus(cmd); stdout = devnull, stderr = devnull))
    @test proc.exitcode != 0   # hang せず nonzero exit すること
end

@testset "JULIA_MVMC_MPI=1 without mpiexec runs as size-1 MPI (A7)" begin
    outdir = mktempdir()
    cmd = addenv(`$(Base.julia_cmd()) --project=$project $worker $outdir`,
                 "JULIA_MVMC_MPI" => "1")
    out = read(cmd, String)
    @test count("worker: root rank ok", out) == 1   # size=1 の MPI context で完走
end
