# R0 smoke gate: serial 実行と `mpiexec -n 2` 実行で rank0 出力が bit 一致すること
# （rank0 の seed は RndSeed + group1 = RndSeed + 0 なので serial と同一 chain）。
# 使い方: JULIA_NUM_THREADS=1 julia --project=. test/mpi/run_mpi_smoke.jl
using MPI: mpiexec
using Test

const worker = joinpath(@__DIR__, "mpi_smoke.jl")
const project = abspath(joinpath(@__DIR__, "..", ".."))

@testset "R0 mpiexec -n 2 smoke (spec §9 R0 gate)" begin
    serial_dir = mktempdir()
    mpi_dir = mktempdir()

    # 1) serial 基準（MPI env なし → SerialContext）。
    run(`$(Base.julia_cmd()) --project=$project $worker $serial_dir`)

    # 2) mpiexec -n 2（両 rank が同じ output_dir を受け取る; rank0 のみ書く）。
    # mpiexec() は Cmd を返す（MPI.jl stable）。do-block 形式は deprecated（plan review F3）。
    run(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $worker $mpi_dir`)

    # 3) rank0-only output: ファイル集合が serial と同一（重複 write なし）。
    @test sort(readdir(serial_dir)) == sort(readdir(mpi_dir))

    # 4) bit 一致（rank0 = serial chain）。
    for f in readdir(serial_dir)
        @test read(joinpath(serial_dir, f)) == read(joinpath(mpi_dir, f))
    end
    println("R0 smoke: serial vs mpiexec -n 2 rank0 output is bit-identical")
end

# 任意追加 smoke（plan Step 3.5）: stdout の rank0 gate と JULIA_MVMC_MPI policy の
# end-to-end 検証。R0 gate の必須条件ではないが、F5/F12/A7 の回帰を防ぐ。
@testset "rank0-only stdout under mpiexec -n 2 (plan review F5)" begin
    outdir = mktempdir()
    out = read(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $worker $outdir`,
               String)
    @test count("Progress of Optimization: 0 %", out) == 1
    @test count("Start: Output opt params.", out) == 1
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
