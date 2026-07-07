# MPI smoke gate:
# - R0: rank0-only output/readback and launch-policy regressions.
# - R1: `mpiexec -n 2` runs independent chains and comm0 allreduce, so MPI output
#   is intentionally not bit-identical to serial output.
# 使い方: JULIA_NUM_THREADS=1 julia --project=. test/mpi/run_mpi_smoke.jl
using MPI: mpiexec
using Test

const worker = joinpath(@__DIR__, "mpi_smoke.jl")
const hubbard_worker = joinpath(@__DIR__, "mpi_hubbard_smoke.jl")
const nsplit_nstore_worker = joinpath(@__DIR__, "mpi_nsplit_nstore_smoke.jl")
const nsplit_standard_projection_worker =
    joinpath(@__DIR__, "mpi_nsplit_standard_projection_smoke.jl")
const physcal_worker = joinpath(@__DIR__, "mpi_physcal_smoke.jl")
const physcal_nsplit_worker = joinpath(@__DIR__, "mpi_physcal_nsplit_smoke.jl")
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
const nsplit_nstore_tol = 1e-8
# NQPFull > 1 changes the floating-point summation order: nsplit=1 sums QP
# sectors locally in one loop, while nsplit>1 sums rank-local partial IP values
# through comm1 allreduce before taking log(IP). Keep this slightly looser than
# the NQPFull=1 smoke where the QP sum has one term.
const nsplit_standard_projection_tol = 5e-8
const physcal_nsplit_tol = 5e-8

function assert_files_present(dir::AbstractString, names)
    for name in names
        @test isfile(joinpath(dir, name))
    end
end

function parse_numeric_file(path::AbstractString)
    return parse.(Float64, split(read(path, String)))
end

function assert_close_vector(label::AbstractString, actual, expected; atol::Float64)
    @test length(actual) == length(expected)
    if length(actual) == length(expected)
        maxdiff = isempty(actual) ? 0.0 : maximum(abs.(actual .- expected))
        @test maxdiff <= atol
        if maxdiff > atol
            @info "$label maxdiff exceeded tolerance" maxdiff atol
        end
    end
end

function run_nsplit_nstore_case(fixture::AbstractString, mode::AbstractString,
                                nsplit::Int, nstore::Int, nranks::Int)
    mpi_dir = mktempdir()
    out = read(
        `$(mpiexec()) -n $nranks $(Base.julia_cmd()) --project=$project $nsplit_nstore_worker $fixture $mode $nsplit $nstore $mpi_dir`,
        String,
    )
    assert_files_present(mpi_dir, para_opt_files)
    @test length(readlines(joinpath(mpi_dir, "zvo_out.dat"))) == 1
    label = "nsplit-nstore worker: $fixture nsplit=$nsplit nstore=$nstore"
    @test count("$label root rank ok", out) == 1
    @test count("$label non-root rank ok", out) == nranks - 1
    return (
        zvo = parse_numeric_file(joinpath(mpi_dir, "zvo_out.dat")),
        zqp = parse_numeric_file(joinpath(mpi_dir, "zqp_opt.dat")),
    )
end

function run_nsplit_standard_projection_case(
    fixture::AbstractString,
    mode::AbstractString,
    nsplit::Int,
    nstore::Int,
    nranks::Int;
    nsp_gauss_leg::Union{Nothing,Int} = nothing,
)
    mpi_dir = mktempdir()
    cmd = `$(mpiexec()) -n $nranks $(Base.julia_cmd()) --project=$project $nsplit_standard_projection_worker $fixture $mode $nsplit $nstore $mpi_dir`
    if nsp_gauss_leg !== nothing
        cmd = addenv(cmd, "JULIA_MVMC_SMOKE_NSPGAUSSLEG" => string(nsp_gauss_leg))
    end
    out = read(cmd, String)
    assert_files_present(mpi_dir, para_opt_files)
    label = "nsplit-standard-projection worker: $fixture nsplit=$nsplit nstore=$nstore"
    @test count("$label root rank ok", out) == 1
    @test count("$label non-root rank ok", out) == nranks - 1
    return (
        zvo = parse_numeric_file(joinpath(mpi_dir, "zvo_out.dat")),
        zqp = parse_numeric_file(joinpath(mpi_dir, "zqp_opt.dat")),
        var = parse_numeric_file(joinpath(mpi_dir, "zvo_var.dat")),
    )
end

function run_physcal_nsplit_case(
    fixture::AbstractString,
    mode::AbstractString,
    nsplit::Int,
    nranks::Int,
)
    mpi_dir = mktempdir()
    out = read(
        `$(mpiexec()) -n $nranks $(Base.julia_cmd()) --project=$project $physcal_nsplit_worker $fixture $mode $nsplit $mpi_dir`,
        String,
    )
    assert_files_present(mpi_dir, physcal_files)
    label = "physcal-nsplit worker: $fixture nsplit=$nsplit"
    @test count("$label root rank ok", out) == 1
    @test count("$label non-root rank ok", out) == nranks - 1
    return (
        zvo = parse_numeric_file(joinpath(mpi_dir, "zvo_out.dat")),
        var = parse_numeric_file(joinpath(mpi_dir, "zvo_var.dat")),
        cisajs = parse_numeric_file(joinpath(mpi_dir, "zvo_cisajs_001.dat")),
        two_body = parse_numeric_file(joinpath(mpi_dir, "zvo_cisajscktalt_001.dat")),
        factored = parse_numeric_file(joinpath(mpi_dir, "zvo_cisajscktaltex_001.dat")),
    )
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

@testset "v0.5 failure mode: NSplitSize > 1 with SR-CG rejects before MPI.Init" begin
    outdir = mktempdir()
    out = read(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $failure_worker nsplit_srcg $outdir`,
               String)
    @test count("failure-mode worker: nsplit_srcg expected rejection ok", out) == 2
    @test isempty(readdir(outdir))
end

@testset "v0.5 failure mode: NSplitSize > 1 with OptTrans rejects before MPI.Init" begin
    outdir = mktempdir()
    out = read(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $failure_worker nsplit_opttrans $outdir`,
               String)
    @test count("failure-mode worker: nsplit_opttrans expected rejection ok", out) == 2
    @test isempty(readdir(outdir))
end

@testset "v0.5 failure mode: NSplitSize > 1 with FSZ projection rejects before MPI.Init" begin
    outdir = mktempdir()
    out = read(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project=$project $failure_worker nsplit_fsz_projection $outdir`,
               String)
    @test count("failure-mode worker: nsplit_fsz_projection expected rejection ok", out) == 2
    @test isempty(readdir(outdir))
end

@testset "v0.5 NSplitSize/NStore direct-SR self-consistency" begin
    cases = (
        (fixture = "hubbard_chain_real", mode = "real"),
        (fixture = "heisenberg_chain_fsz", mode = "fsz"),
    )
    for case in cases
        direct_ref = run_nsplit_nstore_case(case.fixture, case.mode, 1, 0, 2)
        store_ref = run_nsplit_nstore_case(case.fixture, case.mode, 1, 1, 2)
        direct_split = run_nsplit_nstore_case(case.fixture, case.mode, 2, 0, 4)
        store_split = run_nsplit_nstore_case(case.fixture, case.mode, 2, 1, 4)

        for (label, result) in (
            ("store_ref", store_ref),
            ("direct_split", direct_split),
            ("store_split", store_split),
        )
            prefix = "$(case.fixture) $label"
            assert_close_vector("$prefix zvo_out", result.zvo, direct_ref.zvo;
                                atol = nsplit_nstore_tol)
            assert_close_vector("$prefix zqp_opt", result.zqp, direct_ref.zqp;
                                atol = nsplit_nstore_tol)
        end
    end
end

@testset "v0.5 NSplitSize standard-projection NQPFull self-consistency" begin
    cases = (
        (fixture = "heisenberg_chain_real", mode = "real"),
        (fixture = "hubbard_tetragonal_momentum_projection_real", mode = "real"),
        (fixture = "heisenberg_chain_cmp", mode = "cmp"),
    )
    for case in cases
        nsp = haskey(case, :nsp_gauss_leg) ? case.nsp_gauss_leg : nothing
        direct_ref = run_nsplit_standard_projection_case(
            case.fixture,
            case.mode,
            1,
            0,
            2;
            nsp_gauss_leg = nsp,
        )
        store_ref = run_nsplit_standard_projection_case(
            case.fixture,
            case.mode,
            1,
            1,
            2;
            nsp_gauss_leg = nsp,
        )
        direct_split = run_nsplit_standard_projection_case(
            case.fixture,
            case.mode,
            2,
            0,
            4;
            nsp_gauss_leg = nsp,
        )
        store_split = run_nsplit_standard_projection_case(
            case.fixture,
            case.mode,
            2,
            1,
            4;
            nsp_gauss_leg = nsp,
        )

        for (label, result) in (
            ("store_ref", store_ref),
            ("direct_split", direct_split),
            ("store_split", store_split),
        )
            prefix = "$(case.fixture) $label"
            assert_close_vector("$prefix zvo_out", result.zvo, direct_ref.zvo;
                                atol = nsplit_standard_projection_tol)
            assert_close_vector("$prefix zqp_opt", result.zqp, direct_ref.zqp;
                                atol = nsplit_standard_projection_tol)
            assert_close_vector("$prefix zvo_var", result.var, direct_ref.var;
                                atol = nsplit_standard_projection_tol)
        end
    end
end

@testset "v0.5 PhysCal NSplitSize self-consistency" begin
    cases = (
        (fixture = "heisenberg_chain_real", mode = "real"),
        (fixture = "heisenberg_chain_cmp", mode = "cmp"),
    )
    for case in cases
        direct_ref = run_physcal_nsplit_case(case.fixture, case.mode, 1, 2)
        split = run_physcal_nsplit_case(case.fixture, case.mode, 2, 4)

        for field in (:zvo, :var, :cisajs, :two_body, :factored)
            assert_close_vector(
                "$(case.fixture) PhysCal $field",
                getfield(split, field),
                getfield(direct_ref, field);
                atol = physcal_nsplit_tol,
            )
        end
    end
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
