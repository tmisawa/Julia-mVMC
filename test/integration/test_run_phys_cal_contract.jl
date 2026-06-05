using Test
using MVMCOptimizers

# Cheap contract test for run_phys_cal_from_namelist against a committed fixture:
# the namelist must PARSE, then the strict fixed-parameter load must FAIL *before*
# any sampling. This pins the corrected init order's early stages (parse → load)
# without running a full PhysCal sweep — the numeric PhysCal-vs-C comparison is the
# Plan 3b e2e gate's job. A loader error here proves sampling was never reached
# (read_opt_para_file! is step 3; vmc_phys_cal! / mkpath is step 6).

@testset "run_phys_cal_from_namelist: errors at the loader, before sampling" begin
    namelist = joinpath(
        @__DIR__, "reference", "heisenberg_chain_real", "inputs", "namelist.def",
    )
    @test isfile(namelist)

    mktempdir() do dir
        # (a) Missing fixed-parameter file: parse succeeds, loader errors, and no
        #     output directory is created (we never reach mkpath/sampling).
        out_a = joinpath(dir, "out_a")
        err = try
            run_phys_cal_from_namelist(
                namelist; opt_para = joinpath(dir, "nope.dat"), mode = :real,
                output_dir = out_a,
            )
            nothing
        catch e
            e
        end
        @test err !== nothing
        msg = sprint(showerror, err)
        @test occursin("read_opt_para_file!", msg)
        @test occursin("file not found", msg)
        @test !isdir(out_a)

        # (b) Malformed (too-short) fixed-parameter file: parse succeeds, loader
        #     errors before sampling.
        bad = joinpath(dir, "bad_zqp_opt.dat")
        write(bad, "1.0 2.0 3.0\n")  # far fewer than 6 + 3*NPara floats
        out_b = joinpath(dir, "out_b")
        err2 = try
            run_phys_cal_from_namelist(
                namelist; opt_para = bad, mode = :real, output_dir = out_b,
            )
            nothing
        catch e
            e
        end
        @test err2 !== nothing
        msg2 = sprint(showerror, err2)
        @test occursin("read_opt_para_file!", msg2)
        @test occursin("too short", msg2)
        @test !isdir(out_b)
    end
end
