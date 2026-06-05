using Test
using MVMCOptimizers

# Unit-level checks for run_phys_cal_from_namelist argument handling. The
# end-to-end behavior (corrected init order: no double init_parameter! /
# init_qp_weight!, RNG aligned with C, Green output) is validated numerically by
# the PhysCal e2e gate in test/integration (Plan 3b), which needs committed
# fixtures and is too heavy for the subpackage unit suite.

@testset "run_phys_cal_from_namelist: argument validation" begin
    # mode is validated before any file access, so a bad mode errors even with a
    # nonexistent namelist / opt_para.
    err = try
        run_phys_cal_from_namelist("nonexistent_namelist.def"; opt_para = "x.dat", mode = :bogus)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("mode must be", sprint(showerror, err))

    # opt_para is a required keyword (the gate must be explicit about which fixed
    # parameters it consumes) — omitting it throws UndefKeywordError.
    err2 = try
        run_phys_cal_from_namelist("nonexistent_namelist.def"; mode = :real)
        nothing
    catch e
        e
    end
    @test err2 isa UndefKeywordError
end
