using Test
using MVMCOptimizers
using MVMCExpertModeParsers: ExpertModeData, ModParaParameters

# Capture the showerror text of whatever `f()` throws.
# Returns (threw::Bool, message::String). Used so the contract asserts on the
# error *message* (the stable contract) rather than the concrete exception
# type (see design review finding A2).
function capture_error_message(f)
    try
        f()
        return (false, "")
    catch err
        return (true, sprint(showerror, err))
    end
end

@testset "unit/unsupported_inputs: NSplitSize contract" begin
    # NSplitSize = 1 remains accepted for both serial and MPI sample-parallel
    # runs: the global validator is a no-op returning `nothing`.
    @testset "NSplitSize = 1 is accepted" begin
        modpara = ModParaParameters(nsplit_size = 1)
        @test MVMCOptimizers.validate_supported_modpara(modpara) === nothing
    end

    @testset "NSplitSize > 1 is accepted for direct para-opt" begin
        modpara = ModParaParameters(nsplit_size = 2, nsrcg = 0)
        @test MVMCOptimizers.validate_supported_modpara(modpara) === nothing
        @test MVMCOptimizers.validate_supported_para_opt_modpara(modpara) === nothing
        @test MVMCOptimizers.validate_supported_para_opt_parallel_modpara(
            serial_context(),
            modpara,
        ) === nothing
    end

    @testset "NSplitSize > 1 with SR-CG is rejected for para-opt" begin
        modpara = ModParaParameters(nsplit_size = 2, nsrcg = 1)
        @test MVMCOptimizers.validate_supported_modpara(modpara) === nothing
        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_para_opt_modpara(modpara),
        )
        @test threw
        @test occursin("NSplitSize > 1 with SR-CG", msg)
        @test occursin("NSRCG = 1", msg)
    end

    @testset "NSplitSize > 1 is rejected for PhysCal" begin
        modpara = ModParaParameters(nsplit_size = 2)
        @test MVMCOptimizers.validate_supported_modpara(modpara) === nothing
        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_phys_cal_modpara(modpara),
        )
        @test threw
        @test occursin("NSplitSize > 1", msg)
        @test occursin("PhysCal", msg)
    end

    @testset "NSplitSize > 1 with NQPFull > 1 is rejected for para-opt" begin
        data = ExpertModeData()
        data.modpara.nsplit_size = 2
        data.modpara.nmp_trans = 2
        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_para_opt_data(data),
        )
        @test threw
        @test occursin("NSplitSize > 1 with NQPFull > 1", msg)
        @test occursin("NQPFull = 2", msg)
    end

    # NSplitSize < 1 (0 or negative) is an invalid value: a process-split
    # count must be at least 1. Rejected distinctly from the grouped-MPI
    # case (it is an invalid value, not a missing feature).
    @testset "NSplitSize <= 0 is rejected as an invalid value" begin
        for bad in (0, -1)
            modpara = ModParaParameters(nsplit_size = bad)
            threw, msg = capture_error_message(
                () -> MVMCOptimizers.validate_supported_modpara(modpara),
            )
            @test threw
            @test occursin(">= 1", msg)
            # It must not be mislabeled as an unsupported feature combination.
            @test !occursin("NSplitSize > 1 with", msg)
        end
    end

    # Unsupported combinations are enforced at the independent runtime entry
    # points before any optimization or physical calculation proceeds. A default
    # ExpertModeData is enough because the guards run before input-dependent
    # state construction.
    @testset "entry points enforce the contract" begin
        para_data = ExpertModeData()
        para_data.modpara.nsplit_size = 3
        para_data.modpara.nsrcg = 1
        threw, msg = capture_error_message(() -> vmc_para_opt!(para_data))
        @test threw
        @test occursin("NSplitSize > 1 with SR-CG", msg)

        phys_data = ExpertModeData()
        phys_data.modpara.nsplit_size = 3
        threw, msg = capture_error_message(() -> vmc_phys_cal!(phys_data))
        @test threw
        @test occursin("NSplitSize > 1", msg)
        @test occursin("PhysCal", msg)
    end
end

@testset "unit/unsupported_inputs: SR-CG option contract" begin
    @testset "standard direct and CG solvers are accepted" begin
        @test MVMCOptimizers.validate_supported_modpara(
            ModParaParameters(nsrcg = 0),
        ) === nothing
        @test MVMCOptimizers.validate_supported_modpara(
            ModParaParameters(nsrcg = 1),
        ) === nothing
    end

    @testset "NSRCG >= 2 is rejected with a clear message" begin
        modpara = ModParaParameters(nsrcg = 2)
        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_modpara(modpara),
        )
        @test threw
        @test occursin("NSRCG >= 2", msg)
        @test occursin("standard SR-CG solver", msg)
    end

    @testset "unported CG submodes are rejected with clear messages" begin
        for (modpara, label) in (
            (ModParaParameters(use_diag_scale = 1), "useDiagScale != 0"),
            (ModParaParameters(rescale_smat = 1), "RescaleSmat != 0"),
        )
            threw, msg = capture_error_message(
                () -> MVMCOptimizers.validate_supported_modpara(modpara),
            )
            @test threw
            @test occursin(label, msg)
            @test occursin("not supported", msg)
        end
    end

    @testset "entry points enforce the CG option contract" begin
        data = ExpertModeData()
        data.modpara.nsrcg = 2
        threw, msg = capture_error_message(() -> vmc_para_opt!(data))
        @test threw
        @test occursin("NSRCG >= 2", msg)
    end
end

@testset "unit/unsupported_inputs: DoublonHolon complex flag participates in all-complex" begin
    data = ExpertModeData()
    data.doublon_holon_2site_complex = true
    @test MVMCOptimizers.get_all_complex_flag(data)
end

@testset "unit/unsupported_inputs: NLanczosMode contract" begin
    @testset "NLanczosMode = 0/1/2 are globally valid values" begin
        for mode in (0, 1, 2)
            modpara = ModParaParameters(lanczos_mode = mode)
            @test MVMCOptimizers.validate_supported_modpara(modpara) === nothing
        end
    end

    @testset "unknown NLanczosMode values are rejected globally" begin
        for bad in (-1, 3)
            modpara = ModParaParameters(lanczos_mode = bad)
            threw, msg = capture_error_message(
                () -> MVMCOptimizers.validate_supported_modpara(modpara),
            )
            @test threw
            @test occursin("NLanczosMode must be 0, 1, or 2", msg)
        end
    end

    @testset "NLanczosMode = 1 is accepted for PhysCal R1" begin
        modpara = ModParaParameters(lanczos_mode = 1)
        @test MVMCOptimizers.validate_supported_modpara(modpara) === nothing
        @test MVMCOptimizers.validate_supported_phys_cal_modpara(modpara) === nothing
        data = ExpertModeData()
        data.modpara.lanczos_mode = 1
        @test MVMCOptimizers.validate_supported_phys_cal_data(data) === nothing
    end

    @testset "NLanczosMode > 0 is rejected for ParaOpt" begin
        for bad in (1, 2)
            modpara = ModParaParameters(lanczos_mode = bad)
            threw, msg = capture_error_message(
                () -> MVMCOptimizers.validate_supported_para_opt_modpara(modpara),
            )
            @test threw
            @test occursin("NLanczosMode > 0", msg)
            @test occursin("parameter optimization", msg)
        end
    end

    @testset "NLanczosMode > 1 is rejected for PhysCal R1" begin
        modpara = ModParaParameters(lanczos_mode = 2)
        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_phys_cal_modpara(modpara),
        )
        @test threw
        @test occursin("NLanczosMode > 1", msg)
        @test occursin("PhysCal", msg)
    end

    @testset "FSZ/general-orbital Lanczos is rejected for PhysCal R1" begin
        data = ExpertModeData()
        data.modpara.lanczos_mode = 1
        data.i_flg_orbital_general = 1
        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_phys_cal_data(data),
        )
        @test threw
        @test occursin("FSZ / general-orbital", msg)
    end

    @testset "entry points enforce the Lanczos contract" begin
        para_data = ExpertModeData()
        para_data.modpara.lanczos_mode = 1
        threw, msg = capture_error_message(() -> vmc_para_opt!(para_data))
        @test threw
        @test occursin("NLanczosMode > 0", msg)
        @test occursin("parameter optimization", msg)

        phys_data = ExpertModeData()
        phys_data.modpara.lanczos_mode = 2
        threw, msg = capture_error_message(() -> vmc_phys_cal!(phys_data))
        @test threw
        @test occursin("NLanczosMode > 1", msg)
    end
end
