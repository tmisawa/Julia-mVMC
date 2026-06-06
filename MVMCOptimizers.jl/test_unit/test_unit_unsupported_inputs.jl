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
    # NSplitSize = 1 (single process) remains accepted: the validator is a
    # no-op returning `nothing`.
    @testset "NSplitSize = 1 is accepted" begin
        modpara = ModParaParameters(nsplit_size = 1)
        @test MVMCOptimizers.validate_supported_modpara(modpara) === nothing
    end

    # NSplitSize > 1 is rejected. The message must name both the offending
    # input and the missing capability so users know why and what to do.
    @testset "NSplitSize > 1 is rejected with a clear message" begin
        modpara = ModParaParameters(nsplit_size = 2)
        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_modpara(modpara),
        )
        @test threw
        @test occursin("NSplitSize > 1", msg)
        @test occursin("MPI parallelization is not supported", msg)
    end

    # NSplitSize < 1 (0 or negative) is an invalid value: a process-split
    # count must be at least 1. Rejected distinctly from the unsupported-MPI
    # case (it is an invalid value, not a missing feature).
    @testset "NSplitSize <= 0 is rejected as an invalid value" begin
        for bad in (0, -1)
            modpara = ModParaParameters(nsplit_size = bad)
            threw, msg = capture_error_message(
                () -> MVMCOptimizers.validate_supported_modpara(modpara),
            )
            @test threw
            @test occursin(">= 1", msg)
            # It must not be mislabeled as the MPI-unsupported case.
            @test !occursin("MPI parallelization is not supported", msg)
        end
    end

    # The contract is enforced at the independent runtime entry points, before
    # any optimization or physical calculation proceeds. A default
    # ExpertModeData is enough because the guard is the first statement in each
    # entry point, so it fails fast regardless of the rest of the input.
    @testset "entry points enforce the contract" begin
        for entry in (vmc_para_opt!, vmc_phys_cal!)
            data = ExpertModeData()
            data.modpara.nsplit_size = 3
            threw, msg = capture_error_message(() -> entry(data))
            @test threw
            @test occursin("NSplitSize > 1", msg)
            @test occursin("MPI parallelization is not supported", msg)
        end
    end
end

@testset "unit/unsupported_inputs: DoublonHolon complex flag participates in all-complex" begin
    data = ExpertModeData()
    data.doublon_holon_2site_complex = true
    @test MVMCOptimizers.get_all_complex_flag(data)
end

@testset "unit/unsupported_inputs: NLanczosMode contract" begin
    # NLanczosMode = 0 (variational only) is the supported setting: no-op.
    @testset "NLanczosMode = 0 is accepted" begin
        modpara = ModParaParameters(lanczos_mode = 0)
        @test MVMCOptimizers.validate_supported_modpara(modpara) === nothing
    end

    # NLanczosMode > 0 (full Lanczos) is rejected: only step-0 matches C and the
    # C indirect one-body Green list (NLanczosMode > 1) is not reproduced.
    @testset "NLanczosMode > 0 is rejected with a clear message" begin
        for bad in (1, 2)
            modpara = ModParaParameters(lanczos_mode = bad)
            threw, msg = capture_error_message(
                () -> MVMCOptimizers.validate_supported_modpara(modpara),
            )
            @test threw
            @test occursin("NLanczosMode > 0", msg)
            # Must not be mislabeled as the MPI-unsupported case.
            @test !occursin("MPI parallelization is not supported", msg)
        end
    end

    # Enforced at the independent runtime entry points before any work.
    @testset "entry points enforce the contract" begin
        for entry in (vmc_para_opt!, vmc_phys_cal!)
            data = ExpertModeData()
            data.modpara.lanczos_mode = 2
            threw, msg = capture_error_message(() -> entry(data))
            @test threw
            @test occursin("NLanczosMode > 0", msg)
        end
    end
end
