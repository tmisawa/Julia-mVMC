using Test
using MVMCOptimizers
using MVMCExpertModeParsers: ExpertModeData, InterAllTerm, ModParaParameters, TransferTerm

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

    @testset "NSplitSize > 1 is accepted for sz-conserved PhysCal normal Green" begin
        modpara = ModParaParameters(nsplit_size = 2)
        @test MVMCOptimizers.validate_supported_modpara(modpara) === nothing
        @test MVMCOptimizers.validate_supported_phys_cal_modpara(modpara) === nothing

        data = ExpertModeData()
        data.modpara.nsplit_size = 2
        data.modpara.lanczos_mode = 0
        data.i_flg_orbital_general = 0
        @test MVMCOptimizers.validate_supported_phys_cal_data(data) === nothing
    end

    @testset "NSplitSize > 1 accepts standard-projection NQPFull > 1" begin
        spin_data = ExpertModeData()
        spin_data.modpara.nsplit_size = 2
        spin_data.modpara.nsp_gauss_leg = 2
        spin_data.modpara.nmp_trans = 1
        spin_data.n_qp_opt_trans = 1
        spin_data.i_flg_orbital_general = 0
        @test MVMCOptimizers.validate_supported_para_opt_data(spin_data) === nothing

        momentum_data = ExpertModeData()
        momentum_data.modpara.nsplit_size = 2
        momentum_data.modpara.nsp_gauss_leg = 1
        momentum_data.modpara.nmp_trans = 4
        momentum_data.n_qp_opt_trans = 1
        momentum_data.i_flg_orbital_general = 0
        @test MVMCOptimizers.validate_supported_para_opt_data(momentum_data) === nothing

        trivial_opttrans_data = ExpertModeData()
        trivial_opttrans_data.modpara.nsplit_size = 2
        trivial_opttrans_data.modpara.nsp_gauss_leg = 8
        trivial_opttrans_data.modpara.nmp_trans = -1
        trivial_opttrans_data.n_qp_opt_trans = 1
        trivial_opttrans_data.i_flg_orbital_general = 0
        trivial_opttrans_data.opt_trans = ComplexF64[1.0 + 0.0im]
        trivial_opttrans_data.qp_opt_trans = [[0]]
        @test MVMCOptimizers.validate_supported_para_opt_data(trivial_opttrans_data) ===
              nothing
    end

    @testset "NSplitSize > 1 rejects OptTrans-derived NQPFull > 1" begin
        data = ExpertModeData()
        data.modpara.nsplit_size = 2
        data.modpara.nsp_gauss_leg = 1
        data.modpara.nmp_trans = 1
        data.n_qp_opt_trans = 2
        data.opt_trans = ComplexF64[1.0 + 0.0im, 0.5 + 0.0im]
        data.qp_opt_trans = [[0], [0]]

        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_para_opt_data(data),
        )
        @test threw
        @test occursin("NSplitSize > 1 with NQPOptTrans > 1", msg)
        @test occursin("OptTrans", msg)
    end

    @testset "NSplitSize > 1 rejects OptTrans-derived PhysCal split" begin
        trivial = ExpertModeData()
        trivial.modpara.nsplit_size = 2
        trivial.modpara.lanczos_mode = 0
        trivial.i_flg_orbital_general = 0
        trivial.n_qp_opt_trans = 1
        trivial.opt_trans = ComplexF64[1.0 + 0.0im]
        trivial.qp_opt_trans = [[0]]
        @test MVMCOptimizers.validate_supported_phys_cal_data(trivial) === nothing

        data = ExpertModeData()
        data.modpara.nsplit_size = 2
        data.modpara.lanczos_mode = 0
        data.i_flg_orbital_general = 0
        data.n_qp_opt_trans = 2
        data.opt_trans = ComplexF64[1.0 + 0.0im, 0.5 + 0.0im]
        data.qp_opt_trans = [[0], [0]]

        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_phys_cal_data(data),
        )
        @test threw
        @test occursin("NSplitSize > 1 with NQPOptTrans > 1", msg)
        @test occursin("PhysCal", msg)
        @test occursin("OptTrans", msg)
    end

    @testset "NSplitSize > 1 rejects unsupported PhysCal split scopes" begin
        fsz_data = ExpertModeData()
        fsz_data.modpara.nsplit_size = 2
        fsz_data.modpara.lanczos_mode = 0
        fsz_data.i_flg_orbital_general = 1
        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_phys_cal_data(fsz_data),
        )
        @test threw
        @test occursin("NSplitSize > 1", msg)
        @test occursin("FSZ / general-orbital", msg)
        @test occursin("PhysCal", msg)

        lanczos_data = ExpertModeData()
        lanczos_data.modpara.nsplit_size = 2
        lanczos_data.modpara.lanczos_mode = 2
        lanczos_data.i_flg_orbital_general = 0
        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_phys_cal_data(lanczos_data),
        )
        @test threw
        @test occursin("NSplitSize > 1 with NLanczosMode > 0", msg)
        @test occursin("PhysCal", msg)
    end

    @testset "NSplitSize > 1 rejects FSZ standard-projection NQPFull > 1" begin
        for (nsp, nmp) in ((2, 1), (1, 2), (1, -2))
            data = ExpertModeData()
            data.modpara.nsplit_size = 2
            data.modpara.nsp_gauss_leg = nsp
            data.modpara.nmp_trans = nmp
            data.n_qp_opt_trans = 1
            data.i_flg_orbital_general = 1

            threw, msg = capture_error_message(
                () -> MVMCOptimizers.validate_supported_para_opt_data(data),
            )
            @test threw
            @test occursin(
                "NSplitSize > 1 with FSZ standard-projection NQPFull > 1",
                msg,
            )
            @test occursin("NSPGaussLeg = $nsp", msg)
            @test occursin("NMPTrans = $nmp", msg)
        end
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
        phys_data.i_flg_orbital_general = 1
        threw, msg = capture_error_message(() -> vmc_phys_cal!(phys_data))
        @test threw
        @test occursin("NSplitSize > 1", msg)
        @test occursin("FSZ / general-orbital", msg)
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

    @testset "NLanczosMode = 1/2 are accepted for sz-conserved PhysCal" begin
        for mode in (1, 2)
            modpara = ModParaParameters(lanczos_mode = mode)
            @test MVMCOptimizers.validate_supported_modpara(modpara) === nothing
            @test MVMCOptimizers.validate_supported_phys_cal_modpara(modpara) === nothing
            data = ExpertModeData()
            data.modpara.lanczos_mode = mode
            @test MVMCOptimizers.validate_supported_phys_cal_data(data) === nothing
        end
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

    @testset "FSZ/general-orbital Lanczos is rejected for PhysCal" begin
        for mode in (1, 2)
            data = ExpertModeData()
            data.modpara.lanczos_mode = mode
            data.i_flg_orbital_general = 1
            threw, msg = capture_error_message(
                () -> MVMCOptimizers.validate_supported_phys_cal_data(data),
            )
            @test threw
            @test occursin("FSZ / general-orbital", msg)
        end
    end

    @testset "spin-changing Lanczos rejects use mode-independent wording" begin
        transfer_data = ExpertModeData()
        transfer_data.modpara.lanczos_mode = 2
        transfer_data.transfer_terms = [TransferTerm(0, 0, 1, 1, 1.0 + 0.0im)]
        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_phys_cal_data(transfer_data),
        )
        @test threw
        @test occursin("spin-flip Transfer", msg)
        @test !occursin("R1", msg)

        interall_data = ExpertModeData()
        interall_data.modpara.lanczos_mode = 2
        interall_data.inter_all_terms = [InterAllTerm(0, 0, 1, 1, 2, 0, 3, 0, 1.0 + 0.0im, false)]
        threw, msg = capture_error_message(
            () -> MVMCOptimizers.validate_supported_phys_cal_data(interall_data),
        )
        @test threw
        @test occursin("spin-changing InterAll", msg)
        @test !occursin("R1", msg)
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
        phys_data.i_flg_orbital_general = 1
        threw, msg = capture_error_message(() -> vmc_phys_cal!(phys_data))
        @test threw
        @test occursin("FSZ / general-orbital", msg)
    end
end
