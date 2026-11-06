using Test
using MVMCOptimizers
using Printf

# Focused regression tests for the C-compatible PhysCal section timer.
# `run_phys_cal_from_namelist` only writes the report when `MVMC_C_TIMER` is
# set; the heavier end-to-end timer/numerical invariance check lives in
# test/integration/test_run_phys_cal_contract.jl.

# Independent copy of C's OutputTimerPhysCal label/id layout
# (mVMC 1.4.0 src/mVMC/vmcclock.c `OutputTimerPhysCal`, line 140). If the writer reuses the ParaOpt
# labels here, or drops/reorders a line, this list will not match.
const C_PHYS_CAL_TIMER_LINES = [
    "All                         [0] ",
    "Initialization              [1] ",
    "  read options             [10] ",
    "  ReadDefFile              [11] ",
    "  SetMemory                [12] ",
    "  InitParameter            [13] ",
    "VMCPhysCal                  [2] ",
    "  VMCMakeSample             [3] ",
    "    makeInitialSample      [30] ",
    "    make candidate         [31] ",
    "    hopping update         [32] ",
    "      UpdateProjCnt        [60] ",
    "      CalculateNewPfM2     [61] ",
    "      CalculateLogIP       [62] ",
    "      UpdateMAll           [63] ",
    "    exchange update        [33] ",
    "      UpdateProjCnt        [65] ",
    "      CalculateNewPfMTwo2  [66] ",
    "      CalculateLogIP       [67] ",
    "      UpdateMAllTwo        [68] ",
    "    lspinflip update       [36] ",
    "      UpdateProjCnt       [600] ",
    "      CalculateNewPfMTwo2 [601] ",
    "      CalculateLogIP      [602] ",
    "      UpdateMAllTwo       [603] ",
    "    recal PfM and InvM     [34] ",
    "    save electron config   [35] ",
    "  VMCMainCal                [4] ",
    "    CalculateMAll          [40] ",
    "    LocEnergyCal           [41] ",
    "      CalHamiltonian0      [70] ",
    "      CalHamiltonian1      [71] ",
    "      CalHamiltonian2      [72] ",
    "    CalculateGreenFunc     [42] ",
    "      GreenFunc1           [50] ",
    "      GreenFunc2           [51] ",
    "      addPhysCA            [52] ",
    "      addPhysCACA          [53] ",
    "      addPhysTwist         [54] ",
    "    Lanczos1               [43] ",
    "    Lanczos2               [44] ",
    "  UpdateSlaterElm          [20] ",
    "  WeightAverage            [21] ",
    "  outputData               [22] ",
]

@testset "write_ctimer_phys_cal: C OutputTimerPhysCal layout" begin
    @test length(MVMCOptimizers.CTIMER_PHYS_CAL_LINES) == length(C_PHYS_CAL_TIMER_LINES)
    @test [label for (label, _) in MVMCOptimizers.CTIMER_PHYS_CAL_LINES] ==
          C_PHYS_CAL_TIMER_LINES

    timer = MVMCOptimizers.CTimer(true)
    timer.elapsed_ns[0 + 1] = 1_500_000_000            # id 0 -> 1.50000
    timer.elapsed_ns[3 + 1] = 2_250_000_000            # id 3 -> 2.25000
    timer.elapsed_ns[22 + 1] = 123_456                  # id 22 -> 0.00012

    mktempdir() do dir
        path = MVMCOptimizers.write_ctimer_phys_cal(timer, dir)
        @test basename(path) == "zvo_CalcTimer.dat"
        lines = readlines(path)
        @test length(lines) == length(C_PHYS_CAL_TIMER_LINES)

        expected_values = Dict(0 => 1.5, 3 => 2.25, 22 => 0.000123456)
        for (line, (label, id)) in zip(lines, MVMCOptimizers.CTIMER_PHYS_CAL_LINES)
            expected = label * @sprintf("%12.5f", get(expected_values, id, 0.0))
            @test line == expected
        end
    end
end
