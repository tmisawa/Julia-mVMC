"""
C-compatible lightweight timer.

Mirrors C-mVMC's `StartTimer`/`StopTimer`/`OutputTimerParaOpt` (see
`mVMC/src/mVMC/vmcclock.c`) so that Julia runs can emit a `zvo_CalcTimer.dat`
in the same id/label/seconds format, enabling a per-section comparison with C.

Design notes (see docs/plans/2026-05-25-julia-mvmc-c-compatible-timer-plan.md):

- `time_ns()` returns `UInt64`, so the accumulators are `UInt64` too (avoids the
  signed/unsigned mixing and implicit `convert` that `Int64` storage would
  introduce). Seconds conversion happens once, at output time.
- The struct is parameterised on the enabled flag's type `E`
  (`Val{true}`/`Val{false}`). When disabled, `ctimer_start!`/`ctimer_stop!`
  dispatch to no-ops that the compiler inlines away — no array access, no
  branch in the hot path. The `Val(enabled)` -> type conversion is a deliberate
  one-time dynamic dispatch performed at the run entry, followed by a function
  barrier so the optimisation loop specialises on the concrete timer type.
- Timers are inclusive (a parent keeps running while a child runs), matching C.
  Instrument at *call sites*, not inside shared functions: e.g.
  `calculate_m_all_*` is timed as id 30/34 from sampling and id 40 from the main
  calculation, so wrapping its body would conflate them. A given id must not be
  re-entered before it is stopped (single `start_ns[id]` slot, same as C).
- Single shared timer is safe because the sample loops are single-threaded;
  threading is confined to the QP loop inside `calculate_m_all` (joined before
  the surrounding timer stops). If sample-level threading is ever added, move to
  thread-local accumulators.
"""

# C's NTimer (global.h). Sized to cover the maximum id in use, including the
# fsz lspinflip ids 600..603 (so length >= 604). Match C's value exactly.
const CTIMER_N = 1000

struct CTimer{E}
    enabled::E                      # Val{true} or Val{false}
    elapsed_ns::Vector{UInt64}      # per-id accumulated time [ns]; index = id + 1
    start_ns::Vector{UInt64}        # per-id start timestamp [ns]; index = id + 1
end

"""
    CTimer(enabled::Bool) -> CTimer

Construct a timer. `enabled == false` yields a `CTimer{Val{false}}` whose
start/stop calls compile to no-ops. The `Val(enabled)` conversion here is the
intended one-time dynamic dispatch; pass the result through a function barrier.
"""
CTimer(enabled::Bool) =
    CTimer(Val(enabled), zeros(UInt64, CTIMER_N), zeros(UInt64, CTIMER_N))

# Disabled-by-default singleton, used when a caller passes `c_timer = nothing`.
const CTIMER_DISABLED = CTimer(false)

# --- start / stop -----------------------------------------------------------
# Disabled: no-ops (no array access, no branch).
@inline ctimer_start!(::Val{false}, ::CTimer, ::Integer) = nothing
@inline ctimer_stop!(::Val{false}, ::CTimer, ::Integer) = nothing

@inline function ctimer_start!(::Val{true}, timer::CTimer, id::Integer)
    @inbounds timer.start_ns[id + 1] = time_ns()
    return nothing
end

@inline function ctimer_stop!(::Val{true}, timer::CTimer, id::Integer)
    @inbounds timer.elapsed_ns[id + 1] += time_ns() - timer.start_ns[id + 1]
    return nothing
end

# Public entry points: dispatch on the timer's enabled type parameter, resolved
# at compile time for a concretely-typed `CTimer{E}`.
@inline ctimer_start!(timer::CTimer, id::Integer) = ctimer_start!(timer.enabled, timer, id)
@inline ctimer_stop!(timer::CTimer, id::Integer) = ctimer_stop!(timer.enabled, timer, id)

"""
    ctimer_reset!(timer::CTimer)

Zero all accumulators. Call once per run when reusing a timer across repeats so
elapsed times do not accumulate across runs (C runs one process per run).
"""
function ctimer_reset!(timer::CTimer)
    fill!(timer.elapsed_ns, 0)
    fill!(timer.start_ns, 0)
    return timer
end

ctimer_enabled(::CTimer{Val{true}}) = true
ctimer_enabled(::CTimer{Val{false}}) = false

@inline ctimer_seconds(timer::CTimer, id::Integer) =
    @inbounds timer.elapsed_ns[id + 1] / 1.0e9

# --- output ------------------------------------------------------------------
# Fixed label/id prefixes copied verbatim from C's OutputTimerParaOpt
# (mVMC/src/mVMC/vmcclock.c). Each prefix ends with a trailing space so that
# appending an `%12.5f` value reproduces C's `"...[id] %12.5lf\n"` layout
# exactly, letting the existing C zvo_CalcTimer.dat parser read Julia output.
const CTIMER_PARA_OPT_LINES = Tuple{String,Int}[
    ("All                         [0] ", 0),
    ("Initialization              [1] ", 1),
    ("  read options             [10] ", 10),
    ("  ReadDefFile              [11] ", 11),
    ("  SetMemory                [12] ", 12),
    ("  InitParameter            [13] ", 13),
    ("VMCParaOpt                  [2] ", 2),
    ("  VMCMakeSample             [3] ", 3),
    ("    makeInitialSample      [30] ", 30),
    ("    make candidate         [31] ", 31),
    ("    hopping update         [32] ", 32),
    ("      UpdateProjCnt        [60] ", 60),
    ("      CalculateNewPfM2     [61] ", 61),
    ("      CalculateLogIP       [62] ", 62),
    ("      UpdateMAll           [63] ", 63),
    ("    exchange update        [33] ", 33),
    ("      UpdateProjCnt        [65] ", 65),
    ("      CalculateNewPfMTwo2  [66] ", 66),
    ("      CalculateLogIP       [67] ", 67),
    ("      UpdateMAllTwo        [68] ", 68),
    ("    lspinflip update       [36] ", 36),
    ("      UpdateProjCnt       [600] ", 600),
    ("      CalculateNewPfMTwo2 [601] ", 601),
    ("      CalculateLogIP      [602] ", 602),
    ("      UpdateMAllTwo       [603] ", 603),
    ("    recal PfM and InvM     [34] ", 34),
    ("    save electron config   [35] ", 35),
    ("  VMCMainCal                [4] ", 4),
    ("    CalculateMAll          [40] ", 40),
    ("    LocEnergyCal           [41] ", 41),
    ("      CalHamiltonian0      [70] ", 70),
    ("      CalHamiltonian1      [71] ", 71),
    ("      CalHamiltonian2      [72] ", 72),
    ("    ReturnSlaterElmDiff    [42] ", 42),
    ("    calculate OO and HO    [43] ", 43),
    ("    multiply store OO      [45] ", 45),
    ("  StochasticOpt             [5] ", 5),
    ("    preprocess             [50] ", 50),
    ("    stcOptMain             [51] ", 51),
    ("      initBLACS            [55] ", 55),
    ("      calculate S and g    [56] ", 56),
    ("      DPOSV                [57] ", 57),
    ("      gatherParaChange     [58] ", 58),
    ("    postprocess            [52] ", 52),
    ("  UpdateSlaterElm          [20] ", 20),
    ("  WeightAverage            [21] ", 21),
    ("  outputData               [22] ", 22),
    ("  SyncModifiedParameter    [23] ", 23),
    ("  cal                      [24] ", 24),
    ("  SR                       [25] ", 25),
    ("  MAll                     [69] ", 69),
]

"""
    write_ctimer_para_opt(timer::CTimer, output_dir; prefix="zvo")

Write `<prefix>_CalcTimer.dat` under `output_dir` in C's `OutputTimerParaOpt`
format. ids that were never instrumented stay at 0.0, exactly as in C where
unused `Timer[i]` entries print `0.00000`. Single-process: no rank guard (C only
lets rank 0 write; here there is only one process).
"""
function write_ctimer_para_opt(timer::CTimer, output_dir::AbstractString; prefix::AbstractString = "zvo")
    path = joinpath(output_dir, string(prefix, "_CalcTimer.dat"))
    open(path, "w") do fp
        for (label, id) in CTIMER_PARA_OPT_LINES
            print(fp, label, @sprintf("%12.5f\n", ctimer_seconds(timer, id)))
        end
    end
    return path
end
