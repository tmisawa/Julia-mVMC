"""
    run_para_opt_from_namelist(namelist_path;
                               nsteps, mode, nsmp=nothing,
                               output_dir, seed=nothing,
                               initial_def=:auto)

Drive `vmc_para_opt!` for `nsteps` SR steps using the C-mVMC-format
`namelist.def` at `namelist_path`. Outputs are written under `output_dir`
(default: a fresh `tempname()`). Returns a `NamedTuple` summarising the run.

This is the entry point used by `examples/*.jl` and the workspace-level
integration tests in `test/integration/runtests.jl`.

# Arguments

- `namelist_path::AbstractString`: path to `namelist.def`. All other `.def`
  files listed in the namelist are resolved relative to this path.
- `nsteps::Integer`: number of SR steps to execute. The value is written
  into `data.modpara.nsr_opt_itr_step` before invoking `vmc_para_opt!`,
  overriding whatever was in `modpara.def`. Per-step computation is
  deterministic given the RNG state, so capturing the first `nsteps`
  outputs from an `nsteps`-step Julia run reproduces the first `nsteps`
  outputs from a longer C run that uses the same seed.
- `nsmp::Union{Integer,Nothing} = nothing`: number of final optimisation
  samples to average for C ctest-style checks. `nothing` preserves
  `NSROptItrSmp` from `modpara.def`; an integer overrides it. `nsteps` must
  be greater than or equal to the effective `nsmp`. First-N-step tests that
  shorten `nsteps` below the fixture's `NSROptItrSmp` should pass
  `nsmp = nsteps`.
- `mode::Symbol`: sanity label, one of `:real`, `:cmp`, `:fsz`. The actual
  execution mode (complex flag, generalised orbital, etc.) is determined
  from the parsed `.def` files; this argument is validated for the
  documented set but otherwise not passed downstream.
- `output_dir::AbstractString = tempname()`: directory that receives
  `zvo_out.dat`, `zqp_opt.dat`, etc. Created if absent.
- `seed::Union{Integer,Nothing} = nothing`: SFMT19937 seed. `nothing` means
  resolve `modpara.def`'s `RndSeed` with the C-parity rule (`resolve_rnd_seed`,
  v0.4): missing → `11272` (parser default, C `readdef.c:1967`), `< 0` →
  rank0 time seed broadcast over comm0, `== 0` → `0`, `> 0` → the value; the
  per-group `+ group1` offset is added under MPI (C `vmcmain.c:257`). An
  explicit integer overrides the table (still `+ group1` under MPI).
- `initial_def`: starting variational parameters. `:auto` (default) loads
  `initial.def` from the namelist directory if present (mirrors C's
  `vmc.out -s StdFace.def initial.def` test driver). Pass a path to load
  from a specific file, or `nothing` / `:none` to skip even when one is
  present and rely solely on the RNG-seeded defaults.

# Returns

A `NamedTuple` with the following fields:

| Field | Type | Meaning |
|-------|------|---------|
| `status` | `Int` | Exit code from `vmc_para_opt!` (0 = OK). |
| `output_dir` | `String` | Absolute path to the run outputs. |
| `zvo_first_n` | `Vector{String}` | First `nsteps` raw lines of `zvo_out.dat`. |
| `ctest_values` | `Vector{Float64}` | Mean of columns 1 and 2 over the final `nsmp` `zvo_out.dat` rows, matching the two values C's standard ctest compares from `zqp_opt.dat`. |
| `final_energy_per_site` | `Float64` | First column of the last `zvo_out.dat` line divided by `data.modpara.nsite`. |
| `effective_nsteps` | `Int` | `NSROptItrStep` used for the run. |
| `effective_nsmp` | `Int` | `NSROptItrSmp` used for the run and `ctest_values` averaging. |
"""
function run_para_opt_from_namelist(namelist_path::AbstractString;
                                    nsteps::Integer,
                                    mode::Symbol,
                                    nsmp::Union{Integer,Nothing} = nothing,
                                    output_dir::AbstractString = tempname(),
                                    seed::Union{Integer,Nothing} = nothing,
                                    initial_def::Union{AbstractString,Symbol,Nothing} = :auto)
    mode in (:real, :cmp, :fsz) || throw(ArgumentError("mode must be :real, :cmp, or :fsz; got $mode"))
    nsteps > 0 || throw(ArgumentError("nsteps must be positive; got $nsteps"))
    if nsmp !== nothing && nsmp <= 0
        throw(ArgumentError("nsmp must be positive when provided; got $nsmp"))
    end
    if initial_def isa Symbol && !(initial_def in (:auto, :none))
        throw(ArgumentError("initial_def Symbol must be :auto or :none; got :$initial_def"))
    end

    namelist_str = String(namelist_path)

    # C-compatible section timer. Enabled by MVMC_C_TIMER (MVMC_TIMER kept as a
    # deprecated alias during migration off the old TimerOutputs path). This is
    # the single place that reads the env var and constructs the concrete
    # CTimer; it then flows into vmc_para_opt! through a function barrier.
    c_timer_env = get(ENV, "MVMC_C_TIMER", "0") != "0"
    legacy_timer_env = get(ENV, "MVMC_TIMER", "0") != "0"
    if legacy_timer_env && !c_timer_env
        @warn "MVMC_TIMER is deprecated; use MVMC_C_TIMER=1 for the C-compatible zvo_CalcTimer.dat timer."
    end
    timer_enabled = c_timer_env || legacy_timer_env
    c_timer = CTimer(timer_enabled)
    ctimer_reset!(c_timer)   # fresh per run, in case a timer is reused across repeats

    ctimer_start!(c_timer, 0)   # [0] All
    ctimer_start!(c_timer, 1)   # [1] Initialization

    # The phase ordering below mirrors C-mVMC's vmcmain.c:256-281:
    #   init_gen_rand        → seed RNG
    #   InitParameter        → init_parameter!            (random init only)
    #   ReadInitParameter    → read_initial_def!          (initial.def overlay)
    #   ReadInputParameters  → read_input_parameters!     (In*.def overlay)
    #   SyncModifiedParameter→ sync_modified_parameter!   (Slater rescale + DH/GJ shift)
    #   InitQPWeight         → init_qp_weight!
    # NOTE: we deliberately call `init_parameter!`, NOT `initialize_parameters!`
    # — the latter wraps a parser-side sync that would rescale Slater values
    # before the overlays run, deviating from C's order. See step 3 below.
    # Reordering these breaks bit-level reproducibility for fixtures that
    # ship initial.def / In*.def alongside random init.

    # 1. Parse expert-mode .def files (relative paths resolved from namelist_path).
    ctimer_start!(c_timer, 11)   # [11] ReadDefFile
    data = MVMCExpertModeParsers.parse_expert_mode_files(namelist_str)
    ctimer_stop!(c_timer, 11)

    # v0.4: unsupported input の検証は MPI context 構築より前（plan review F4）。
    # 現行は vmc_para_opt!（vmc_para_opt.jl:83）内でのみ呼ばれるが、R0 で
    # build_parallel_context を parse 直後へ移すため、invalid な NSplitSize
    #（< 1、または R0 では > 1 も未 support）で MPI context を作らないよう
    # ここで先に検証する。vmc_para_opt! 側の既存呼び出しは defense-in-depth
    # としてそのまま残す（重複呼び出しは無害）。
    validate_supported_modpara(data.modpara)

    # v0.4: MPI context は seed / init_parameter! より前に作る（spec §4.2, F3）。
    ctx = build_parallel_context(data.modpara.nsplit_size)
    validate_supported_para_opt_parallel_modpara(ctx, data.modpara)

    effective_nsteps = Int(nsteps)
    effective_nsmp = nsmp === nothing ? data.modpara.nsr_opt_itr_smp : Int(nsmp)
    effective_nsmp > 0 || throw(ArgumentError("effective nsmp must be positive; got $effective_nsmp"))
    if effective_nsteps < effective_nsmp
        throw(ArgumentError("nsteps ($effective_nsteps) must be >= nsmp ($effective_nsmp); smaller nsteps would zero-pad optimisation averages"))
    end

    # 2. Construct and seed the SFMT19937 RNG with the C-compatible convention.
    rng = SFMT19937RNG()
    # C parity seed 解決 + group1 offset（spec §5-1; C vmcmain.c:257）。
    actual_seed = resolve_rnd_seed(ctx, data.modpara.rnd_seed, seed)
    Random.seed!(rng, actual_seed)

    # 3. Random-seeded variational parameter initialisation. We call
    #    `init_parameter!` directly rather than `initialize_parameters!`,
    #    because the latter wraps `init_parameter!` + a parser-side
    #    `sync_modified_parameter!` (Slater rescale only, parameter_init.jl).
    #    That inner sync would rescale random-init Slater values *before*
    #    the initial.def / In*.def overlays, which is not what C does — C's
    #    `InitParameter()` performs no sync, and the only sync happens at
    #    vmcmain.c:276 after both overlays. Step 6 below is the matching
    #    sync; using `init_parameter!` here keeps that the only sync.
    # [13] InitParameter: init_parameter! + overlays + sync (C vmcmain.c:261-277)
    ctimer_start!(c_timer, 13)
    init_parameter!(data; rng = rng)

    # 4. Optional initial.def overlay (C's ReadInitParameter, vmcmain.c:265).
    #    Resolution rules:
    #      - `:auto` (default): if `inputs/initial.def` is absent → silently
    #        skip; if present but read fails → ERROR. Random-init fallback
    #        for a present-but-broken fixture would mask reproducibility
    #        regressions, so we fail loud.
    #      - explicit path: ERROR on any failure (missing or malformed).
    #      - `:none` / `nothing`: skip entirely, even if a sibling
    #        initial.def is present (e.g. for ablation runs).
    #    Required for bit-level match against C references generated via
    #    `vmc.out -s StdFace.def initial.def` (HubbardChain, *_fsz fixtures).
    initial_path, present_for_auto = if initial_def === :auto
        candidate = joinpath(dirname(namelist_str), "initial.def")
        if isfile(candidate)
            (candidate, true)
        else
            (nothing, false)
        end
    elseif initial_def === :none || initial_def === nothing
        (nothing, false)
    else
        (String(initial_def), false)
    end
    if initial_path !== nothing
        ok = read_initial_def!(data, initial_path)
        if !ok
            if present_for_auto
                throw(ArgumentError("read_initial_def! failed on auto-detected initial.def at $initial_path; pass initial_def=:none to skip explicitly"))
            else
                throw(ArgumentError("read_initial_def! failed for explicitly requested path: $initial_path"))
            end
        end
    end

    # 5. Selective In*.def overlay (C's ReadInputParameters, vmcmain.c:268).
    read_input_parameters!(data, namelist_str)

    # 6. Sync / rescale parameters before InitQPWeight, mirroring
    #    C's SyncModifiedParameter call at vmcmain.c:276. This shifts
    #    correlation factors, rescales Slater, and normalizes OptTrans
    #    when the OptTrans mode is active.
    sync_modified_parameter!(ctx, data)
    ctimer_stop!(c_timer, 13)

    # 7. Quantum projection weights (C's InitQPWeight, vmcmain.c:281).
    init_qp_weight!(data)
    ctimer_stop!(c_timer, 1)   # end [1] Initialization

    # Cap the optimisation length at nsteps so that test/example runs stay
    # fast. Each SR step is deterministic given the RNG state and parsed
    # input, so the first `nsteps` outputs from an `nsteps`-step run match
    # the first `nsteps` outputs from a longer run that uses the same seed.
    data.modpara.nsr_opt_itr_step = effective_nsteps
    data.modpara.nsr_opt_itr_smp = effective_nsmp

    # 4. Run optimisation. `mode` is intentionally not passed downstream:
    #    real/cmp/fsz behaviour is encoded in the parsed input files.
    mkpath(output_dir)
    status = vmc_para_opt!(
        data;
        rng = rng,
        output_dir = String(output_dir),
        c_timer = c_timer,
        ctx = ctx,
    )

    ctimer_stop!(c_timer, 0)   # end [0] All (post-run zvo readback below is bookkeeping, not timed)
    if timer_enabled && is_output_rank(ctx)
        write_ctimer_para_opt(c_timer, String(output_dir))
    end

    # rank0 の write 完了を待ってから読み返す（spec §5-9、F11）。
    barrier(ctx)
    if !is_output_rank(ctx)
        # 非 rank0 は readback しない。minimal result を返す。
        return (
            status = status,
            output_dir = abspath(output_dir),
            zvo_first_n = String[],
            ctest_values = Float64[],
            final_energy_per_site = NaN,
            effective_nsteps = effective_nsteps,
            effective_nsmp = effective_nsmp,
        )
    end

    # 5. Read back outputs for caller convenience. Julia writes zvo_out.dat
    #    (not the C-style zvo_out_001.dat).
    zvo_path = joinpath(output_dir, "zvo_out.dat")
    zvo_lines = strip.(readlines(zvo_path))
    if length(zvo_lines) < effective_nsteps
        throw(ArgumentError("expected at least $effective_nsteps zvo_out.dat rows, found $(length(zvo_lines)) at $zvo_path"))
    end
    zvo_first_n = zvo_lines[1:effective_nsteps]
    zvo_rows = [parse.(Float64, split(line)) for line in zvo_first_n]
    if any(row -> length(row) < 2, zvo_rows)
        throw(ArgumentError("zvo_out.dat must contain at least two columns for C ctest comparison: $zvo_path"))
    end
    ctest_start = effective_nsteps - effective_nsmp + 1
    ctest_rows = zvo_rows[ctest_start:effective_nsteps]
    ctest_values = [
        sum(row[1] for row in ctest_rows) / effective_nsmp,
        sum(row[2] for row in ctest_rows) / effective_nsmp,
    ]
    final_e = zvo_rows[end][1]
    nsite = data.modpara.nsite
    nsite > 0 || throw(ArgumentError("modpara.nsite must be positive to compute energy per site"))

    return (
        status = status,
        output_dir = abspath(output_dir),
        zvo_first_n = zvo_first_n,
        ctest_values = ctest_values,
        final_energy_per_site = final_e / nsite,
        effective_nsteps = effective_nsteps,
        effective_nsmp = effective_nsmp,
    )
end
