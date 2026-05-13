"""
    run_para_opt_from_namelist(namelist_path;
                               nsteps, mode,
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
- `mode::Symbol`: sanity label, one of `:real`, `:cmp`, `:fsz`. The actual
  execution mode (complex flag, generalised orbital, etc.) is determined
  from the parsed `.def` files; this argument is validated for the
  documented set but otherwise not passed downstream.
- `output_dir::AbstractString = tempname()`: directory that receives
  `zvo_out.dat`, `zqp_opt.dat`, etc. Created if absent.
- `seed::Union{Integer,Nothing} = nothing`: SFMT19937 seed. `nothing` means
  use `modpara.def`'s `RndSeed` field, falling back to the C-compatible
  default `11272` when `RndSeed` is non-positive.
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
| `final_energy_per_site` | `Float64` | First column of the last `zvo_out.dat` line divided by `data.modpara.nsite`. |
"""
function run_para_opt_from_namelist(namelist_path::AbstractString;
                                    nsteps::Integer,
                                    mode::Symbol,
                                    output_dir::AbstractString = tempname(),
                                    seed::Union{Integer,Nothing} = nothing,
                                    initial_def::Union{AbstractString,Symbol,Nothing} = :auto)
    mode in (:real, :cmp, :fsz) || throw(ArgumentError("mode must be :real, :cmp, or :fsz; got $mode"))
    nsteps > 0 || throw(ArgumentError("nsteps must be positive; got $nsteps"))
    if initial_def isa Symbol && !(initial_def in (:auto, :none))
        throw(ArgumentError("initial_def Symbol must be :auto or :none; got :$initial_def"))
    end

    namelist_str = String(namelist_path)

    # The phase ordering below mirrors C-mVMC's vmcmain.c:256-281:
    #   init_gen_rand        → seed RNG
    #   InitParameter        → init_parameter!            (random init only)
    #   ReadInitParameter    → read_initial_def!          (initial.def overlay)
    #   ReadInputParameters  → read_input_parameters!     (In*.def overlay)
    #   SyncModifiedParameter→ sync_modified_parameter!   (Slater rescale + GJ shift)
    #   InitQPWeight         → init_qp_weight!
    # NOTE: we deliberately call `init_parameter!`, NOT `initialize_parameters!`
    # — the latter wraps a parser-side sync that would rescale Slater values
    # before the overlays run, deviating from C's order. See step 3 below.
    # Reordering these breaks bit-level reproducibility for fixtures that
    # ship initial.def / In*.def alongside random init.

    # 1. Parse expert-mode .def files (relative paths resolved from namelist_path).
    data = MVMCExpertModeParsers.parse_expert_mode_files(namelist_str)

    # 2. Construct and seed the SFMT19937 RNG with the C-compatible convention.
    rng = SFMT19937RNG()
    actual_seed = if seed === nothing
        data.modpara.rnd_seed > 0 ? data.modpara.rnd_seed : 11272
    else
        Int(seed)
    end
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
    #    C's SyncModifiedParameter call at vmcmain.c:276. With OptTrans
    #    not yet implemented in Julia this currently rescales Slater and
    #    shifts Gutzwiller/Jastrow only; see parameter_sync.jl for the
    #    explicit OptTrans carve-out.
    sync_modified_parameter!(data)

    # 7. Quantum projection weights (C's InitQPWeight, vmcmain.c:281).
    init_qp_weight!(data)

    # Cap the optimisation length at nsteps so that test/example runs stay
    # fast. Each SR step is deterministic given the RNG state and parsed
    # input, so the first `nsteps` outputs from an `nsteps`-step run match
    # the first `nsteps` outputs from a longer run that uses the same seed.
    data.modpara.nsr_opt_itr_step = Int(nsteps)

    # 4. Run optimisation. `mode` is intentionally not passed downstream:
    #    real/cmp/fsz behaviour is encoded in the parsed input files.
    mkpath(output_dir)
    status = vmc_para_opt!(
        data;
        rng = rng,
        output_dir = String(output_dir),
    )

    # 5. Read back outputs for caller convenience. Julia writes zvo_out.dat
    #    (not the C-style zvo_out_001.dat).
    zvo_path = joinpath(output_dir, "zvo_out.dat")
    zvo_lines = open(zvo_path) do io
        [strip(readline(io)) for _ in 1:nsteps]
    end
    final_e = parse(Float64, split(zvo_lines[end])[1])
    nsite = data.modpara.nsite
    nsite > 0 || throw(ArgumentError("modpara.nsite must be positive to compute energy per site"))

    return (
        status = status,
        output_dir = abspath(output_dir),
        zvo_first_n = zvo_lines,
        final_energy_per_site = final_e / nsite,
    )
end
