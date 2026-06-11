"""
    run_phys_cal_from_namelist(namelist_path; opt_para, mode, seed=nothing,
                               output_dir=tempname()) -> NamedTuple

Drive `vmc_phys_cal!` (VMCPhysCal mode) from a C-mVMC `namelist.def`, loading a
committed *fixed* optimized-parameter file (`zqp_opt.dat`) so the Green-function
output is deterministic and decoupled from optimization. This is the entry point
used by the PhysCal reference gate.

# Init order (PhysCal — note the difference from `run_para_opt_from_namelist`)

    parse → seed RNG → read_opt_para_file! → read_input_parameters!
          → sync_modified_parameter! → vmc_phys_cal!(data; rng, output_dir)

Unlike `run_para_opt_from_namelist`, this runner deliberately does **not** call
`init_parameter!` or `init_qp_weight!`. `vmc_phys_cal!` performs both internally:
its save → `init_parameter!` → restore dance consumes the SFMT RNG exactly once
to match C's single `InitParameter()` (`vmc_phys_cal.jl`), and it calls
`init_qp_weight!` itself. Calling either here as well would advance the RNG a
second time before sampling and desynchronise the Monte-Carlo stream from the C
reference — a bit-level failure that looks like a BLAS/platform issue but is not.
The runner only loads + overlays + syncs the parameters; `vmc_phys_cal!`'s
internal save/restore preserves those values across its single internal
`init_parameter!`.

`read_opt_para_file!` runs right after parsing: it has no `init_parameter!`
precondition (the parser already created the Gutzwiller/Jastrow/orbital term
arrays and set `modpara.n_orbital_idx`; the loader only sets each term's
`.value`).

# Arguments
- `namelist_path::AbstractString`: path to `namelist.def`; sibling `.def` files
  resolve relative to it.
- `opt_para::AbstractString` (**required**): path to the committed fixed
  `zqp_opt.dat`. Explicit and required — the gate must be deterministic about
  which parameters it consumed.
- `mode::Symbol`: `:real` / `:cmp` / `:fsz` sanity label (validated for the
  documented set; the real execution path is determined by the parsed `.def`
  files, as in `run_para_opt_from_namelist`, so it is not passed downstream).
- `seed::Union{Integer,Nothing} = nothing`: SFMT19937 seed; `nothing` uses
  `modpara.RndSeed` with the **legacy rule** (`<= 0` → `11272`). Note: unlike
  `run_para_opt_from_namelist`, which uses the C-parity `resolve_rnd_seed`
  (`0` → `0`, `< 0` → time seed) as of v0.4 R0, this runner keeps the legacy
  rule until its R1 MPI-aware rewrite (review 2026-06-11 F6-(2)) — `RndSeed 0`
  therefore seeds `11272` here but `0` in the para-opt runner.
- `output_dir::AbstractString = tempname()`: directory for the `zvo_*` outputs
  (created if absent).

Sampling counts (`NDataQtySmp`, `NVMCSample`) come from `modpara`; there is no
`nsmp` override — in PhysCal those counts are distinct and easy to confuse, so
`modpara` (from the namelist) is the single source of truth.

# Returns
A `NamedTuple` `(; status, output_dir, n_para_consumed)`:
- `status::Int`: exit code from `vmc_phys_cal!` (0 = OK).
- `output_dir::String`: absolute path to the run outputs.
- `n_para_consumed::Int`: number of fixed parameters applied
  (`n_proj + n_slater`); the gate can assert this is `> 0`.
"""
function run_phys_cal_from_namelist(
    namelist_path::AbstractString;
    opt_para::AbstractString,
    mode::Symbol,
    seed::Union{Integer,Nothing} = nothing,
    output_dir::AbstractString = tempname(),
)
    mode in (:real, :cmp, :fsz) ||
        throw(ArgumentError("mode must be :real, :cmp, or :fsz; got :$mode"))

    # v0.4 R0: PhysCal runner は MPI 未対応（R1 で対応予定）。mpiexec 配下で
    # そのまま走らせると全 rank が同一 output_dir の zvo_* へ並行書き込みして
    # silent corruption になるため fail-fast する（review 2026-06-11 F7）。
    # JULIA_MVMC_MPI=0 の明示時（:mpi_guarded_serial）は build_parallel_context が
    # size>1 なら MPI.Abort、size==1 なら serial で続行する（F12 と同じ guard）。
    mpi_mode = resolve_mpi_mode()
    if mpi_mode === :mpi
        error("run_phys_cal_from_namelist is not MPI-aware in v0.4 R0 " *
              "(planned for R1). Run it without mpiexec, or set JULIA_MVMC_MPI=0 " *
              "to force a single-rank serial run.")
    elseif mpi_mode === :mpi_guarded_serial
        build_parallel_context(1)
    end

    namelist_str = String(namelist_path)

    # 1. Parse expert-mode .def files (relative paths resolved from namelist_path).
    data = MVMCExpertModeParsers.parse_expert_mode_files(namelist_str)

    # 2. Seed the SFMT19937 RNG (C-compatible convention) and pass it to
    #    vmc_phys_cal! so its single internal RNG consumption is reproducible.
    rng = SFMT19937RNG()
    actual_seed = if seed === nothing
        data.modpara.rnd_seed > 0 ? data.modpara.rnd_seed : 11272
    else
        Int(seed)
    end
    Random.seed!(rng, actual_seed)

    # 3. Load the fixed optimized parameters (strict; returns the consumed count).
    n_para_consumed = read_opt_para_file!(data, String(opt_para))

    # 4. In*.def overlays (C ReadInputParameters), after the fixed load so any
    #    overlay wins — matching C's order.
    read_input_parameters!(data, namelist_str)

    # 5. Slater rescale before PhysCal. C enables DH/GJ shift flags only in Opt mode.
    sync_modified_parameter!(data; shift_correlations = false)

    # 6. PhysCal. vmc_phys_cal! owns the single init_parameter! (RNG match) and
    #    init_qp_weight!; its save/restore preserves the params loaded above.
    out_dir = abspath(String(output_dir))  # absolute, matching run_para_opt_from_namelist
    mkpath(out_dir)
    status = vmc_phys_cal!(data; rng = rng, output_dir = out_dir)

    return (; status = status, output_dir = out_dir, n_para_consumed = n_para_consumed)
end
