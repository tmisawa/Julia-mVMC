using MVMCExpertModeParsers: count_rbm_parameters

"""
    _load_para_triples!(data::ExpertModeData, text::AbstractString)
        -> (ok::Bool, n_consumed::Int, reason::String)

Shared core for the optimized-parameter file layout used by both
`read_initial_def!` (warn+false on problems) and `read_opt_para_file!`
(error on problems). Parses mVMC's `zqp_opt.dat` / `initial.def` record
(`NSROptItrSmp > 1` layout):

- 6 leading floats (energy diagnostics) skipped;
- `NProj` triples `(real, imag, gradient)` for Gutzwiller then Jastrow;
- `NSlater` triples for orbital parameters, scattered into `orbital_terms`
  via each term's `idx`.

Scope: **Gutzwiller + Jastrow + Slater only.** RBM and DoublonHolon are
rejected up front: the C file carries extra triples for them (between NProj
and NSlater) that this loader does not place, so silently reading them would
corrupt Slater (design-review LOADER-1). SpinJastrow is not parsed by the
toolchain at all; if a future model adds it, the float-count check below
fails loud rather than mis-attributing.

Validate-before-commit: on any problem returns `(false, 0, reason)` and
leaves `data` untouched. On success commits and returns
`(true, n_proj + n_slater, "")`. Non-numeric tokens are reported as a stable
reason (not a raw `ArgumentError`) so strict callers can test for them.
"""
function _load_para_triples!(data::ExpertModeData, text::AbstractString)
    n_rbm = count_rbm_parameters(data)
    if n_rbm > 0
        return (false, 0, "RBM-bearing models are not supported by this loader (n_rbm = $n_rbm)")
    end
    n_dh =
        length(data.doublon_holon_2site_terms) + length(data.doublon_holon_4site_terms)
    if n_dh > 0
        return (
            false,
            0,
            "DoublonHolon parameters are not supported by this loader (n_dh = $n_dh)",
        )
    end

    tokens = split(strip(text))
    values = Vector{Float64}(undef, length(tokens))
    for (i, tok) in enumerate(tokens)
        v = tryparse(Float64, tok)
        v === nothing && return (false, 0, "non-numeric token '$(tok)' at field $i")
        # A reference-gate parameter file must carry finite values; `tryparse`
        # accepts "NaN"/"Inf", so reject them here rather than feeding a poisoned
        # parameter into the deterministic comparison.
        isfinite(v) || return (false, 0, "non-finite token '$(tok)' at field $i")
        values[i] = v
    end

    n_gutzwiller = length(data.gutzwiller_terms)
    n_jastrow = length(data.jastrow_terms)
    n_proj = n_gutzwiller + n_jastrow
    n_slater = data.modpara.n_orbital_idx

    expected_floats = 6 + 3 * (n_proj + n_slater)  # NRBM = 0, NDH = 0 verified above
    if length(values) < expected_floats
        return (
            false,
            0,
            "too short: got $(length(values)) floats, expected $expected_floats " *
            "(6 + 3*(NProj=$n_proj + NSlater=$n_slater))",
        )
    end
    # Trailing floats beyond the expected count indicate an OptTrans block
    # (unsupported) or a malformed file. Compare floats, not triples, so 1–2
    # stray tokens cannot round to zero remaining triples and slip through.
    extra_floats = length(values) - expected_floats
    if extra_floats > 0
        n_extra_triples, rem = divrem(extra_floats, 3)
        reason =
            rem == 0 ?
            "OptTrans-style block of $n_extra_triples triples (not supported)" :
            "$extra_floats trailing floats (not a whole number of triples; file likely malformed)"
        return (false, 0, reason)
    end

    # ── Commit phase: validation passed, now mutate `data`. ──────────────
    idx = 7  # skip 6 leading floats (1-based: start at index 7)
    @inbounds for i = 1:n_gutzwiller
        data.gutzwiller_terms[i].value = ComplexF64(values[idx], values[idx+1])
        idx += 3
    end
    @inbounds for i = 1:n_jastrow
        data.jastrow_terms[i].value = ComplexF64(values[idx], values[idx+1])
        idx += 3
    end
    slater_values = Vector{ComplexF64}(undef, n_slater)
    @inbounds for i = 1:n_slater
        slater_values[i] = ComplexF64(values[idx], values[idx+1])
        idx += 3
    end
    @inbounds for term in data.orbital_terms
        if term.idx >= 0 && term.idx < n_slater
            term.value = slater_values[term.idx+1]
        end
    end

    return (true, n_proj + n_slater, "")
end

"""
    read_initial_def!(data::ExpertModeData, initial_path::AbstractString) -> Bool

Overwrite the variational parameters in `data` with values stored in
`initial_path` (mVMC's `initial.def` format).

The C-mVMC test driver invokes the optimiser with
`vmc.out -s StdFace.def initial.def`, where `initial.def` is generated from a
prior `zqp_opt.dat` plus a Python `random.uniform(-0.01, 0.01)` perturbation
(see `mVMC/test/python/data/<Model>/AddRand.py`). To reproduce that reference
output bit-for-bit from Julia, the same starting parameters must be loaded
*after* `init_parameter!` (which seeds parameters from the SFMT RNG without
the parser-side sync wrap) and *before* `read_input_parameters!` overlays
any `In*.def` files, matching C's order in `vmcmain.c:264-272`
(`InitParameter` → `ReadInitParameter` → `ReadInputParameters`).

Layout (whitespace-separated, all on a single record):
1. 6 leading floats are skipped (energy diagnostics written by C).
2. NProj triples `(real, imag, gradient)` for Gutzwiller + Jastrow.
3. NRBM triples for RBM (only when FlagRBM > 0). **Not yet supported.**
4. NSlater triples for orbital parameters; values are scattered into
   `data.orbital_terms` via each term's `idx`.
5. NOptTrans triples for OptTrans (`FlagOptTrans > 0`). **Not yet supported.**

Returns `true` on success, `false` on a recoverable problem (missing file,
file too short for the parameter count, or unsupported block present such
as RBM or OptTrans). On `false` an `@warn` is emitted **and `data` is left
untouched** — all validation runs before any mutation, so a refused load
never leaves a partial overlay behind. The caller (e.g.
`run_para_opt_from_namelist`) decides whether to abort or continue.
"""
function read_initial_def!(data::ExpertModeData, initial_path::AbstractString)
    if !isfile(initial_path)
        @warn "initial.def not found" path=initial_path
        return false
    end
    ok, _, reason = _load_para_triples!(data, read(initial_path, String))
    if !ok
        @warn "read_initial_def!: $reason; refusing to apply initial.def" path=initial_path
        return false
    end
    return true
end

"""
    read_opt_para_file!(data::ExpertModeData, opt_para_path::AbstractString) -> Int

Strict, non-perturbing loader for a committed C `zqp_opt.dat` (the optimized
variational parameters), for the PhysCal reference gate. Same 6+triples layout
as [`read_initial_def!`](@ref) (see [`_load_para_triples!`](@ref)), but built for
a deterministic gate rather than an optional warm-start overlay:

- **Fails hard** (`error`) instead of `@warn`+`false` on a missing file, a
  non-numeric or non-finite (`NaN`/`Inf`) token, a too-short or trailing-garbage
  record, or an unsupported block (RBM / DoublonHolon / OptTrans).
- **Returns the number of parameters consumed** (`n_proj + n_slater`) so the
  caller can assert the fixed parameters were actually applied; errors if zero.
- Does **not** perturb (the C test driver's `random.uniform` perturbation is
  external), so a committed unperturbed `zqp_opt.dat` loads verbatim.

Scope is the `NSROptItrSmp > 1` layout for Gutzwiller + Jastrow + Slater models
(the Plan 3 PhysCal fixtures); the `NSROptItrSmp == 1` "pairs" layout and
RBM/DoublonHolon/OptTrans models are rejected.
"""
function read_opt_para_file!(data::ExpertModeData, opt_para_path::AbstractString)::Int
    isfile(opt_para_path) ||
        error("read_opt_para_file!: file not found: $opt_para_path")
    ok, n_consumed, reason = _load_para_triples!(data, read(opt_para_path, String))
    ok || error("read_opt_para_file!: $reason (path: $opt_para_path)")
    n_consumed > 0 ||
        error("read_opt_para_file!: no parameters consumed (no Gutzwiller/Jastrow/Slater terms) from $opt_para_path")
    return n_consumed
end
