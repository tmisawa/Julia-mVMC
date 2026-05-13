using MVMCExpertModeParsers: count_rbm_parameters

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

    n_rbm = count_rbm_parameters(data)
    if n_rbm > 0
        # An initial.def for an RBM model would carry NRBM triples between
        # the NProj and NSlater blocks. Reading them as Slater would silently
        # corrupt every Slater parameter, so refuse rather than guess.
        @warn "read_initial_def!: RBM-bearing models are not yet supported by the loader; refusing to apply initial.def" n_rbm path=initial_path
        return false
    end

    values = parse.(Float64, split(strip(read(initial_path, String))))

    n_gutzwiller = length(data.gutzwiller_terms)
    n_jastrow = length(data.jastrow_terms)
    n_proj = n_gutzwiller + n_jastrow
    n_slater = data.modpara.n_orbital_idx

    # ── Validation phase: every check must run *before* any mutation, so a
    # rejected file leaves the caller's `data` unchanged. ────────────────
    expected_floats = 6 + 3 * (n_proj + n_slater)  # NRBM = 0 verified above
    if length(values) < expected_floats
        @warn "initial.def too short" got=length(values) expected_min=expected_floats path=initial_path
        return false
    end

    # Any trailing floats beyond `expected_floats` indicate either an
    # OptTrans block (C reads NOptTrans triples when FlagOptTrans > 0;
    # Julia has no consumer for them — no OptTrans[] storage, no
    # FlagOptTrans gate, no qp_weight_update consumption) or a malformed
    # / truncated file. Compare floats, *not* triples: with `÷ 3`, 1 or 2
    # stray tokens silently round to zero remaining triples and would let
    # a malformed file pass. Refuse on any non-zero excess.
    extra_floats = length(values) - expected_floats
    if extra_floats > 0
        n_extra_triples, rem = divrem(extra_floats, 3)
        msg = if rem == 0
            "OptTrans-style block of $n_extra_triples triples (not supported)"
        else
            "$extra_floats trailing floats (not a whole number of triples; file likely malformed)"
        end
        @warn "read_initial_def!: $msg; refusing to apply initial.def" expected_floats got=length(values) path=initial_path
        return false
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

    return true
end
