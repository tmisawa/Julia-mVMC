# Green-function file comparison helpers for the PhysCal end-to-end gate (Plan 3).
#
# Test-only, no package dependency. The PhysCal reference gate compares Julia's
# `zvo_cisajs` / `zvo_cisajscktalt` (direct) / `zvo_cisajscktaltex` (factored)
# against committed C references.
#
# Two-stage comparison (spec Finding 8):
#   1. Raw bytes: read(path, String) == read(path, String). This must NOT strip
#      or tokenize — C's exact format (trailing spaces, a trailing blank line for
#      the indexed files) is part of the bit-for-bit completion condition, and
#      the existing harnesses' strip()+split() would silently mask a regression.
#   2. Numeric fallback when bytes differ: integer index columns must still match
#      exactly; float values compared with `isapprox(rtol, atol)`.
#
# Per-quantity tolerances mirror test/integration/runtests.jl's split between
# linear accumulators and product/squared accumulators (which pick up amplified
# C/Julia BLAS summation-order noise):
#   - one-body  zvo_cisajs        : linear           → rtol 1e-10
#   - factored  zvo_cisajscktaltex: Σ w·g·conj(g)    → rtol 1e-9  (like <H^2>)
#   - direct DC zvo_cisajscktalt  : 4-operator prod  → rtol 1e-9
# atol = 1e-12 floors near-zero Green values so they don't fail on relative error.

const GREEN_TOL_LINEAR = 1e-10   # cf. runtests.jl TOL_DEFAULT
const GREEN_TOL_PRODUCT = 1e-9   # cf. runtests.jl TOL_LOOSE
const GREEN_ATOL = 1e-12

"Result of comparing one Green file against its reference."
struct GreenCompareResult
    ok::Bool             # passed: exact, or indices match and values within tol
    exact::Bool          # raw-byte identical
    fallback_used::Bool  # numeric (non-exact) comparison was needed
    n_values::Int        # number of float values compared
    max_abs_err::Float64
    max_rel_err::Float64
    rtol::Float64
    atol::Float64
    detail::String
end

function _green_fail(detail; rtol, atol)
    return GreenCompareResult(false, false, true, 0, NaN, NaN, rtol, atol, detail)
end

# Non-empty, whitespace-split lines.
_green_lines(text::AbstractString) =
    [split(strip(l)) for l in split(text, '\n') if !isempty(strip(l))]

# Columnar comparison: each data row is `n_int` integer columns followed by
# exactly 2 float columns (real, imag). Used for one-body (n_int=4) and direct
# DC (n_int=8). Integer columns must match exactly; floats within tolerance.
function _compare_columnar(jl_path, c_path, n_int::Int; rtol, atol)
    isfile(jl_path) || return _green_fail("Julia output missing: $jl_path"; rtol, atol)
    isfile(c_path) || return _green_fail("reference missing: $c_path"; rtol, atol)
    raw_jl = read(jl_path, String)
    raw_c = read(c_path, String)
    if raw_jl == raw_c
        nfloat = 2 * length(_green_lines(raw_c))
        return GreenCompareResult(true, true, false, nfloat, 0.0, 0.0, rtol, atol, "exact")
    end

    jl_rows = _green_lines(raw_jl)
    c_rows = _green_lines(raw_c)
    if length(jl_rows) != length(c_rows)
        return _green_fail(
            "row count mismatch: julia=$(length(jl_rows)) ref=$(length(c_rows))"; rtol, atol,
        )
    end

    width = n_int + 2
    max_abs = 0.0
    max_rel = 0.0
    nvals = 0
    for (r, (jr, cr)) in enumerate(zip(jl_rows, c_rows))
        if length(jr) != width || length(cr) != width
            return _green_fail(
                "row $r column count: julia=$(length(jr)) ref=$(length(cr)) (expected $width)";
                rtol, atol,
            )
        end
        for k = 1:n_int
            if jr[k] != cr[k]
                return _green_fail(
                    "index mismatch row $r col $k: julia='$(jr[k])' ref='$(cr[k])'"; rtol, atol,
                )
            end
        end
        for k = (n_int+1):width
            a = parse(Float64, jr[k])
            b = parse(Float64, cr[k])
            ae = abs(a - b)
            max_abs = max(max_abs, ae)
            max_rel = max(max_rel, ae / max(abs(b), atol))
            nvals += 1
            if !isapprox(a, b; rtol = rtol, atol = atol)
                return GreenCompareResult(
                    false, false, true, nvals, max_abs, max_rel, rtol, atol,
                    "value out of tolerance at row $r col $k: julia=$a ref=$b " *
                    "(abs=$ae, rtol=$rtol, atol=$atol)",
                )
            end
        end
    end
    return GreenCompareResult(
        true, false, true, nvals, max_abs, max_rel, rtol, atol,
        "within tolerance (rtol=$rtol, atol=$atol; max_abs=$max_abs)",
    )
end

# Value-only comparison: the factored file is all `real imag` floats on a single
# line (no index columns). Compare element-wise.
function _compare_value_only(jl_path, c_path; rtol, atol)
    isfile(jl_path) || return _green_fail("Julia output missing: $jl_path"; rtol, atol)
    isfile(c_path) || return _green_fail("reference missing: $c_path"; rtol, atol)
    raw_jl = read(jl_path, String)
    raw_c = read(c_path, String)

    # The factored file is value-only with ALL values on a single line (C's
    # vmcmain.c writes "% .18e  % .18e " in a loop then one newline). Enforce that
    # invariant *before* the raw-byte exact check, so even a byte-identical
    # multi-line pair is rejected — the helper's contract is "factored is one
    # line", and the numeric fallback flattens whitespace so it must not be a way
    # in for a layout regression.
    jl_nlines = length(_green_lines(raw_jl))
    c_nlines = length(_green_lines(raw_c))
    if jl_nlines != 1 || c_nlines != 1
        return _green_fail(
            "factored file must be a single line of values; julia=$jl_nlines ref=$c_nlines non-empty lines";
            rtol, atol,
        )
    end

    if raw_jl == raw_c
        return GreenCompareResult(
            true, true, false, length(split(strip(raw_c))), 0.0, 0.0, rtol, atol, "exact",
        )
    end

    jl_vals = split(strip(raw_jl))
    c_vals = split(strip(raw_c))
    if length(jl_vals) != length(c_vals)
        return _green_fail(
            "value count mismatch: julia=$(length(jl_vals)) ref=$(length(c_vals))"; rtol, atol,
        )
    end

    max_abs = 0.0
    max_rel = 0.0
    for k in eachindex(c_vals)
        a = parse(Float64, jl_vals[k])
        b = parse(Float64, c_vals[k])
        ae = abs(a - b)
        max_abs = max(max_abs, ae)
        max_rel = max(max_rel, ae / max(abs(b), atol))
        if !isapprox(a, b; rtol = rtol, atol = atol)
            return GreenCompareResult(
                false, false, true, k, max_abs, max_rel, rtol, atol,
                "value out of tolerance at field $k: julia=$a ref=$b " *
                "(abs=$ae, rtol=$rtol, atol=$atol)",
            )
        end
    end
    return GreenCompareResult(
        true, false, true, length(c_vals), max_abs, max_rel, rtol, atol,
        "within tolerance (rtol=$rtol, atol=$atol; max_abs=$max_abs)",
    )
end

"""
    compare_green_one(julia_path, c_path; rtol=GREEN_TOL_LINEAR, atol=GREEN_ATOL)

Compare a `zvo_cisajs` one-body Green file (6 columns `ri si rj sj real imag`,
linear accumulator → tight tolerance).
"""
compare_green_one(julia_path, c_path; rtol = GREEN_TOL_LINEAR, atol = GREEN_ATOL) =
    _compare_columnar(julia_path, c_path, 4; rtol = rtol, atol = atol)

"""
    compare_green_two_dc(julia_path, c_path; rtol=GREEN_TOL_PRODUCT, atol=GREEN_ATOL)

Compare a `zvo_cisajscktalt` direct two-body Green file (10 columns
`ri si rj sj rk sk rl sl real imag`, 4-operator product → looser tolerance).
"""
compare_green_two_dc(julia_path, c_path; rtol = GREEN_TOL_PRODUCT, atol = GREEN_ATOL) =
    _compare_columnar(julia_path, c_path, 8; rtol = rtol, atol = atol)

"""
    compare_green_factored(julia_path, c_path; rtol=GREEN_TOL_PRODUCT, atol=GREEN_ATOL)

Compare a `zvo_cisajscktaltex` factored two-body Green file (value-only, all
`real imag` floats on one line; product accumulator → looser tolerance).
"""
compare_green_factored(julia_path, c_path; rtol = GREEN_TOL_PRODUCT, atol = GREEN_ATOL) =
    _compare_value_only(julia_path, c_path; rtol = rtol, atol = atol)
