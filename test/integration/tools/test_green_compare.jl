using Test

include(joinpath(@__DIR__, "green_compare.jl"))

# Write `content` to a fresh file under `dir` and return its path.
_wf(dir, name, content) = (p = joinpath(dir, name); write(p, content); p)

const ONE_A = "0 0 1 0 1.0 0.0 \n1 1 0 1 2.0 0.0 \n\n"        # C-like: trailing space + blank line
const DC_A = "0 0 1 0 2 0 3 0 1.5 -0.5\n\n"
const FAC_A = "2.5  -1.0  0.3  0.1 \n"                          # value-only, one line

@testset "green_compare: exact (raw-byte) match" begin
    mktempdir() do dir
        a = _wf(dir, "a_cisajs.dat", ONE_A)
        b = _wf(dir, "b_cisajs.dat", ONE_A)
        r = compare_green_one(a, b)
        @test r.ok && r.exact && !r.fallback_used
        @test r.n_values == 4  # 2 rows × 2 floats

        fa = _wf(dir, "a_ex.dat", FAC_A)
        fb = _wf(dir, "b_ex.dat", FAC_A)
        rf = compare_green_factored(fa, fb)
        @test rf.ok && rf.exact && rf.n_values == 4
    end
end

@testset "green_compare: numeric fallback within tolerance" begin
    mktempdir() do dir
        a = _wf(dir, "a.dat", ONE_A)
        # Perturb one value by 5e-11 (< one-body rtol 1e-10) → not exact, but ok.
        b = _wf(dir, "b.dat", replace(ONE_A, "1.0 0.0" => "1.00000000005 0.0"))
        r = compare_green_one(a, b)
        @test r.ok && !r.exact && r.fallback_used
        @test r.max_abs_err < 1e-10
    end
end

@testset "green_compare: beyond tolerance fails" begin
    mktempdir() do dir
        a = _wf(dir, "a.dat", ONE_A)
        b = _wf(dir, "b.dat", replace(ONE_A, "1.0 0.0" => "1.000000005 0.0"))  # 5e-9 > 1e-10
        r = compare_green_one(a, b)
        @test !r.ok && occursin("out of tolerance", r.detail)
    end
end

@testset "green_compare: per-quantity tolerance differs (one-body tight, factored loose)" begin
    mktempdir() do dir
        # A 5e-10 difference: fails the one-body 1e-10 bound, passes factored 1e-9.
        one_a = _wf(dir, "1a.dat", "0 0 1 0 1.0 0.0\n")
        one_b = _wf(dir, "1b.dat", "0 0 1 0 1.0000000005 0.0\n")
        @test !compare_green_one(one_a, one_b).ok

        fac_a = _wf(dir, "fa.dat", "1.0 0.0\n")
        fac_b = _wf(dir, "fb.dat", "1.0000000005 0.0\n")
        @test compare_green_factored(fac_a, fac_b).ok
    end
end

@testset "green_compare: structural mismatches hard-fail regardless of tol" begin
    mktempdir() do dir
        a = _wf(dir, "a.dat", ONE_A)
        # Row count mismatch.
        b_rows = _wf(dir, "b_rows.dat", "0 0 1 0 1.0 0.0 \n")
        @test occursin("row count mismatch", compare_green_one(a, b_rows).detail)
        # Index (integer) column mismatch — even with identical float values.
        b_idx = _wf(dir, "b_idx.dat", replace(ONE_A, "1 1 0 1" => "1 1 0 0"))
        r_idx = compare_green_one(a, b_idx)
        @test !r_idx.ok && occursin("index mismatch", r_idx.detail)
        # Column-count mismatch (wrong number of fields per row).
        b_cols = _wf(dir, "b_cols.dat", "0 0 1 0 1.0\n1 1 0 1 2.0 0.0\n")
        @test occursin("column count", compare_green_one(a, b_cols).detail)
        # Factored value-count mismatch.
        fa = _wf(dir, "fa.dat", FAC_A)
        fb = _wf(dir, "fb.dat", "2.5 -1.0\n")
        @test occursin("value count mismatch", compare_green_factored(fa, fb).detail)
    end
end

@testset "green_compare: raw-byte check is whitespace-sensitive but numeric fallback passes" begin
    mktempdir() do dir
        # Same numeric values, but the reference drops C's trailing space / blank
        # line. Raw-byte must see them as different (exact=false), proving the
        # check is not silently strip()ed; the numeric fallback still passes.
        a = _wf(dir, "a.dat", "0 0 1 0 1.0 0.0 \n\n")   # trailing space + blank line
        b = _wf(dir, "b.dat", "0 0 1 0 1.0 0.0\n")      # no trailing space, no blank line
        r = compare_green_one(a, b)
        @test !r.exact          # bytes differ
        @test r.ok              # values identical → numeric fallback passes
        @test r.fallback_used
    end
end

@testset "green_compare: missing files fail cleanly" begin
    mktempdir() do dir
        a = _wf(dir, "a.dat", ONE_A)
        @test occursin("missing", compare_green_one(a, joinpath(dir, "nope.dat")).detail)
        @test occursin("missing", compare_green_one(joinpath(dir, "nope.dat"), a).detail)
    end
end
