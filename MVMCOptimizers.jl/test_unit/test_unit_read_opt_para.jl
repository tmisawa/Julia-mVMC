using Test
using MVMCOptimizers
using MVMCExpertModeParsers:
    ExpertModeData, GutzwillerTerm, JastrowTerm, OrbitalTerm, DoublonHolon2SiteTerm

# Capture the showerror text of whatever `f()` throws (message is the stable
# contract; assert on substrings, not the concrete exception type).
function _capture_msg(f)
    try
        f()
        return (false, "")
    catch err
        return (true, sprint(showerror, err))
    end
end

# Build an ExpertModeData with 2 Gutzwiller, 1 Jastrow, 2 Slater (orbital idx
# 0,1) terms — the layout read_opt_para_file! supports.
function _make_data()
    data = ExpertModeData()
    data.modpara.n_orbital_idx = 2
    data.gutzwiller_terms = [GutzwillerTerm(0, 0.0 + 0im, true), GutzwillerTerm(1, 0.0 + 0im, true)]
    data.jastrow_terms = [JastrowTerm(0, 1, 0.0 + 0im, true)]
    data.orbital_terms = [
        OrbitalTerm(0, 1, 0, 0.0 + 0im, true),
        OrbitalTerm(1, 2, 1, 0.0 + 0im, true),
    ]
    return data
end

# A valid zqp_opt.dat record: 6 leading floats + 5 triples (re, im, grad) for
# gutz1, gutz2, jastrow1, slater0, slater1.
const _GOLDEN_OPT = join(
    [
        "1.0", "2.0", "3.0", "4.0", "5.0", "6.0",          # 6 energy diagnostics (skipped)
        "0.10", "0.0", "9.9",                              # gutz1 -> 0.10 + 0im
        "0.20", "0.0", "9.9",                              # gutz2 -> 0.20 + 0im
        "0.30", "0.0", "9.9",                              # jastrow1 -> 0.30 + 0im
        "0.40", "-0.10", "9.9",                            # slater idx0 -> 0.40 - 0.10im
        "0.50", "-0.20", "9.9",                            # slater idx1 -> 0.50 - 0.20im
    ],
    " ",
) * " \n"

@testset "read_opt_para_file!: golden load + consumed count" begin
    mktempdir() do dir
        path = joinpath(dir, "zqp_opt.dat")
        write(path, _GOLDEN_OPT)
        data = _make_data()
        n = MVMCOptimizers.read_opt_para_file!(data, path)

        @test n == 5  # n_proj (3) + n_slater (2)
        @test data.gutzwiller_terms[1].value ≈ 0.10 + 0im
        @test data.gutzwiller_terms[2].value ≈ 0.20 + 0im
        @test data.jastrow_terms[1].value ≈ 0.30 + 0im
        # Slater values scattered by orbital idx (idx0 -> first triple, idx1 -> second).
        @test data.orbital_terms[1].value ≈ 0.40 - 0.10im
        @test data.orbital_terms[2].value ≈ 0.50 - 0.20im
    end
end

@testset "read_opt_para_file!: non-perturbing (idempotent) load" begin
    mktempdir() do dir
        path = joinpath(dir, "zqp_opt.dat")
        write(path, _GOLDEN_OPT)
        d1 = _make_data()
        MVMCOptimizers.read_opt_para_file!(d1, path)
        d2 = _make_data()
        MVMCOptimizers.read_opt_para_file!(d2, path)
        MVMCOptimizers.read_opt_para_file!(d2, path)  # load twice
        @test [t.value for t in d1.gutzwiller_terms] == [t.value for t in d2.gutzwiller_terms]
        @test [t.value for t in d1.orbital_terms] == [t.value for t in d2.orbital_terms]
    end
end

@testset "read_opt_para_file!: strict failures error with a clear message" begin
    mktempdir() do dir
        # Missing file.
        threw, msg = _capture_msg(
            () -> MVMCOptimizers.read_opt_para_file!(_make_data(), joinpath(dir, "nope.dat")),
        )
        @test threw && occursin("file not found", msg)

        # Empty file -> too short.
        p = joinpath(dir, "empty.dat")
        write(p, "")
        threw, msg = _capture_msg(() -> MVMCOptimizers.read_opt_para_file!(_make_data(), p))
        @test threw && occursin("too short", msg)

        # Short (drop the last triple).
        p = joinpath(dir, "short.dat")
        write(p, join(split(_GOLDEN_OPT)[1:end-3], " "))
        threw, msg = _capture_msg(() -> MVMCOptimizers.read_opt_para_file!(_make_data(), p))
        @test threw && occursin("too short", msg)

        # One trailing stray float (not a whole triple).
        p = joinpath(dir, "trail1.dat")
        write(p, strip(_GOLDEN_OPT) * " 0.123")
        threw, msg = _capture_msg(() -> MVMCOptimizers.read_opt_para_file!(_make_data(), p))
        @test threw && occursin("trailing floats", msg)

        # A whole extra triple -> OptTrans-style block (unsupported).
        p = joinpath(dir, "trail3.dat")
        write(p, strip(_GOLDEN_OPT) * " 0.1 0.2 0.3")
        threw, msg = _capture_msg(() -> MVMCOptimizers.read_opt_para_file!(_make_data(), p))
        @test threw && occursin("OptTrans", msg)

        # Non-numeric token -> stable error (not a raw ArgumentError leak).
        p = joinpath(dir, "garbled.dat")
        write(p, replace(_GOLDEN_OPT, "0.30" => "NaNsense"))
        threw, msg = _capture_msg(() -> MVMCOptimizers.read_opt_para_file!(_make_data(), p))
        @test threw && occursin("non-numeric token", msg)
    end
end

@testset "read_opt_para_file!: DoublonHolon models are rejected (scope guard)" begin
    mktempdir() do dir
        path = joinpath(dir, "zqp_opt.dat")
        write(path, _GOLDEN_OPT)
        data = _make_data()
        data.doublon_holon_2site_terms = [DoublonHolon2SiteTerm(0, 1, 0.0 + 0im, true)]
        threw, msg = _capture_msg(() -> MVMCOptimizers.read_opt_para_file!(data, path))
        @test threw && occursin("DoublonHolon", msg)
    end
end

@testset "read_initial_def! still loads the same layout (delegation regression)" begin
    mktempdir() do dir
        path = joinpath(dir, "initial.def")
        write(path, _GOLDEN_OPT)
        data = _make_data()
        @test MVMCOptimizers.read_initial_def!(data, path) === true
        @test data.gutzwiller_terms[1].value ≈ 0.10 + 0im
        @test data.orbital_terms[2].value ≈ 0.50 - 0.20im
        # Missing file is a recoverable (warn+false) problem, not an error.
        @test MVMCOptimizers.read_initial_def!(_make_data(), joinpath(dir, "nope.def")) === false
    end
end
