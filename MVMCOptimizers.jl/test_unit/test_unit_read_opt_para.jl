using Test
using MVMCOptimizers
using MVMCExpertModeParsers:
    ExpertModeData,
    GutzwillerTerm,
    JastrowTerm,
    OrbitalTerm,
    DoublonHolon2SiteIndex,
    ChargeRBMPhysLayerTerm

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

const _GOLDEN_OPT_DH = join(
    [
        "1.0", "2.0", "3.0", "4.0", "5.0", "6.0",          # diagnostics
        "0.10", "0.0", "9.9",                              # gutz1
        "0.20", "0.0", "9.9",                              # gutz2
        "0.30", "0.0", "9.9",                              # jastrow1
        "1.10", "-0.10", "9.9",                            # DH2 local 0
        "1.20", "-0.20", "9.9",
        "1.30", "-0.30", "9.9",
        "1.40", "-0.40", "9.9",
        "1.50", "-0.50", "9.9",
        "1.60", "-0.60", "9.9",                            # DH2 local 5
        "0.40", "-0.10", "9.9",                            # slater idx0
        "0.50", "-0.20", "9.9",                            # slater idx1
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

        # NaN / Inf parse as Float64 but must be rejected for a reference gate.
        for tok in ("NaN", "Inf", "-Inf")
            p = joinpath(dir, "nonfinite.dat")
            write(p, replace(_GOLDEN_OPT, "0.40" => tok))
            threw, msg = _capture_msg(() -> MVMCOptimizers.read_opt_para_file!(_make_data(), p))
            @test threw && occursin("non-finite token", msg)
        end
    end
end

@testset "read_opt_para_file!: DH projection block loads before Slater" begin
    mktempdir() do dir
        path = joinpath(dir, "zqp_opt.dat")
        write(path, _GOLDEN_OPT_DH)
        data = _make_data()
        data.doublon_holon_2site_indices = [DoublonHolon2SiteIndex([1 0; 0 1])]
        data.doublon_holon_2site_params = fill(0.0 + 0.0im, 6)
        data.doublon_holon_2site_opt_flags = fill(true, 6)

        n = MVMCOptimizers.read_opt_para_file!(data, path)

        @test n == 11  # NProj (2 + 1 + 6) + NSlater (2)
        @test data.doublon_holon_2site_params == ComplexF64[
            1.10 - 0.10im,
            1.20 - 0.20im,
            1.30 - 0.30im,
            1.40 - 0.40im,
            1.50 - 0.50im,
            1.60 - 0.60im,
        ]
        @test data.orbital_terms[1].value ≈ 0.40 - 0.10im
        @test data.orbital_terms[2].value ≈ 0.50 - 0.20im
    end
end

@testset "read_opt_para_file!: OptTrans tail loads when active" begin
    mktempdir() do dir
        path = joinpath(dir, "zqp_opt.dat")
        write(path, strip(_GOLDEN_OPT) * " 0.70 -0.80 9.9\n")
        data = _make_data()
        data.para_qp_opt_trans = [1.0 + 0.0im]
        data.opt_trans = [1.0 + 0.0im]

        n = MVMCOptimizers.read_opt_para_file!(data, path)

        @test n == 6  # NProj (3) + NSlater (2) + NOptTrans (1)
        @test data.opt_trans == ComplexF64[0.70 - 0.80im]
        @test data.orbital_terms[2].value ≈ 0.50 - 0.20im
    end
end

@testset "read_opt_para_file!: RBM-bearing models still fail loud" begin
    mktempdir() do dir
        path = joinpath(dir, "zqp_opt.dat")
        write(path, _GOLDEN_OPT)
        data = _make_data()
        data.charge_rbm_phys_layer_terms = [
            ChargeRBMPhysLayerTerm(0, 0.0 + 0.0im, false, 0),
        ]
        threw, msg = _capture_msg(() -> MVMCOptimizers.read_opt_para_file!(data, path))
        @test threw && occursin("RBM-bearing", msg)
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
