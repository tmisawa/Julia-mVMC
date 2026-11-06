using Test
using MVMCExpertModeParsers

# C mVMC 1.4.0 (v1.4.0, e83aa6cb) reader contracts adopted for Julia-mVMC#95.

@testset "ModPara SR-CG keys and C 1.4.0 normalization" begin
    base = """
    --------------------
    Model_Parameters   0
    --------------------
    Nsite          4
    """
    parse(extra) = MVMCExpertModeParsers.parse_modpara_content(base * extra).data
    # SetDefaultValuesModPara (readdef.c:2130-2136).
    m = parse("")
    @test (m.nsrcg, m.use_diag_scale, m.rescale_smat) == (0, 0, 0)
    @test (m.nsrcg_fallback, m.nsrcg_abort_on_fail, m.reweight) == (0, 1, 0)
    # readdef.c:927-944: NSRCG 2 -> DiagScale CG; flags become 0/1.
    m = parse("NSRCG 2\nNSRCGFallback 5\nNSRCGAbortOnFail 0\nreweight 1\n")
    @test (m.nsrcg, m.use_diag_scale) == (1, 1)
    @test (m.nsrcg_fallback, m.nsrcg_abort_on_fail, m.reweight) == (1, 0, 1)
    m = parse("NSRCG 1\nuseDiagScale 1\nRescaleSmat 1\nNSRCGAbortOnFail 7\n")
    @test (m.nsrcg, m.use_diag_scale, m.rescale_smat, m.nsrcg_abort_on_fail) == (1, 1, 1, 1)
    # useDiagScale with NSRCG != 1 is cleared ("Use direct method instead").
    @test parse("NSRCG 0\nuseDiagScale 1\n").use_diag_scale == 0
    @test parse("NSRCG 3\nuseDiagScale 1\n").use_diag_scale == 0
end

function jastrow_def(n, rows, opts)
    header = """
    =============================================
    NJastrowIdx          $n
    ComplexType          0
    =============================================
    =============================================
    """
    return header * join(("$i $j $c" for (i, j, c) in rows), "\n") * "\n" *
           join(("$k $f" for (k, f) in opts), "\n") * "\n"
end

@testset "strict Jastrow table (C 1.4.0 GetInfoJastrow)" begin
    nsite = 3
    full(classes) = [(i, j, classes[minmax(i, j)]) for i = 0:2 for j = 0:2 if i != j]
    classes = Dict((0, 1) => 0, (0, 2) => 1, (1, 2) => 0)
    ok = jastrow_def(2, full(classes), [(0, 1), (1, 1)])
    table = MVMCExpertModeParsers.jastrow_idx_matrix_c_contract(ok, nsite, 2)
    @test table == [-1 0 1; 0 -1 0; 1 0 -1]
    bad(content, label) = begin
        err = try
            MVMCExpertModeParsers.jastrow_idx_matrix_c_contract(content, nsite, 2)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        err isa ErrorException && @test occursin(label, err.msg)
    end
    rows = full(classes)
    # (1,0) with a different class from (0,1).
    r = copy(rows); r[3] = (1, 0, 1)
    bad(jastrow_def(2, r, [(0, 1), (1, 1)]), "inconsistent class for pair (0,1)")
    # Duplicated directed pair (and therefore one missing).
    r = copy(rows); r[2] = (0, 1, 0)
    bad(jastrow_def(2, r, [(0, 1), (1, 1)]), "duplicated pair (0,1)")
    r = copy(rows); r[1] = (1, 1, 0)
    bad(jastrow_def(2, r, [(0, 1), (1, 1)]), "i neq j")
    r = copy(rows); r[1] = (0, 3, 0)
    bad(jastrow_def(2, r, [(0, 1), (1, 1)]), "Site index is incorrect")
    r = copy(rows); r[1] = (0, 1, 2)
    bad(jastrow_def(2, r, [(0, 1), (1, 1)]), "class index 2 out of range")
    bad(jastrow_def(2, rows[1:5], []), "failed to read Jastrow row")
    bad(jastrow_def(2, rows, [(0, 1)]), "1 opt-flag rows, expected NJastrowIdx = 2")
end

@testset "parse_expert_mode_files rejects a duplicated Jastrow pair" begin
    src = joinpath(@__DIR__, "samples", "HeisenbergChain")
    mktempdir() do dir
        for f in readdir(src)
            cp(joinpath(src, f), joinpath(dir, f))
        end
        @test parse_expert_mode_files(joinpath(dir, "namelist.def")) isa MVMCExpertModeParsers.ExpertModeData
        path = joinpath(dir, "jastrowidx.def")
        lines = readlines(path)
        k = findfirst(l -> split(l) == ["0", "1", "0"], lines)
        @test k !== nothing
        if k !== nothing
            lines[k] = "    0      2      0"  # (0,2) twice, (0,1) missing
            write(path, join(lines, "\n") * "\n")
            @test_throws ErrorException parse_expert_mode_files(joinpath(dir, "namelist.def"))
        end
    end
end
