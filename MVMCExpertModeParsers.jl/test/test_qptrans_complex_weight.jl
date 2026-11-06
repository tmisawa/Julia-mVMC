using Test
using MVMCExpertModeParsers

# Independent C GetInfoTransSym input contract. Actual native observation:
# 0 1.00000 0.50000 -> 1+0.5im under both APFlag values. Other dyadic
# coefficients are literal input controls, not generated Julia expectations.
# C authority: mVMC 1.4.0 readdef.c 2646–2679 (GetInfoTransSym, byte-identical
# to the original 1.3.0 acquisition at readdef.c 2232–2265, source SHA256
# 6c53cb832f93d6cbfd7cea955fbb693738af5536b913d36af32b98eed38c32d9);
# expert.rst TransSym specifies idx/real/imag. Standalone GCC13.3.0
# -std=c11 -O0 -ffp-contract=off acquired STATUS0 and WEIGHT
# 0x1p+0 0x1p-1 for both AP0/AP1; BLAS/RNG/MPI not used.
# QPTransTerm.phase is a separate scalar representation, NOT tested as a
# complex-coefficient API here. Historical RED used unchanged production loader.
@testset "public TransSym complex coefficient loader" begin
    for nmp in (1, -1)
        for (columns, expected) in (
            ("1.00000 0.50000", ComplexF64(1.0, 0.5)),
            ("-1.25000 -0.50000", ComplexF64(-1.25, -0.5)),
            ("1.00000", ComplexF64(1.0, 0.0)),
            ("1.00000 0.00000", ComplexF64(1.0, 0.0)),
        )
            @testset "NMP=$nmp coefficient=$columns" begin
                mktempdir() do dir
                    write(joinpath(dir, "modpara.def"),
                        "===\nModel_Parameters\n===\n===\n===\n" *
                        "Nsite 4\nNe 1\nNLocSpin 0\nNCond -1\n" *
                        "NMPTrans $nmp\nNSPGaussLeg 1\n")
                    write(joinpath(dir, "qptransidx.def"),
                        "===\nNQPTrans 1\n===\nTrIdx_TrWeight_and_TrIdx_i_xi\n===\n" *
                        "0 $columns\n0 0 1 -1\n0 1 2 1\n0 2 3 -1\n0 3 0 1\n")
                    namelist = joinpath(dir, "namelist.def")
                    write(namelist, "ModPara modpara.def\nTransSym qptransidx.def\n")
                    data = MVMCExpertModeParsers.parse_expert_mode_files(namelist)
                    @test data.modpara.nsite == 4
                    @test data.modpara.nmp_trans == nmp
                    @test data.n_qp_trans == 1
                    @test data.qp_trans == [[1, 2, 3, 0]]
                    @test data.qp_trans_sgn == [nmp < 0 ? [-1, 1, -1, 1] : [1, 1, 1, 1]]
                    @test length(data.para_qp_trans) == 1
                    @test data.para_qp_trans == [expected]
                end
            end
        end
    end
end
