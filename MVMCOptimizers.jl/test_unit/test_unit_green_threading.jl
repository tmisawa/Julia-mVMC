using Test, MVMCOptimizers, LinearAlgebra, Random
using MVMCExpertModeParsers: MVMCExpertModeParsers, ExpertModeData, ModParaParameters, GreenOneTerm, GreenTwoTerm, GreenTwoExTerm

@testset "Cyclic direct Green workers preserve each operator and ordered sample accumulation" begin
    M=MVMCOptimizers
    sites,nelec,qps=64,32,8
    data=ExpertModeData()
    data.modpara=ModParaParameters(nsite=sites,nelec=nelec,ncond=2nelec,nsp_gauss_leg=qps,nmp_trans=1)
    data.green_one_terms=[GreenOneTerm(nelec,0,:up,:up),GreenOneTerm(0,0,:down,:down)]
    terms=[GreenTwoTerm(0,0,1,1,:up,:up,:down,:down),
           GreenTwoTerm(nelec,0,nelec+1,1,:up,:up,:down,:down),
           GreenTwoTerm(nelec,0,1,1,:up,:up,:down,:down),
           GreenTwoTerm(nelec,0,0,nelec,:up,:up,:up,:up)]
    data.green_two_terms=repeat(terms,32)
    data.green_two_ex_terms=[GreenTwoExTerm(nelec,0,0,0,0,1,0,1)]
    MVMCExpertModeParsers.init_qp_weight!(data)
    # This fixed-state operator test supplies positive QP weights directly;
    # no variational OptTrans parameters are initialized by a public runner.
    fill!(data.qp_weights.qp_full_weight,ComplexF64(1/qps))
    state=M.VMCOptimizationState(sites,nelec,0,0,qps,1,false,false)
    generator=MersenneTwister(19)
    raw=randn(generator,2sites,2sites); matrix=raw-transpose(raw)
    for q in 1:qps
        state.slater_matrix.slater_elm_real[(q-1)*(2sites)^2+1:q*(2sites)^2] .= vec(permutedims(matrix)).*(1+q/16)
    end
    indices=vcat(collect(0:nelec-1),collect(0:nelec-1))
    cfg=fill(-1,2sites)
    for spin in 0:1, electron in 0:nelec-1
        cfg[spin*sites+electron+1]=electron
    end
    num=Int.(cfg.>=0)
    M.initialize_phys_quantities!(state,data)
    @test M.calculate_m_all_real!(indices,1,qps+1,data,state)==0
    ip=ComplexF64(M.calculate_ip_real(state.slater_matrix.pf_m_real,1,qps+1,data))
    @test isfinite(ip) && !iszero(ip)
    reference=deepcopy(state)
    withenv("JULIA_MVMC_INNER_THREADS"=>"0") do
        for _ in 1:2
            M.calculate_green_func!(data,reference,1.0,ip,indices,cfg,num,Int[])
        end
    end
    withenv("JULIA_MVMC_INNER_THREADS"=>"1") do
        for repeat in 1:4
            actual=deepcopy(state)
            for _ in 1:2
                M.calculate_green_func!(data,actual,1.0,ip,indices,cfg,num,Int[])
            end
            for field in (:local_cis_ajs,:phys_cis_ajs,:local_cis_ajs_ckt_alt_dc,:phys_cis_ajs_ckt_alt_dc,:phys_cis_ajs_ckt_alt)
                @test isapprox(getfield(actual.phys_quantities,field),getfield(reference.phys_quantities,field);atol=1e-11,rtol=1e-11)
            end
        end
    end
    @test indices==vcat(collect(0:nelec-1),collect(0:nelec-1))
    @test num==Int.(cfg.>=0)
end
