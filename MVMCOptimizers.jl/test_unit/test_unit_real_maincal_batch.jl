using Test, MVMCOptimizers, LinearAlgebra, Random
using MVMCExpertModeParsers: ExpertModeData, ModParaParameters

@testset "Saved-sample real batches preserve independent inverses, tails and failure publication" begin
    M=MVMCOptimizers
    for n in (4,32,128), qps in (2,8)
        sites=n÷2+1
        data=ExpertModeData()
        data.modpara=ModParaParameters(nsite=sites,nelec=n÷2)
        state=M.VMCOptimizationState(sites,n÷2,0,0,qps,3,false,false)
        rng=MersenneTwister(19)
        raw=randn(rng,2sites,2sites)
        matrix=raw-transpose(raw)
        for q in 1:qps
            state.slater_matrix.slater_elm_real[(q-1)*(2sites)^2+1:q*(2sites)^2] .= vec(permutedims(matrix))
        end
        for sample in 0:2
            indices=vcat(collect(0:n÷2-1),collect(0:n÷2-1))
            sample>0 && (indices[sample]=sites-1)
            state.electron_config.ele_idx[sample*n+1:(sample+1)*n] .= indices
        end
        batch=M.VMCRealMainCalBatch(n,qps,2)
        for first_sample in (0,2,0)
            count=min(2,3-first_sample)
            M._prepare_real_maincal_batch!(batch,first_sample,count,data,state)
            for slot in 1:count
                @test M._publish_real_maincal_batch!(batch,slot,state)==0
                sample=first_sample+slot-1
                indices=state.electron_config.ele_idx[sample*n+1:(sample+1)*n]
                selected=vcat(indices[1:n÷2].+1,indices[n÷2+1:n].+sites.+1)
                operator=matrix[selected,selected]
                for q in 1:qps
                    computed=transpose(reshape(state.slater_matrix.inv_m_real[(q-1)*n^2+1:q*n^2],n,n))
                    residual=opnorm(operator*computed-I,Inf)/(opnorm(operator,Inf)*opnorm(computed,Inf)+1)
                    @test residual<=256n*eps(Float64)
                end
            end
        end
        fill!(state.slater_matrix.slater_elm_real,0)
        M._prepare_real_maincal_batch!(batch,0,2,data,state)
        fill!(state.slater_matrix.inv_m_real,17)
        fill!(state.slater_matrix.pf_m_real,19)
        @test M._publish_real_maincal_batch!(batch,1,state)!=0
        @test all(==(17),state.slater_matrix.inv_m_real)
        @test all(==(19),state.slater_matrix.pf_m_real)
    end
end
