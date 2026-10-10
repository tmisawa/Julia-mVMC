using Test, MVMCOptimizers, LinearAlgebra
@testset "Ordinary real QP threading preserves independent matrix operations" begin
    BLAS.set_num_threads(1)
    for n in (4, 16, 32, 64), qps in (0, 1, 2, 3, 4, 7, 8, 17)
        slater = Float64[]
        operators = Matrix{Float64}[]
        for q in 1:qps
            A = zeros(n,n)
            for j in 2:n, i in 1:j-1
                A[i,j] = isodd(i) && j == i+1 ? 2.0 + q/16 : (i-j+q)/1024
                A[j,i] = -A[i,j]
            end
            push!(operators,A)
            append!(slater,vec(permutedims(A)))
        end
        indices = vcat(collect(0:n÷2-1),collect(0:n÷2-1))
        serial = zeros(n,n,qps); parallel = similar(serial)
        ps = zeros(qps); pp = similar(ps)
        sw = MVMCOptimizers.PfaPackWorkspace(n; real_only=true)
        tw = MVMCOptimizers.ThreadedPfaPackWorkspace(n; real_only=true)
        @test length(tw.workspaces) >= Base.Threads.maxthreadid()
        for reuse in 1:20
            fill!(parallel,17.0); fill!(pp,17.0)
            @test MVMCOptimizers.calculate_m_all_real!(indices,slater,serial,ps,1,qps+1,n÷2,n,sw) == 0
            @test MVMCOptimizers.calculate_m_all_real!(indices,slater,parallel,pp,1,qps+1,n÷2,n,tw) == 0
            @test all(abs(a-b) <= 256eps(Float64)*(1+abs(b)) for (a,b) in zip(parallel,serial))
            @test all(abs(a-b) <= 256eps(Float64)*(1+abs(b)) for (a,b) in zip(pp,ps))
            for q in 1:qps
                A = operators[q]; inverse = transpose(parallel[:,:,q])
                residual = opnorm(A*inverse-I,Inf)/(opnorm(A,Inf)*opnorm(inverse,Inf)+1)
                @test residual <= 256n*eps(Float64)
            end
        end
    end
end
@testset "Complex QP threaded workspace remains private" begin
    for n in (4,16,32), qps in (1,4,8)
        slater = ComplexF64[]; operators = Matrix{ComplexF64}[]
        for q in 1:qps
            A = zeros(ComplexF64,n,n)
            for j in 2:n, i in 1:j-1
                A[i,j] = (isodd(i) && j==i+1 ? 2.0+q/16 : (i-j+q)/1024)*(1+im/8)
                A[j,i] = -A[i,j]
            end
            push!(operators,A); append!(slater,vec(permutedims(A)))
        end
        indices = vcat(collect(0:n÷2-1),collect(0:n÷2-1))
        serial=zeros(ComplexF64,n,n,qps); parallel=similar(serial)
        ps=zeros(ComplexF64,qps); pp=similar(ps)
        sw=MVMCOptimizers.PfaPackWorkspace(n;complex_only=true)
        tw=MVMCOptimizers.ThreadedPfaPackWorkspace(n;complex_only=true)
        for reuse in 1:20
            fill!(parallel,17); fill!(pp,17)
            @test MVMCOptimizers.calculate_m_all_fcmp!(indices,slater,serial,ps,1,qps+1,n÷2,n,sw)==0
            @test MVMCOptimizers.calculate_m_all_fcmp!(indices,slater,parallel,pp,1,qps+1,n÷2,n,tw)==0
            @test all(abs(a-b)<=256eps(Float64)*(1+abs(b)) for (a,b) in zip(parallel,serial))
            @test all(abs(a-b)<=256eps(Float64)*(1+abs(b)) for (a,b) in zip(pp,ps))
            for q in 1:qps
                A=operators[q]; inverse=transpose(parallel[:,:,q])
                @test opnorm(A*inverse-I,Inf)/(opnorm(A,Inf)*opnorm(inverse,Inf)+1)<=256n*eps(Float64)
            end
        end
    end
end
using MVMCExpertModeParsers: ExpertModeData, ModParaParameters
@testset "SIMD real accepted hop satisfies independent inverse residual" begin
    for n in (4,16,32,64), reuse in 1:3
        sites = n÷2+1
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite=sites,nelec=n÷2)
        state = MVMCOptimizers.VMCOptimizationState(sites,n÷2,0,0,1,1,false,false)
        full = zeros(2sites,2sites)
        for j in 2:2sites, i in 1:j-1
            full[i,j] = j==i+sites && i<=sites ? 3.0+(i+reuse)/16 : (i-j+reuse)/1024
            full[j,i] = -full[i,j]
        end
        full[sites,sites+1] = 2.0+reuse/8
        full[sites+1,sites] = -full[sites,sites+1]
        # Rank-2 skew matrix updates are valid for arbitrary invertible skew
        # operators; site/spin offsets select the physical electron submatrix.
        indices = vcat(collect(0:n÷2-1),collect(0:n÷2-1))
        selected = vcat(indices[1:n÷2].+1,indices[n÷2+1:n].+sites.+1)
        old = full[selected,selected]
        # A tiny independent solve provides the initial inverse. Its residual
        # is checked, then the changed operator is checked independently again.
        oldinverse = inv(old)
        state.slater_matrix.inv_m_real[1:n*n] .= vec(permutedims(oldinverse))
        state.slater_matrix.slater_elm_real .= vec(permutedims(full))
        state.slater_matrix.pf_m_real[1] = 1.0
        indices[1] = sites-1
        selected[1] = sites
        MVMCOptimizers.update_m_all_real!(0,0,indices,1,2,data,state)
        actual = transpose(reshape(state.slater_matrix.inv_m_real[1:n*n],n,n))
        changed = full[selected,selected]
        residual = opnorm(changed*actual-I,Inf)/(opnorm(changed,Inf)*opnorm(actual,Inf)+1)
        @test residual <= 256n*eps(Float64)
        @test isfinite(state.slater_matrix.pf_m_real[1])
    end
end
@testset "Worker and nested QP callers retain serial scratch ownership" begin
    n=4; qps=8; indices=[0,1,0,1]
    slater=Float64[]
    for q in 1:qps
        A=[0.0 2+q/16 .125 0; -(2+q/16) 0 0 -.125; -.125 0 0 3+q/16; 0 .125 -(3+q/16) 0]
        append!(slater,vec(permutedims(A)))
    end
    reference=zeros(n,n,qps); pref=zeros(qps)
    @test MVMCOptimizers.calculate_m_all_real!(indices,slater,reference,pref,1,qps+1,2,n,MVMCOptimizers.PfaPackWorkspace(n;real_only=true))==0
    tw=MVMCOptimizers.ThreadedPfaPackWorkspace(n;real_only=true)
    errors=Vector{Tuple{Int,Float64}}(undef,Threads.nthreads())
    Threads.@threads :static for slot in 1:Threads.nthreads()
        actual=zeros(n,n,qps); pf=zeros(qps)
        status=MVMCOptimizers.calculate_m_all_real!(indices,slater,actual,pf,1,qps+1,2,n,tw)
        errors[slot]=(status,maximum(abs.(actual-reference)))
    end
    @test all(status==0 && error<=256eps(Float64) for (status,error) in errors)
    status,error=fetch(Threads.@spawn begin
        actual=zeros(n,n,qps); pf=zeros(qps)
        status=MVMCOptimizers.calculate_m_all_real!(indices,slater,actual,pf,1,qps+1,2,n,tw)
        (status,maximum(abs.(actual-reference)))
    end)
    @test status==0 && error<=256eps(Float64)
end

@testset "Unit-upper inverse preserves ignored storage and satisfies residual" begin
    for n in (1, 3, 15, 31, 63), reuse in 1:3
        parent = fill(-91.0, n+2, n+2)
        A = view(parent, 2:n+1, 2:n+1)
        operator = Matrix{Float64}(I, n, n)
        for j in 1:n, i in 1:n
            A[i,j] = i < j ? (i-j+reuse)/1024 : NaN
            if i < j
                operator[i,j] = A[i,j]
            end
        end
        expected = inv(operator)
        MVMCOptimizers._ordinary_unit_upper_trtri_real!(A)
        for j in 1:n, i in 1:n
            if i < j
                @test abs(A[i,j]-expected[i,j]) <= 256eps(Float64)*(1+abs(expected[i,j]))
            else
                @test isnan(A[i,j])
            end
        end
        inverse = Matrix{Float64}(I, n, n)
        for j in 2:n, i in 1:j-1
            inverse[i,j] = A[i,j]
        end
        residual = opnorm(operator*inverse-I, Inf)/(opnorm(operator,Inf)*opnorm(inverse,Inf)+1)
        @test residual <= 256n*eps(Float64)
        @test all(==(-91.0), parent[[1,n+2],:]) && all(==(-91.0), parent[:,[1,n+2]])
    end
end

@testset "Real QP partial pools require sufficient matrix work" begin
    # Explicit pool-size boundaries for these matrix/QP counts. Integer pivot
    # sentinels distinguish actual worker execution from a requested capacity.
    cases=((32,2,2),(32,3,3),(32,8,8),
           (48,2,6),(48,3,10),(48,8,27),
           (64,2,16),(64,3,24),(64,8,64))
    for (n,qps,last_parallel_pool_size) in cases
        indices=vcat(collect(0:n÷2-1),collect(0:n÷2-1))
        slater=Float64[]
        for q in 1:qps
            A=zeros(n,n)
            # Dense matrix avoids zero intermediate columns: C DSKTF2 normal
            # mode reports INFO>0 for those even if a block Pfaffian is nonzero.
            for j in 2:n, i in 1:j-1
                A[i,j]=isodd(i) && j==i+1 ? 2+q/16 : (i-j+q)/1024
                A[j,i]=-A[i,j]
            end
            append!(slater,vec(permutedims(A)))
        end
        workspace=MVMCOptimizers.ThreadedPfaPackWorkspace(n;real_only=true)
        for scratch in workspace.workspaces
            fill!(scratch.iwork,-1)
        end
        inverse=zeros(n,n,qps); pfaffian=zeros(qps)
        @test MVMCOptimizers.calculate_m_all_real!(
            indices,slater,inverse,pfaffian,1,qps+1,n÷2,n,workspace)==0
        used=count(scratch->all(>(0),scratch.iwork),workspace.workspaces)
        expected=Threads.nthreads()<=last_parallel_pool_size ? min(qps,Threads.nthreads()) : 1
        @test used==expected
    end
end
