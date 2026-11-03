using Test, MVMCOptimizers, MVMCExpertModeParsers, LinearAlgebra
using MVMCExpertModeParsers: ExpertModeData, ModParaParameters

@testset "Accepted real hops use private workers and retain the serial matrix operation" begin
    BLAS.set_num_threads(1)
    M = MVMCOptimizers
    for n in (16, 128), qps in (1, 8)
        sites = n ÷ 2 + 1
        data = ExpertModeData()
        data.modpara = ModParaParameters(nsite=sites, nelec=n÷2)
        state = M.VMCOptimizationState(sites,n÷2,0,0,qps,1,false,false)
        full = zeros(2sites,2sites)
        for j in 2:2sites, i in 1:j-1
            full[i,j] = j == i+sites ? 3.0+i/128 : (i-j)/8192
            full[j,i] = -full[i,j]
        end
        full[sites,sites+1] = 2.5
        full[sites+1,sites] = -2.5
        original = vcat(collect(0:n÷2-1),collect(0:n÷2-1))
        selected = vcat(original[1:n÷2].+1,original[n÷2+1:n].+sites.+1)
        inverse = inv(full[selected,selected])
        for q in 1:qps
            state.slater_matrix.slater_elm_real[(q-1)*(2sites)^2+1:q*(2sites)^2] .= vec(permutedims(full))
            state.slater_matrix.inv_m_real[(q-1)*n^2+1:q*n^2] .= vec(permutedims(inverse))
            state.slater_matrix.pf_m_real[q] = 1+q/16
        end
        moved = copy(original)
        moved[1] = sites-1
        selected[1] = sites
        reference = deepcopy(state)
        withenv("JULIA_MVMC_INNER_THREADS"=>"0") do
            M.update_m_all_real!(0,0,moved,1,qps+1,data,reference)
        end
        withenv("JULIA_MVMC_INNER_THREADS"=>"1") do
            for repeat in 1:4
                actual = deepcopy(state)
                for scratch in actual.workspace.hop_thread_vec1_real
                    fill!(scratch,NaN)
                end
                M.with_vmc_hop_team(actual.workspace) do
                    M.update_m_all_real!(0,0,moved,1,qps+1,data,actual)
                end
                @test isapprox(actual.slater_matrix.inv_m_real,reference.slater_matrix.inv_m_real;atol=1e-13,rtol=1e-13)
                @test isapprox(actual.slater_matrix.pf_m_real,reference.slater_matrix.pf_m_real;atol=1e-13,rtol=1e-13)
                used = count(s -> all(isfinite,s), actual.workspace.hop_thread_vec1_real)
                expected = Threads.nthreads()>1 && qps>1 ? min(qps,Threads.nthreads()) : 0
                @test used == expected
                for q in 1:qps
                    computed = transpose(reshape(actual.slater_matrix.inv_m_real[(q-1)*n^2+1:q*n^2],n,n))
                    operator = full[selected,selected]
                    residual = opnorm(operator*computed-I,Inf)/(opnorm(operator,Inf)*opnorm(computed,Inf)+1)
                    @test residual <= 256n*eps(Float64)
                end
            end
            # An already threaded caller must retain a private serial route.
            tasks = [Threads.@spawn begin
                actual=deepcopy(state)
                M.update_m_all_real!(0,0,moved,1,qps+1,data,actual)
                actual.slater_matrix.inv_m_real
            end for _ in 1:2]
            @test all(isapprox(fetch(task),reference.slater_matrix.inv_m_real;atol=1e-13,rtol=1e-13) for task in tasks)
        end
    end
end

using Random, SFMT

mutable struct RealHopTape <: AbstractRNG
    words::Vector{UInt32}
    cursor::Int
    owner::Int
    state::MVMCOptimizers.VMCOptimizationState
    events::Vector{Vector{Int}}
end
function Random.rand(rng::RealHopTape, ::Type{UInt32})
    Threads.threadid() == rng.owner || error("a worker consumed sampling RNG")
    push!(rng.events,vcat(rng.state.electron_config.tmp_ele_idx,rng.state.electron_config.counter))
    rng.cursor += 1
    return rng.words[rng.cursor]
end
Random.rand(rng::RealHopTape) = Float64(rand(rng,UInt32)) * (1.0/4294967296.0)

@testset "Real sampler preserves supplied tape, decisions and saved configurations" begin
    M = MVMCOptimizers
    sites, nelec, qps, samples = 128, 64, 8, 4
    data = ExpertModeData()
    data.modpara = ModParaParameters(nsite=sites,nelec=nelec,ncond=2nelec,
        nsp_gauss_leg=qps,nmp_trans=1,nvmc_sample=samples,nvmc_warmup=1,nvmc_interval=1)
    MVMCExpertModeParsers.init_qp_weight!(data)
    fill!(data.qp_weights.qp_full_weight,ComplexF64(1/qps))
    state = M.VMCOptimizationState(sites,nelec,0,0,qps,samples,false,false)
    generator = MersenneTwister(19)
    raw = randn(generator,2sites,2sites)
    matrix = raw-transpose(raw)
    for q in 1:qps
        range = (q-1)*(2sites)^2+1:q*(2sites)^2
        state.slater_matrix.slater_elm[range] .= vec(permutedims(matrix)) .* (1+q/16)
    end
    Random.seed!(SFMT19937RNG(),79)
    tape = SFMT.sfmt_dump_rand32(20000)
    reference = deepcopy(state)
    serial_rng = RealHopTape(tape,0,Threads.threadid(),reference,Vector{Int}[])
    withenv("JULIA_MVMC_INNER_THREADS"=>"0") do
        M.vmc_make_sample_real!(data,reference,serial_rng)
    end
    @test serial_rng.cursor > 1000
    @test reference.electron_config.counter[2] > 0
    withenv("JULIA_MVMC_INNER_THREADS"=>"1") do
        for repeat in 1:2
            actual = deepcopy(state)
            rng = RealHopTape(tape,0,Threads.threadid(),actual,Vector{Int}[])
            M.with_vmc_real_sampler_team(data,actual,M.serial_context()) do
                M.vmc_make_sample_real!(data,actual,rng)
            end
            @test rng.cursor == serial_rng.cursor
            @test rng.events == serial_rng.events
            for field in fieldnames(M.ElectronConfiguration)
                @test getfield(actual.electron_config,field) == getfield(reference.electron_config,field)
            end
            @test isapprox(actual.slater_matrix.inv_m_real,reference.slater_matrix.inv_m_real;atol=1e-11,rtol=1e-11)
            @test isapprox(actual.slater_matrix.pf_m_real,reference.slater_matrix.pf_m_real;atol=1e-11,rtol=1e-11)
            @test actual.workspace.hop_team === nothing
        end
    end
end

@testset "Persistent hop workers are joined on worker and caller errors" begin
    M=MVMCOptimizers
    state=M.VMCOptimizationState(4,2,0,0,8,1,false,false)
    if Threads.nthreads()>1
        @test_throws CapturedException M.with_vmc_hop_team(state.workspace) do
            M.vmc_hop_dispatch!(state.workspace.hop_team) do slot
                slot==1 && error("worker failure")
            end
        end
        @test state.workspace.hop_team === nothing
        @test_throws CapturedException M.with_vmc_hop_team(state.workspace) do
            M.vmc_hop_dispatch!(state.workspace.hop_team) do slot
                slot==0 && error("caller failure")
            end
        end
        @test state.workspace.hop_team === nothing
    end
end

@testset "Persistent worker admission preserves other proposal families" begin
    M=MVMCOptimizers
    state=M.VMCOptimizationState(128,64,0,0,8,1,false,false)
    withenv("JULIA_MVMC_INNER_THREADS"=>"1") do
        for path in (1,2,3)
            data=ExpertModeData()
            data.modpara=ModParaParameters(nsite=128,nelec=64,nex_update_path=path)
            M.with_vmc_real_sampler_team(data,state,M.serial_context()) do
                @test state.workspace.hop_team === nothing
            end
        end
    end
end
