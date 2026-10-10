# Standalone dispatch crossover; not a whole-VMC performance result.
using MVMCOptimizers, LinearAlgebra, PfaPack, Libdl
const O=MVMCOptimizers
BLAS.set_num_threads(1)
@assert Threads.nthreads()==16
native=Libdl.dlopen(PfaPack.libltl2inv)
@assert ccall(Libdl.dlsym(native,:openblas_get_num_threads),Cint,())==1
@eval O const real_qp_budget_candidate=Ref{Bool}(false)
source=read(joinpath(dirname(pathof(O)),"calculate_m_all.jl"),String)
method=strip(source[first(findlast("function calculate_m_all_real!(",source)):end])
@assert count("qp_num < 2",method)==1
Core.eval(O,Meta.parse(replace(method,"qp_num < 2"=>
    "qp_num < (real_qp_budget_candidate[] ? 2 : n_threads)")))

function case(n,qps,gap_ns)
    slater=Float64[]; operators=Matrix{Float64}[]
    for q in 1:qps
        A=zeros(n,n)
        for j in 2:n, i in 1:j-1
            A[i,j]=isodd(i) && j==i+1 ? 2+q/16 : (i-j+q)/1024
            A[j,i]=-A[i,j]
        end
        push!(operators,A)
        append!(slater,vec(permutedims(A)))
    end
    indices=vcat(collect(0:n÷2-1),collect(0:n÷2-1))
    inverse=zeros(n,n,qps); pfaffian=zeros(qps)
    workspace=O.ThreadedPfaPackWorkspace(n;real_only=true)
    production()=O.calculate_m_all_real!(indices,slater,inverse,pfaffian,1,qps+1,n÷2,n,workspace)
    for candidate in (false,true)
        O.real_qp_budget_candidate[]=candidate
        for _ in 1:100
            @assert production()==0
        end
    end
    for rep in 1:3
        for candidate in (isodd(rep) ? (false,true) : (true,false))
            O.real_qp_budget_candidate[]=candidate
            total=UInt64(0)
            for _ in 1:100
                deadline=time_ns()+gap_ns
                while time_ns()<deadline
                end
                start=time_ns()
                info=production()
                total+=time_ns()-start
                @assert info==0
            end
            residual=0.0
            for q in 1:qps
                A=operators[q]; X=transpose(inverse[:,:,q])
                residual=max(residual,opnorm(A*X-I,Inf)/(opnorm(A,Inf)*opnorm(X,Inf)+1))
            end
            @assert residual<=256n*eps(Float64)
            println("CROSSOVER $n $qps $gap_ns $candidate $rep $(total/100) $residual")
            flush(stdout)
        end
    end
end

for gap in (0,100_000,1_000_000), n in (24,32,40,48,64,96)
    case(n,8,gap)
end
