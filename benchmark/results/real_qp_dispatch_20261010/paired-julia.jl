# Optional diagnostic: install a typed real-QP budget selector before all warmups.
# Both variants use the production API and identical arithmetic; no timed eval/JIT.
using MPI, MVMCOptimizers, MVMCExpertModeParsers, LinearAlgebra, Libdl, PfaPack, SHA
const O=MVMCOptimizers
BLAS.set_num_threads(1)
namelist, steps_s, reps_s, output, ranks_s = ARGS
steps,reps,ranks=parse.(Int,(steps_s,reps_s,ranks_s))
@assert steps in (20,300) && reps==3
modpara=MVMCExpertModeParsers.parse_expert_mode_files(namelist).modpara
@assert modpara.nsplit_size==1 && modpara.nvmc_sample*ranks==320
provided=MPI.Init_thread(MPI.THREAD_FUNNELED)
@assert provided>=MPI.THREAD_FUNNELED && MPI.Is_thread_main()
comm=MPI.COMM_WORLD
rank=MPI.Comm_rank(comm)
@assert MPI.Comm_size(comm)==ranks
@assert BLAS.get_num_threads()==1
native=Libdl.dlopen(PfaPack.libltl2inv)
@assert ccall(Libdl.dlsym(native,:openblas_get_num_threads),Cint,())==1
@eval O const real_qp_budget_candidate=Ref{Bool}(false)
kernel_source=read(joinpath(dirname(pathof(O)),"calculate_m_all.jl"),String)
method_start=first(findlast("function calculate_m_all_real!(",kernel_source))
method=strip(kernel_source[method_start:end])
@assert occursin("threaded_workspace::ThreadedPfaPackWorkspace",method)
@assert count("qp_num < 2",method)==1
Core.eval(O, Meta.parse(replace(method,"qp_num < 2"=>
    "qp_num < (real_qp_budget_candidate[] ? 2 : n_threads)")))
if rank==0
    mkpath(output)
    write(joinpath(output,"gate-override.txt"),"Typed Ref selector installed before warmups; baseline requires QP count >= thread count, candidate requires at least two QPs; other work gates unchanged.\n")
end
println("WORLD $rank $ranks\nTHREADS $rank $(Threads.nthreads())\nBLAS_THREADS $rank $(BLAS.get_num_threads())\nNATIVE_BLAS_THREADS $rank 1")
flush(stdout)
function production(label)
    @assert MPI.Is_thread_main()
    O.run_para_opt_from_namelist(namelist;nsteps=steps,nsmp=steps,mode=:real,
                                output_dir=joinpath(output,label))
end
try
    for candidate in (false,true)
        O.real_qp_budget_candidate[]=candidate
        MPI.Barrier(comm)
        result=production(candidate ? "warm-candidate" : "warm-baseline")
        @assert result.status==0 && result.effective_nsteps==steps
        rank==0 && @assert length(result.zvo_first_n)==steps
        MPI.Barrier(comm)
    end
    for rep in 1:reps
        for candidate in (isodd(rep) ? (false,true) : (true,false))
            O.real_qp_budget_candidate[]=candidate
            name=candidate ? "candidate" : "baseline"
            MPI.Barrier(comm)
            start=time_ns()
            result=production("$name-$rep")
            MPI.Barrier(comm)
            elapsed=(time_ns()-start)/1e9
            seconds=MPI.Allreduce(elapsed,max,comm)
            @assert result.status==0 && result.effective_nsteps==steps
            if rank==0
                @assert length(result.zvo_first_n)==steps && isfinite(result.final_energy_per_site)
                println("PAIRED $name $rep $seconds $(result.final_energy_per_site)")
                flush(stdout)
            end
        end
    end
catch err
    showerror(stderr,err,catch_backtrace())
    MPI.Abort(comm,1)
finally
    MPI.Finalize()
end
