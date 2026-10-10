# Optional diagnostic: install a typed SR-store gate selector before all warmups.
# Both variants use the production API and identical arithmetic; no timed eval/JIT.
using MPI, MVMCOptimizers, LinearAlgebra, Libdl, PfaPack, SHA
const O=MVMCOptimizers
BLAS.set_num_threads(1)
namelist, steps_s, reps_s, output, ranks_s = ARGS
steps,reps,ranks=parse.(Int,(steps_s,reps_s,ranks_s))
@assert steps==300 && reps==3
provided=MPI.Init_thread(MPI.THREAD_FUNNELED)
@assert provided>=MPI.THREAD_FUNNELED && MPI.Is_thread_main()
comm=MPI.COMM_WORLD
rank=MPI.Comm_rank(comm)
@assert MPI.Comm_size(comm)==ranks
@assert BLAS.get_num_threads()==1
native=Libdl.dlopen(PfaPack.libltl2inv)
@assert ccall(Libdl.dlsym(native,:openblas_get_num_threads),Cint,())==1
@eval O begin
    const sr_store_gate_candidate=Ref{Bool}(false)
    @inline function vmc_sr_store_threading_enabled(work_items::Integer,threaded::Bool,::Type{T}) where T
        threshold = sr_store_gate_candidate[] ?
            (T <: Complex ? VMC_SR_STORE_MIN_COMPLEX_ITEMS : VMC_SR_STORE_MIN_REAL_ITEMS) : 64
        return vmc_inner_threading_enabled(work_items,threaded;min_work_per_thread=threshold)
    end
end
if rank==0
    mkpath(output)
    write(joinpath(output,"gate-override.txt"),"Typed Ref selector installed before warmups; baseline threshold=64, candidate thresholds=65536 real/32768 complex.\n")
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
        O.sr_store_gate_candidate[]=candidate
        MPI.Barrier(comm)
        result=production(candidate ? "warm-candidate" : "warm-baseline")
        @assert result.status==0 && result.effective_nsteps==steps
        rank==0 && @assert length(result.zvo_first_n)==steps
        MPI.Barrier(comm)
    end
    for rep in 1:reps
        for candidate in (isodd(rep) ? (false,true) : (true,false))
            O.sr_store_gate_candidate[]=candidate
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
