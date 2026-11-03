# julia +1.13 --project=. --threads=N benchmark/cpu_l128/measure.jl INPUTS NEW_OUTPUT [STEPS REPS opt|phys|both]
using MVMCOptimizers, LinearAlgebra, SHA, TOML, Statistics, InteractiveUtils, Libdl

input, output = abspath.(ARGS[1:2])
steps = length(ARGS)>=3 ? parse(Int,ARGS[3]) : 20
reps = length(ARGS)>=4 ? parse(Int,ARGS[4]) : 3
workloads = length(ARGS)>=5 && ARGS[5]!="both" ? [ARGS[5]] : ["opt","phys"]
@assert steps>0 && reps>0 && all(in(("opt","phys")),workloads)
ispath(output) && error("output must be a fresh directory: $output")
BLAS.set_num_threads(1)
mkpath(output)
cp(input,joinpath(output,"inputs"))
repo = dirname(pkgdir(MVMCOptimizers))
hashes = Dict{String,String}()
for (root,_,files) in walkdir(joinpath(output,"inputs")), file in files
    path=joinpath(root,file)
    hashes[relpath(path,joinpath(output,"inputs"))]=bytes2hex(sha256(read(path)))
end
metadata=Dict("julia_channel"=>"1.13","julia_version"=>string(VERSION),
    "driver_sha256"=>bytes2hex(sha256(read(@__FILE__))),
    "source"=>readchomp(`git -C $repo rev-parse HEAD`),
    "source_status"=>readchomp(`git -C $repo status --porcelain`),
    "source_diff"=>readchomp(`git -C $repo diff --binary`),
    "manifest_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Manifest-v1.13.toml")))),
    "threads"=>Threads.nthreads(),"interactive_threads"=>Threads.nthreads(:interactive),
    "blas"=>string(BLAS.get_config()),"blas_threads"=>BLAS.get_num_threads(),
    "libraries"=>filter(x->occursin(r"blas|lapack|sfmt|pfapack"i,x),Libdl.dllist()),
    "environment"=>Dict(k=>v for (k,v) in ENV if any(p->startswith(k,p),("JULIA_","OPENBLAS_","OMP_","MKL_","BLIS_","MVMC_"))),
    "input_sha256"=>hashes,"steps"=>steps,"repetitions"=>reps,
    "scope"=>"Warm public CPU API including parsing and file output; process/JIT/explicit GC excluded")
open(io->TOML.print(io,metadata),joinpath(output,"metadata.toml"),"w")
open(joinpath(output,"environment.txt"),"w") do io
    versioninfo(io;verbose=true)
    if Sys.islinux()
        println(io,read("/proc/self/status",String))
        println(io,read("/proc/loadavg",String))
        println(io,read(`lscpu`,String))
    end
end
function production(workload, destination)
    root=joinpath(output,"inputs",workload)
    result=if workload=="opt"
        run_para_opt_from_namelist(joinpath(root,"namelist.def");nsteps=steps,nsmp=1,
            mode=:real,seed=79,output_dir=destination)
    else
        run_phys_cal_from_namelist(joinpath(root,"namelist.def");
            opt_para=joinpath(root,"zqp_opt.dat"),mode=:real,seed=79,output_dir=destination)
    end
    @assert result.status==0
    return nothing
end
open(joinpath(output,"timings.csv"),"w") do io
    println(io,"workload,rep,seconds,bytes,gc_seconds")
    for workload in workloads
        times=Float64[]
        for rep in 0:reps
            GC.gc()
            stats=@timed production(workload,joinpath(output,"$workload-rep$rep"))
            println(io,"$workload,$rep,$(stats.time),$(stats.bytes),$(stats.gctime)")
            flush(io)
            rep>0 && push!(times,stats.time)
            println("BENCH workload=$workload rep=$rep seconds=$(stats.time)")
            flush(stdout)
        end
        println("MEDIAN workload=$workload seconds=$(median(times))")
    end
end
