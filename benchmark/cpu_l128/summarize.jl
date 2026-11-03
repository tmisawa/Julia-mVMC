# Recompute published medians and verify metadata from an extracted raw archive.
using TOML, Statistics
root=abspath(ARGS[1])
metadata_reference=nothing
cases=Dict{String,Any}[]
for threads in (1,4,8,16)
    rows=Dict{String,Any}()
    for side in ("before","after")
        directory=joinpath(root,"$side-t$threads")
        metadata=TOML.parsefile(joinpath(directory,"metadata.toml"))
        @assert metadata["threads"]==threads
        @assert metadata["steps"]==20 && metadata["repetitions"]==3
        @assert metadata["blas_threads"]==1
        @assert isempty(metadata["source_status"]) && isempty(metadata["source_diff"])
        @assert !any(k->occursin(r"token|secret|password|credential"i,k),keys(metadata["environment"]))
        if metadata_reference===nothing
            global metadata_reference=metadata
        else
            for key in ("input_sha256","manifest_sha256","driver_sha256","julia_version")
                @assert metadata[key]==metadata_reference[key] (directory,key)
            end
        end
        lines=readlines(joinpath(directory,"timings.csv"))
        @assert first(lines)=="workload,rep,seconds,bytes,gc_seconds"
        parsed=[split(line,',') for line in lines[2:end]]
        @assert length(parsed)==8
        for workload in ("opt","phys")
            group=filter(row->row[1]==workload,parsed)
            @assert sort(parse.(Int,[row[2] for row in group]))==collect(0:3)
            times=[parse(Float64,row[3]) for row in group if row[2]!="0"]
            @assert length(times)==3 && all(x->isfinite(x)&&x>0,times)
            rows["$side-$workload"]=Dict("median"=>median(times),"min"=>minimum(times),
                "max"=>maximum(times),"seconds"=>times)
        end
    end
    for workload in ("opt","phys")
        before,after=rows["before-$workload"],rows["after-$workload"]
        push!(cases,Dict("threads"=>threads,"workload"=>workload,"before"=>before,
            "after"=>after,"speedup"=>before["median"]/after["median"]))
    end
end
TOML.print(stdout,Dict("julia_version"=>metadata_reference["julia_version"],"cases"=>cases))
