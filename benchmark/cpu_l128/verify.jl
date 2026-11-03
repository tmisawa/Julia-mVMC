# Numerical output comparison; raw RNG/tape checks are a separate contract.
using TOML
before, after = abspath.(ARGS[1:2])
function output_names(path)
    sort(filter(name->endswith(name,".dat") &&
        any(prefix->startswith(name,prefix),("zvo_out","zvo_var","zqp_opt","zvo_cisajs")),readdir(path)))
end
max_error=0.0
max_ratio=0.0
pairs=0
for name in readdir(before)
    occursin(r"^(opt|phys)-rep\d+$",name) || continue
    a,b=joinpath(before,name),joinpath(after,name)
    @assert isdir(b) name
    files=output_names(a)
    @assert !isempty(files) && files==output_names(b) name
    for file in files
        lines_a=filter(line->!isempty(strip(line)) && !startswith(strip(line),"#"),readlines(joinpath(a,file)))
        lines_b=filter(line->!isempty(strip(line)) && !startswith(strip(line),"#"),readlines(joinpath(b,file)))
        @assert length(lines_a)==length(lines_b) (name,file)
        indices=startswith(file,"zvo_cisajscktalt") ? 8 : startswith(file,"zvo_cisajs") ? 4 : 0
        for (la,lb) in zip(lines_a,lines_b)
            xs,ys=parse.(Float64,split(la)),parse.(Float64,split(lb))
            @assert length(xs)==length(ys) (name,file)
            for i in eachindex(xs,ys)
                x,y=xs[i],ys[i]
                @assert isfinite(x) && isfinite(y) (name,file,i,x,y)
                if i<=indices
                    @assert x==y && isinteger(x) (name,file,i,x,y)
                else
                    bound=1e-11+1e-11*max(abs(x),abs(y))
                    err=abs(x-y)
                    @assert err<=bound (name,file,i,x,y,err,bound)
                    global max_error=max(max_error,err)
                    global max_ratio=max(max_ratio,err/bound)
                end
            end
        end
        global pairs+=1
    end
end
@assert pairs>0
println("PASS file_pairs=$pairs max_abs=$max_error max_error_over_bound=$max_ratio")
