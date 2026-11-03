# Optional, untimed RNG diagnostic: julia +1.13 --project=. audit.jl INPUTS NEW_OUTPUT opt|phys
using MVMCOptimizers, SFMT, LinearAlgebra, SHA, TOML
const AUDIT_BYTES=IOBuffer()
const AUDIT_U32=Ref(0)
const AUDIT_REAL2=Ref(0)
const AUDIT_SEEDS=UInt32[]
function SFMT.C_API.gen_rand32()
    word=ccall((:gen_rand32,SFMT.C_API.libsfmt),UInt32,())
    write(AUDIT_BYTES,UInt8(1),word)
    AUDIT_U32[]+=1
    return word
end
function SFMT.C_API.genrand_real2()
    word=only(SFMT.sfmt_dump_rand32(1))
    value=ccall((:genrand_real2,SFMT.C_API.libsfmt),Cdouble,())
    expected=ccall((:to_real2,SFMT.C_API.libsfmt),Cdouble,(UInt32,),word)
    @assert value==expected # Exact RNG conversion, independently of numerical outputs.
    write(AUDIT_BYTES,UInt8(2),word)
    AUDIT_REAL2[]+=1
    return value
end
function SFMT.C_API.init_gen_rand(seed)
    word=UInt32(seed)
    push!(AUDIT_SEEDS,word)
    write(AUDIT_BYTES,UInt8(3),word)
    ccall((:init_gen_rand,SFMT.C_API.libsfmt),Cvoid,(UInt32,),word)
end
input,output,workload=ARGS
@assert workload in ("opt","phys")
ispath(output) && error("output must be fresh")
mkpath(output)
BLAS.set_num_threads(1)
root=joinpath(abspath(input),workload)
result=if workload=="opt"
    run_para_opt_from_namelist(joinpath(root,"namelist.def");nsteps=20,nsmp=1,mode=:real,seed=79,output_dir=output)
else
    run_phys_cal_from_namelist(joinpath(root,"namelist.def");opt_para=joinpath(root,"zqp_opt.dat"),mode=:real,seed=79,output_dir=output)
end
@assert result.status==0
bytes=take!(AUDIT_BYTES)
@assert length(bytes)==5*(AUDIT_U32[]+AUDIT_REAL2[]+length(AUDIT_SEEDS))
open(io->write(io,bytes),joinpath(output,"rng.bin"),"w")
open(joinpath(output,"rng.toml"),"w") do io
    TOML.print(io,Dict("julia_version"=>string(VERSION),"u32"=>AUDIT_U32[],
        "real2"=>AUDIT_REAL2[],"seeds"=>Int.(AUDIT_SEEDS),
        "trace_sha256"=>bytes2hex(sha256(bytes)),"next_624_words"=>Int.(SFMT.sfmt_dump_rand32(624))))
end
println("AUDIT u32=$(AUDIT_U32[]) real2=$(AUDIT_REAL2[]) seeds=$AUDIT_SEEDS sha256=$(bytes2hex(sha256(bytes)))")
