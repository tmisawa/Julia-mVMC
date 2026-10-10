# Optional SR sample-store crossover probe. Run from the Julia workspace:
# UCX_ERROR_SIGNALS=SIGILL,SIGBUS,SIGFPE UCX_MEMTYPE_CACHE=no \
# OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=16,0 JULIA_NUM_GC_THREADS=1 \
# JULIA_MVMC_INNER_THREADS=1 julia +1.13.1 --project=. benchmark/sr_store_dispatch.jl
using MVMCOptimizers, LinearAlgebra, Statistics

BLAS.set_num_threads(1)
const OPT = MVMCOptimizers
const UNUSED_OO = ComplexF64[]

# The original NStore dispatch starts a static region for vectors >=64 items.
# Keep its arithmetic as a benchmark control, independent of the new gate.
function legacy_store!(ho, store, values, w, e, n)
    we = w * e
    sqrtw = sqrt(w)
    Base.Threads.@threads :static for i in eachindex(values)
        @inbounds begin
            store[i] = sqrtw * values[i]
            ho[i] += we * values[i]
        end
    end
end

function production_store!(ho::Vector{Float64}, store, values, w, e, n, threaded)
    OPT.calculate_oo_store_real!(ho, store, values, w, e, 0, n; threaded)
end

function production_store!(ho::Vector{ComplexF64}, store, values, w, e, n, threaded)
    OPT.calculate_oo_store!(UNUSED_OO, ho, store, values, w, e, 0, n; threaded)
end

function batch!(f::F, ho, store, values, w, e, n, repeats) where {F}
    start = time_ns()
    for _ in 1:repeats
        f(ho, store, values, w, e, n)
    end
    return (time_ns() - start) / repeats
end

serial_store!(args...) = production_store!(args..., false)
candidate_store!(args...) = production_store!(args..., true)

function main()
    OPT.vmc_inner_threading_requested(true) || error("Use at least two default-pool threads and JULIA_MVMC_INNER_THREADS=1")
    println("Julia=$VERSION threads=$(Threads.nthreads()) BLAS=$(BLAS.get_config())")
    println("kind,size,serial_ns,legacy_threaded_ns,candidate_ns")
    for T in (Float64, ComplexF64), n in (200, 400, 1024, 4096, 16384, 65536, 262144)
        len = T <: Complex ? 2n : n
        values = T[sin(i / 17) for i in 1:len]
        energy = T <: Complex ? T(-0.625, 0.125) : T(-0.625)
        arrays = [(zeros(T, len), zeros(T, len)) for _ in 1:3]
        functions = (serial_store!, legacy_store!, candidate_store!)
        timings = [Float64[] for _ in 1:3]
        repeats = clamp(2^21 ÷ len, 16, 1000)
        for (f, (ho, store)) in zip(functions, arrays)
            for _ in 1:20
                f(ho, store, values, 0.75, energy, n)
            end
        end
        for trial in 0:6, offset in 0:2
            index = mod(trial + offset, 3) + 1
            ho, store = arrays[index]
            push!(timings[index], batch!(functions[index], ho, store, values, 0.75, energy, n, repeats))
        end
        for (ho, store) in arrays[2:3]
            # Same independent elementwise operations. Four eps cover lowering
            # roundoff; a missed write is far outside this budget.
            @assert isapprox(ho, arrays[1][1]; atol=4eps(Float64), rtol=4eps(Float64))
            @assert isapprox(store, arrays[1][2]; atol=4eps(Float64), rtol=4eps(Float64))
        end
        println("$T,$n,$(median(timings[1])),$(median(timings[2])),$(median(timings[3]))")
        flush(stdout)
    end
end

main()
