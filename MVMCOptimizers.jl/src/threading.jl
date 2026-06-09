"""
Threading support for VMCMainCal sample-level parallelism.

Provides deterministic per-worker accumulators, timer reduction, and chunk
selection used by `vmc_main_cal!` and `vmc_main_cal_fsz!`.
"""

struct VMCThreadConfig
    requested_threads::Int
    effective_threads::Int
    work_items::Int
    min_work_per_thread::Int
end

function VMCThreadConfig(
    work_items::Integer;
    requested_threads::Integer = Base.Threads.nthreads(),
    min_work_per_thread::Integer = 1,
)
    requested_threads >= 1 ||
        error("requested_threads must be >= 1, got $requested_threads")
    work_items >= 0 || error("work_items must be >= 0, got $work_items")
    min_work_per_thread >= 1 ||
        error("min_work_per_thread must be >= 1, got $min_work_per_thread")

    max_threads_by_work = work_items == 0 ? 1 : cld(work_items, min_work_per_thread)
    effective = min(Int(requested_threads), Int(max_threads_by_work), max(1, Int(work_items)))
    return VMCThreadConfig(
        Int(requested_threads),
        effective,
        Int(work_items),
        Int(min_work_per_thread),
    )
end

@inline vmc_threading_enabled(config::VMCThreadConfig) = config.effective_threads > 1
@inline effective_thread_count(config::VMCThreadConfig) = config.effective_threads

function vmc_main_cal_requested_threads()
    raw = strip(get(ENV, "JULIA_MVMC_MAINCAL_THREADS", ""))
    isempty(raw) && return 1

    parsed = tryparse(Int, raw)
    parsed !== nothing && parsed >= 1 ||
        error("JULIA_MVMC_MAINCAL_THREADS must be a positive integer, got '$raw'")
    return min(parsed, Base.Threads.nthreads())
end

@inline function vmc_inner_threading_requested(threaded::Bool)
    raw = strip(get(ENV, "JULIA_MVMC_INNER_THREADS", "0"))
    if isempty(raw) || raw == "0"
        return false
    end
    raw == "1" ||
        error("JULIA_MVMC_INNER_THREADS must be 0 or 1, got '$raw'")
    return threaded &&
        Base.Threads.nthreads() > 1
end

@inline function vmc_pfapack_threading_requested(threaded::Bool)
    raw = strip(get(ENV, "JULIA_MVMC_PFAPACK_THREADS", "0"))
    if isempty(raw) || raw == "0"
        return false
    end
    raw == "1" ||
        error("JULIA_MVMC_PFAPACK_THREADS must be 0 or 1, got '$raw'")
    return vmc_inner_threading_requested(threaded)
end

const PFAPACK_CALL_LOCK = ReentrantLock()

function with_pfapack_call_lock(f)
    if get(ENV, "JULIA_MVMC_PFAPACK_LOCK", "1") == "0"
        return f()
    end

    lock(PFAPACK_CALL_LOCK)
    try
        return f()
    finally
        unlock(PFAPACK_CALL_LOCK)
    end
end

@inline function vmc_inner_threading_enabled(
    work_items::Integer,
    threaded::Bool;
    min_work_per_thread::Integer = 64,
)
    if !vmc_inner_threading_requested(threaded) || work_items <= 0
        return false
    end
    threshold = max(Int(min_work_per_thread), Base.Threads.nthreads())
    return Int(work_items) >= threshold
end

function copy_real_to_complex!(
    dst::AbstractVector{ComplexF64},
    src::AbstractVector{Float64},
    n::Integer = min(length(dst), length(src));
    threaded::Bool = false,
)
    n_copy = min(Int(n), length(dst), length(src))
    n_copy <= 0 && return dst
    if vmc_inner_threading_enabled(n_copy, threaded)
        Base.Threads.@threads :static for i = 1:n_copy
            @inbounds dst[i] = ComplexF64(src[i], 0.0)
        end
    else
        @inbounds @simd for i = 1:n_copy
            dst[i] = ComplexF64(src[i], 0.0)
        end
    end
    return dst
end

function copy_complex_realpart!(
    dst::AbstractVector{Float64},
    src::AbstractVector{ComplexF64},
    n::Integer = min(length(dst), length(src));
    threaded::Bool = false,
)
    n_copy = min(Int(n), length(dst), length(src))
    n_copy <= 0 && return dst
    if vmc_inner_threading_enabled(n_copy, threaded)
        Base.Threads.@threads :static for i = 1:n_copy
            @inbounds dst[i] = real(src[i])
        end
    else
        @inbounds @simd for i = 1:n_copy
            dst[i] = real(src[i])
        end
    end
    return dst
end

mutable struct VMCEnergyAccumulator
    wc::ComplexF64
    etot::ComplexF64
    etot2::ComplexF64
    sztot::ComplexF64
    sztot2::ComplexF64
end

VMCEnergyAccumulator() =
    VMCEnergyAccumulator(0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im)

function clear_energy_accumulator!(acc::VMCEnergyAccumulator)
    acc.wc = 0.0 + 0.0im
    acc.etot = 0.0 + 0.0im
    acc.etot2 = 0.0 + 0.0im
    acc.sztot = 0.0 + 0.0im
    acc.sztot2 = 0.0 + 0.0im
    return acc
end

function accumulate_energy!(
    acc::VMCEnergyAccumulator,
    w::Real,
    e::ComplexF64;
    sz::Union{Nothing,Real} = nothing,
)
    wc = ComplexF64(w, 0.0)
    acc.wc += wc
    acc.etot += w * e
    acc.etot2 += w * conj(e) * e
    if sz !== nothing
        acc.sztot += w * sz
        acc.sztot2 += w * sz * sz
    end
    return acc
end

function merge_energy_accumulator!(dst::VMCEnergyAccumulator, src::VMCEnergyAccumulator)
    dst.wc += src.wc
    dst.etot += src.etot
    dst.etot2 += src.etot2
    dst.sztot += src.sztot
    dst.sztot2 += src.sztot2
    return dst
end

function merge_energy_accumulator!(dst::EnergyData, src::VMCEnergyAccumulator)
    dst.wc += src.wc
    dst.etot += src.etot
    dst.etot2 += src.etot2
    dst.sztot += src.sztot
    dst.sztot2 += src.sztot2
    return dst
end

function merge_energy_accumulators!(dst::EnergyData, locals)
    for acc in locals
        merge_energy_accumulator!(dst, acc)
    end
    return dst
end

mutable struct VMCSROptAccumulator
    sr_opt_oo::Vector{ComplexF64}
    sr_opt_ho::Vector{ComplexF64}
    sr_opt_o::Vector{ComplexF64}
    sr_opt_o_store::Vector{ComplexF64}
    sr_opt_oo_real::Vector{Float64}
    sr_opt_ho_real::Vector{Float64}
    sr_opt_o_real::Vector{Float64}
    sr_opt_o_store_real::Vector{Float64}
end

function VMCSROptAccumulator(sr::SROptData)
    return VMCSROptAccumulator(
        zeros(ComplexF64, length(sr.sr_opt_oo)),
        zeros(ComplexF64, length(sr.sr_opt_ho)),
        zeros(ComplexF64, length(sr.sr_opt_o)),
        zeros(ComplexF64, length(sr.sr_opt_o_store)),
        zeros(Float64, length(sr.sr_opt_oo_real)),
        zeros(Float64, length(sr.sr_opt_ho_real)),
        zeros(Float64, length(sr.sr_opt_o_real)),
        zeros(Float64, length(sr.sr_opt_o_store_real)),
    )
end

function clear_sropt_accumulator!(acc::VMCSROptAccumulator)
    fill!(acc.sr_opt_oo, 0.0 + 0.0im)
    fill!(acc.sr_opt_ho, 0.0 + 0.0im)
    fill!(acc.sr_opt_o, 0.0 + 0.0im)
    fill!(acc.sr_opt_o_store, 0.0 + 0.0im)
    fill!(acc.sr_opt_oo_real, 0.0)
    fill!(acc.sr_opt_ho_real, 0.0)
    fill!(acc.sr_opt_o_real, 0.0)
    fill!(acc.sr_opt_o_store_real, 0.0)
    return acc
end

function clear_sropt_store!(sr::SROptData)
    fill!(sr.sr_opt_o_store, 0.0 + 0.0im)
    fill!(sr.sr_opt_o_store_real, 0.0)
    return sr
end

function merge_sropt_accumulator!(dst::VMCSROptAccumulator, src::VMCSROptAccumulator)
    dst.sr_opt_oo .+= src.sr_opt_oo
    dst.sr_opt_ho .+= src.sr_opt_ho
    dst.sr_opt_o_store .+= src.sr_opt_o_store
    dst.sr_opt_oo_real .+= src.sr_opt_oo_real
    dst.sr_opt_ho_real .+= src.sr_opt_ho_real
    dst.sr_opt_o_store_real .+= src.sr_opt_o_store_real
    return dst
end

function merge_sropt_accumulator!(dst::SROptData, src::VMCSROptAccumulator)
    dst.sr_opt_oo .+= src.sr_opt_oo
    dst.sr_opt_ho .+= src.sr_opt_ho
    dst.sr_opt_o_store .+= src.sr_opt_o_store
    dst.sr_opt_oo_real .+= src.sr_opt_oo_real
    dst.sr_opt_ho_real .+= src.sr_opt_ho_real
    dst.sr_opt_o_store_real .+= src.sr_opt_o_store_real
    return dst
end

function merge_sropt_accumulators!(dst::SROptData, locals)
    for acc in locals
        merge_sropt_accumulator!(dst, acc)
    end
    return dst
end

mutable struct VMCPhysAccumulator
    local_cis_ajs::Vector{ComplexF64}
    phys_cis_ajs::Vector{ComplexF64}
    phys_cis_ajs_ckt_alt::Vector{ComplexF64}
    local_cis_ajs_ckt_alt_dc::Vector{ComplexF64}
    phys_cis_ajs_ckt_alt_dc::Vector{ComplexF64}
end

VMCPhysAccumulator(::Nothing) =
    VMCPhysAccumulator(ComplexF64[], ComplexF64[], ComplexF64[], ComplexF64[], ComplexF64[])

function VMCPhysAccumulator(phys::PhysicalQuantities)
    return VMCPhysAccumulator(
        zeros(ComplexF64, length(phys.local_cis_ajs)),
        zeros(ComplexF64, length(phys.phys_cis_ajs)),
        zeros(ComplexF64, length(phys.phys_cis_ajs_ckt_alt)),
        zeros(ComplexF64, length(phys.local_cis_ajs_ckt_alt_dc)),
        zeros(ComplexF64, length(phys.phys_cis_ajs_ckt_alt_dc)),
    )
end

function clear_phys_accumulator!(acc::VMCPhysAccumulator)
    fill!(acc.local_cis_ajs, 0.0 + 0.0im)
    fill!(acc.phys_cis_ajs, 0.0 + 0.0im)
    fill!(acc.phys_cis_ajs_ckt_alt, 0.0 + 0.0im)
    fill!(acc.local_cis_ajs_ckt_alt_dc, 0.0 + 0.0im)
    fill!(acc.phys_cis_ajs_ckt_alt_dc, 0.0 + 0.0im)
    return acc
end

function merge_phys_accumulator!(dst::VMCPhysAccumulator, src::VMCPhysAccumulator)
    dst.phys_cis_ajs .+= src.phys_cis_ajs
    dst.phys_cis_ajs_ckt_alt .+= src.phys_cis_ajs_ckt_alt
    dst.phys_cis_ajs_ckt_alt_dc .+= src.phys_cis_ajs_ckt_alt_dc
    return dst
end

function merge_phys_accumulator!(dst::PhysicalQuantities, src::VMCPhysAccumulator)
    dst.phys_cis_ajs .+= src.phys_cis_ajs
    dst.phys_cis_ajs_ckt_alt .+= src.phys_cis_ajs_ckt_alt
    dst.phys_cis_ajs_ckt_alt_dc .+= src.phys_cis_ajs_ckt_alt_dc
    return dst
end

function merge_phys_accumulators!(dst::Union{PhysicalQuantities,Nothing}, locals)
    dst === nothing && return nothing
    for acc in locals
        merge_phys_accumulator!(dst, acc)
    end
    return dst
end

mutable struct VMCCounterAccumulator
    counter::Vector{Int}
end

VMCCounterAccumulator(n_counter::Integer) = VMCCounterAccumulator(zeros(Int, Int(n_counter)))

function record_counter!(acc::VMCCounterAccumulator, index::Integer, amount::Integer = 1)
    1 <= index <= length(acc.counter) ||
        error("counter index $index out of range 1:$(length(acc.counter))")
    acc.counter[Int(index)] += Int(amount)
    return acc
end

function merge_counter_accumulator!(dst::VMCCounterAccumulator, src::VMCCounterAccumulator)
    if length(dst.counter) < length(src.counter)
        old_len = length(dst.counter)
        resize!(dst.counter, length(src.counter))
        fill!(@view(dst.counter[(old_len+1):end]), 0)
    end
    @inbounds for i in eachindex(src.counter)
        dst.counter[i] += src.counter[i]
    end
    return dst
end

function merge_counter_accumulator!(dst::Vector{Int}, src::VMCCounterAccumulator)
    if length(dst) < length(src.counter)
        old_len = length(dst)
        resize!(dst, length(src.counter))
        fill!(@view(dst[(old_len+1):end]), 0)
    end
    @inbounds for i in eachindex(src.counter)
        dst[i] += src.counter[i]
    end
    return dst
end

function merge_counter_accumulators!(dst::Vector{Int}, locals)
    for acc in locals
        merge_counter_accumulator!(dst, acc)
    end
    return dst
end

@inline function _check_ctimer_id(id::Integer)
    0 <= id < CTIMER_N || error("timer id $id out of range [0, $(CTIMER_N - 1)]")
    return Int(id)
end

@inline ctimer_add_elapsed!(::CTimer{Val{false}}, ::Integer, ::UInt64) = nothing

function ctimer_add_elapsed!(timer::CTimer{Val{true}}, id::Integer, elapsed_ns::UInt64)
    idx = _check_ctimer_id(id) + 1
    @inbounds timer.elapsed_ns[idx] += elapsed_ns
    return timer
end

ctimer_merge!(dst::CTimer{Val{false}}, ::CTimer) = dst
ctimer_merge!(dst::CTimer{Val{true}}, ::CTimer{Val{false}}) = dst

function ctimer_merge!(dst::CTimer{Val{true}}, src::CTimer{Val{true}})
    @inbounds for i in eachindex(dst.elapsed_ns, src.elapsed_ns)
        dst.elapsed_ns[i] += src.elapsed_ns[i]
    end
    return dst
end

function ctimer_merge_all!(dst::CTimer, locals)
    for timer in locals
        ctimer_merge!(dst, timer)
    end
    return dst
end

mutable struct VMCThreadAccumulator
    energy::VMCEnergyAccumulator
    sr_opt::VMCSROptAccumulator
    phys::VMCPhysAccumulator
    counter::VMCCounterAccumulator
    timer::CTimer
end

function VMCThreadAccumulator(state::VMCOptimizationState, parent_timer::CTimer = CTIMER_DISABLED)
    return VMCThreadAccumulator(
        VMCEnergyAccumulator(),
        VMCSROptAccumulator(state.sr_opt),
        VMCPhysAccumulator(state.phys_quantities),
        VMCCounterAccumulator(length(state.electron_config.counter)),
        CTimer(ctimer_enabled(parent_timer)),
    )
end

function make_thread_accumulators(
    state::VMCOptimizationState,
    config::VMCThreadConfig,
    parent_timer::CTimer = CTIMER_DISABLED,
)
    return [VMCThreadAccumulator(state, parent_timer) for _ = 1:config.effective_threads]
end

function merge_thread_accumulator!(
    state::VMCOptimizationState,
    parent_timer::CTimer,
    local_acc::VMCThreadAccumulator,
)
    merge_energy_accumulator!(state.energy, local_acc.energy)
    merge_sropt_accumulator!(state.sr_opt, local_acc.sr_opt)
    merge_phys_accumulators!(state.phys_quantities, (local_acc.phys,))
    merge_counter_accumulator!(state.electron_config.counter, local_acc.counter)
    ctimer_merge!(parent_timer, local_acc.timer)
    return state
end

function merge_thread_accumulators!(
    state::VMCOptimizationState,
    parent_timer::CTimer,
    locals,
)
    for local_acc in locals
        merge_thread_accumulator!(state, parent_timer, local_acc)
    end
    return state
end
