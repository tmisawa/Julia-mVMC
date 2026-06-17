"""
Threading and local-accumulator support for C-mVMC OpenMP-equivalent inner loops.

`VMCMainCal` keeps the C rank-local sample loop sequential. Threading in this
file is limited to explicitly gated inner loops with disjoint writes and to
serial local accumulation before merging into the parent state.
"""

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

mutable struct VMCMainCalScratch
    proj_cnt_new::Vector{Int}
    pf_m_new_real::Vector{Float64}
    sample_ele_idx::Vector{Int}
    sample_ele_cfg::Vector{Int}
    sample_ele_num::Vector{Int}
    sample_ele_proj_cnt::Vector{Int}
    sample_ele_spn::Vector{Int}
    slater_buffer::Vector{ComplexF64}
    slater_trans_orb_idx::Vector{Int}
    slater_trans_orb_sgn::Vector{Int}
    calh1_transfer_ri::Vector{Int}
    calh1_transfer_rj::Vector{Int}
    calh1_transfer_spin::Vector{Int}
    calh1_transfer_value::Vector{Float64}
    calh1_transfer_source_len::Int
    calh1_transfer_n_site::Int
    calh1_gutz_site_value::Vector{Float64}
    calh1_jastrow_pair_value::Vector{Float64}
    calh1_projection_n_site::Int
    calh1_projection_n_gutzwiller::Int
    calh1_projection_n_jastrow::Int
    calh1_thread_ele_idx::Vector{Vector{Int}}
    calh1_thread_ele_num::Vector{Vector{Int}}
    calh1_thread_proj_cnt_new::Vector{Vector{Int}}
    calh1_thread_pf_m_new_real::Vector{Vector{Float64}}
    calh1_thread_acc::Vector{Float64}
end

function VMCMainCalScratch(state::VMCOptimizationState)
    return VMCMainCalScratch(
        zeros(Int, length(state.electron_config.tmp_ele_proj_cnt)),
        zeros(Float64, length(state.slater_matrix.pf_m_real)),
        zeros(Int, length(state.electron_config.tmp_ele_idx)),
        zeros(Int, length(state.electron_config.tmp_ele_cfg)),
        zeros(Int, length(state.electron_config.tmp_ele_num)),
        zeros(Int, length(state.electron_config.tmp_ele_proj_cnt)),
        zeros(Int, length(state.electron_config.tmp_ele_spn)),
        ComplexF64[],
        Int[],
        Int[],
        Int[],
        Int[],
        Int[],
        Float64[],
        -1,
        -1,
        Float64[],
        Float64[],
        -1,
        -1,
        -1,
        Vector{Int}[],
        Vector{Int}[],
        Vector{Int}[],
        Vector{Float64}[],
        Float64[],
    )
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
    main_cal_scratch::VMCMainCalScratch
    timer::CTimer
    all_complex::Bool
    use_sr_store::Bool
    nsrcg::Bool
end

function VMCThreadAccumulator(
    state::VMCOptimizationState,
    parent_timer::CTimer = CTIMER_DISABLED;
    all_complex::Bool = isempty(state.sr_opt.sr_opt_oo_real),
    use_sr_store::Bool = false,
    nsrcg::Bool = false,
)
    return VMCThreadAccumulator(
        VMCEnergyAccumulator(),
        VMCSROptAccumulator(state.sr_opt),
        VMCPhysAccumulator(state.phys_quantities),
        VMCCounterAccumulator(length(state.electron_config.counter)),
        VMCMainCalScratch(state),
        CTimer(ctimer_enabled(parent_timer)),
        all_complex,
        use_sr_store,
        nsrcg,
    )
end

function _resize_only!(v::Vector, n::Integer)
    n_int = Int(n)
    if length(v) == n_int
        return false
    end
    resize!(v, n_int)
    return true
end

function _resize_fill!(v::Vector{T}, n::Integer, value::T) where {T}
    _resize_only!(v, n)
    fill!(v, value)
    return v
end

@inline function _zero_tail!(v::Vector{T}, first_tail::Integer, value::T) where {T}
    first = Int(first_tail)
    if first <= length(v)
        fill!(@view(v[first:end]), value)
    end
    return v
end

@inline function _maincal_active_sr_opt_oo_length(
    sr::SROptData;
    all_complex::Bool,
    nsrcg::Bool,
)
    if all_complex
        size_2 = 2 * sr.sr_opt_size
        return nsrcg ? 2 * size_2 : size_2 * size_2
    end
    return nsrcg ? 2 * sr.sr_opt_size : sr.sr_opt_size * sr.sr_opt_size
end

function _resize_sropt_accumulator!(acc::VMCSROptAccumulator, sr::SROptData)
    changed = false
    changed |= _resize_only!(acc.sr_opt_oo, length(sr.sr_opt_oo))
    changed |= _resize_only!(acc.sr_opt_ho, length(sr.sr_opt_ho))
    changed |= _resize_only!(acc.sr_opt_o, length(sr.sr_opt_o))
    changed |= _resize_only!(acc.sr_opt_o_store, length(sr.sr_opt_o_store))
    changed |= _resize_only!(acc.sr_opt_oo_real, length(sr.sr_opt_oo_real))
    changed |= _resize_only!(acc.sr_opt_ho_real, length(sr.sr_opt_ho_real))
    changed |= _resize_only!(acc.sr_opt_o_real, length(sr.sr_opt_o_real))
    changed |= _resize_only!(acc.sr_opt_o_store_real, length(sr.sr_opt_o_store_real))
    return changed
end

function reset_sropt_accumulator_for_maincal!(
    acc::VMCSROptAccumulator,
    sr::SROptData;
    all_complex::Bool,
    use_sr_store::Bool,
    nsrcg::Bool = false,
    use_sr_opt::Bool = true,
)
    if !use_sr_opt
        return clear_sropt_accumulator!(acc)
    end

    # Store mode overwrites the active OO range in finalize_oo_store*!.
    # Keep the store buffer zeroed for skipped samples and keep the inactive
    # tail zero so real->complex SR conversion never observes stale values.
    if all_complex
        if use_sr_store
            fill!(acc.sr_opt_ho, 0.0 + 0.0im)
            fill!(acc.sr_opt_o_store, 0.0 + 0.0im)
            active = min(
                _maincal_active_sr_opt_oo_length(sr; all_complex = true, nsrcg = nsrcg),
                length(acc.sr_opt_oo),
            )
            _zero_tail!(acc.sr_opt_oo, active + 1, 0.0 + 0.0im)
        else
            fill!(acc.sr_opt_oo, 0.0 + 0.0im)
            fill!(acc.sr_opt_ho, 0.0 + 0.0im)
        end
    else
        if use_sr_store
            fill!(acc.sr_opt_ho_real, 0.0)
            fill!(acc.sr_opt_o_store_real, 0.0)
            active = min(
                _maincal_active_sr_opt_oo_length(sr; all_complex = false, nsrcg = nsrcg),
                length(acc.sr_opt_oo_real),
            )
            _zero_tail!(acc.sr_opt_oo_real, active + 1, 0.0)
        else
            fill!(acc.sr_opt_oo_real, 0.0)
            fill!(acc.sr_opt_ho_real, 0.0)
        end
    end
    return acc
end

function reset_phys_accumulator!(acc::VMCPhysAccumulator, ::Nothing)
    _resize_fill!(acc.local_cis_ajs, 0, 0.0 + 0.0im)
    _resize_fill!(acc.phys_cis_ajs, 0, 0.0 + 0.0im)
    _resize_fill!(acc.phys_cis_ajs_ckt_alt, 0, 0.0 + 0.0im)
    _resize_fill!(acc.local_cis_ajs_ckt_alt_dc, 0, 0.0 + 0.0im)
    _resize_fill!(acc.phys_cis_ajs_ckt_alt_dc, 0, 0.0 + 0.0im)
    return acc
end

function reset_phys_accumulator!(acc::VMCPhysAccumulator, phys::PhysicalQuantities)
    _resize_fill!(acc.local_cis_ajs, length(phys.local_cis_ajs), 0.0 + 0.0im)
    _resize_fill!(acc.phys_cis_ajs, length(phys.phys_cis_ajs), 0.0 + 0.0im)
    _resize_fill!(
        acc.phys_cis_ajs_ckt_alt,
        length(phys.phys_cis_ajs_ckt_alt),
        0.0 + 0.0im,
    )
    _resize_fill!(
        acc.local_cis_ajs_ckt_alt_dc,
        length(phys.local_cis_ajs_ckt_alt_dc),
        0.0 + 0.0im,
    )
    _resize_fill!(
        acc.phys_cis_ajs_ckt_alt_dc,
        length(phys.phys_cis_ajs_ckt_alt_dc),
        0.0 + 0.0im,
    )
    return acc
end

function reset_counter_accumulator!(acc::VMCCounterAccumulator, n_counter::Integer)
    _resize_fill!(acc.counter, Int(n_counter), 0)
    return acc
end

function reset_main_cal_scratch!(scratch::VMCMainCalScratch, state::VMCOptimizationState)
    _resize_only!(scratch.proj_cnt_new, length(state.electron_config.tmp_ele_proj_cnt))
    _resize_only!(scratch.pf_m_new_real, length(state.slater_matrix.pf_m_real))
    _resize_only!(scratch.sample_ele_idx, length(state.electron_config.tmp_ele_idx))
    _resize_only!(scratch.sample_ele_cfg, length(state.electron_config.tmp_ele_cfg))
    _resize_only!(scratch.sample_ele_num, length(state.electron_config.tmp_ele_num))
    _resize_only!(scratch.sample_ele_proj_cnt, length(state.electron_config.tmp_ele_proj_cnt))
    _resize_only!(scratch.sample_ele_spn, length(state.electron_config.tmp_ele_spn))

    # Rebuild the transfer cache each MainCal pass. The Hamiltonian is normally
    # fixed, but invalidating here preserves the pre-reuse semantics if terms are
    # edited between calls on the same state.
    scratch.calh1_transfer_source_len = -1
    scratch.calh1_transfer_n_site = -1
    return scratch
end

function reset_thread_accumulator!(
    acc::VMCThreadAccumulator,
    state::VMCOptimizationState;
    all_complex::Bool,
    use_sr_store::Bool,
    nsrcg::Bool = false,
    use_sr_opt::Bool = true,
)
    same_mode =
        acc.all_complex == all_complex &&
        acc.use_sr_store == use_sr_store &&
        acc.nsrcg == nsrcg
    resized = _resize_sropt_accumulator!(acc.sr_opt, state.sr_opt)

    clear_energy_accumulator!(acc.energy)
    if same_mode && !resized
        reset_sropt_accumulator_for_maincal!(
            acc.sr_opt,
            state.sr_opt;
            all_complex = all_complex,
            use_sr_store = use_sr_store,
            nsrcg = nsrcg,
            use_sr_opt = use_sr_opt,
        )
    else
        clear_sropt_accumulator!(acc.sr_opt)
    end
    reset_phys_accumulator!(acc.phys, state.phys_quantities)
    reset_counter_accumulator!(acc.counter, length(state.electron_config.counter))
    reset_main_cal_scratch!(acc.main_cal_scratch, state)
    ctimer_reset!(acc.timer)

    acc.all_complex = all_complex
    acc.use_sr_store = use_sr_store
    acc.nsrcg = nsrcg
    return acc
end

function main_cal_accumulator!(
    state::VMCOptimizationState,
    parent_timer::CTimer = CTIMER_DISABLED;
    all_complex::Bool = isempty(state.sr_opt.sr_opt_oo_real),
    use_sr_store::Bool = false,
    nsrcg::Bool = false,
    use_sr_opt::Bool = true,
)
    cached = state.workspace.main_cal_accumulator
    if cached isa VMCThreadAccumulator &&
       ctimer_enabled(cached.timer) == ctimer_enabled(parent_timer)
        return reset_thread_accumulator!(
            cached,
            state;
            all_complex = all_complex,
            use_sr_store = use_sr_store,
            nsrcg = nsrcg,
            use_sr_opt = use_sr_opt,
        )
    end

    acc = VMCThreadAccumulator(
        state,
        parent_timer;
        all_complex = all_complex,
        use_sr_store = use_sr_store,
        nsrcg = nsrcg,
    )
    state.workspace.main_cal_accumulator = acc
    return acc
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
