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

# --- Optional Linux CPU affinity for the persistent hop team -----------------
# Worker threads that migrate between physical cores lose their private L2
# copy of the accepted-hop inverse planes. On the benchmark hosts this alone
# costs about 1.3x in the rank-1 update. Pin the caller and each worker to a
# distinct CPU from the process affinity set for the lifetime of one sampling
# call, then restore the original set. Linux-only and best-effort: an
# unsupported call or a smaller CPU set leaves scheduling untouched.
const _CPU_SET_BYTES = 128

@inline _cpu_pinning_enabled() =
    Sys.islinux() && get(ENV, "JULIA_MVMC_PIN_CPUS", "1") != "0"

function _get_affinity_mask()
    mask = zeros(UInt8, _CPU_SET_BYTES)
    rc = ccall(:sched_getaffinity, Cint, (Cint, Csize_t, Ptr{UInt8}),
               0, _CPU_SET_BYTES, mask)
    rc == 0 || return nothing
    return mask
end

function _set_affinity_mask(mask::Vector{UInt8})
    rc = ccall(:sched_setaffinity, Cint, (Cint, Csize_t, Ptr{UInt8}),
               0, length(mask), mask)
    return rc == 0
end

function _affinity_cpus(mask::Vector{UInt8})
    cpus = Int[]
    for byte in eachindex(mask)
        bits = mask[byte]
        bits == 0x00 && continue
        for bit in 0:7
            (bits >> bit) & 0x01 == 0x01 && push!(cpus, (byte - 1) * 8 + bit)
        end
    end
    return cpus
end

function _cpu_core_key(cpu::Int)
    top = "/sys/devices/system/cpu/cpu$cpu/topology"
    pkg_file = joinpath(top, "physical_package_id")
    core_file = joinpath(top, "core_id")
    (isfile(pkg_file) && isfile(core_file)) || return (-1, -1)
    pk = tryparse(Int, strip(read(pkg_file, String)))
    co = tryparse(Int, strip(read(core_file, String)))
    (pk === nothing || co === nothing) && return (-1, -1)
    return (pk, co)
end

function _pin_team_cpus(workers::Int)
    mask = _get_affinity_mask()
    mask === nothing && return Int[], UInt8[]
    cpus = _affinity_cpus(mask)
    length(cpus) >= workers + 1 || return Int[], UInt8[]
    # Prefer one CPU per physical core so SMT siblings do not share a worker
    # slot's L1/L2. Fall back to the allowed order when topology is unavailable.
    selected = Int[]
    seen = Set{Tuple{Int,Int}}()
    for cpu in cpus
        key = _cpu_core_key(cpu)
        if key !== (-1, -1) && !(key in seen)
            push!(seen, key)
            push!(selected, cpu)
        end
        length(selected) >= workers + 1 && break
    end
    if length(selected) < workers + 1
        for cpu in cpus
            cpu in selected || push!(selected, cpu)
            length(selected) >= workers + 1 && break
        end
    end
    length(selected) >= workers + 1 || return Int[], UInt8[]
    return selected[1:workers+1], mask
end

@inline function _pin_current_thread(cpu::Int)
    mask = zeros(UInt8, _CPU_SET_BYTES)
    mask[cpu ÷ 8 + 1] = UInt8(1) << (cpu % 8)
    return _set_affinity_mask(mask)
end

# Accepted hops are short and frequent: creating static tasks for every hop
# costs more than the rank-2 update. Keep N-1 workers for one sampling call;
# the calling thread performs the remaining partition and owns the RNG.
mutable struct VMCHopTeam
    epoch::Threads.Atomic{Int}
    completed::Threads.Atomic{Int}
    started::Threads.Atomic{Int}
    stop::Threads.Atomic{Bool}
    wake::Base.Event
    sleeping::Bool
    job::Any
    tasks::Vector{Task}
    errors::Vector{Any}
    workers::Int
    pin_cpus::Vector{Int}
    original_affinity::Vector{UInt8}
end

@inline function _vmc_spin_pause()
    ccall(:jl_cpu_pause, Cvoid, ())
end

@inline function _vmc_spin_safepoint(counter::Int)
    (counter & 0xff) == 0 && GC.safepoint()
    return nothing
end

function _vmc_hop_worker!(team::VMCHopTeam, slot::Int)
    # A spinning worker cannot migrate; keep it sticky for its whole lifetime.
    current_task().sticky = true
    if slot < length(team.pin_cpus)
        _pin_current_thread(team.pin_cpus[slot+1])
    end
    epoch = 0
    spins = 0
    Threads.atomic_add!(team.started, 1)
    try
        while !team.stop[]
            next_epoch = team.epoch[]
            if next_epoch == epoch
                if team.sleeping
                    wait(team.wake)
                else
                    _vmc_spin_pause()
                    spins += 1
                    _vmc_spin_safepoint(spins)
                end
                continue
            end
            spins = 0
            epoch = next_epoch
            try
                team.job(slot)
            catch err
                team.errors[slot] = CapturedException(err, catch_backtrace())
            finally
                Threads.atomic_add!(team.completed, 1)
            end
        end
    finally
        isempty(team.original_affinity) || _set_affinity_mask(team.original_affinity)
    end
    return nothing
end

function with_vmc_hop_team(f, workspace::SamplingWorkspace;
                           qp_count=length(workspace.pf_m_real_temp),
                           sleep_workers::Bool=false)
    workspace.hop_team === nothing || return f()
    # Use only default-pool threads that are not the caller's own thread. With
    # `--threads=N,0` the caller itself is a default thread; pinning a
    # spinning helper there would starve it.
    caller_tid = Threads.threadid()
    available_tids = filter(tid -> tid != caller_tid, Threads.threadpooltids(:default))
    workers = min(min(Threads.nthreads(), qp_count) - 1, length(available_tids))
    workers > 0 || return f()
    pin_cpus, original_affinity = _cpu_pinning_enabled() ?
        _pin_team_cpus(workers) : (Int[], UInt8[])
    team = VMCHopTeam(Threads.Atomic{Int}(0), Threads.Atomic{Int}(0),
                      Threads.Atomic{Int}(0), Threads.Atomic{Bool}(false),
                      Base.Event(), sleep_workers,
                      nothing, Task[], Any[nothing for _ in 1:workers], workers,
                      pin_cpus, original_affinity)
    caller_pinned = !isempty(pin_cpus) && _pin_current_thread(pin_cpus[1])
    storage = task_local_storage()
    previous = get(storage,:mvmc_hop_team_active,false)
    try
        # Pin each worker to its own default-pool thread. `Threads.@spawn`
        # does not guarantee distinct threads, and co-located workers both
        # halve the hop-update parallelism and make the thread-indexed scratch
        # vectors alias. `jl_set_task_tid` uses 0-based ids.
        for slot in 1:workers
            task = Task(() -> _vmc_hop_worker!(team, slot))
            task.sticky = true
            ccall(:jl_set_task_tid, Cint, (Any, Cint), task, available_tids[slot] - 1)
            push!(team.tasks, task)
            schedule(task)
        end
        # Do not yield the caller onto a spinning task on its own thread.
        spins = 0
        while team.started[] != workers
            _vmc_spin_pause()
            spins += 1
            _vmc_spin_safepoint(spins)
        end
        workspace.hop_team = team
        storage[:mvmc_hop_team_active] = true
        return f()
    finally
        workspace.hop_team = nothing
        storage[:mvmc_hop_team_active] = previous
        team.stop[] = true
        team.sleeping && notify(team.wake)
        foreach(wait, team.tasks)
        caller_pinned && !isempty(original_affinity) && _set_affinity_mask(original_affinity)
    end
end

function with_vmc_real_sampler_team(f, data::ExpertModeData,
                                    state::VMCOptimizationState, ctx::ParallelContext)
    n_size = 2data.modpara.nelec
    first_qp,last_qp = qp_split_range(length(state.slater_matrix.pf_m_real),ctx)
    n_qp = last_qp-first_qp
    # Other proposal families retain their existing scheduling until their
    # accepted updates have a measured parallel route of their own.
    if data.modpara.nex_update_path == 0 && vmc_inner_threading_requested(true) &&
       n_qp>=2 && n_qp*n_size^2>=65536 &&
       Threads.threadid()==1 && ccall(:jl_in_threaded_region,Cint,())==0
        return with_vmc_hop_team(f,state.workspace;qp_count=n_qp)
    end
    return f()
end

function vmc_hop_dispatch!(job, team::VMCHopTeam)
    team.completed[] = 0
    fill!(team.errors, nothing)
    team.job = job
    # Sequentially consistent atomic publication also publishes the closure
    # and all matrix/configuration writes before workers read them.
    Threads.atomic_add!(team.epoch, 1)
    team.sleeping && notify(team.wake)
    caller_error = nothing
    try
        job(0)
    catch err
        caller_error = CapturedException(err, catch_backtrace())
    finally
        spins = 0
        while team.completed[] != team.workers
            _vmc_spin_pause()
            spins += 1
            _vmc_spin_safepoint(spins)
        end
        team.job = nothing
        team.sleeping && reset(team.wake)
    end
    caller_error === nothing || throw(caller_error)
    for err in team.errors
        err === nothing || throw(err)
    end
    return nothing
end

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

# Array copies carry much less work per item than numerical kernels. Starting
# static thread tasks for 64 values costs more than the serial SIMD copy.
# The standalone copy benchmark in benchmark/inner_copy_dispatch.jl records
# the crossover; keep the generic kernel gate unchanged.
const VMC_COPY_MIN_WORK_ITEMS = 65536

@inline function vmc_copy_threading_enabled(work_items::Integer, threaded::Bool)
    # Persistent sampling workers own their threads until the sampling call
    # returns. Its few bulk copies stay serial instead of spawning other tasks.
    get(task_local_storage(),:mvmc_hop_team_active,false) && return false
    return vmc_inner_threading_enabled(
        work_items, threaded; min_work_per_thread = VMC_COPY_MIN_WORK_ITEMS,
    )
end

# SR sample storage is two independent vector updates, not a matrix kernel.
# On the reference Ryzen host, launching static tasks for a few hundred values
# costs several microseconds while the serial SIMD loop takes under one.
# The standalone SR-store probe records profitable large-vector crossovers.
const VMC_SR_STORE_MIN_REAL_ITEMS = 65536
const VMC_SR_STORE_MIN_COMPLEX_ITEMS = 32768

@inline function vmc_sr_store_threading_enabled(
    work_items::Integer, threaded::Bool, ::Type{T},
) where {T}
    threshold = T <: Complex ? VMC_SR_STORE_MIN_COMPLEX_ITEMS : VMC_SR_STORE_MIN_REAL_ITEMS
    return vmc_inner_threading_enabled(
        work_items, threaded; min_work_per_thread = threshold,
    )
end

function copy_real_to_complex!(
    dst::AbstractVector{ComplexF64},
    src::AbstractVector{Float64},
    n::Integer = min(length(dst), length(src));
    threaded::Bool = false,
)
    n_copy = min(Int(n), length(dst), length(src))
    n_copy <= 0 && return dst
    if vmc_copy_threading_enabled(n_copy, threaded)
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
    if vmc_copy_threading_enabled(n_copy, threaded)
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

function VMCSROptAccumulator(sr::SROptData; owned::Bool = false)
    if owned
        # MainCal cleared aggregates before acquiring this private accumulator.
        # Derivative scratch remains private: it was never published by merge.
        clear_sropt_store!(sr)
        return VMCSROptAccumulator(
            sr.sr_opt_oo, sr.sr_opt_ho, zeros(ComplexF64, length(sr.sr_opt_o)), sr.sr_opt_o_store,
            sr.sr_opt_oo_real, sr.sr_opt_ho_real, zeros(Float64, length(sr.sr_opt_o_real)), sr.sr_opt_o_store_real,
        )
    end
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

# Only a private MainCal accumulator may alias every published SR array.
function _owns_sropt_accumulator(sr::SROptData, acc::VMCSROptAccumulator)
    return all(name -> getproperty(sr, name) === getproperty(acc, name),
        (:sr_opt_oo, :sr_opt_ho, :sr_opt_o_store, :sr_opt_oo_real, :sr_opt_ho_real, :sr_opt_o_store_real))
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
    phys_lanczos_qqqq::Vector{ComplexF64}
    phys_lanczos_qcisajsq::Vector{ComplexF64}
    phys_lanczos_qcisajscktaltq::Vector{ComplexF64}
    phys_lanczos_qcisajscktaltq_dc::Vector{ComplexF64}
end

VMCPhysAccumulator(::Nothing) =
    VMCPhysAccumulator(
        ComplexF64[],
        ComplexF64[],
        ComplexF64[],
        ComplexF64[],
        ComplexF64[],
        ComplexF64[],
        ComplexF64[],
        ComplexF64[],
        ComplexF64[],
    )

function VMCPhysAccumulator(phys::PhysicalQuantities)
    return VMCPhysAccumulator(
        zeros(ComplexF64, length(phys.local_cis_ajs)),
        zeros(ComplexF64, length(phys.phys_cis_ajs)),
        zeros(ComplexF64, length(phys.phys_cis_ajs_ckt_alt)),
        zeros(ComplexF64, length(phys.local_cis_ajs_ckt_alt_dc)),
        zeros(ComplexF64, length(phys.phys_cis_ajs_ckt_alt_dc)),
        zeros(ComplexF64, length(phys.phys_lanczos_qqqq)),
        zeros(ComplexF64, length(phys.phys_lanczos_qcisajsq)),
        zeros(ComplexF64, length(phys.phys_lanczos_qcisajscktaltq)),
        zeros(ComplexF64, length(phys.phys_lanczos_qcisajscktaltq_dc)),
    )
end

function clear_phys_accumulator!(acc::VMCPhysAccumulator)
    fill!(acc.local_cis_ajs, 0.0 + 0.0im)
    fill!(acc.phys_cis_ajs, 0.0 + 0.0im)
    fill!(acc.phys_cis_ajs_ckt_alt, 0.0 + 0.0im)
    fill!(acc.local_cis_ajs_ckt_alt_dc, 0.0 + 0.0im)
    fill!(acc.phys_cis_ajs_ckt_alt_dc, 0.0 + 0.0im)
    fill!(acc.phys_lanczos_qqqq, 0.0 + 0.0im)
    fill!(acc.phys_lanczos_qcisajsq, 0.0 + 0.0im)
    fill!(acc.phys_lanczos_qcisajscktaltq, 0.0 + 0.0im)
    fill!(acc.phys_lanczos_qcisajscktaltq_dc, 0.0 + 0.0im)
    return acc
end

function merge_phys_accumulator!(dst::VMCPhysAccumulator, src::VMCPhysAccumulator)
    dst.phys_cis_ajs .+= src.phys_cis_ajs
    dst.phys_cis_ajs_ckt_alt .+= src.phys_cis_ajs_ckt_alt
    dst.phys_cis_ajs_ckt_alt_dc .+= src.phys_cis_ajs_ckt_alt_dc
    dst.phys_lanczos_qqqq .+= src.phys_lanczos_qqqq
    dst.phys_lanczos_qcisajsq .+= src.phys_lanczos_qcisajsq
    dst.phys_lanczos_qcisajscktaltq .+= src.phys_lanczos_qcisajscktaltq
    dst.phys_lanczos_qcisajscktaltq_dc .+= src.phys_lanczos_qcisajscktaltq_dc
    return dst
end

function merge_phys_accumulator!(dst::PhysicalQuantities, src::VMCPhysAccumulator)
    dst.phys_cis_ajs .+= src.phys_cis_ajs
    dst.phys_cis_ajs_ckt_alt .+= src.phys_cis_ajs_ckt_alt
    dst.phys_cis_ajs_ckt_alt_dc .+= src.phys_cis_ajs_ckt_alt_dc
    dst.phys_lanczos_qqqq .+= src.phys_lanczos_qqqq
    dst.phys_lanczos_qcisajsq .+= src.phys_lanczos_qcisajsq
    dst.phys_lanczos_qcisajscktaltq .+= src.phys_lanczos_qcisajscktaltq
    dst.phys_lanczos_qcisajscktaltq_dc .+= src.phys_lanczos_qcisajscktaltq_dc
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

struct VMCRealMainCalBatch
    indices::Vector{Vector{Int}}
    inverse::Vector{Vector{Float64}}
    pfaffian::Matrix{Float64}
    status::Matrix{Int}
end

function VMCRealMainCalBatch(n_size::Int, n_qp::Int, capacity::Int)
    # Keep the same storage size as `SlaterMatrixData.inv_m_real` so a slot can
    # be swapped into the state without copying: `n_qp` planes of `n_size^2`
    # entries plus the historical trailing pad entry per plane.
    VMCRealMainCalBatch([zeros(Int,n_size) for _ in 1:capacity],
        [zeros(Float64,n_qp*(n_size*n_size+1)) for _ in 1:capacity],
        zeros(n_qp,capacity),zeros(Int,n_qp,capacity))
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
    slater_buffer_real::Vector{Float64}
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
    green_thread_idx::Vector{Vector{Int}}
    green_thread_num::Vector{Vector{Int}}
    green_thread_proj::Vector{Vector{Int}}
    green_thread_pf::Vector{Vector{Float64}}
    real_batch::Union{Nothing,VMCRealMainCalBatch}
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
        Float64[],
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
        Vector{Int}[],
        Vector{Int}[],
        Vector{Int}[],
        Vector{Float64}[],
        nothing,
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
    sr_opt::Union{Nothing,VMCSROptAccumulator}
    phys::VMCPhysAccumulator
    counter::VMCCounterAccumulator
    main_cal_scratch::VMCMainCalScratch
    timer::CTimer
    all_complex::Bool
    use_sr_store::Bool
    nsrcg::Bool
    owned_sr::Bool
    use_sr_opt::Bool
end

function VMCThreadAccumulator(
    state::VMCOptimizationState,
    parent_timer::CTimer = CTIMER_DISABLED;
    all_complex::Bool = isempty(state.sr_opt.sr_opt_oo_real),
    use_sr_store::Bool = false,
    nsrcg::Bool = false,
    owned_sr::Bool = false,
    use_sr_opt::Bool = true,
)
    sr_acc = use_sr_opt ? VMCSROptAccumulator(state.sr_opt; owned = owned_sr) : nothing
    return VMCThreadAccumulator(
        VMCEnergyAccumulator(),
        sr_acc,
        VMCPhysAccumulator(state.phys_quantities),
        VMCCounterAccumulator(length(state.electron_config.counter)),
        VMCMainCalScratch(state),
        CTimer(ctimer_enabled(parent_timer)),
        all_complex,
        use_sr_store,
        nsrcg,
        owned_sr,
        use_sr_opt,
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
    if _owns_sropt_accumulator(sr, acc)
        clear_sropt_store!(sr)
        return acc
    end
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
    _resize_fill!(acc.phys_lanczos_qqqq, 0, 0.0 + 0.0im)
    _resize_fill!(acc.phys_lanczos_qcisajsq, 0, 0.0 + 0.0im)
    _resize_fill!(acc.phys_lanczos_qcisajscktaltq, 0, 0.0 + 0.0im)
    _resize_fill!(acc.phys_lanczos_qcisajscktaltq_dc, 0, 0.0 + 0.0im)
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
    _resize_fill!(
        acc.phys_lanczos_qqqq,
        length(phys.phys_lanczos_qqqq),
        0.0 + 0.0im,
    )
    _resize_fill!(
        acc.phys_lanczos_qcisajsq,
        length(phys.phys_lanczos_qcisajsq),
        0.0 + 0.0im,
    )
    _resize_fill!(
        acc.phys_lanczos_qcisajscktaltq,
        length(phys.phys_lanczos_qcisajscktaltq),
        0.0 + 0.0im,
    )
    _resize_fill!(
        acc.phys_lanczos_qcisajscktaltq_dc,
        length(phys.phys_lanczos_qcisajscktaltq_dc),
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
    clear_energy_accumulator!(acc.energy)
    if use_sr_opt
        if acc.sr_opt === nothing
            acc.sr_opt = VMCSROptAccumulator(state.sr_opt; owned = acc.owned_sr)
        end
        sr_acc = acc.sr_opt::VMCSROptAccumulator
        resized = _resize_sropt_accumulator!(sr_acc, state.sr_opt)
        if same_mode && !resized
            reset_sropt_accumulator_for_maincal!(
                sr_acc,
                state.sr_opt;
                all_complex = all_complex,
                use_sr_store = use_sr_store,
                nsrcg = nsrcg,
                use_sr_opt = use_sr_opt,
            )
        else
            clear_sropt_accumulator!(sr_acc)
        end
    end
    reset_phys_accumulator!(acc.phys, state.phys_quantities)
    reset_counter_accumulator!(acc.counter, length(state.electron_config.counter))
    reset_main_cal_scratch!(acc.main_cal_scratch, state)
    ctimer_reset!(acc.timer)

    acc.all_complex = all_complex
    acc.use_sr_store = use_sr_store
    acc.nsrcg = nsrcg
    acc.use_sr_opt = use_sr_opt
    return acc
end

function main_cal_accumulator!(
    state::VMCOptimizationState,
    parent_timer::CTimer = CTIMER_DISABLED;
    all_complex::Bool = isempty(state.sr_opt.sr_opt_oo_real),
    use_sr_store::Bool = false,
    nsrcg::Bool = false,
    use_sr_opt::Bool = true,
    owned_sr::Bool = false,
)
    owned_sr = owned_sr && use_sr_opt && !all_complex
    cached = state.workspace.main_cal_accumulator
    if cached isa VMCThreadAccumulator &&
       ctimer_enabled(cached.timer) == ctimer_enabled(parent_timer) &&
       cached.owned_sr == owned_sr &&
       (!owned_sr || _owns_sropt_accumulator(state.sr_opt, cached.sr_opt))
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
        owned_sr = owned_sr,
        use_sr_opt = use_sr_opt,
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
    if local_acc.use_sr_opt
        sr_acc = local_acc.sr_opt::VMCSROptAccumulator
        if _owns_sropt_accumulator(state.sr_opt, sr_acc)
            # Preserve the original zero-plus-local publication, including signed zero.
            for name in (:sr_opt_oo, :sr_opt_ho, :sr_opt_o_store,
                         :sr_opt_oo_real, :sr_opt_ho_real, :sr_opt_o_store_real)
                values = getproperty(state.sr_opt, name)
                values .+= zero(eltype(values))
            end
        else
            merge_sropt_accumulator!(state.sr_opt, sr_acc)
        end
    end

    merge_phys_accumulators!(state.phys_quantities, (local_acc.phys,))
    merge_counter_accumulator!(state.electron_config.counter, local_acc.counter)
    ctimer_merge!(parent_timer, local_acc.timer)
    return state
end
