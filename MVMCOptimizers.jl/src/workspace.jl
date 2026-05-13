"""
Workspace structures for Pfaffian calculations.

Pre-allocated workspace arrays to avoid repeated memory allocations
in hot loops during VMC sampling.

Includes thread-safe workspaces for parallel execution.
"""

using Base.Threads

"""
    PfaPackWorkspace

Pre-allocated workspace for Pfaffian and inverse matrix calculations.
This significantly reduces memory allocation overhead in hot loops.

# Fields
- `buf_m_real::Matrix{Float64}`: Work matrix for real calculations (n_size × n_size)
- `buf_m_complex::Matrix{ComplexF64}`: Work matrix for complex calculations (n_size × n_size)
- `iwork::Vector{Cint}`: Integer work array for pivot indices (n_size)
- `v_t_real::Vector{Float64}`: Real tridiagonal work vector (n_size - 1)
- `v_t_complex::Vector{ComplexF64}`: Complex tridiagonal work vector (n_size - 1)
- `m_work_real::Matrix{Float64}`: Real work matrix for utu2inv! (n_size × n_size)
- `m_work_complex::Matrix{ComplexF64}`: Complex work matrix for utu2inv! (n_size × n_size)
"""
mutable struct PfaPackWorkspace
    # For real calculations
    buf_m_real::Matrix{Float64}
    v_t_real::Vector{Float64}
    m_work_real::Matrix{Float64}

    # For complex calculations
    buf_m_complex::Matrix{ComplexF64}
    v_t_complex::Vector{ComplexF64}
    m_work_complex::Matrix{ComplexF64}

    # Shared (integer work array)
    iwork::Vector{Cint}

    """
        PfaPackWorkspace(n_size::Int; complex_only::Bool=false, real_only::Bool=false)

    Create a pre-allocated workspace for Pfaffian calculations.

    # Arguments
    - `n_size`: Size of the matrices (2 * number of electrons)
    - `complex_only`: If true, only allocate complex arrays
    - `real_only`: If true, only allocate real arrays
    """
    function PfaPackWorkspace(n_size::Int; complex_only::Bool=false, real_only::Bool=false)
        if complex_only && real_only
            throw(ArgumentError("Cannot specify both complex_only and real_only"))
        end

        if complex_only
            new(
                Matrix{Float64}(undef, 0, 0),
                Float64[],
                Matrix{Float64}(undef, 0, 0),
                zeros(ComplexF64, n_size, n_size),
                zeros(ComplexF64, n_size - 1),
                zeros(ComplexF64, n_size, n_size),
                zeros(Cint, n_size)
            )
        elseif real_only
            new(
                zeros(Float64, n_size, n_size),
                zeros(Float64, n_size - 1),
                zeros(Float64, n_size, n_size),
                Matrix{ComplexF64}(undef, 0, 0),
                ComplexF64[],
                Matrix{ComplexF64}(undef, 0, 0),
                zeros(Cint, n_size)
            )
        else
            new(
                zeros(Float64, n_size, n_size),
                zeros(Float64, n_size - 1),
                zeros(Float64, n_size, n_size),
                zeros(ComplexF64, n_size, n_size),
                zeros(ComplexF64, n_size - 1),
                zeros(ComplexF64, n_size, n_size),
                zeros(Cint, n_size)
            )
        end
    end
end

"""
    ThreadedPfaPackWorkspace

Thread-safe workspace container for parallel Pfaffian calculations.
Holds one PfaPackWorkspace per thread to avoid race conditions.

# Fields
- `workspaces::Vector{PfaPackWorkspace}`: One workspace per thread
- `n_size::Int`: Size of matrices (for reference)
- `complex_only::Bool`: Only allocate complex arrays
- `real_only::Bool`: Only allocate real arrays
- `lock::ReentrantLock`: Lock for thread-safe workspace expansion
"""
mutable struct ThreadedPfaPackWorkspace
    workspaces::Vector{PfaPackWorkspace}
    n_size::Int
    complex_only::Bool
    real_only::Bool
    lock::ReentrantLock

    """
        ThreadedPfaPackWorkspace(n_size::Int; complex_only::Bool=false, real_only::Bool=false)

    Create thread-local workspaces for parallel Pfaffian calculations.
    Workspaces are allocated lazily as threads access them.

    # Arguments
    - `n_size`: Size of the matrices (2 * number of electrons)
    - `complex_only`: If true, only allocate complex arrays
    - `real_only`: If true, only allocate real arrays
    """
    function ThreadedPfaPackWorkspace(n_size::Int; complex_only::Bool=false, real_only::Bool=false)
        # Pre-allocate for current number of threads, but can grow dynamically
        n_threads = max(nthreads(), 1)
        workspaces = [PfaPackWorkspace(n_size; complex_only=complex_only, real_only=real_only) for _ in 1:n_threads]
        new(workspaces, n_size, complex_only, real_only, ReentrantLock())
    end
end

"""
    get_thread_workspace(tws::ThreadedPfaPackWorkspace) -> PfaPackWorkspace

Get the workspace for the current thread. Dynamically allocates new workspaces
if the current thread ID exceeds the pre-allocated count.
"""
@inline function get_thread_workspace(tws::ThreadedPfaPackWorkspace)
    tid = threadid()
    if tid <= length(tws.workspaces)
        return @inbounds tws.workspaces[tid]
    else
        # Need to expand the workspaces vector (rare case)
        lock(tws.lock) do
            # Double-check after acquiring lock
            while tid > length(tws.workspaces)
                push!(tws.workspaces, PfaPackWorkspace(tws.n_size;
                    complex_only=tws.complex_only, real_only=tws.real_only))
            end
        end
        return @inbounds tws.workspaces[tid]
    end
end

"""
    ensure_thread_capacity!(tws::ThreadedPfaPackWorkspace)

Ensure the workspace has capacity for all current threads.
Call this at the beginning of parallel regions to avoid lock contention.
"""
function ensure_thread_capacity!(tws::ThreadedPfaPackWorkspace)
    n_threads = nthreads()
    if length(tws.workspaces) < n_threads
        lock(tws.lock) do
            while length(tws.workspaces) < n_threads
                push!(tws.workspaces, PfaPackWorkspace(tws.n_size;
                    complex_only=tws.complex_only, real_only=tws.real_only))
            end
        end
    end
end
