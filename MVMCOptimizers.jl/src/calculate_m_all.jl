"""
Calculate Pfaffian and inverse matrix from Slater elements.

Implements calculateMAll_fcmp / calculate_m_all_real / FSZ variants which
compute Pfaffian and inverse matrix for a given electron configuration.
Uses PfaPack (external dep) for LTL decomposition + utu2 routines, and
local workspace.jl (same module) for thread-local buffers.

Reference implementation:
- C code: mVMC/src/mVMC/matrix.c:285-387
"""

using LinearAlgebra
using PfaPack: cimpl_utu2inv!, utu2pfa, utu2inv!
using PfaPack: julia_zsktf2!, julia_dsktf2!, julia_zsktf2_turbo!
using Base.Threads

"""
    calculate_m_all_child_fcmp!(
        ele_idx::Vector{Int},
        slater_elm::AbstractVector{ComplexF64},
        inv_m_original::AbstractMatrix{ComplexF64},
        inv_m::Matrix{ComplexF64},
        pf_m::Ref{ComplexF64},
        buf_m::Matrix{ComplexF64},
        iwork::Vector{Int},
        work::Vector{ComplexF64},
        v_t::Vector{ComplexF64},
        m_work::Matrix{ComplexF64},
        n_site::Int,
        n_elec::Int
    )::Int

Calculate Pfaffian and inverse matrix for a single QP index (child function).

This function performs:
1. Construct invM from SlaterElm using ele_idx
2. LTL decomposition (upper triangular)
3. Pfaffian calculation
4. Inverse matrix calculation

# Arguments
- `ele_idx`: Electron indices [mi+si*Ne] (1-based)
- `slater_elm`: Slater matrix elements [ri+si*Nsite][rj+sj*Nsite] (2D array)
- `inv_m`: Output inverse matrix [mi+si*Ne][mj+sj*Ne] (will be overwritten)
- `pf_m`: Output Pfaffian value (Ref for in-place modification)
- `buf_m`: Work matrix (n_elec x n_elec)
- `iwork`: Work integer array (length n_elec)
- `work`: Work complex array (for utu2inv!)
- `v_t`: Work vector for tridiagonal elements (length n_elec-1)
- `m_work`: Work matrix for utu2inv! (n_elec x n_elec)
- `n_site`: Number of sites
- `n_elec`: Number of electrons (total, including both spins)

# Returns
- `info::Int`: Return code (0 = success, >0 = error)

# Reference
- C implementation: mVMC/src/mVMC/matrix.c:332-387 (calculateMAll_child_fcmp)
"""
function calculate_m_all_child_fcmp!(
    ele_idx::Vector{Int},
    slater_elm::AbstractVector{ComplexF64},  # 1D array with C row-major layout
    inv_m::AbstractMatrix{ComplexF64},
    pf_m::Ref{ComplexF64},
    buf_m::Matrix{ComplexF64},
    iwork::Vector{<:Integer},  # Cint for Fortran compatibility
    v_t::Vector{ComplexF64},
    m_work::Matrix{ComplexF64},
    n_site::Int,
    n_elec::Int
)::Int
    n_size = n_elec  # Nsize = 2*Ne (total electrons)
    n_site2 = 2 * n_site

    # Construct invM from SlaterElm
    # C code:
    #   for(msi=0;msi<nsize;msi++) {
    #     rsi = eleIdx[msi] + (msi/Ne)*Nsite;
    #     invM_i = invM + msi*Nsize;
    #     sltE_i = sltE + rsi*Nsite2;
    #     for(msj=0;msj<nsize;msj++) {
    #       rsj = eleIdx[msj] + (msj/Ne)*Nsite;
    #       invM_i[msj] = -sltE_i[rsj];
    #     }
    #   }
    # This means: invM[msi*Nsize + msj] = -sltE[rsi*Nsite2 + rsj]
    # In Julia with 1-based indexing: invM[msi, msj] = -slater_elm[rsi*n_site2 + rsj + 1]

    # Determine Ne (electrons per spin) from n_elec
    # n_elec should be even (2*Ne)
    ne = div(n_elec, 2)

    for msi in 1:n_size
        # Calculate rsi: rsi = ele_idx[msi] + (msi-1)/ne * n_site
        # ele_idx[msi] contains 0-based site index (ri)
        # msi is 1-based Julia index: msi = mi + si*ne + 1
        # si = div(msi - 1, ne) gives spin index (0 or 1)
        si = div(msi - 1, ne)  # spin index (0 or 1)
        ri = ele_idx[msi]  # 0-based site index
        rsi = ri + si * n_site  # 0-based rsi (ri + si*Nsite)

        # Bounds check for rsi
        if rsi < 0 || rsi >= n_site2
            @error "rsi out of bounds: rsi=$rsi, msi=$msi, ri=$ri, si=$si, n_site=$n_site"
            return 2
        end

        for msj in 1:n_size
            # Calculate rsj: rsj = ele_idx[msj] + (msj-1)/ne * n_site
            sj = div(msj - 1, ne)  # spin index (0 or 1)
            rj = ele_idx[msj]  # 0-based site index
            rsj = rj + sj * n_site  # 0-based rsj (rj + sj*Nsite)

            # Bounds check for rsj
            if rsj < 0 || rsj >= n_site2
                @error "rsj out of bounds: rsj=$rsj, msj=$msj, rj=$rj, sj=$sj, n_site=$n_site"
                return 2
            end

            # C: invM[msi*Nsize + msj] = -sltE[rsi*Nsite2 + rsj]
            # Julia: inv_m[msi, msj] = -slater_elm[rsi * n_site2 + rsj + 1]
            linear_idx = rsi * n_site2 + rsj + 1  # 1-based index

            if linear_idx > length(slater_elm)
                @error "linear_idx out of bounds: linear_idx=$linear_idx, rsi=$rsi, rsj=$rsj, n_site2=$n_site2, length(slater_elm)=$(length(slater_elm))"
                return 2
            end

            # C stores: invM[msi*Nsize + msj] = -sltE[rsi*Nsite2 + rsj]
            # C実装では、メモリ上では行優先として格納: invM[msi*Nsize + msj]
            # しかし、LAPACK/Fortranの規約では列優先として扱われる
            # LTL分解（M_ZSKTRF）は列優先を期待するため、Julia実装では列優先として格納: inv_m[msj, msi]
            # これにより、LTL分解が正しく動作する
            # その後、M_ZSCALで-1を掛け、使用時には行優先レイアウトでアクセスする
            inv_m[msj, msi] = -slater_elm[linear_idx]
        end
    end

    # Check if inv_m is all zeros or has issues
    max_abs2 = typemin(Float64)
    for i in 1:n_size
        for j in 1:n_size
            @inbounds max_abs2 = max(max_abs2, abs2(inv_m[i, j]))
        end
    end
    if max_abs2 < 1e-14 ^ 2
        return 2
    end

    # Check skew-symmetry: inv_m should be skew-symmetric (A = -A^T)
    # For a skew-symmetric matrix, diagonal elements should be zero
    # and A[i,j] = -A[j,i]
    #=
    max_skew_error = 0.0
    for i in 1:n_size
        for j in 1:n_size
            if i == j
                # Diagonal should be zero for skew-symmetric
                if abs(inv_m[i, j]) > 1e-10
                    max_skew_error = max(max_skew_error, abs(inv_m[i, j]))
                end
            else
                # Off-diagonal: A[i,j] should equal -A[j,i]
                skew_error = abs(inv_m[i, j] + inv_m[j, i])
                max_skew_error = max(max_skew_error, skew_error)
            end
        end
    end
    =#

    # LTL decomposition (upper triangular)
    # M_ZSKTRF("U", "P", &n, invM, &lda, iwork, bufM, &nsq, &info)
    # Use julia_zsktf2_turbo! for optimized SIMD vectorization (faster than Fortran for n >= 128)
    info = julia_zsktf2_turbo!(inv_m, iwork)
    if info != 0
        return info
    end

    # Copy inv_m to buf_m for LTL decomposition result
    buf_m .= inv_m

    # Calculate Pfaffian using utu2pfa from LTL-decomposed matrix
    # This avoids copying the original matrix and is more efficient
    # utu2pfa computes Pfaffian directly from the LTL-decomposed matrix
    pfaff = utu2pfa(n_size, inv_m, n_size, iwork)

    # Check if Pfaffian is finite
    if !isfinite(real(pfaff)) || !isfinite(imag(pfaff))
        return 1  # Error: non-finite Pfaffian
    end

    pf_m[] = pfaff

    # Calculate inverse matrix using utu2inv! on the LTL-decomposed matrix
    # C++: utu2inv(n, invM, lda, iPiv, vT, bufM, lda)
    # Input: inv_m contains LTL-decomposed matrix, iwork contains pivot info
    # Output: inv_m contains the inverse matrix

    # Copy LTL-decomposed matrix back to inv_m (it was already modified by julia_zsktf2_turbo!)
    # buf_m already has the LTL result, but we'll use inv_m directly
    utu2inv!(n_size, inv_m, n_size, iwork, v_t, m_work, n_size)
    #cimpl_utu2inv!(n_size, inv_m, n_size, iwork, v_t, m_work, n_size)

    # C implementation applies M_ZSCAL(&nsq, &minus_one, invM, &one) for row-major
    # to column-major conversion: InvM -> -InvM
    #
    # In Julia, we also need this sign flip because:
    # 1. Julia's inv_m is stored column-major (inv_m[msj, msi] = -sltE[...])
    # 2. Later, calculate_new_pf_m_two2! accesses inv_m as row-major using linear indexing
    #    (inv_m_arr[inv_offset + msa * n_size + msi + 1])
    # 3. This transpose effectively negates the antisymmetric matrix values
    # 4. To match C's behavior, we need to apply the sign flip here
    # Use rmul! for efficient BLAS-based scaling
    rmul!(inv_m, -1.0)

    # Debug: Print first few elements of inv_m after M_ZSCAL
    # Disabled for now as it causes issues with output
    # debug_inv = get(ENV, "MVMC_DEBUG_INV", "0") != "0"
    # if debug_inv
    #     println(stderr, "DEBUG-INV Julia: pfaff=$(real(pfaff))+$(imag(pfaff))i")
    #     for i in 1:min(4, n_size)
    #         for j in 1:min(4, n_size)
    #             println(stderr, "DEBUG-INV Julia: inv_m[$(i-1),$(j-1)] = $(real(inv_m[i, j])) + $(imag(inv_m[i, j]))i")
    #         end
    #     end
    #     flush(stderr)
    # end

    return 0
end

"""
    calculate_m_all_fcmp!(
        ele_idx::Vector{Int},
        slater_elm_1d::AbstractVector{ComplexF64},
        inv_m_array::AbstractArray{ComplexF64, 3},
        pf_m_array::AbstractVector{ComplexF64},
        qp_start::Int,
        qp_end::Int,
        n_site::Int,
        n_elec::Int,
        workspace::PfaPackWorkspace
    )::Int

Calculate Pfaffian and inverse matrix for multiple QP indices.

This version uses pre-allocated workspace arrays to avoid memory allocation overhead.

# Arguments
- `ele_idx`: Electron indices [mi+si*Ne] (1-based)
- `slater_elm_1d`: Slater matrix elements as 1D array with C row-major layout
                   [qp_idx * Nsite2 * Nsite2 + rsi * Nsite2 + rsj]
- `inv_m_array`: Output inverse matrices [qp_idx][mi+si*Ne][mj+sj*Ne] (3D array, will be overwritten)
- `pf_m_array`: Output Pfaffian values [qp_idx] (will be overwritten)
- `qp_start`, `qp_end`: Define the QP count via `qp_num = qp_end - qp_start`.
                       `slater_elm_1d`, `inv_m_array`, and `pf_m_array` are
                       expected to be already sliced to cover the
                       `[qp_start, qp_end)` range, so this function reads
                       and writes them from index 1 regardless of the
                       `qp_start` value. (The C reference
                       `CalculateMAll_*` addresses arrays by absolute
                       `qp_start`/`qp_end`; the Julia port shifts that
                       responsibility to the caller.)
- `n_site`: Number of sites
- `n_elec`: Number of electrons (total, including both spins)
- `workspace`: Pre-allocated workspace for avoiding allocations

# Returns
- `info::Int`: Return code (0 = success, >0 = error)

# Reference
- C implementation: mVMC/src/mVMC/matrix.c:285-330 (CalculateMAll_fcmp)
"""
function calculate_m_all_fcmp!(
    ele_idx::Vector{Int},
    slater_elm_1d::AbstractVector{ComplexF64},
    inv_m_array::AbstractArray{ComplexF64, 3},
    pf_m_array::AbstractVector{ComplexF64},
    qp_start::Int,
    qp_end::Int,
    n_site::Int,
    n_elec::Int,
    workspace::PfaPackWorkspace
)::Int
    qp_num = qp_end - qp_start
    n_size = n_elec
    n_site2 = 2 * n_site

    # Use pre-allocated work arrays from workspace
    buf_m = workspace.buf_m_complex
    iwork = workspace.iwork
    v_t = workspace.v_t_complex
    m_work = workspace.m_work_complex

    info = 0

    # Process each QP index
    # slater_elm_1d is already the subset for the relevant QP range
    # So we iterate from 0 to qp_num-1 (0-based offset into the subset)
    for qpidx in 1:qp_num
        if info != 0
            continue
        end

        # Calculate the offset for this QP index in the 1D slater_elm array subset
        # Since slater_elm_1d is already the subset, we use (qpidx - 1) as 0-based offset
        qp_idx_0based = qpidx - 1
        slater_elm_offset = qp_idx_0based * n_site2 * n_site2

        # Create a view into the 1D slater_elm array for this QP
        slater_elm = view(slater_elm_1d, (slater_elm_offset + 1):(slater_elm_offset + n_site2 * n_site2))

        # Get inverse matrix for this QP index
        # inv_m_array has shape (n_size, n_size, qp_num) - Julia column-major
        inv_m = view(inv_m_array, :, :, qpidx)
        # Pfaffian value (as Ref for in-place modification)
        pf_m = Ref(pf_m_array[qpidx])

        # Call child function
        my_info = calculate_m_all_child_fcmp!(
            ele_idx,
            slater_elm,
            inv_m,
            pf_m,
            buf_m,
            iwork,
            v_t,
            m_work,
            n_site,
            n_elec
        )

        if my_info != 0
            info = my_info
        end

        # Update pf_m_array
        pf_m_array[qpidx] = pf_m[]
    end

    return info
end

"""
    calculate_m_all_fcmp!(
        ele_idx::Vector{Int},
        slater_elm_1d::AbstractVector{ComplexF64},
        inv_m_array::AbstractArray{ComplexF64, 3},
        pf_m_array::Vector{ComplexF64},
        qp_start::Int,
        qp_end::Int,
        n_site::Int,
        n_elec::Int
    )::Int

Calculate Pfaffian and inverse matrix for multiple QP indices.

This version allocates work arrays internally (for backward compatibility).
For better performance in hot loops, use the version that accepts a PfaPackWorkspace.

# Arguments
- `ele_idx`: Electron indices [mi+si*Ne] (1-based)
- `slater_elm_1d`: Slater matrix elements as 1D array with C row-major layout
                   [qp_idx * Nsite2 * Nsite2 + rsi * Nsite2 + rsj]
- `inv_m_array`: Output inverse matrices [qp_idx][mi+si*Ne][mj+sj*Ne] (3D array, will be overwritten)
- `pf_m_array`: Output Pfaffian values [qp_idx] (will be overwritten)
- `qp_start`, `qp_end`: Define the QP count via `qp_num = qp_end - qp_start`.
                       `slater_elm_1d`, `inv_m_array`, and `pf_m_array` are
                       expected to be already sliced to cover the
                       `[qp_start, qp_end)` range, so this function reads
                       and writes them from index 1 regardless of the
                       `qp_start` value. (The C reference
                       `CalculateMAll_*` addresses arrays by absolute
                       `qp_start`/`qp_end`; the Julia port shifts that
                       responsibility to the caller.)
- `n_site`: Number of sites
- `n_elec`: Number of electrons (total, including both spins)

# Returns
- `info::Int`: Return code (0 = success, >0 = error)
"""
function calculate_m_all_fcmp!(
    ele_idx::Vector{Int},
    slater_elm_1d::AbstractVector{ComplexF64},
    inv_m_array::AbstractArray{ComplexF64, 3},
    pf_m_array::AbstractVector{ComplexF64},
    qp_start::Int,
    qp_end::Int,
    n_site::Int,
    n_elec::Int
)::Int
    # Create temporary workspace (allocates memory)
    workspace = PfaPackWorkspace(n_elec; complex_only=true)
    return calculate_m_all_fcmp!(
        ele_idx, slater_elm_1d, inv_m_array, pf_m_array,
        qp_start, qp_end, n_site, n_elec, workspace
    )
end

"""
    calculate_m_all_fcmp!(... , threaded_workspace::ThreadedPfaPackWorkspace)

Calculate Pfaffian and inverse matrix with multi-threaded parallelization (complex version).
"""
function calculate_m_all_fcmp!(
    ele_idx::Vector{Int},
    slater_elm_1d::AbstractVector{ComplexF64},
    inv_m_array::AbstractArray{ComplexF64, 3},
    pf_m_array::AbstractVector{ComplexF64},
    qp_start::Int,
    qp_end::Int,
    n_site::Int,
    n_elec::Int,
    threaded_workspace::ThreadedPfaPackWorkspace
)::Int
    qp_num = qp_end - qp_start
    n_size = n_elec
    n_site2 = 2 * n_site

    # Fall back to sequential execution if only 1 thread or workload is small
    # This avoids @threads overhead when parallelization won't help
    n_threads = nthreads()
    if n_threads == 1 || qp_num < n_threads
        # Use sequential execution with first workspace
        workspace = threaded_workspace.workspaces[1]
        return calculate_m_all_fcmp!(
            ele_idx, slater_elm_1d, inv_m_array, pf_m_array,
            qp_start, qp_end, n_site, n_elec, workspace
        )
    end

    # Ensure workspace has capacity for all threads (avoids lock contention in loop)
    ensure_thread_capacity!(threaded_workspace)

    # Use atomic for thread-safe error handling
    info = Atomic{Int}(0)

    # Process each QP index in parallel
    @threads for qpidx in 1:qp_num
        # Skip if error already occurred
        if info[] != 0
            continue
        end

        # Get thread-local workspace
        workspace = get_thread_workspace(threaded_workspace)
        buf_m = workspace.buf_m_complex
        iwork = workspace.iwork
        v_t = workspace.v_t_complex
        m_work = workspace.m_work_complex

        # Calculate the offset for this QP index in the 1D slater_elm array subset
        qp_idx_0based = qpidx - 1
        slater_elm_offset = qp_idx_0based * n_site2 * n_site2

        # Create a view into the 1D slater_elm array for this QP
        slater_elm = view(slater_elm_1d, (slater_elm_offset + 1):(slater_elm_offset + n_site2 * n_site2))

        # Get inverse matrix for this QP index
        inv_m = view(inv_m_array, :, :, qpidx)
        # Pfaffian value (as Ref for in-place modification)
        pf_m = Ref(pf_m_array[qpidx])

        # Call child function
        my_info = calculate_m_all_child_fcmp!(
            ele_idx,
            slater_elm,
            inv_m,
            pf_m,
            buf_m,
            iwork,
            v_t,
            m_work,
            n_site,
            n_elec
        )

        if my_info != 0
            atomic_cas!(info, 0, my_info)
        end

        # Update pf_m_array (thread-safe since each thread writes to different index)
        pf_m_array[qpidx] = pf_m[]
    end

    return info[]
end

"""
    calculate_m_all_child_real!(
        ele_idx::Vector{Int},
        slater_elm::AbstractVector{Float64},  # 1D array with C row-major layout
        inv_m::AbstractMatrix{Float64},
        pf_m::Ref{Float64},
        buf_m::Matrix{Float64},
        iwork::Vector{<:Integer},  # Cint for Fortran compatibility
        v_t::Vector{Float64},
        m_work::Matrix{Float64},
        n_site::Int,
        n_elec::Int
    )::Int

Calculate Pfaffian and inverse matrix for a single QP index using real arithmetic (child function).

This function performs:
1. Construct invM from SlaterElm_real using ele_idx
2. LTL decomposition (upper triangular)
3. Pfaffian calculation
4. Inverse matrix calculation

# Arguments
- `ele_idx`: Electron indices [mi+si*Ne] (1-based)
- `slater_elm`: Slater matrix elements [ri+si*Nsite][rj+sj*Nsite] (2D array, real)
- `inv_m`: Output inverse matrix [mi+si*Ne][mj+sj*Ne] (will be overwritten, real)
- `pf_m`: Output Pfaffian value (Ref for in-place modification, real)
- `buf_m`: Work matrix (n_elec x n_elec, real)
- `iwork`: Work integer array (length n_elec)
- `v_t`: Work vector for tridiagonal elements (length n_elec-1, real)
- `m_work`: Work matrix for utu2inv! (n_elec x n_elec, real)
- `n_site`: Number of sites
- `n_elec`: Number of electrons (total, including both spins)

# Returns
- `info::Int`: Return code (0 = success, >0 = error)

# Reference
- C implementation: mVMC/src/mVMC/matrix.c:567-624 (calculateMAll_child_real)
"""
function calculate_m_all_child_real!(
    ele_idx::Vector{Int},
    slater_elm::AbstractVector{Float64},  # 1D array with C row-major layout
    inv_m::AbstractMatrix{Float64},
    pf_m::Ref{Float64},
    buf_m::Matrix{Float64},
    iwork::Vector{<:Integer},  # Cint for Fortran compatibility
    v_t::Vector{Float64},
    m_work::Matrix{Float64},
    n_site::Int,
    n_elec::Int
)::Int
    n_size = n_elec  # Nsize = 2*Ne (total electrons)
    n_site2 = 2 * n_site

    # Construct invM from SlaterElm_real
    # C code:
    #   for(msi=0;msi<nsize;msi++) {
    #     rsi = eleIdx[msi] + (msi/Ne)*Nsite;
    #     invM_i = invM + msi*Nsize;
    #     sltE_i = sltE + rsi*Nsite2;
    #     for(msj=0;msj<nsize;msj++) {
    #       rsj = eleIdx[msj] + (msj/Ne)*Nsite;
    #       invM_i[msj] = -sltE_i[rsj];
    #     }
    #   }

    # Determine Ne (electrons per spin) from n_elec
    # n_elec should be even (2*Ne)
    ne = div(n_elec, 2)

    for msi in 1:n_size
        # Calculate rsi: rsi = ele_idx[msi] + (msi-1)/ne * n_site
        si = div(msi - 1, ne)  # spin index (0 or 1)
        ri = ele_idx[msi]  # 0-based site index
        rsi = ri + si * n_site  # 0-based rsi (ri + si*Nsite)

        # Bounds check for rsi
        if rsi < 0 || rsi >= n_site2
            @error "rsi out of bounds: rsi=$rsi, msi=$msi, ri=$ri, si=$si, n_site=$n_site"
            return 2
        end

        for msj in 1:n_size
            # Calculate rsj: rsj = ele_idx[msj] + (msj-1)/ne * n_site
            sj = div(msj - 1, ne)  # spin index (0 or 1)
            rj = ele_idx[msj]  # 0-based site index
            rsj = rj + sj * n_site  # 0-based rsj (rj + sj*Nsite)

            # Bounds check for rsj
            if rsj < 0 || rsj >= n_site2
                @error "rsj out of bounds: rsj=$rsj, msj=$msj, rj=$rj, sj=$sj, n_site=$n_site"
                return 2
            end

            # C: invM[msi*Nsize + msj] = -sltE[rsi*Nsite2 + rsj]
            # Julia: inv_m[msi, msj] = -slater_elm[rsi * n_site2 + rsj + 1]
            linear_idx = rsi * n_site2 + rsj + 1  # 1-based index

            if linear_idx > length(slater_elm)
                @error "linear_idx out of bounds: linear_idx=$linear_idx, rsi=$rsi, rsj=$rsj, n_site2=$n_site2, length(slater_elm)=$(length(slater_elm))"
                return 2
            end

            # Store in column-major format for LTL decomposition
            inv_m[msj, msi] = -slater_elm[linear_idx]
        end
    end

    # Check if inv_m is all zeros or has issues
    max_abs2 = typemin(Float64)
    for i in 1:n_size
        for j in 1:n_size
            @inbounds max_abs2 = max(max_abs2, abs2(inv_m[i, j]))
        end
    end
    if max_abs2 < 1e-14 ^ 2
        return 2
    end

    # LTL decomposition (upper triangular)
    # M_DSKTRF("U", "N", &n, invM, &lda, iwork, bufM, &nsq, &info)
    # Use optimized Julia DSKTF2 (same performance as Fortran)
    info = julia_dsktf2!(inv_m, iwork)
    if info != 0
        return info
    end

    # Copy inv_m to buf_m for LTL decomposition result
    buf_m .= inv_m

    # Calculate Pfaffian using utu2pfa from LTL-decomposed matrix
    pfaff = utu2pfa(n_size, inv_m, n_size, iwork)

    # Check if Pfaffian is finite
    if !isfinite(pfaff)
        return 1  # Error: non-finite Pfaffian
    end

    pf_m[] = pfaff

    # Calculate inverse matrix using utu2inv! on the LTL-decomposed matrix
    utu2inv!(n_size, inv_m, n_size, iwork, v_t, m_work, n_size)

    # C implementation applies M_DSCAL(&nsq, &minus_one, invM, &one)
    # InvM -> InvM' = -InvM
    # Use rmul! for efficient BLAS-based scaling
    rmul!(inv_m, -1.0)

    return 0
end

"""
    calculate_m_all_real!(
        ele_idx::Vector{Int},
        slater_elm_1d::AbstractVector{Float64},
        inv_m_array::AbstractArray{Float64, 3},
        pf_m_array::AbstractVector{Float64},
        qp_start::Int,
        qp_end::Int,
        n_site::Int,
        n_elec::Int,
        workspace::PfaPackWorkspace
    )::Int

Calculate Pfaffian and inverse matrix for multiple QP indices using real arithmetic.

This version uses pre-allocated workspace arrays to avoid memory allocation overhead.

# Arguments
- `ele_idx`: Electron indices [mi+si*Ne] (1-based)
- `slater_elm_1d`: Slater matrix elements as 1D array with C row-major layout (real)
                   [qp_idx * Nsite2 * Nsite2 + rsi * Nsite2 + rsj]
- `inv_m_array`: Output inverse matrices [qp_idx][mi+si*Ne][mj+sj*Ne] (3D array, will be overwritten, real)
- `pf_m_array`: Output Pfaffian values [qp_idx] (will be overwritten, real)
- `qp_start`, `qp_end`: Define the QP count via `qp_num = qp_end - qp_start`.
                       `slater_elm_1d`, `inv_m_array`, and `pf_m_array` are
                       expected to be already sliced to cover the
                       `[qp_start, qp_end)` range, so this function reads
                       and writes them from index 1 regardless of the
                       `qp_start` value. (The C reference
                       `CalculateMAll_*` addresses arrays by absolute
                       `qp_start`/`qp_end`; the Julia port shifts that
                       responsibility to the caller.)
- `n_site`: Number of sites
- `n_elec`: Number of electrons (total, including both spins)
- `workspace`: Pre-allocated workspace for avoiding allocations

# Returns
- `info::Int`: Return code (0 = success, >0 = error)

# Reference
- C implementation: mVMC/src/mVMC/matrix.c:526-565 (CalculateMAll_real)
"""
function calculate_m_all_real!(
    ele_idx::Vector{Int},
    slater_elm_1d::AbstractVector{Float64},
    inv_m_array::AbstractArray{Float64, 3},
    pf_m_array::AbstractVector{Float64},
    qp_start::Int,
    qp_end::Int,
    n_site::Int,
    n_elec::Int,
    workspace::PfaPackWorkspace
)::Int
    qp_num = qp_end - qp_start
    n_size = n_elec
    n_site2 = 2 * n_site

    # Use pre-allocated work arrays from workspace
    buf_m = workspace.buf_m_real
    iwork = workspace.iwork
    v_t = workspace.v_t_real
    m_work = workspace.m_work_real

    info = 0

    # Process each QP index
    for qpidx in 1:qp_num
        if info != 0
            continue
        end

        # Calculate the offset for this QP index in the 1D slater_elm array subset
        qp_idx_0based = qpidx - 1
        slater_elm_offset = qp_idx_0based * n_site2 * n_site2

        # Create a view into the 1D slater_elm array for this QP
        slater_elm = view(slater_elm_1d, (slater_elm_offset + 1):(slater_elm_offset + n_site2 * n_site2))

        # Get inverse matrix for this QP index
        inv_m = view(inv_m_array, :, :, qpidx)
        # Pfaffian value (as Ref for in-place modification)
        pf_m = Ref(pf_m_array[qpidx])

        # Call child function
        my_info = calculate_m_all_child_real!(
            ele_idx,
            slater_elm,
            inv_m,
            pf_m,
            buf_m,
            iwork,
            v_t,
            m_work,
            n_site,
            n_elec
        )

        if my_info != 0
            info = my_info
        end

        # Update pf_m_array
        pf_m_array[qpidx] = pf_m[]
    end

    return info
end

"""
    calculate_m_all_real!(
        ele_idx::Vector{Int},
        slater_elm_1d::AbstractVector{Float64},
        inv_m_array::AbstractArray{Float64, 3},
        pf_m_array::Vector{Float64},
        qp_start::Int,
        qp_end::Int,
        n_site::Int,
        n_elec::Int
    )::Int

Calculate Pfaffian and inverse matrix for multiple QP indices using real arithmetic.

This version allocates work arrays internally (for backward compatibility).
For better performance in hot loops, use the version that accepts a PfaPackWorkspace.

# Arguments
- `ele_idx`: Electron indices [mi+si*Ne] (1-based)
- `slater_elm_1d`: Slater matrix elements as 1D array with C row-major layout (real)
                   [qp_idx * Nsite2 * Nsite2 + rsi * Nsite2 + rsj]
- `inv_m_array`: Output inverse matrices [qp_idx][mi+si*Ne][mj+sj*Ne] (3D array, will be overwritten, real)
- `pf_m_array`: Output Pfaffian values [qp_idx] (will be overwritten, real)
- `qp_start`, `qp_end`: Define the QP count via `qp_num = qp_end - qp_start`.
                       `slater_elm_1d`, `inv_m_array`, and `pf_m_array` are
                       expected to be already sliced to cover the
                       `[qp_start, qp_end)` range, so this function reads
                       and writes them from index 1 regardless of the
                       `qp_start` value. (The C reference
                       `CalculateMAll_*` addresses arrays by absolute
                       `qp_start`/`qp_end`; the Julia port shifts that
                       responsibility to the caller.)
- `n_site`: Number of sites
- `n_elec`: Number of electrons (total, including both spins)

# Returns
- `info::Int`: Return code (0 = success, >0 = error)
"""
function calculate_m_all_real!(
    ele_idx::Vector{Int},
    slater_elm_1d::AbstractVector{Float64},
    inv_m_array::AbstractArray{Float64, 3},
    pf_m_array::AbstractVector{Float64},
    qp_start::Int,
    qp_end::Int,
    n_site::Int,
    n_elec::Int
)::Int
    # Create temporary workspace (allocates memory)
    workspace = PfaPackWorkspace(n_elec; real_only=true)
    return calculate_m_all_real!(
        ele_idx, slater_elm_1d, inv_m_array, pf_m_array,
        qp_start, qp_end, n_site, n_elec, workspace
    )
end

"""
    calculate_m_all_real!(
        ele_idx::Vector{Int},
        slater_elm_1d::AbstractVector{Float64},
        inv_m_array::AbstractArray{Float64, 3},
        pf_m_array::AbstractVector{Float64},
        qp_start::Int,
        qp_end::Int,
        n_site::Int,
        n_elec::Int,
        threaded_workspace::ThreadedPfaPackWorkspace
    )::Int

Calculate Pfaffian and inverse matrix for multiple QP indices using real arithmetic
with multi-threaded parallelization.

This version uses thread-local workspaces for parallel execution with Threads.@threads.

# Arguments
- `ele_idx`: Electron indices [mi+si*Ne] (1-based)
- `slater_elm_1d`: Slater matrix elements as 1D array with C row-major layout (real)
- `inv_m_array`: Output inverse matrices (3D array, will be overwritten, real)
- `pf_m_array`: Output Pfaffian values (will be overwritten, real)
- `qp_start`, `qp_end`: Define the QP count via `qp_num = qp_end - qp_start`.
                       `slater_elm_1d`, `inv_m_array`, and `pf_m_array` are
                       expected to be already sliced to cover the
                       `[qp_start, qp_end)` range, so this function reads
                       and writes them from index 1 regardless of the
                       `qp_start` value. (The C reference
                       `CalculateMAll_*` addresses arrays by absolute
                       `qp_start`/`qp_end`; the Julia port shifts that
                       responsibility to the caller.)
- `n_site`: Number of sites
- `n_elec`: Number of electrons (total, including both spins)
- `threaded_workspace`: Thread-local workspaces for parallel execution

# Returns
- `info::Int`: Return code (0 = success, >0 = error)
"""
function calculate_m_all_real!(
    ele_idx::Vector{Int},
    slater_elm_1d::AbstractVector{Float64},
    inv_m_array::AbstractArray{Float64, 3},
    pf_m_array::AbstractVector{Float64},
    qp_start::Int,
    qp_end::Int,
    n_site::Int,
    n_elec::Int,
    threaded_workspace::ThreadedPfaPackWorkspace
)::Int
    qp_num = qp_end - qp_start
    n_size = n_elec
    n_site2 = 2 * n_site

    # Fall back to sequential execution if only 1 thread or workload is small
    # This avoids @threads overhead when parallelization won't help
    n_threads = nthreads()
    if n_threads == 1 || qp_num < n_threads
        # Use sequential execution with first workspace
        workspace = threaded_workspace.workspaces[1]
        return calculate_m_all_real!(
            ele_idx, slater_elm_1d, inv_m_array, pf_m_array,
            qp_start, qp_end, n_site, n_elec, workspace
        )
    end

    # Ensure workspace has capacity for all threads
    ensure_thread_capacity!(threaded_workspace)

    # Use atomic for thread-safe error handling
    info = Atomic{Int}(0)

    # Process each QP index in parallel
    @threads for qpidx in 1:qp_num
        # Skip if error already occurred
        if info[] != 0
            continue
        end

        # Get thread-local workspace
        workspace = get_thread_workspace(threaded_workspace)
        buf_m = workspace.buf_m_real
        iwork = workspace.iwork
        v_t = workspace.v_t_real
        m_work = workspace.m_work_real

        # Calculate the offset for this QP index in the 1D slater_elm array subset
        qp_idx_0based = qpidx - 1
        slater_elm_offset = qp_idx_0based * n_site2 * n_site2

        # Create a view into the 1D slater_elm array for this QP
        slater_elm = view(slater_elm_1d, (slater_elm_offset + 1):(slater_elm_offset + n_site2 * n_site2))

        # Get inverse matrix for this QP index
        inv_m = view(inv_m_array, :, :, qpidx)
        # Pfaffian value (as Ref for in-place modification)
        pf_m = Ref(pf_m_array[qpidx])

        # Call child function
        my_info = calculate_m_all_child_real!(
            ele_idx,
            slater_elm,
            inv_m,
            pf_m,
            buf_m,
            iwork,
            v_t,
            m_work,
            n_site,
            n_elec
        )

        if my_info != 0
            atomic_cas!(info, 0, my_info)
        end

        # Update pf_m_array (thread-safe since each thread writes to different index)
        pf_m_array[qpidx] = pf_m[]
    end

    return info[]
end

