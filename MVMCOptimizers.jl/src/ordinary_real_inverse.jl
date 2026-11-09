# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Private real translation of RuQing Xu's invert.tcc, via PfaPack Utu2.
# See source-adjacent ordinary_real_inverse.NOTICE and .LICENSE.
# Unit-upper DTRTI2 column algorithm, independently expressed in Julia.
# LAPACK DTRTI2 calls DTRMV(U,N,U) then DSCAL(-1) for each column.
# The diagonal and lower triangle must remain untouched: they store other
# parts of the skew factorization. Only independent i entries are vectorized.
function _ordinary_unit_upper_trtri_real!(A::AbstractMatrix{Float64})
    Base.require_one_based_indexing(A)
    n = size(A, 1)
    size(A, 2) == n || throw(DimensionMismatch("triangular inverse must be square"))
    @inbounds for j in 2:n
        for k in 1:j-1
            temp = A[k,j]
            if temp != 0.0
                @simd ivdep for i in 1:k-1
                    A[i,j] += temp * A[i,k]
                end
            end
        end
        @simd ivdep for i in 1:j-1
            A[i,j] = -A[i,j]
        end
    end
    return A
end

function _ordinary_sktdsmx_real_c_order!(n::Int, vt::Vector{Float64},
                                        B::AbstractMatrix{Float64}, C::AbstractMatrix{Float64})
    Base.require_one_based_indexing(vt, B, C)
    n >= 2 && iseven(n) || throw(ArgumentError("solve requires a positive even dimension"))
    size(B) == (n,n) && size(C) == (n,n) || throw(DimensionMismatch("solve workspace dimensions"))
    length(vt) >= n-1 || throw(DimensionMismatch("solve vector workspace"))
    # Process 16-column blocks. Each block's three row slices stay in L1/L2
    # (the full 128-column pass streams ~1.5 MB of strided row data per plane);
    # this is 1.6x faster serially and 1.3x at 8 workers on L=128.
    BLK = 16
    for j0 in 1:BLK:n
        j1 = min(j0 + BLK - 1, n)
        @inbounds @simd ivdep for j in j0:j1
            C[2,j] = B[1,j] / -vt[1]
        end
        @inbounds for i in 2:2:n-2
            @simd ivdep for j in j0:j1
                C[i+2,j] = (B[i+1,j] - C[i,j] * vt[i]) / -vt[i+1]
            end
        end
        @inbounds @simd ivdep for j in j0:j1
            C[n-1,j] = B[n,j] / vt[n-1]
        end
        @inbounds for i in n-3:-2:1
            @simd ivdep for j in j0:j1
                C[i,j] = (B[i+1,j] + C[i+2,j] * vt[i+1]) / vt[i]
            end
        end
    end
    return nothing
end

const _ORDINARY_TRMMT_PANEL = 64

# Panelized upper-triangle form of `A ← Mᵀ A` for a unit upper-triangular M.
# Mirrors the Rust port's `trmmt_upper_lutu` (deps/trmmt.tcc): each 64-column
# panel only needs the leading `rows × rows` block of M and the leading `rows`
# rows of A, so the product stays cache resident and touches half the entries.
function _ordinary_trmmt_upper_real!(M::AbstractMatrix{Float64},
                                     A::AbstractMatrix{Float64}, n::Int)
    @inbounds for panel_col in 0:_ORDINARY_TRMMT_PANEL:(n-1)
        nloc = min(_ORDINARY_TRMMT_PANEL, n - panel_col)
        rows = panel_col + nloc
        Mblock = @view M[1:rows, 1:rows]
        Ablock = @view A[1:rows, (panel_col+1):(panel_col+nloc)]
        LinearAlgebra.BLAS.trmm!('L', 'U', 'T', 'U', 1.0, Mblock, Ablock)
    end
    return nothing
end

# Restore the strict lower triangle from skew-symmetry after the upper-only
# panel product, row by row in 8-column blocks (the Rust `fill_lower_from_upper_skew`
# order). The diagonal becomes exactly zero.
function _ordinary_fill_lower_from_upper_skew!(A::AbstractMatrix{Float64}, n::Int)
    @inbounds for j0 in 0:8:(n-1)
        j1 = min(j0 + 8, n)
        for j in j0:(j1-1)
            A[j+1, j+1] = 0.0
        end
        for i in (j0+1):(n-1)
            jmax = min(j1 - 1, i - 1)
            for j in j0:jmax
                A[i+1, j+1] = -A[j+1, i+1]
            end
        end
    end
    return nothing
end

function _ordinary_utu2inv_real_c_order!(n::Int, A::AbstractMatrix{Float64}, ldA::Int,
                                        ipiv::Vector{<:Integer}, vt::Vector{Float64},
                                        M::AbstractMatrix{Float64}, ldM::Int)
    Base.require_one_based_indexing(A, ipiv, vt, M)
    n >= 2 && iseven(n) || throw(ArgumentError("inverse requires a positive even dimension"))
    size(A) == (n,n) && size(M) == (n,n) || throw(DimensionMismatch("inverse workspace dimensions"))
    ldA >= n && ldM >= n || throw(DimensionMismatch("inverse leading dimensions"))
    length(ipiv) >= n && length(vt) >= n-1 || throw(DimensionMismatch("inverse vector workspace"))
    all(i -> 1 <= ipiv[i] <= n, 1:n) || throw(ArgumentError("inverse pivot outside matrix"))
    # Same existing Julia/C full pipeline; only the scalar solver differs.
    # `sktdsmx` reads the lower triangle and diagonal of `M` (they are not
    # written by the superdiagonal copy below), so the identity fill is live.
    fill!(M, 0.0)
    @inbounds for i in 1:n
        M[i,i] = 1.0
    end
    # Optimize the measured matrix sizes. Tiny SR-CG systems amplify provider
    # roundoff: on macOS Julia 1.11, a 2.8e-17 inverse difference at n=6 changed
    # the first parameter update. Retain their established LAPACK provider.
    if 32 <= n <= 64
        _ordinary_unit_upper_trtri_real!(@view(A[1:n-1,2:n]))
    else
        LinearAlgebra.LAPACK.trtri!('U', 'U', @view(A[1:n-1,2:n]))
    end
    @inbounds for j in 1:n-2, i in 1:j
        M[i,j+1] = A[i,j+2]
    end
    @inbounds for i in 1:n-1
        vt[i] = -A[i,i+1]
    end
    _ordinary_sktdsmx_real_c_order!(n, vt, M, A)
    # For n > 64, compute only the upper triangle of Mᵀ A in 64-column panels
    # and restore the lower triangle from skew-symmetry (Rust `utu2inv`). The
    # column permutation below acts on the right and commutes with this product,
    # so applying it before or after leaves every column's arithmetic unchanged.
    panel_trmmt = n > _ORDINARY_TRMMT_PANEL
    if panel_trmmt
        _ordinary_trmmt_upper_real!(M, A, n)
        _ordinary_fill_lower_from_upper_skew!(A, n)
    end
    @inbounds for j in 1:n
        target = ipiv[j]
        if target != j
            for i in 1:n
                A[i,j], A[i,target] = A[i,target], A[i,j]
            end
        end
    end
    if !panel_trmmt
        LinearAlgebra.BLAS.trmm!('L', 'U', 'T', 'U', 1.0, M, A)
    end
    # Row interchanges, column-major: keep the inner loop on one contiguous
    # column instead of striding through rows (the Rust `permute_rows_forward`
    # order). Columns are independent and receive the same swap sequence.
    @inbounds for j in 1:n
        for i in 1:n
            target = ipiv[i]
            if target != i
                A[i,j], A[target,j] = A[target,j], A[i,j]
            end
        end
    end
    return nothing
end
