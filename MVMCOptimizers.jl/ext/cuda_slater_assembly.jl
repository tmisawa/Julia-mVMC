# Resident Slater assembly follows the C four-spin-block expressions directly,
# avoiding basis GEMV reassociation. Source: mVMC 1.4.0 src/mVMC/slater.c:38-91 (unchanged from 1.3.0);
# CPU counterpart: slater_update.jl:_update_slater_elm_fcmp_kernel!.
struct SlaterAssembly
    parameter_a::CuArray{Int32,3}
    parameter_b::CuArray{Int32,3}
    sign_a::CuArray{Int32,3}
    sign_b::CuArray{Int32,3}
    cs::CuVector{Float64}
    cc::CuVector{Float64}
    ss::CuVector{Float64}
    nsite::Int
    nqp::Int
end

function SlaterAssembly(data::ExpertModeData,nproj::Integer,ns::Integer,nq::Integer)
    ns == data.modpara.nsite || throw(DimensionMismatch("Slater assembly site count"))
    data.i_flg_orbital_general == 0 || throw(ArgumentError("Slater assembly requires normal orbitals"))
    nsp = max(1,data.modpara.nsp_gauss_leg)
    nmp = max(1,abs(data.modpara.nmp_trans))
    nopt = max(1,data.n_qp_opt_trans)
    nq == nsp*nmp*nopt || throw(DimensionMismatch("Slater assembly projection count"))
    oi,os,_ = CPU.build_orbital_idx_sgn_matrices(data)
    qt,_,qsg,qo,qosg = CPU.build_qp_trans_matrices(data)
    length(qt) >= nmp && length(qo) >= nopt || throw(DimensionMismatch("Slater assembly translations"))
    params = CPU.pack_parameters(data)
    all(isreal,params) || throw(ArgumentError("Slater assembly requires real parameters"))
    all(i -> 0 <= i < length(params)-nproj,oi) || throw(ArgumentError("incomplete Slater orbital mapping"))
    weights = data.qp_weights
    a = Array{Int32}(undef,ns,ns,nq)
    b = similar(a); sa = similar(a); sb = similar(a)
    cs = Vector{Float64}(undef,nq); cc = similar(cs); ss = similar(cs)
    for q in 1:nq
        z = q-1; opt = z÷(nsp*nmp)+1; mp = (z%(nsp*nmp))÷nsp+1; sp = z%nsp+1
        cs[q] = weights.spgl_cos_sin[sp]
        cc[q] = weights.spgl_cos_cos[sp]
        ss[q] = weights.spgl_sin_sin[sp]
        for ri in 1:ns, rj in 1:ns
            ori = qo[opt][ri]+1; orj = qo[opt][rj]+1
            tri = qt[mp][ori]+1; trj = qt[mp][orj]+1
            sgni = qsg[mp][ori]*qosg[opt][ri]
            sgnj = qsg[mp][orj]*qosg[opt][rj]
            a[ri,rj,q] = nproj+oi[tri,trj]+1
            b[ri,rj,q] = nproj+oi[trj,tri]+1
            sa[ri,rj,q] = os[tri,trj]*sgni*sgnj
            sb[ri,rj,q] = os[trj,tri]*sgni*sgnj
        end
    end
    return SlaterAssembly(CuArray(a),CuArray(b),CuArray(sa),CuArray(sb),
                          CuArray(cs),CuArray(cc),CuArray(ss),Int(ns),Int(nq))
end

function direct_slater_kernel!(out,params,pa,pb,sa,sb,cs,cc,ss,ns,nq)
    index = (blockIdx().x-1)*blockDim().x+threadIdx().x-1
    index >= ns*ns*nq && return nothing
    ri = index%ns+1
    rj = (index÷ns)%ns+1
    q = index÷(ns*ns)+1
    @inbounds begin
        a = CUDA.mul_rn(params[pa[ri,rj,q]],Float64(sa[ri,rj,q]))
        b = CUDA.mul_rn(params[pb[ri,rj,q]],Float64(sb[ri,rj,q]))
        difference = CUDA.add_rn(a,-b)
        # Explicit RN intrinsics prevent multiplication/addition contraction.
        # CuArray column-major indexing is the transpose of C's flat layout.
        out[rj,ri,q] = CUDA.mul_rn(-difference,cs[q])
        out[rj+ns,ri,q] = CUDA.add_rn(CUDA.mul_rn(a,cc[q]),CUDA.mul_rn(b,ss[q]))
        out[rj,ri+ns,q] = CUDA.add_rn(CUDA.mul_rn(-a,ss[q]),-CUDA.mul_rn(b,cc[q]))
        out[rj+ns,ri+ns,q] = CUDA.mul_rn(difference,cs[q])
    end
    return nothing
end

function assemble_slater!(s,assembly::SlaterAssembly)
    ns,nq = assembly.nsite,assembly.nqp
    size(s.slater) == (2ns,2ns,nq) || throw(DimensionMismatch("resident Slater output"))
    entries = ns*ns*nq
    @cuda threads=256 blocks=cld(entries,256) fastmath=false direct_slater_kernel!(
        s.slater,s.parameters,assembly.parameter_a,assembly.parameter_b,
        assembly.sign_a,assembly.sign_b,assembly.cs,assembly.cc,assembly.ss,ns,nq)
    return s.slater
end

# Unit-parameter coefficients from the same C four-block expressions. Coincident
# a/b parameter indices are combined before evaluation, exactly as CPU unit
# assembly does; no sparse duplicate-reduction arithmetic is introduced.
@inline function slater_unit_coefficient(a,b,cs,cc,ss,rowspin,columnspin)
    if rowspin == 0
        return columnspin == 0 ? -(a-b)*cs : a*cc+b*ss
    end
    return columnspin == 0 ? -a*ss-b*cc : (a-b)*cs
end

function slater_sparse_pass!(positions,rows,values,pa,pb,sa,sb,cs,cc,ss,ns,nq,::Val{write}) where {write}
    row = 0
    # This is C's flattened row order. Each column therefore receives sorted
    # entries directly, without tuples, duplicate summation or sorting.
    @inbounds for q in 1:nq, rowspin in 0:1, ri in 1:ns, columnspin in 0:1, rj in 1:ns
        row += 1
        ia,ib = Int(pa[ri,rj,q]),Int(pb[ri,rj,q])
        a,b = Float64(sa[ri,rj,q]),Float64(sb[ri,rj,q])
        x = slater_unit_coefficient(a,ia == ib ? b : 0.0,cs[q],cc[q],ss[q],rowspin,columnspin)
        if x != 0
            if write
                slot = positions[ia]
                rows[slot] = Int32(row); values[slot] = x
            end
            positions[ia] += 1
        end
        if ia != ib
            y = slater_unit_coefficient(0.0,b,cs[q],cc[q],ss[q],rowspin,columnspin)
            if y != 0
                if write
                    slot = positions[ib]
                    rows[slot] = Int32(row); values[slot] = y
                end
                positions[ib] += 1
            end
        end
    end
    return nothing
end

function sparse_slater_basis(assembly::SlaterAssembly,np::Integer)
    ns,nq = assembly.nsite,assembly.nqp
    pa,pb = Array(assembly.parameter_a),Array(assembly.parameter_b)
    sa,sb = Array(assembly.sign_a),Array(assembly.sign_b)
    cs,cc,ss = Array(assembly.cs),Array(assembly.cc),Array(assembly.ss)
    counts = zeros(Int32,np)
    slater_sparse_pass!(counts,nothing,nothing,pa,pb,sa,sb,cs,cc,ss,ns,nq,Val(false))
    pointers = Vector{Int32}(undef,np+1)
    pointers[1] = 1
    for p in 1:np
        pointers[p+1] = pointers[p]+counts[p]
    end
    rows = Vector{Int32}(undef,Int(pointers[end])-1)
    values = Vector{Float64}(undef,length(rows))
    positions = pointers[1:np]
    slater_sparse_pass!(positions,rows,values,pa,pb,sa,sb,cs,cc,ss,ns,nq,Val(true))
    return CUDA.cuSPARSE.CuSparseMatrixCSC(CuArray(pointers),CuArray(rows),CuArray(values),(4ns*ns*nq,Int(np)))
end
