"""
Quantum Projection Weight Utilities

Functions for initializing quantum projection weights.
Equivalent to C's InitQPWeight() and UpdateQPWeight().
Based on mVMC/src/mVMC/qp.c
"""

# Constants
const GAUSSLEG_EPS = 5.0e-14

"""
    gauss_legendre(x1::Float64, x2::Float64, n::Int) -> (Vector{Float64}, Vector{Float64})

Calculate n points and weights for Gauss-Legendre quadrature integration from x1 to x2.
Equivalent to C's GaussLeg() function.

# Arguments
- `x1::Float64`: Lower bound of integration
- `x2::Float64`: Upper bound of integration
- `n::Int`: Number of quadrature points

# Returns
- `(x::Vector{Float64}, w::Vector{Float64})`: Quadrature points and weights

# Example
```julia
x, w = gauss_legendre(0.0, π, 4)
```
"""
function gauss_legendre(
    x1::Float64,
    x2::Float64,
    n::Int,
)::Tuple{Vector{Float64},Vector{Float64}}
    if n <= 0
        return Float64[], Float64[]
    end

    x = zeros(Float64, n)
    w = zeros(Float64, n)

    m = (n + 1) ÷ 2
    xm = 0.5 * (x2 + x1)
    xl = 0.5 * (x2 - x1)

    for i = 0:(m-1)
        z = cos(π * (i + 0.75) / (n + 0.5))

        # Newton-Raphson iteration
        z1 = z
        pp = 0.0
        while true
            p1 = 1.0
            p2 = 0.0

            for j = 1:n
                p3 = p2
                p2 = p1
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j
            end

            pp = n * (z * p1 - p2) / (z * z - 1.0)
            z1 = z
            z = z1 - p1 / pp

            if abs(z - z1) <= GAUSSLEG_EPS
                break
            end
        end

        x[i+1] = xm - xl * z
        x[n-i] = xm + xl * z
        w[i+1] = 2.0 * xl / ((1.0 - z * z) * pp * pp)
        w[n-i] = w[i+1]
    end

    return x, w
end

"""
    legendre_poly(x::Float64, n::Int) -> Float64

Calculate Legendre polynomial P_n(x) using recurrence relation.
Equivalent to C's LegendrePoly() function.

# Arguments
- `x::Float64`: Point at which to evaluate the polynomial
- `n::Int`: Order of the Legendre polynomial

# Returns
- `Float64`: Value of P_n(x)

# Example
```julia
p3 = legendre_poly(0.5, 3)  # P_3(0.5)
```
"""
function legendre_poly(x::Float64, n::Int)::Float64
    if n <= 0
        return 1.0
    elseif n == 1
        return x
    else
        P01 = 1.0
        P02 = x
        P03 = x  # Initialize P03

        for i = 2:n
            P03 = (1.0 / i) * ((2.0 * i - 1.0) * x * P02 - (i - 1.0) * P01)
            P01 = P02
            P02 = P03
        end

        return P03
    end
end

"""
    QuantumProjectionWeights

Structure to store quantum projection weights.
Equivalent to C's QPFullWeight, QPFixWeight, SPGLCos, SPGLSin arrays.
"""
mutable struct QuantumProjectionWeights
    # Full weights: QPFullWeight[NQPFull]
    qp_full_weight::Vector{ComplexF64}

    # Fixed weights: QPFixWeight[NQPFix]
    qp_fix_weight::Vector{ComplexF64}

    # Spin projection trigonometric values: [NSPGaussLeg]
    spgl_cos::Vector{ComplexF64}      # cos(beta/2)
    spgl_sin::Vector{ComplexF64}      # sin(beta/2)
    spgl_cos_sin::Vector{ComplexF64}   # cos(beta/2) * sin(beta/2)
    spgl_cos_cos::Vector{ComplexF64}   # cos(beta/2)^2
    spgl_sin_sin::Vector{ComplexF64}   # sin(beta/2)^2

    function QuantumProjectionWeights()
        new(
            ComplexF64[],
            ComplexF64[],
            ComplexF64[],
            ComplexF64[],
            ComplexF64[],
            ComplexF64[],
            ComplexF64[],
        )
    end
end

"""
    init_qp_weight!(weights::QuantumProjectionWeights, modpara::ModParaParameters, para_qp_trans::Vector{ComplexF64}, opt_trans::Vector{ComplexF64}=ComplexF64[])

Initialize quantum projection weights.
Equivalent to C's InitQPWeight() function.

# Arguments
- `weights::QuantumProjectionWeights`: Structure to store weights
- `modpara::ModParaParameters`: Main simulation parameters (contains NSPGaussLeg, NSPStot, NMPTrans)
- `para_qp_trans::Vector{ComplexF64}`: ParaQPTrans values from qptransidx.def
- `opt_trans::Vector{ComplexF64}`: OptTrans values (optional, for FlagOptTrans > 0)

# Notes
- NQPFix = NSPGaussLeg * NMPTrans
- NQPFull = NQPFix * (NOptTrans > 0 ? NOptTrans : 1)
"""
function init_qp_weight!(
    weights::QuantumProjectionWeights,
    modpara::ModParaParameters,
    para_qp_trans::Vector{ComplexF64},
    opt_trans::Vector{ComplexF64} = ComplexF64[],
)
    nsp_gauss_leg = modpara.nsp_gauss_leg
    nsp_stot = modpara.nsp_stot
    # C実装では、NMPTrans < 0の場合、APFlag = 1を設定し、NMPTrans *= -1で正の値に変換する
    # (readdef.c:742-747行目を参照)
    # したがって、ここでは絶対値を使用する
    nmp_trans = abs(modpara.nmp_trans)

    if nsp_gauss_leg <= 0 || nmp_trans <= 0
        @warn "Invalid quantum projection parameters: NSPGaussLeg=$nsp_gauss_leg, NMPTrans=$(modpara.nmp_trans) (abs=$nmp_trans)"
        return weights
    end

    nqp_fix = nsp_gauss_leg * nmp_trans
    nqp_opt_trans = length(opt_trans) > 0 ? length(opt_trans) : 1
    nqp_full = nqp_fix * nqp_opt_trans

    # Allocate arrays
    weights.qp_full_weight = zeros(ComplexF64, nqp_full)
    weights.qp_fix_weight = zeros(ComplexF64, nqp_fix)
    weights.spgl_cos = zeros(ComplexF64, nsp_gauss_leg)
    weights.spgl_sin = zeros(ComplexF64, nsp_gauss_leg)
    weights.spgl_cos_sin = zeros(ComplexF64, nsp_gauss_leg)
    weights.spgl_cos_cos = zeros(ComplexF64, nsp_gauss_leg)
    weights.spgl_sin_sin = zeros(ComplexF64, nsp_gauss_leg)

    if nsp_gauss_leg == 1
        # Special case: NSPGaussLeg == 1
        weights.spgl_cos[1] = ComplexF64(1.0, 0.0)
        weights.spgl_sin[1] = ComplexF64(0.0, 0.0)
        weights.spgl_cos_sin[1] = ComplexF64(0.0, 0.0)
        weights.spgl_cos_cos[1] = ComplexF64(1.0, 0.0)
        weights.spgl_sin_sin[1] = ComplexF64(0.0, 0.0)

        # QPFixWeight = ParaQPTrans
        for j = 1:nmp_trans
            if j <= length(para_qp_trans)
                weights.qp_fix_weight[j] = para_qp_trans[j]
            end
        end
    else
        # Calculate Gauss-Legendre quadrature points and weights
        beta, weight_gl = gauss_legendre(0.0, Float64(π), nsp_gauss_leg)

        # Calculate spin projection trigonometric values and fixed weights
        for i = 1:nsp_gauss_leg
            beta_i = beta[i]
            weights.spgl_cos[i] = ComplexF64(cos(0.5 * beta_i), 0.0)
            weights.spgl_sin[i] = ComplexF64(sin(0.5 * beta_i), 0.0)
            weights.spgl_cos_sin[i] = weights.spgl_cos[i] * weights.spgl_sin[i]
            weights.spgl_cos_cos[i] = weights.spgl_cos[i] * weights.spgl_cos[i]
            weights.spgl_sin_sin[i] = weights.spgl_sin[i] * weights.spgl_sin[i]

            # Calculate weight: w = 0.5*sin(beta[i])*weight[i]*LegendrePoly(cos(beta[i]), NSPStot)
            cos_beta = cos(beta_i)
            w = 0.5 * sin(beta_i) * weight_gl[i] * legendre_poly(cos_beta, nsp_stot)

            # QPFixWeight[idx] = w * ParaQPTrans[j]
            # idx = i + j*NSPGaussLeg (0-based in C, 1-based in Julia)
            for j = 1:nmp_trans
                idx = i + (j - 1) * nsp_gauss_leg  # Convert to 1-based indexing
                if j <= length(para_qp_trans)
                    weights.qp_fix_weight[idx] = w * para_qp_trans[j]
                end
            end
        end
    end

    # Update full weights
    update_qp_weight!(weights, opt_trans)

    return weights
end

"""
    update_qp_weight!(weights::QuantumProjectionWeights, opt_trans::Vector{ComplexF64}=ComplexF64[])

Update quantum projection full weights from fixed weights.
Equivalent to C's UpdateQPWeight() function.

# Arguments
- `weights::QuantumProjectionWeights`: Structure containing weights
- `opt_trans::Vector{ComplexF64}`: OptTrans values (optional)

# Notes
- If opt_trans is provided and non-empty, QPFullWeight = OptTrans[i] * QPFixWeight[j]
- Otherwise, QPFullWeight = QPFixWeight
"""
function update_qp_weight!(
    weights::QuantumProjectionWeights,
    opt_trans::Vector{ComplexF64} = ComplexF64[],
)
    nqp_fix = length(weights.qp_fix_weight)

    if length(opt_trans) > 0
        # FlagOptTrans > 0: QPFullWeight[offset+j] = OptTrans[i] * QPFixWeight[j]
        nqp_opt_trans = length(opt_trans)
        nqp_full = nqp_fix * nqp_opt_trans

        # Resize if needed
        if length(weights.qp_full_weight) != nqp_full
            weights.qp_full_weight = zeros(ComplexF64, nqp_full)
        end

        for i = 1:nqp_opt_trans
            offset = (i - 1) * nqp_fix
            tmp = opt_trans[i]
            for j = 1:nqp_fix
                weights.qp_full_weight[offset+j] = tmp * weights.qp_fix_weight[j]
            end
        end
    else
        # FlagOptTrans == 0: QPFullWeight = QPFixWeight
        if length(weights.qp_full_weight) != nqp_fix
            weights.qp_full_weight = zeros(ComplexF64, nqp_fix)
        end

        for j = 1:nqp_fix
            weights.qp_full_weight[j] = weights.qp_fix_weight[j]
        end
    end

    return weights
end

"""
    init_qp_weight!(data::ExpertModeData)

Initialize quantum projection weights for ExpertModeData.
Convenience wrapper that uses data.modpara and data.para_qp_trans.

# Arguments
- `data::ExpertModeData`: Expert Mode data structure

# Notes
- Creates QuantumProjectionWeights if not already created
- Uses data.modpara.nsp_gauss_leg, data.modpara.nsp_stot, data.modpara.nmp_trans
- Uses data.para_qp_trans for ParaQPTrans values
- Uses data.opt_trans for OptTrans values when FlagOptTrans-equivalent inputs
  are active
"""
function init_qp_weight!(data::ExpertModeData)
    if data.qp_weights === nothing
        data.qp_weights = QuantumProjectionWeights()
    end

    init_qp_weight!(data.qp_weights, data.modpara, data.para_qp_trans, data.opt_trans)

    return data
end
