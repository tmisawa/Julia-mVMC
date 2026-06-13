"""
Expert Mode Data Types

Core data structures for mVMC Expert Mode file parsing and generation.
These types represent the various components of mVMC Expert Mode files.
"""

"""
    ParseResult{T}

Generic result container for parsing operations.
"""
struct ParseResult{T}
    success::Bool
    data::Union{T,Nothing}
    error_message::String
    line_number::Int
end

"""
    ValidationResult

Result of parameter validation.
"""
struct ValidationResult
    is_valid::Bool
    errors::Vector{String}
    warnings::Vector{String}
end

"""
    FileInfo

Information about a file.
"""
struct FileInfo
    filename::String
    filepath::String
    exists::Bool
    size_bytes::Int
    last_modified::Dates.DateTime
end

"""
    ParsingContext

Context for parsing operations, tracking errors, warnings, and line numbers.
"""
mutable struct ParsingContext
    filename::String
    line_number::Int
    errors::Vector{String}
    warnings::Vector{String}

    function ParsingContext(filename::String)
        new(filename, 0, String[], String[])
    end
end

"""
    ModParaParameters

Parameters for modpara.def file.
Contains all the main simulation parameters.
"""
mutable struct ModParaParameters
    # Basic system parameters
    nsite::Int
    nelec::Int
    nlocspin::Int
    ncond::Int

    # Calculation modes
    vmc_calc_mode::Int  # 0: optimization, 1: physics calculation
    lanczos_mode::Int   # 0: none, 1: energy only, 2: Green functions

    # VMC parameters
    nsr_opt_itr_step::Int
    nsr_opt_itr_smp::Int
    nsr_opt_fix_smp::Int
    nvmc_warmup::Int
    nvmc_interval::Int
    nvmc_sample::Int

    # SR parameters
    dsr_opt_red_cut::Float64
    dsr_opt_sta_del::Float64
    dsr_opt_step_dt::Float64
    dsr_opt_cg_tol::Float64
    nsr_opt_cg_max_iter::Int

    # SR solver selection
    nsrcg::Int          # 0: direct solver (LAPACK), !=0: CG solver
    nstore_o::Int       # 0: normal, !=0: store O samples

    # Random number generation
    rnd_seed::Int
    nsplit_size::Int

    # Quantum projection
    nsp_gauss_leg::Int
    nsp_stot::Int
    nmp_trans::Int
    two_sz::Int

    # Data output
    n_data_idx_start::Int
    n_data_qty_smp::Int
    c_data_file_head::String
    c_para_file_head::String

    # File control
    n_file_flush_interval::Int

    # Complex flag
    complex_flag::Int

    # RBM parameters
    nneuron::Int
    nneuron_general::Int
    nneuron_charge::Int
    nneuron_spin::Int
    nblock_size_rbm_ratio::Int

    # Lanczos parameters
    n_one_body_g::Int
    n_two_body_g::Int
    n_two_body_g_ex::Int

    # Exchange update
    nex_update_path::Int

    # Orbital parameters
    n_orbital_idx::Int  # Number of unique orbital parameters (from NOrbitalIdx in orbitalidx.def)

    function ModParaParameters(;
        nsite::Int = 0,
        nelec::Int = 0,
        nlocspin::Int = 0,
        ncond::Int = -1,
        vmc_calc_mode::Int = 0,
        lanczos_mode::Int = 0,
        nsr_opt_itr_step::Int = 1000,
        nsr_opt_itr_smp::Int = 1000,
        nsr_opt_fix_smp::Int = 0,
        nvmc_warmup::Int = 1000,
        nvmc_interval::Int = 1,
        nvmc_sample::Int = 10000,
        dsr_opt_red_cut::Float64 = 1e-6,
        dsr_opt_sta_del::Float64 = 0.0,
        dsr_opt_step_dt::Float64 = 0.01,
        dsr_opt_cg_tol::Float64 = 1e-10,
        nsr_opt_cg_max_iter::Int = 0,
        nsrcg::Int = 0,
        nstore_o::Int = 1,
        rnd_seed::Int = 11272,  # C parity: RndSeed 行欠落時の default (readdef.c:1967)
        nsplit_size::Int = 1,
        nsp_gauss_leg::Int = 1,
        nsp_stot::Int = 0,
        nmp_trans::Int = 0,
        two_sz::Int = -1,  # -1 means Sz not conserved (FSZ mode), matching C default
        n_data_idx_start::Int = 0,
        n_data_qty_smp::Int = 1,
        c_data_file_head::String = "zvo",
        c_para_file_head::String = "zqp",
        n_file_flush_interval::Int = 1,
        complex_flag::Int = 0,
        nneuron::Int = 0,
        nneuron_general::Int = 0,
        nneuron_charge::Int = 0,
        nneuron_spin::Int = 0,
        # C default (readdef.c): NBlockSize_RBMRatio = 200
        nblock_size_rbm_ratio::Int = 200,
        n_one_body_g::Int = 0,
        n_two_body_g::Int = 0,
        n_two_body_g_ex::Int = 0,
        nex_update_path::Int = 1,
        n_orbital_idx::Int = 0,
    )
        new(
            nsite,
            nelec,
            nlocspin,
            ncond,
            vmc_calc_mode,
            lanczos_mode,
            nsr_opt_itr_step,
            nsr_opt_itr_smp,
            nsr_opt_fix_smp,
            nvmc_warmup,
            nvmc_interval,
            nvmc_sample,
            dsr_opt_red_cut,
            dsr_opt_sta_del,
            dsr_opt_step_dt,
            dsr_opt_cg_tol,
            nsr_opt_cg_max_iter,
            nsrcg,
            nstore_o,
            rnd_seed,
            nsplit_size,
            nsp_gauss_leg,
            nsp_stot,
            nmp_trans,
            two_sz,
            n_data_idx_start,
            n_data_qty_smp,
            c_data_file_head,
            c_para_file_head,
            n_file_flush_interval,
            complex_flag,
            nneuron,
            nneuron_general,
            nneuron_charge,
            nneuron_spin,
            nblock_size_rbm_ratio,
            n_one_body_g,
            n_two_body_g,
            n_two_body_g_ex,
            nex_update_path,
            n_orbital_idx,
        )
    end
end

"""
    TransferTerm

Single transfer (hopping) term.
"""
struct TransferTerm
    site1::Int
    spin1::Int
    site2::Int
    spin2::Int
    value::ComplexF64
    spin::Symbol  # :up, :down, :both
end

TransferTerm(site1::Int, spin1::Int, site2::Int, spin2::Int, value::ComplexF64) = begin
    spin = spin1 == 0 ? :up : :down
    TransferTerm(site1, spin1, site2, spin2, value, spin)
end

TransferTerm(site1::Int, site2::Int, value::ComplexF64, spin::Symbol) = begin
    spin1 = spin == :down ? 1 : 0
    spin2 = spin1
    TransferTerm(site1, spin1, site2, spin2, value, spin)
end

"""
    CoulombIntraTerm

On-site Coulomb interaction term.
"""
struct CoulombIntraTerm
    site::Int
    value::Float64
end

"""
    CoulombInterTerm

Inter-site Coulomb interaction term.
"""
struct CoulombInterTerm
    site1::Int
    site2::Int
    value::Float64
end

"""
    HundTerm

Hund coupling term.
"""
struct HundTerm
    site1::Int
    site2::Int
    value::Float64
end

"""
    ExchangeTerm

Exchange coupling term.
"""
struct ExchangeTerm
    site1::Int
    site2::Int
    value::Float64
end

"""
    PairHopTerm

Pair hopping term.
"""
struct PairHopTerm
    site1::Int
    site2::Int
    value::Float64
end

"""
    GutzwillerTerm

Gutzwiller projection term.
"""
mutable struct GutzwillerTerm
    site::Int
    value::ComplexF64
    is_complex::Bool
end

"""
    JastrowTerm

Jastrow factor term.
"""
mutable struct JastrowTerm
    site1::Int
    site2::Int
    value::ComplexF64
    is_complex::Bool
end

"""
    OrbitalTerm

Orbital parameter term.
"""
mutable struct OrbitalTerm
    site1::Int
    site2::Int
    idx::Int      # Orbital parameter index (0 to NOrbitalIdx-1)
    value::ComplexF64
    is_complex::Bool
    sign::Int  # OrbitalSgn: +1 or -1 (default: 1)
end

# Constructor with default sign = 1 for backward compatibility
OrbitalTerm(site1::Int, site2::Int, value::ComplexF64, is_complex::Bool) =
    OrbitalTerm(site1, site2, 0, value, is_complex, 1)

# Constructor with explicit sign (backward compatibility for older test/helpers)
OrbitalTerm(site1::Int, site2::Int, value::ComplexF64, is_complex::Bool, sign::Int) =
    OrbitalTerm(site1, site2, 0, value, is_complex, sign)

# Constructor with idx for orbitalidx.def format
OrbitalTerm(site1::Int, site2::Int, idx::Int, value::ComplexF64, is_complex::Bool) =
    OrbitalTerm(site1, site2, idx, value, is_complex, 1)

"""
    GreenOneTerm

One-body Green's function measurement term.
"""
struct GreenOneTerm
    site1::Int
    site2::Int
    spin1::Symbol
    spin2::Symbol
end

"""
    GreenTwoTerm

Two-body Green's function measurement term.
"""
struct GreenTwoTerm
    site1::Int
    site2::Int
    site3::Int
    site4::Int
    spin1::Symbol
    spin2::Symbol
    spin3::Symbol
    spin4::Symbol
end

"""
    GreenTwoExTerm

Factored (product-side) two-body Green's function term, `TwoBodyGEx` /
`greentwoex.def`. Each term names two one-body Green's functions whose product
`⟨c†_{i1} c_{j1}⟩ · conj(⟨c†_{i2} c_{j2}⟩)` is accumulated by PhysCal.

The 8 input columns `x0 x1 x2 x3 x4 x5 x6 x7` map with C's reorder absorbed
(see `GetInfoTwoBodyGEx` in mVMC `readdef.c`):
- first one-body Green  `⟨c†_{x0,x1} c_{x2,x3}⟩`  → `(site_i1,spin_i1,site_j1,spin_j1)`
- second one-body Green `⟨c†_{x6,x7} c_{x4,x5}⟩`  → `(site_i2,spin_i2,site_j2,spin_j2)`

Spins are kept as integers (0 = up, 1 = down) to match the integer lookup used
when resolving these to one-body Green indices (Plan 2).
"""
struct GreenTwoExTerm
    site_i1::Int
    spin_i1::Int
    site_j1::Int
    spin_j1::Int
    site_i2::Int
    spin_i2::Int
    site_j2::Int
    spin_j2::Int
end

"""
    QPTransTerm

Quantum projection translation term.
"""
struct QPTransTerm
    site::Int
    momentum::Vector{Float64}
    phase::Float64
end

"""
    InterAllTerm

General interaction term.

The format stores 4 site-spin pairs:
- site0, spin0: first creation operator (c†_{site0, spin0})
- site1, spin1: first annihilation operator (c_{site1, spin1})
- site2, spin2: second creation operator (c†_{site2, spin2})
- site3, spin3: second annihilation operator (c_{site3, spin3})

The term represents: value * <c†_{site0,spin0} c_{site1,spin1} c†_{site2,spin2} c_{site3,spin3}>
"""
struct InterAllTerm
    site0::Int
    spin0::Int
    site1::Int
    spin1::Int
    site2::Int
    spin2::Int
    site3::Int
    spin3::Int
    value::ComplexF64
    is_complex::Bool
end

"""
    LocSpinTerm

Local spin term.
"""
struct LocSpinTerm
    site::Int
    spin_value::Int
end

"""
    RBMTerm

Base RBM term structure.
"""
abstract type RBMTerm end

"""
    ChargeRBMPhysLayerTerm

Charge RBM physical layer term.
"""
mutable struct ChargeRBMPhysLayerTerm <: RBMTerm
    site::Int
    value::ComplexF64
    is_complex::Bool
    idx::Int
end
ChargeRBMPhysLayerTerm(site::Int, value::ComplexF64, is_complex::Bool) =
    ChargeRBMPhysLayerTerm(site, value, is_complex, site)

"""
    SpinRBMPhysLayerTerm

Spin RBM physical layer term.
"""
mutable struct SpinRBMPhysLayerTerm <: RBMTerm
    site::Int
    value::ComplexF64
    is_complex::Bool
    idx::Int
end
SpinRBMPhysLayerTerm(site::Int, value::ComplexF64, is_complex::Bool) =
    SpinRBMPhysLayerTerm(site, value, is_complex, site)

"""
    GeneralRBMPhysLayerTerm

General RBM physical layer term.
"""
mutable struct GeneralRBMPhysLayerTerm <: RBMTerm
    site::Int
    spin::Int
    value::ComplexF64
    is_complex::Bool
    idx::Int
end
GeneralRBMPhysLayerTerm(site::Int, value::ComplexF64, is_complex::Bool) =
    GeneralRBMPhysLayerTerm(site, 0, value, is_complex, site)

"""
    ChargeRBMHiddenLayerTerm

Charge RBM hidden layer term.
"""
mutable struct ChargeRBMHiddenLayerTerm <: RBMTerm
    site::Int
    value::ComplexF64
    is_complex::Bool
    idx::Int
end
ChargeRBMHiddenLayerTerm(site::Int, value::ComplexF64, is_complex::Bool) =
    ChargeRBMHiddenLayerTerm(site, value, is_complex, site)

"""
    SpinRBMHiddenLayerTerm

Spin RBM hidden layer term.
"""
mutable struct SpinRBMHiddenLayerTerm <: RBMTerm
    site::Int
    value::ComplexF64
    is_complex::Bool
    idx::Int
end
SpinRBMHiddenLayerTerm(site::Int, value::ComplexF64, is_complex::Bool) =
    SpinRBMHiddenLayerTerm(site, value, is_complex, site)

"""
    GeneralRBMHiddenLayerTerm

General RBM hidden layer term.
"""
mutable struct GeneralRBMHiddenLayerTerm <: RBMTerm
    site::Int
    value::ComplexF64
    is_complex::Bool
    idx::Int
end
GeneralRBMHiddenLayerTerm(site::Int, value::ComplexF64, is_complex::Bool) =
    GeneralRBMHiddenLayerTerm(site, value, is_complex, site)

"""
    ChargeRBMPhysHiddenTerm

Charge RBM physical-hidden connection term.
"""
mutable struct ChargeRBMPhysHiddenTerm <: RBMTerm
    site1::Int
    site2::Int
    value::ComplexF64
    is_complex::Bool
    idx::Int
end
ChargeRBMPhysHiddenTerm(site1::Int, site2::Int, value::ComplexF64, is_complex::Bool) =
    ChargeRBMPhysHiddenTerm(site1, site2, value, is_complex, 0)

"""
    SpinRBMPhysHiddenTerm

Spin RBM physical-hidden connection term.
"""
mutable struct SpinRBMPhysHiddenTerm <: RBMTerm
    site1::Int
    site2::Int
    value::ComplexF64
    is_complex::Bool
    idx::Int
end
SpinRBMPhysHiddenTerm(site1::Int, site2::Int, value::ComplexF64, is_complex::Bool) =
    SpinRBMPhysHiddenTerm(site1, site2, value, is_complex, 0)

"""
    GeneralRBMPhysHiddenTerm

General RBM physical-hidden connection term.
"""
mutable struct GeneralRBMPhysHiddenTerm <: RBMTerm
    site1::Int
    spin::Int
    site2::Int
    value::ComplexF64
    is_complex::Bool
    idx::Int
end
GeneralRBMPhysHiddenTerm(site1::Int, site2::Int, value::ComplexF64, is_complex::Bool) =
    GeneralRBMPhysHiddenTerm(site1, 0, site2, value, is_complex, 0)

"""
    DoublonHolon2SiteTerm

Deprecated compatibility shim for the pre-DH-1 value-bearing DH2 term API.
Runtime data uses `DoublonHolon2SiteIndex` plus projection parameters instead.
"""
struct DoublonHolon2SiteTerm
    site1::Int
    site2::Int
    value::ComplexF64
    is_complex::Bool
end

"""
    DoublonHolon4SiteTerm

Deprecated compatibility shim for the pre-DH-1 value-bearing DH4 term API.
Runtime data uses `DoublonHolon4SiteIndex` plus projection parameters instead.
"""
struct DoublonHolon4SiteTerm
    site1::Int
    site2::Int
    site3::Int
    site4::Int
    value::ComplexF64
    is_complex::Bool
end

"""
    DoublonHolon2SiteIndex

C-compatible DH2 neighbor table for one `NDoublonHolon2siteIdx` index.
`neighbors[i+1, :]` stores the two zero-based neighbor site ids for center
site `i`, matching C's `DoublonHolon2siteIdx[idx][2*i : 2*i+1]`.
"""
struct DoublonHolon2SiteIndex
    neighbors::Matrix{Int}

    function DoublonHolon2SiteIndex(neighbors::Matrix{Int})
        size(neighbors, 2) == 2 ||
            throw(ArgumentError("DoublonHolon2SiteIndex.neighbors must be Nsite x 2"))
        new(neighbors)
    end
end

"""
    DoublonHolon4SiteIndex

C-compatible DH4 neighbor table for one `NDoublonHolon4siteIdx` index.
`neighbors[i+1, :]` stores the four zero-based neighbor site ids for center
site `i`, matching C's `DoublonHolon4siteIdx[idx][4*i : 4*i+3]`.
"""
struct DoublonHolon4SiteIndex
    neighbors::Matrix{Int}

    function DoublonHolon4SiteIndex(neighbors::Matrix{Int})
        size(neighbors, 2) == 4 ||
            throw(ArgumentError("DoublonHolon4SiteIndex.neighbors must be Nsite x 4"))
        new(neighbors)
    end
end

"""
    DoublonHolon2SiteDefinition / DoublonHolon4SiteDefinition

Parsed C-compatible DH definition files. `opt_flags` stores the real-part
optimization flags in C file order. Imaginary flags are derived from
`is_complex` when folded into `ExpertModeData.optimization_flags`.
"""
struct DoublonHolon2SiteDefinition
    indices::Vector{DoublonHolon2SiteIndex}
    opt_flags::Vector{Bool}
    is_complex::Bool
end

struct DoublonHolon4SiteDefinition
    indices::Vector{DoublonHolon4SiteIndex}
    opt_flags::Vector{Bool}
    is_complex::Bool
end

"""
    ProjectionLayout

C-compatible projection-factor layout:
Gutzwiller | Jastrow | SpinJastrow | DH2 | DH4.
Offsets are zero-based C offsets into the `Proj` / `projCnt` slices.
"""
struct ProjectionLayout
    n_gutzwiller::Int
    n_jastrow::Int
    n_spinjastrow::Int
    n_dh2::Int
    n_dh4::Int
    gutzwiller_offset::Int
    jastrow_offset::Int
    spinjastrow_offset::Int
    dh2_offset::Int
    dh4_offset::Int
    n_proj::Int
end

"""
    ExpertModeData

Container for all Expert Mode data.
"""
mutable struct ExpertModeData
    # Main parameters
    modpara::ModParaParameters

    # Transfer terms
    transfer_terms::Vector{TransferTerm}

    # Interaction terms
    coulomb_intra_terms::Vector{CoulombIntraTerm}
    coulomb_inter_terms::Vector{CoulombInterTerm}
    hund_terms::Vector{HundTerm}
    exchange_terms::Vector{ExchangeTerm}
    pair_hop_terms::Vector{PairHopTerm}
    inter_all_terms::Vector{InterAllTerm}

    # Variational parameters
    gutzwiller_terms::Vector{GutzwillerTerm}
    jastrow_terms::Vector{JastrowTerm}
    orbital_terms::Vector{OrbitalTerm}

    # Index arrays for projection factors (C implementation equivalents)
    gutzwiller_idx::Vector{Int}  # GutzwillerIdx[ri] = idx for site ri (0-based indexing)
    jastrow_idx::Matrix{Int}  # JastrowIdx[ri, rj] = idx for site pair (ri, rj) (0-based indexing)
    n_gutzwiller_idx::Int  # NGutzwillerIdx
    n_jastrow_idx::Int  # NJastrowIdx

    # Green's function measurements
    green_one_terms::Vector{GreenOneTerm}
    green_two_terms::Vector{GreenTwoTerm}
    green_two_ex_terms::Vector{GreenTwoExTerm}  # TwoBodyGEx / greentwoex.def (factored)

    # Quantum projection
    qptrans_terms::Vector{QPTransTerm}
    para_qp_trans::Vector{ComplexF64}  # ParaQPTrans values from qptransidx.def
    para_qp_opt_trans::Vector{ComplexF64}  # ParaQPOptTrans values from opttrans.def
    opt_trans::Vector{ComplexF64}  # Runtime OptTrans values (empty when FlagOptTrans-equivalent is off)
    n_qp_trans::Int  # NQPTrans from qptransidx.def
    n_qp_opt_trans::Int  # NQPOptTrans (default: 1)
    qp_weights::Union{Nothing,Any}  # QuantumProjectionWeights (initialized by init_qp_weight!)

    # QPTrans mappings (from qptransidx.def)
    qp_trans::Vector{Vector{Int}}  # QPTrans[mpidx][ori] = trj (translated site)
    qp_trans_inv::Vector{Vector{Int}}  # QPTransInv[mpidx][trj] = ori (inverse mapping)
    qp_trans_sgn::Vector{Vector{Int}}  # QPTransSgn[mpidx][ori] = sign (+1 or -1)
    qp_opt_trans::Vector{Vector{Int}}  # QPOptTrans[optidx][ri] = ori
    qp_opt_trans_sgn::Vector{Vector{Int}}  # QPOptTransSgn[optidx][ri] = sign (+1 or -1)

    # OrbitalIdx and OrbitalSgn matrices (from orbitalidx.def)
    # OrbitalIdx[i+1, j+1] = orbital parameter index for site pair (i, j)
    orbital_idx_matrix::Union{Nothing,Matrix{Int}}  # Like C's OrbitalIdx[ri][rj]
    orbital_sgn::Union{Nothing,Matrix{Int}}  # OrbitalSgn[i+1, j+1] = sign (+1 or -1) for site pair (i, j)

    # Local spins
    locspin_terms::Vector{LocSpinTerm}

    # RBM terms
    charge_rbm_phys_layer_terms::Vector{ChargeRBMPhysLayerTerm}
    spin_rbm_phys_layer_terms::Vector{SpinRBMPhysLayerTerm}
    general_rbm_phys_layer_terms::Vector{GeneralRBMPhysLayerTerm}
    charge_rbm_hidden_layer_terms::Vector{ChargeRBMHiddenLayerTerm}
    spin_rbm_hidden_layer_terms::Vector{SpinRBMHiddenLayerTerm}
    general_rbm_hidden_layer_terms::Vector{GeneralRBMHiddenLayerTerm}
    charge_rbm_phys_hidden_terms::Vector{ChargeRBMPhysHiddenTerm}
    spin_rbm_phys_hidden_terms::Vector{SpinRBMPhysHiddenTerm}
    general_rbm_phys_hidden_terms::Vector{GeneralRBMPhysHiddenTerm}

    # Doublon-Holon index tables and projection parameters
    doublon_holon_2site_indices::Vector{DoublonHolon2SiteIndex}
    doublon_holon_4site_indices::Vector{DoublonHolon4SiteIndex}
    doublon_holon_2site_params::Vector{ComplexF64}
    doublon_holon_4site_params::Vector{ComplexF64}
    doublon_holon_2site_opt_flags::Vector{Bool}
    doublon_holon_4site_opt_flags::Vector{Bool}
    doublon_holon_2site_complex::Bool
    doublon_holon_4site_complex::Bool

    # Optimization flags (corresponds to C implementation's OptFlag)
    optimization_flags::Vector{Bool}
    complex_flags::Vector{Int}  # Complex number flags for each term type

    # Orbital mode flags (corresponds to C implementation's iFlgOrbitalGeneral)
    i_flg_orbital_general::Int  # 0: sz conserved, 1: general (fsz)
    i_flg_orbital_anti_parallel::Int  # 0: not used, 1: OrbitalAntiParallel exists
    i_flg_orbital_parallel::Int  # 0: not used, 1: OrbitalParallel exists
    # Number of anti-parallel orbital parameters (C's iNOrbitalAntiParallel / NArrayAP).
    # Recorded exactly at parse time and used as the offset where the parallel
    # orbital block begins, mirroring readdef.c (Slater[iNOrbitalAntiParallel + idx]).
    n_orbital_anti_parallel::Int

    function ExpertModeData()
        new(
            ModParaParameters(),
            TransferTerm[],
            CoulombIntraTerm[],
            CoulombInterTerm[],
            HundTerm[],
            ExchangeTerm[],
            PairHopTerm[],
            InterAllTerm[],
            GutzwillerTerm[],
            JastrowTerm[],
            OrbitalTerm[],
            Int[],
            zeros(Int, 0, 0),
            0,
            0,  # gutzwiller_idx, jastrow_idx, n_gutzwiller_idx, n_jastrow_idx
            GreenOneTerm[],
            GreenTwoTerm[],
            GreenTwoExTerm[],
            QPTransTerm[],
            ComplexF64[],
            ComplexF64[],
            ComplexF64[],
            0,
            1,
            nothing,  # para_qp_trans, para_qp_opt_trans, opt_trans, n_qp_trans, n_qp_opt_trans, qp_weights
            Vector{Vector{Int}}[],  # qp_trans
            Vector{Vector{Int}}[],  # qp_trans_inv
            Vector{Vector{Int}}[],  # qp_trans_sgn
            Vector{Vector{Int}}[],  # qp_opt_trans
            Vector{Vector{Int}}[],  # qp_opt_trans_sgn
            nothing,  # orbital_idx_matrix
            nothing,  # orbital_sgn
            LocSpinTerm[],
            ChargeRBMPhysLayerTerm[],
            SpinRBMPhysLayerTerm[],
            GeneralRBMPhysLayerTerm[],
            ChargeRBMHiddenLayerTerm[],
            SpinRBMHiddenLayerTerm[],
            GeneralRBMHiddenLayerTerm[],
            ChargeRBMPhysHiddenTerm[],
            SpinRBMPhysHiddenTerm[],
            GeneralRBMPhysHiddenTerm[],
            # Doublon-Holon index tables and projection parameters
            DoublonHolon2SiteIndex[],
            DoublonHolon4SiteIndex[],
            ComplexF64[],
            ComplexF64[],
            Bool[],
            Bool[],
            false,
            false,
            Bool[],
            Int[],  # Initialize empty optimization and complex flags
            0,
            0,
            0,  # i_flg_orbital_general, i_flg_orbital_anti_parallel, i_flg_orbital_parallel (default: 0)
            0,  # n_orbital_anti_parallel (NArrayAP, default: 0)
        )
    end
end

@inline function _projection_count_from_header_or_terms(header_count::Int, term_count::Int)
    return header_count > 0 ? header_count : term_count
end

"""
    projection_layout(data::ExpertModeData) -> ProjectionLayout

Return C-compatible projection offsets for the currently parsed data. SpinJastrow
is deliberately zero until the parser implements it; `parse_expert_mode_files`
hard-fails present SpinJastrow inputs rather than silently accepting them.
"""
function projection_layout(data::ExpertModeData)::ProjectionLayout
    n_gutz = _projection_count_from_header_or_terms(
        data.n_gutzwiller_idx,
        length(data.gutzwiller_terms),
    )
    n_jast = _projection_count_from_header_or_terms(
        data.n_jastrow_idx,
        length(data.jastrow_terms),
    )
    n_spinjast = 0
    n_dh2 = length(data.doublon_holon_2site_indices)
    n_dh4 = length(data.doublon_holon_4site_indices)

    gutz_offset = 0
    jast_offset = gutz_offset + n_gutz
    spin_offset = jast_offset + n_jast
    dh2_offset = spin_offset + n_spinjast
    dh4_offset = dh2_offset + 6 * n_dh2
    n_proj = dh4_offset + 10 * n_dh4

    return ProjectionLayout(
        n_gutz,
        n_jast,
        n_spinjast,
        n_dh2,
        n_dh4,
        gutz_offset,
        jast_offset,
        spin_offset,
        dh2_offset,
        dh4_offset,
        n_proj,
    )
end

n_projection_parameters(data::ExpertModeData)::Int = projection_layout(data).n_proj

function projection_parameters(
    data::ExpertModeData,
    layout::ProjectionLayout = projection_layout(data),
)::Vector{ComplexF64}
    params = Vector{ComplexF64}(undef, layout.n_proj)
    fill!(params, 0.0 + 0.0im)

    for i = 1:min(layout.n_gutzwiller, length(data.gutzwiller_terms))
        params[layout.gutzwiller_offset+i] = data.gutzwiller_terms[i].value
    end
    for i = 1:min(layout.n_jastrow, length(data.jastrow_terms))
        params[layout.jastrow_offset+i] = data.jastrow_terms[i].value
    end
    for i = 1:min(6 * layout.n_dh2, length(data.doublon_holon_2site_params))
        params[layout.dh2_offset+i] = data.doublon_holon_2site_params[i]
    end
    for i = 1:min(10 * layout.n_dh4, length(data.doublon_holon_4site_params))
        params[layout.dh4_offset+i] = data.doublon_holon_4site_params[i]
    end

    return params
end
