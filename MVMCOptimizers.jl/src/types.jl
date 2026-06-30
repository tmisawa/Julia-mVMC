"""
Data types for VMC optimization
"""

"""
    EnergyData

Energy-related data for VMC calculation.
"""
mutable struct EnergyData
    wc::ComplexF64      # Weight for correlation sampling
    etot::ComplexF64    # <H>
    etot2::ComplexF64   # <H^2>
    sztot::ComplexF64   # <Sz>
    sztot2::ComplexF64  # <Sz^2>

    function EnergyData()
        new(0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im)
    end
end

"""
    SROptData

Stochastic Reconfiguration optimization data.
"""
mutable struct SROptData
    # Size
    sr_opt_size::Int  # 1 + NPara

    # Complex version
    sr_opt_oo::Vector{ComplexF64}      # <O^† O> matrix
    sr_opt_ho::Vector{ComplexF64}      # <H O> vector
    sr_opt_o::Vector{ComplexF64}       # Calculation buffer
    sr_opt_o_store::Vector{ComplexF64} # Sample storage buffer

    # Real version
    sr_opt_oo_real::Vector{Float64}
    sr_opt_ho_real::Vector{Float64}
    sr_opt_o_real::Vector{Float64}
    sr_opt_o_store_real::Vector{Float64}

    function SROptData(sr_opt_size::Int, n_vmc_sample::Int, all_complex::Bool)
        if all_complex
            new(
                sr_opt_size,
                zeros(ComplexF64, 2 * sr_opt_size * (2 * sr_opt_size + 2)),  # SROptOO: C: (2*SROptSize)*(2*SROptSize+2)
                zeros(ComplexF64, 2 * sr_opt_size),                      # SROptHO
                zeros(ComplexF64, 2 * sr_opt_size),                      # SROptO
                zeros(ComplexF64, 2 * sr_opt_size * n_vmc_sample),       # SROptO_Store
                Float64[],  # Real versions not used
                Float64[],
                Float64[],
                Float64[],
            )
        else
            # Real mode: also allocate complex arrays for intermediate calculations
            # Complex arrays are used as working buffers even in real mode
            new(
                sr_opt_size,
                zeros(ComplexF64, 2 * sr_opt_size * (2 * sr_opt_size + 2)),  # SROptOO: C: (2*SROptSize)*(2*SROptSize+2)
                zeros(ComplexF64, 2 * sr_opt_size),                      # SROptHO (complex working buffer)
                zeros(ComplexF64, 2 * sr_opt_size),                      # SROptO (complex working buffer)
                zeros(ComplexF64, 2 * sr_opt_size * n_vmc_sample),       # SROptO_Store (complex working buffer)
                zeros(Float64, sr_opt_size * (sr_opt_size + 2)),  # SROptOO_real: C: SROptSize*(SROptSize+2)
                zeros(Float64, sr_opt_size),                      # SROptHO_real
                zeros(Float64, sr_opt_size),                      # SROptO_real
                zeros(Float64, sr_opt_size * n_vmc_sample),        # SROptO_Store_real
            )
        end
    end
end

"""
    OptDataPoint

Single point of optimization data for averaging.
"""
struct OptDataPoint
    energy::ComplexF64
    parameters::Vector{ComplexF64}
end

"""
    ElectronConfiguration

Electron configuration for VMC sampling.
"""
mutable struct ElectronConfiguration
    ele_idx::Vector{Int}      # Electron indices [sample][mi+si*Ne]
    ele_cfg::Vector{Int}      # Electron configuration [sample][ri+si*Nsite]
    ele_num::Vector{Int}      # Electron number [sample][ri+si*Nsite]
    ele_proj_cnt::Vector{Int} # Projection count [sample][proj]
    ele_spn::Vector{Int}      # Electron spin [sample][mi+si*Ne] (for fsz)

    # Temporary arrays for sampling (single sample)
    tmp_ele_idx::Vector{Int}      # Temporary electron indices [mi+si*Ne]
    tmp_ele_cfg::Vector{Int}      # Temporary electron configuration [ri+si*Nsite]
    tmp_ele_num::Vector{Int}      # Temporary electron number [ri+si*Nsite]
    tmp_ele_proj_cnt::Vector{Int} # Temporary projection count [proj]
    tmp_ele_spn::Vector{Int}      # Temporary electron spin [mi+si*Ne] (for fsz)

    # Burn-in sample storage
    burn_ele_idx::Vector{Int}
    burn_ele_cfg::Vector{Int}
    burn_ele_num::Vector{Int}
    burn_ele_proj_cnt::Vector{Int}
    burn_ele_spn::Vector{Int}

    # Counters for statistics
    counter::Vector{Int}  # Various counters (hopping attempts, accepts, etc.)

    function ElectronConfiguration(
        n_sample::Int,
        n_site::Int,
        n_elec::Int,
        n_proj::Int,
        use_fsz::Bool,
    )
        n_size = 2 * n_elec
        n_site2 = 2 * n_site

        if use_fsz
            new(
                zeros(Int, n_sample * n_size),
                zeros(Int, n_sample * n_site2),
                zeros(Int, n_sample * n_site2),
                zeros(Int, n_sample * n_proj),
                zeros(Int, n_sample * n_size),
                # Temporary arrays
                zeros(Int, n_size),
                zeros(Int, n_site2),
                zeros(Int, n_site2),
                zeros(Int, n_proj),
                zeros(Int, n_size),
                # Burn-in arrays
                zeros(Int, n_size + n_site2 + n_site2 + n_proj + n_size),  # Combined storage
                zeros(Int, n_site2),
                zeros(Int, n_site2),
                zeros(Int, n_proj),
                zeros(Int, n_size),
                # Counters (max 10 counters)
                zeros(Int, 10),
            )
        else
            new(
                zeros(Int, n_sample * n_size),
                zeros(Int, n_sample * n_site2),
                zeros(Int, n_sample * n_site2),
                zeros(Int, n_sample * n_proj),
                Int[],  # Not used for sz-conserved case
                # Temporary arrays
                zeros(Int, n_size),
                zeros(Int, n_site2),
                zeros(Int, n_site2),
                zeros(Int, n_proj),
                Int[],  # Not used for sz-conserved case
                # Burn-in arrays
                zeros(Int, n_size + n_site2 + n_site2 + n_proj),  # Combined storage
                zeros(Int, n_site2),
                zeros(Int, n_site2),
                zeros(Int, n_proj),
                Int[],  # Not used for sz-conserved case
                # Counters
                zeros(Int, 10),
            )
        end
    end
end

"""
    SlaterMatrixData

Slater matrix and inverse matrix data.
"""
mutable struct SlaterMatrixData
    slater_elm::Vector{ComplexF64}      # Slater matrix elements [QPidx][ri+si*Nsite][rj+sj*Nsite]
    inv_m::Vector{ComplexF64}           # Inverse matrix [QPidx][mi+si*Ne][mj+sj*Ne]
    pf_m::Vector{ComplexF64}            # Pfaffian [QPidx]

    # Real versions (for real TBC)
    slater_elm_real::Vector{Float64}
    inv_m_real::Vector{Float64}
    pf_m_real::Vector{Float64}

    function SlaterMatrixData(n_qp_full::Int, n_site::Int, n_elec::Int, all_complex::Bool)
        # Validate inputs
        n_qp_full = max(1, n_qp_full)
        n_site = max(1, n_site)
        n_elec = max(1, n_elec)

        n_site2 = 2 * n_site
        n_size = 2 * n_elec

        if all_complex
            new(
                zeros(ComplexF64, n_qp_full * n_site2 * n_site2),
                zeros(ComplexF64, n_qp_full * (n_size * n_size + 1)),
                zeros(ComplexF64, n_qp_full),
                Float64[],  # Real versions not used
                Float64[],
                Float64[],
            )
        else
            new(
                zeros(ComplexF64, n_qp_full * n_site2 * n_site2),
                zeros(ComplexF64, n_qp_full * (n_size * n_size + 1)),
                zeros(ComplexF64, n_qp_full),
                zeros(Float64, n_qp_full * n_site2 * n_site2),
                zeros(Float64, n_qp_full * (n_size * n_size + 1)),
                zeros(Float64, n_qp_full),
            )
        end
    end
end

"""
    SamplingWorkspace

Pre-allocated workspace arrays for VMC sampling to avoid repeated allocations.
This significantly reduces memory allocation overhead in hot loops.
"""
mutable struct SamplingWorkspace
    # For calculate_m_all_real!
    inv_m_real_temp::Array{Float64,3}
    pf_m_real_temp::Vector{Float64}

    # For calculate_m_all! (complex version)
    inv_m_temp::Array{ComplexF64,3}
    pf_m_temp::Vector{ComplexF64}

    # For vmc_make_sample_real! / vmc_make_sample!
    proj_cnt_new::Vector{Int}
    pf_m_new_real::Vector{Float64}
    pf_m_new::Vector{ComplexF64}

    # Cached arrays (computed once)
    loc_spn::Vector{Int}

    # Workspace for PfaPack (Pfaffian calculations)
    # Use ThreadedPfaPackWorkspace for parallel execution (one workspace per thread)
    pfapack_workspace::ThreadedPfaPackWorkspace

    # Cached VMCMainCal local accumulator. Kept as Any because VMCThreadAccumulator
    # is defined later in threading.jl.
    main_cal_accumulator::Any

    function SamplingWorkspace(n_size::Int, n_qp_full::Int, n_proj::Int, n_site::Int)
        new(
            zeros(Float64, n_size, n_size, n_qp_full),
            zeros(Float64, n_qp_full),
            zeros(ComplexF64, n_size, n_size, n_qp_full),
            zeros(ComplexF64, n_qp_full),
            zeros(Int, n_proj),
            zeros(Float64, n_qp_full),
            zeros(ComplexF64, n_qp_full),
            zeros(Int, n_site),  # loc_spn will be initialized later
            ThreadedPfaPackWorkspace(n_size),  # Thread-local workspaces for parallel Pfaffian calculations
            nothing,
        )
    end
end

"""
    PhysicalQuantities

Physical quantities for VMC measurement mode (NVMCCalMode=1).
Stores Green's functions and other measured quantities.
"""
mutable struct PhysicalQuantities
    # 1-body Green's function: <c†_{ri,s} c_{rj,s}>
    # Size: NCisAjs
    local_cis_ajs::Vector{ComplexF64}      # Local value for each sample
    phys_cis_ajs::Vector{ComplexF64}        # Weighted average accumulation

    # 2-body Green's function (product): <c†_i c_j> × <c†_k c_l>
    # Size: NCisAjsCktAlt
    phys_cis_ajs_ckt_alt::Vector{ComplexF64}

    # 2-body Green's function (direct): <c†_i c_j c†_k c_l>
    # Size: NCisAjsCktAltDC
    local_cis_ajs_ckt_alt_dc::Vector{ComplexF64}
    phys_cis_ajs_ckt_alt_dc::Vector{ComplexF64}

    # Canonical one-body Green list (ri, si, rj, sj), C-compatible order.
    # When TwoBodyGEx is present this includes appended factored constituents
    # (C's iOneBodyGIdx); otherwise it equals greenone.def order.
    cis_ajs_idx::Vector{NTuple{4,Int}}
    # Factored two-body pairs: 1-based indices into cis_ajs_idx / local_cis_ajs.
    cis_ajs_ckt_alt_idx::Vector{Tuple{Int,Int}}

    # Full Lanczos R1 accumulator: QQQQ tensor flattened in C order.
    phys_lanczos_qqqq::Vector{ComplexF64}

    function PhysicalQuantities(
        n_cis_ajs::Int,
        n_cis_ajs_ckt_alt::Int,
        n_cis_ajs_ckt_alt_dc::Int,
    )
        new(
            zeros(ComplexF64, n_cis_ajs),
            zeros(ComplexF64, n_cis_ajs),
            zeros(ComplexF64, n_cis_ajs_ckt_alt),
            zeros(ComplexF64, n_cis_ajs_ckt_alt_dc),
            zeros(ComplexF64, n_cis_ajs_ckt_alt_dc),
            NTuple{4,Int}[],
            Tuple{Int,Int}[],
            zeros(ComplexF64, 16),
        )
    end
end

"""
    VMCOptimizationState

State data for VMC optimization.
"""
mutable struct VMCOptimizationState
    energy::EnergyData
    slater_matrix::SlaterMatrixData
    electron_config::ElectronConfiguration
    sr_opt::SROptData
    opt_data::Vector{OptDataPoint}
    workspace::SamplingWorkspace
    phys_quantities::Union{PhysicalQuantities,Nothing}  # For VMCPhysCal mode

    function VMCOptimizationState(
        n_site::Int,
        n_elec::Int,
        n_proj::Int,
        n_para::Int,
        n_qp_full::Int,
        n_vmc_sample::Int,
        all_complex::Bool,
        use_fsz::Bool,
    )
        n_size = 2 * n_elec
        new(
            EnergyData(),
            SlaterMatrixData(n_qp_full, n_site, n_elec, all_complex),
            ElectronConfiguration(n_vmc_sample, n_site, n_elec, n_proj, use_fsz),
            SROptData(1 + n_para, n_vmc_sample, all_complex),
            OptDataPoint[],
            SamplingWorkspace(n_size, n_qp_full, n_proj, n_site),
            nothing,  # Initialize as nothing, will be set when needed
        )
    end
end
