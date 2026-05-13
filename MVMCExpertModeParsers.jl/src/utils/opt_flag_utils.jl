"""
Optimization Flag Utilities

Functions for checking if parameters are optimization targets.
Corresponds to C implementation's OptFlag array.
"""

"""
    is_slater_optimized(data::ExpertModeData, slater_idx::Int) -> Bool

Check if a Slater parameter (OrbitalTerm) is an optimization target.

# Arguments
- `data::ExpertModeData`: Expert Mode data structure
- `slater_idx::Int`: Slater parameter index (1-based)

# Returns
- `Bool`: true if optimized, false if fixed

# C Implementation Correspondence
- C: OptFlag[2*i+2*NProj + 2*FlagRBM*NRBM] > 0
- Julia: optimization_flags[2*(i-1) + 2*n_proj + 2*flag_rbm*n_rbm + 1]
"""
function is_slater_optimized(data::ExpertModeData, slater_idx::Int)::Bool
    n_proj = length(data.gutzwiller_terms) + length(data.jastrow_terms)

    # FlagRBMの判定
    flag_rbm = (
        length(data.charge_rbm_phys_layer_terms) > 0 ||
        length(data.spin_rbm_phys_layer_terms) > 0 ||
        length(data.general_rbm_phys_layer_terms) > 0 ||
        length(data.charge_rbm_hidden_layer_terms) > 0 ||
        length(data.spin_rbm_hidden_layer_terms) > 0 ||
        length(data.general_rbm_hidden_layer_terms) > 0 ||
        length(data.charge_rbm_phys_hidden_terms) > 0 ||
        length(data.spin_rbm_phys_hidden_terms) > 0 ||
        length(data.general_rbm_phys_hidden_terms) > 0
    )

    n_rbm = (flag_rbm ? count_rbm_parameters(data) : 0)

    # OptFlagのインデックス計算（Juliaは1-based）
    opt_flag_idx = 2 * (slater_idx - 1) + 2 * n_proj + 2 * flag_rbm * n_rbm + 1

    # 範囲チェック
    if opt_flag_idx > length(data.optimization_flags)
        return false
    end

    return data.optimization_flags[opt_flag_idx]
end

"""
    is_gutzwiller_optimized(data::ExpertModeData, gutzwiller_idx::Int) -> Bool

Check if a Gutzwiller parameter is an optimization target.

# C Implementation Correspondence
- C: OptFlag[2*i] > 0
- Julia: optimization_flags[2*(i-1) + 1]
"""
function is_gutzwiller_optimized(data::ExpertModeData, gutzwiller_idx::Int)::Bool
    opt_flag_idx = 2 * (gutzwiller_idx - 1) + 1

    if opt_flag_idx > length(data.optimization_flags)
        return false
    end

    return data.optimization_flags[opt_flag_idx]
end

"""
    is_jastrow_optimized(data::ExpertModeData, jastrow_idx::Int) -> Bool

Check if a Jastrow parameter is an optimization target.

# C Implementation Correspondence
- C: OptFlag[2*(NGutzwillerIdx + i)] > 0
- Julia: optimization_flags[2*(NGutzwillerIdx + i - 1) + 1]
"""
function is_jastrow_optimized(data::ExpertModeData, jastrow_idx::Int)::Bool
    n_gutzwiller = length(data.gutzwiller_terms)
    opt_flag_idx = 2 * (n_gutzwiller + jastrow_idx - 1) + 1

    if opt_flag_idx > length(data.optimization_flags)
        return false
    end

    return data.optimization_flags[opt_flag_idx]
end

"""
    get_slater_opt_flag_index(data::ExpertModeData, slater_idx::Int) -> Int

Get the OptFlag index for a Slater parameter.

# Returns
- `Int`: Index in optimization_flags array (1-based)
"""
function get_slater_opt_flag_index(data::ExpertModeData, slater_idx::Int)::Int
    n_proj = length(data.gutzwiller_terms) + length(data.jastrow_terms)

    flag_rbm = (
        length(data.charge_rbm_phys_layer_terms) > 0 ||
        length(data.spin_rbm_phys_layer_terms) > 0 ||
        length(data.general_rbm_phys_layer_terms) > 0 ||
        length(data.charge_rbm_hidden_layer_terms) > 0 ||
        length(data.spin_rbm_hidden_layer_terms) > 0 ||
        length(data.general_rbm_hidden_layer_terms) > 0 ||
        length(data.charge_rbm_phys_hidden_terms) > 0 ||
        length(data.spin_rbm_phys_hidden_terms) > 0 ||
        length(data.general_rbm_phys_hidden_terms) > 0
    )

    n_rbm = (flag_rbm ? count_rbm_parameters(data) : 0)

    return 2 * (slater_idx - 1) + 2 * n_proj + 2 * flag_rbm * n_rbm + 1
end

"""
    print_optimization_status(data::ExpertModeData)

Print optimization status for all parameters.
"""
function print_optimization_status(data::ExpertModeData)
    println("=== パラメータの最適化状態 ===")

    # Gutzwiller
    if !isempty(data.gutzwiller_terms)
        println("\nGutzwillerパラメータ:")
        for (i, term) in enumerate(data.gutzwiller_terms)
            is_opt = is_gutzwiller_optimized(data, i)
            println("  Gutzwiller[$i] (site $term.site): $(is_opt ? "最適化対象" : "固定")")
        end
    end

    # Jastrow
    if !isempty(data.jastrow_terms)
        println("\nJastrowパラメータ:")
        for (i, term) in enumerate(data.jastrow_terms)
            is_opt = is_jastrow_optimized(data, i)
            println(
                "  Jastrow[$i] (site $term.site1-$term.site2): $(is_opt ? "最適化対象" : "固定")",
            )
        end
    end

    # Slater
    if !isempty(data.orbital_terms)
        println("\nSlaterパラメータ:")
        optimized_slater = 0
        for (i, term) in enumerate(data.orbital_terms)
            is_opt = is_slater_optimized(data, i)
            if is_opt
                optimized_slater += 1
            end
            if i <= 5 || i > length(data.orbital_terms) - 5
                println(
                    "  Slater[$i] (site $term.site1-$term.site2): $(is_opt ? "最適化対象" : "固定")",
                )
            elseif i == 6
                println("  ...")
            end
        end
        println(
            "  合計: 最適化対象=$optimized_slater, 固定=$(length(data.orbital_terms) - optimized_slater)",
        )
    end
end
