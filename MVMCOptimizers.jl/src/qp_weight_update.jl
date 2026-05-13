"""
Quantum Projection Weight Update Functions

Update quantum projection weights.
"""

"""
    update_qp_weight!(data::ExpertModeData)

Update quantum projection full weights.
Equivalent to C's `UpdateQPWeight()`.

Updates QPFullWeight from QPFixWeight using OptTrans if enabled.

This function delegates to MVMCExpertModeParsers.update_qp_weight!.
"""
function update_qp_weight!(data::ExpertModeData)
    # Check if qp_weights is initialized
    if data.qp_weights === nothing
        @warn "Quantum projection weights not initialized. Call init_qp_weight! first."
        return
    end

    # Get OptTrans values if available
    # OptTrans is stored in data.orbital_terms or as a separate field
    # For now, use empty vector (no OptTrans)
    opt_trans = ComplexF64[]

    # Use MVMCExpertModeParsers function
    MVMCExpertModeParsers.update_qp_weight!(data.qp_weights, opt_trans)
end
