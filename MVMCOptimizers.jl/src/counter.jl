"""
Counter Reduction Functions

Reduce counters across MPI processes (for compatibility with C implementation).
"""

"""
    reduce_counter!(state::VMCOptimizationState)

Reduce counters across processes.
Equivalent to C's `ReduceCounter()`.

In the C implementation, this function uses MPI_Allreduce to sum counters
across all MPI processes. In the Julia implementation (single process),
this is a no-op but kept for compatibility with the C code structure.

# Arguments
- `state::VMCOptimizationState`: Optimization state containing counters

# Note
- This function is called after `weight_average_sr_opt!()` in the optimization loop
- In a single-process implementation, counters are already local, so no reduction is needed
"""
function reduce_counter!(state::VMCOptimizationState)
    # In single-process implementation, no reduction is needed
    # Counters are already local to the process
    # This function is kept for compatibility with C code structure
    # where MPI_Allreduce is used to sum counters across processes

    # If MPI support is added in the future, implement reduction here
    # For now, this is a no-op
    return
end

"""
    reduce_counter!(ctx, state)

C `ReduceCounter(comm_child2)` の state 版。MPI では `Counter_max=6` 相当だけを
`ctx.comm2` で allreduce し、`ctx.rank2 == 0` の rank にだけ書き戻す。
serial context では既存 `reduce_counter!(state)` と同じ no-op。
"""
function reduce_counter!(ctx::ParallelContext, state::VMCOptimizationState)
    ctx.is_mpi || return reduce_counter!(state)
    reduce_counter!(ctx, state.electron_config.counter)
    return
end
