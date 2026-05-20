# examples/heisenberg_chain_real.jl
# Run VMCParaOpt on a 16-site Heisenberg chain (real mode).
# Inputs are in examples/inputs/heisenberg_chain_real/.
# Step count is overridable via JULIA_MVMC_EXAMPLE_STEPS (default 50).

using MVMCOptimizers

const INPUT_DIR = joinpath(@__DIR__, "inputs", "heisenberg_chain_real")
const NAMELIST  = joinpath(INPUT_DIR, "namelist.def")
const NSTEPS    = parse(Int, get(ENV, "JULIA_MVMC_EXAMPLE_STEPS", "50"))

println("=== Heisenberg chain (real) — $(NSTEPS) SR steps ===")
result = MVMCOptimizers.run_para_opt_from_namelist(
    NAMELIST;
    nsteps = NSTEPS,
    nsmp = NSTEPS,
    mode = :real,
)
println("Final energy / site = ", result.final_energy_per_site)
