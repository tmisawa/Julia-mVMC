# C-mVMC ctest models supported by the Julia ctest-equivalent runner.
#
# `fixture` is the public Julia-mVMC fixture directory name.
# `c_model` is the C-mVMC `runtest.py <ModelName>` argument.
# The 2 deferred models complete the private 15-model set, but need additional
# Julia support before they can use this runner:
#   - SpinChainLanczos, HubbardChainLanczos: mode1/Lanczos output path.

const CTEST_STANDARD_MODELS = [
    (
        fixture = "heisenberg_chain_real",
        c_model = "HeisenbergChain",
        mode = :real,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
    (
        fixture = "hubbard_chain_real",
        c_model = "HubbardChain",
        mode = :real,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
    (
        fixture = "hubbard_tetragonal_real",
        c_model = "HubbardTetragonal",
        mode = :real,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
    (
        fixture = "hubbard_tetragonal_momentum_projection_real",
        c_model = "HubbardTetragonal_MomentumProjection",
        mode = :real,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
    (
        fixture = "kondo_chain_real",
        c_model = "KondoChain",
        mode = :real,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
    (
        fixture = "heisenberg_chain_cmp",
        c_model = "HeisenbergChain_cmp",
        mode = :cmp,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
    (
        fixture = "hubbard_chain_cmp",
        c_model = "HubbardChain_cmp",
        mode = :cmp,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
    (
        fixture = "kondo_chain_cmp",
        c_model = "KondoChain_cmp",
        mode = :cmp,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
    (
        fixture = "kondo_chain_stot1_cmp",
        c_model = "KondoChain_Stot1_cmp",
        mode = :cmp,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
    (
        fixture = "general_rbm_cmp",
        c_model = "GeneralRBM_cmp",
        mode = :cmp,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
    (
        fixture = "heisenberg_chain_fsz",
        c_model = "HeisenbergChain_fsz",
        mode = :fsz,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
    (
        fixture = "hubbard_chain_fsz",
        c_model = "HubbardChain_fsz",
        mode = :fsz,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
    (
        fixture = "kondo_chain_fsz",
        c_model = "KondoChain_fsz",
        mode = :fsz,
        nsteps_override = nothing,
        nsmp_override = nothing,
    ),
]

const CTEST_DEFERRED_MODELS = [
    (
        fixture = "spin_chain_lanczos",
        c_model = "SpinChainLanczos",
        mode = :real,
        reason = "NVMCCalMode=1/Lanczos output path is not wired to run_para_opt_from_namelist",
    ),
    (
        fixture = "hubbard_chain_lanczos",
        c_model = "HubbardChainLanczos",
        mode = :real,
        reason = "NVMCCalMode=1/Lanczos output path is not wired to run_para_opt_from_namelist",
    ),
]
