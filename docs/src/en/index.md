# Overview

[日本語](../ja/index.md)

Julia-mVMC is a Julia implementation of the mVMC many-variable variational
Monte Carlo workflow. It reads C-mVMC expert-mode `.def` input files and focuses
on reproducible compatibility with the C reference implementation.

The v0.5.0 manual describes the public, verified scope of the current release.
It is not a full replacement for every C-mVMC feature. Unsupported combinations
are documented explicitly and should fail early rather than silently producing
incompatible results.

## v0.5.0 status

| Component | Status | Notes |
|-----------|--------|-------|
| `VMCParaOpt` | Verified | Strict first-10-step C-reference checks, C ctest-equivalent gates, and serial / MPI `NSRCG = 1` first-step tolerance gates. |
| `VMCPhysCal` | Verified for supported paths | One-body, direct two-body, and factored/product two-body Green functions are C-referenced for supported fixtures. |
| Full Lanczos PhysCal | Verified for supported paths | `NLanczosMode = 1/2` output is C-referenced on the sz-conserved `NSplitSize = 1` path. |
| MPI parallelization | Partially verified / smoke-tested | Direct-SR split, standard-projection split, SR-CG `NSplitSize = 1`, rank0 output, and PhysCal split paths have targeted smoke gates. |
| Shared-memory threading | Experimental opt-in | Conservative inner-loop opt-ins only. Sample-level Markov-chain threading is intentionally disabled for C parity. |
| BackFlow | Not supported | Planned for a future release. |

## Wave-function convention

Julia-mVMC follows the same high-level variational form as C-mVMC:

```math
|\Psi\rangle = \mathcal{P} |\Phi\rangle,
```

where ``|\Phi\rangle`` is the Slater or generalized orbital reference state and
``\mathcal{P}`` collects projection and correlation factors such as Gutzwiller,
Jastrow, doublon-holon, and quantum-number projections.

The sampled variational energy is

```math
E = \frac{\langle \Psi | H | \Psi \rangle}
         {\langle \Psi | \Psi \rangle}.
```

The implementation goal for v0.5.0 is not to introduce new stochastic
semantics. Verified paths should match the C reference within the documented
floating-point tolerances.

## Manual map

- [Installation](installation.md): native prerequisites and clone-based setup.
- [Tutorial](tutorial.md): run a bundled Heisenberg-chain optimization.
- [Input files](input_files.md): supported expert-mode `.def` files.
- [Optimization](optimization.md): `VMCParaOpt` usage and MPI scope.
- [Physics calculation](physics_calc.md): `VMCPhysCal` Green and Lanczos output.
- [Output files](output_files.md): main `zvo_*` and `zqp_*` files.
- [Compatibility](compatibility.md): C-reference tests and unsupported features.
