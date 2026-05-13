# MVMCExpertModeParsers.jl

> Part of [**Julia-mVMC**](https://github.com/tmisawa/Julia-mVMC) — the Julia port of [issp-center-dev/mVMC](https://github.com/issp-center-dev/mVMC). See the top-level repository for installation guidance and the full user manual.

A Julia package for parsing mVMC Expert Mode definition files (`*.def`) and initializing variational parameters / quantum-projection weights from them. Used internally by [`MVMCOptimizers.jl`](../MVMCOptimizers.jl).

## Scope

- **Parsing**: read `.def` files in mVMC Expert Mode format (`modpara.def`, `trans.def`, `coulombintra.def`, `gutzwilleridx.def`, `jastrowidx.def`, `orbitalidx.def`, `qptransidx.def`, `namelist.def`, …).
- **Parameter initialization** (mirrors C's `InitParameter()` + `SyncModifiedParameter()`).
- **Input parameter overlay** from `In*.def` files (mirrors C's `ReadInputParameters()`).
- **Quantum-projection weight initialization** (Gauss-Legendre + Legendre polynomials, mirrors C's `InitQPWeight()` / `UpdateQPWeight()`).
- **Validation helpers** for individual term types.

Note: this v0.1 release does **not** include writers / `.def` file generation. Reading and consuming existing `.def` files is the supported path; producing new ones from scratch is a v0.2+ concern.

## Installation

This subpackage is **not** published as a standalone registered package and `Pkg.add("MVMCExpertModeParsers")` will not work in v0.1. It is intended to be used as part of the Julia-mVMC workspace via the root [`Project.toml`](../Project.toml) `[sources]` block. From a clone of [tmisawa/Julia-mVMC](https://github.com/tmisawa/Julia-mVMC):

```bash
git clone --recurse-submodules https://github.com/tmisawa/Julia-mVMC
cd Julia-mVMC
julia --project=@. -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
```

After that, `using MVMCExpertModeParsers` works inside any code activated against the root workspace.

## Quick start

```julia
using MVMCExpertModeParsers

# Parse all .def files referenced by namelist.def
data = parse_expert_mode_files("namelist.def")

println("NSite: $(data.modpara.nsite)")
println("NElec: $(data.modpara.nelec)")
println("Transfer terms: $(length(data.transfer_terms))")
```

`parse_expert_mode_files` returns a fully-populated `ExpertModeData` value. The type itself is not exported; if you need to refer to it by name, qualify it: `MVMCExpertModeParsers.ExpertModeData`.

### Parameter initialization

Most callers will go through `MVMCOptimizers.run_para_opt_from_namelist`, which wraps the steps below in the correct C-compatible order. If you need to do it by hand, the symbols live inside the module but are not exported, so call them qualified:

```julia
using Random
using MVMCExpertModeParsers
using SFMT  # provides SFMT19937RNG, the C-compatible RNG

data = parse_expert_mode_files("namelist.def")

# Optimization flag layout (which entries of Proj / Slater are optimized).
n_proj = length(data.gutzwiller_terms) + length(data.jastrow_terms)
n_slater = length(data.orbital_terms)
data.optimization_flags = vcat(
    repeat([true, false], n_proj),
    repeat([true, false], n_slater),
)

data.modpara.rnd_seed = 123456789
rng = SFMT19937RNG()
Random.seed!(rng, data.modpara.rnd_seed)

# C: InitParameter + SyncModifiedParameter
MVMCExpertModeParsers.initialize_parameters!(data; rng=rng)

# C: ReadInputParameters (signature: data, namelist_path)
MVMCExpertModeParsers.read_input_parameters!(data, "namelist.def")

# C: InitQPWeight
MVMCExpertModeParsers.init_qp_weight!(data)
```

## Supported file types

| File              | Read |
|-------------------|------|
| `modpara.def`     | ✓ |
| `trans.def`       | ✓ |
| `coulombintra.def`| ✓ |
| `coulombinter.def`| ✓ |
| `hund.def`        | ✓ |
| `exchange.def`    | ✓ |
| `pairhop.def`     | ✓ |
| `gutzwilleridx.def` | ✓ |
| `jastrowidx.def`  | ✓ |
| `orbitalidx.def`  | ✓ |
| `greenone.def`    | ✓ |
| `greentwo.def`    | ✓ |
| `qptransidx.def`  | ✓ |
| `namelist.def`    | ✓ |

Writing back out to `.def` files is not implemented in v0.1.

## Parameter initialization details

The `initialize_parameters!()` function follows the same logic as the C implementation:

- **Proj parameters** (Gutzwiller + Jastrow): set to `0.0`.
- **RBM parameters**: randomly initialized if `FlagRBM > 0` and `OptFlag > 0`.
  - Real: `RBM[i] = 0.01*(rand() - 0.5)` (uniform).
  - Complex: `RBM[i] = 1e-2*rand()*exp(2π*I*rand())`.
- **Slater parameters** (Orbital): randomly initialized if `OptFlag > 0`.
  - Real: `Slater[i] = 2*(rand() - 0.5)` (uniform on `[-1, 1)`).
  - Complex: `Slater[i] = (2*(rand()-0.5) + 2*I*(rand()-0.5)) / sqrt(2.0)`.
- **Scaling**: Slater parameters are rescaled so `max(|Slater[i]|) ≤ 4.0`.

### `read_input_parameters!`

Reads from `In*.def` files referenced by `namelist.def` and overlays the values onto the already-initialized `data`:

- `InGutzwiller.def` → `data.gutzwiller_terms`.
- `InJastrow.def` → `data.jastrow_terms`.
- `InOrbital.def` → `data.orbital_terms`.
- `InChargeRBM_PhysLayer.def`, etc. → RBM term vectors.

Signature: `read_input_parameters!(data::ExpertModeData, namelist_path::String)`. The base directory is taken from the namelist path itself; pass a relative or absolute path to `namelist.def`.

### `init_qp_weight!`

Initializes `data.qp_weights` (Gauss-Legendre quadrature for spin projection, Legendre polynomials, `QPFixWeight`, `QPFullWeight`). Equivalent to C's `InitQPWeight()`.

```julia
data.modpara.nsp_gauss_leg = 4
data.modpara.nsp_stot = 1
data.modpara.nmp_trans = 1
data.para_qp_trans = [1.0+0.0im]

MVMCExpertModeParsers.init_qp_weight!(data)
```

## C implementation compatibility

| C function | Julia function |
|------------|----------------|
| `InitParameter` + `SyncModifiedParameter` | `MVMCExpertModeParsers.initialize_parameters!` |
| `ReadInputParameters` | `MVMCExpertModeParsers.read_input_parameters!` |
| `InitQPWeight` | `MVMCExpertModeParsers.init_qp_weight!` |
| `UpdateQPWeight` | `MVMCExpertModeParsers.update_qp_weight!` |
| `GaussLeg` | `MVMCExpertModeParsers.gauss_legendre` |
| `LegendrePoly` | `MVMCExpertModeParsers.legendre_poly` |

File format compatibility:

- 18-digit precision for floating-point inputs.
- Parses every `.def` file produced by C mVMC's Expert Mode.
- Uses the same SFMT19937 random stream as the C implementation (via [SFMT.jl](../SFMT.jl)).

## Validation

Individual validators are provided per term type:

```julia
result = MVMCExpertModeParsers.validate_modpara_params(data.modpara)
if !result.is_valid
    for err in result.errors
        println("Error: $err")
    end
end
```

See `src/utils/validation.jl` for the full list (`validate_transfer_terms`, `validate_coulomb_intra_terms`, `validate_gutzwiller_terms`, …).

## Testing

```bash
cd MVMCExpertModeParsers.jl
julia --project=@. -e 'using Pkg; Pkg.test()'
```

The integration tests in the root workspace (`test/integration/runtests.jl` from the repository root) additionally exercise this parser in combination with the rest of the Julia-mVMC stack.

## License

This subpackage is part of [Julia-mVMC](https://github.com/tmisawa/Julia-mVMC) and is licensed under **GPL-3.0-or-later**. See the subpackage [`LICENSE`](LICENSE) (a verbatim copy of GPLv3) and the workspace-level [`THIRD_PARTY_LICENSES.md`](../THIRD_PARTY_LICENSES.md) for bundled-source notices.

## References

- [mVMC](https://github.com/issp-center-dev/mVMC) — original C reference implementation.
