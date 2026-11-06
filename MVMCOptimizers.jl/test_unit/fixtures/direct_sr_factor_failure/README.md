# Direct SR factorization failure

Source base: Julia-mVMC PR54 commit
`62b0f97f076fb55c71c3ab0caa041a9adff94e04`. This regression changes no
sampling, input, retained-parameter layout or CG algorithm.

C authority is mVMC-1.3.0 revision
`d73d06bd529d3b2573f38eb5817c4a5f52971006`:

- `src/mVMC/stcopt_dposv.c:32–49` calls `M_DPOSV` and returns its status.
  SHA256: `2bd48d880dcbd95ea1b1b92931178c07fd08f981b8ebf04e7db1709e909e57ca`.
- `src/mVMC/stcopt.c:142–174` propagates failure and updates parameters only
  when status is zero.
  SHA256: `43ed8790cff2715284849f0f8906e4645179dcbb100b51f819918b279d6a36f2`.

DPOSV stops after positive Cholesky `info`, without substituting into the
right-hand side. Julia `potrf!` returns this positive status rather than
throwing; ignoring it previously called `potrs!` on a partial factor.
Julia's existing public success/failure convention is retained: positive
factorization failures map to status1, not the raw leading-minor index.
The scratch matrix may be modified by factorization, as in C.

The immutable matrices require no native oracle to regenerate:

| Matrix | Positive POTRF info | Reason |
| --- | --- | --- |
| `[0 0; 0 0]` | 1 | First pivot is zero |
| `[1 1; 1 1]` | 2 | Second Schur pivot is zero |
| `[1 2; 2 1]` | 2 | Second Schur pivot is -3 |

The private production solve phase checks RHS `[-.5,-.5]` remains unchanged
for all three failures. A diagonal SPD control has the analytical solution
`[-.5,.5]`; its only numerical budget is `2eps(Float64)` absolute, zero
relative, for the small Cholesky/substitution sequence.

Six actual `stochastic_opt!` cases cover real input conversion and complex
input storage separately. Assertions verify the exact assembled matrices
and RHS, status1, unchanged packed parameters, flags and source HO arrays.
The default Julia RNG's next624 UInt32 outputs are compared non-consumingly;
this is a no-draw check, not a C/SFMT sampling-trajectory claim.

Unpatched production at 62b0f97, run38214: 26 assertions passed, four failed.
Both indefinite cases returned success and changed parameters finitely.
This is not a test that merely detects NaNs after a singular solve.

Reproduce in Julia1.13.1 using the workspace Manifest-v1.13.toml (not part
of this patch), Linux x86_64, ILP64 OpenBLAS, BLAS threads1:

```sh
julia +1.13.1 --project=. -e 'using LinearAlgebra; BLAS.set_num_threads(1); println(VERSION); println(BLAS.get_config()); include("MVMCOptimizers.jl/test_unit/test_unit_direct_sr_factor_failure.jl")'
```

Normal Julia tests read literal inputs only; no C compiler, native probe,
Rust result or fixture regeneration dependency. C authority above is source
and LAPACK control-flow evidence, not a full native C/MPI execution claim.
The focused test contains 59 assertions. Broader MPI/runner recovery is not
established by this regression.

## mVMC 1.4.0 re-verification (Julia-mVMC#95, 2026-11-06)

The C reference moved to mVMC 1.4.0 (`v1.4.0`, `e83aa6cb`).
`src/mVMC/stcopt_dposv.c` and `src/mVMC/stcopt.c` are byte-identical in 1.4.0
(same SHA256 as above, same line numbers), so this contract is unchanged.
