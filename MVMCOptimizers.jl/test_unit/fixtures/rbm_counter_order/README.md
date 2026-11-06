# RBM hidden-counter arithmetic order

Authoritative source: mVMC 1.3.0 `src/mVMC/rbm.c`, `MakeRBMCnt`, lines
187–288. Upstream revision `d73d06bd529d3b2573f38eb5817c4a5f52971006`;
whole-file SHA-256
`e35af053ce06e127309b1636dd2776650b4f1b2abd1ea0eaca10cd5257c3d742`.
`probe.c` contains the verbatim function plus minimal parameter/mapping inputs;
its upstream GPL notice is preserved. It is an optional standalone kernel
probe, not a C executable/MPI/sampling-trajectory validation.

For each Charge/Spin/General hidden neuron, C starts the coupling sum at zero,
then adds that sum to the already initialized hidden bias once. Julia previously
added each coupling directly to the bias. Bias `1+2i` and coupling terms
`+(2^54+2^55 i)` and `-(2^54+2^55 i)` distinguish these operations exactly:
the C result is the bias, whereas the old Julia result is zero. Real-only cases
use zero imaginary parts. Every mapping is complete for the two-site input;
General uses spin-major order and zero-valued down-spin couplings.

Reproduce the optional probe from the fork root:

```sh
cc -std=c11 -O0 -ffp-contract=off \
  MVMCOptimizers.jl/test_unit/fixtures/rbm_counter_order/probe.c \
  -o /tmp/rbm-counter-order-probe
/tmp/rbm-counter-order-probe
```

Verified on Linux x86_64 with Ubuntu GCC 13.3.0: all six cases return `1+0i`
or `1+2i`. Normal Julia tests use immutable analytical expectations only;
they do not compile or invoke this probe. No floating-point tolerance is used.
This patch preserves existing term traversal and incremental hopping updates;
it changes only full-counter bias/coupling grouping. It does not by itself
establish full-model sampling parity or justify any numerical-budget change.

## mVMC 1.4.0 re-verification (Julia-mVMC#95, 2026-11-06)

The C reference moved to mVMC 1.4.0 (`v1.4.0`, `e83aa6cb`). `MakeRBMCnt` is
byte-identical and still at `rbm.c` lines 187–288 (1.4.0 whole-file SHA256
`bbe4a1ab7d6f3ab0cbaa41d6ffabe10f729ce5bc560a22802a6939e287fbd870`; the only
change elsewhere in `rbm.c` is `logSqPfFullSlater += 2*creal(x)` instead of
`2*cabs(x)`, used only by 1.4.0 reweighting, which Julia does not port). The probe
and its expectations are unchanged.
