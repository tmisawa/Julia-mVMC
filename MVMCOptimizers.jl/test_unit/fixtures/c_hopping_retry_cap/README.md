# C candidate retry-cap fixture

`expected.txt` comes from a standalone probe of the authoritative C candidate
functions in mVMC 1.4.0 (`v1.4.0`, `e83aa6cb`): `makeCandidate_hopping`
(`src/mVMC/vmcmake.c` lines 620-650) and `makeCandidate_hopping_fsz`,
`makeCandidate_hopping_csz`, `makeCandidate_LocalSpinFlip_conduction`
(`src/mVMC/vmcmake_fsz.c` lines 633-703 and 866-898), copied verbatim into the
`.inc` files. These bodies are byte-identical to 1.3.0 (`d73d06bd`, lines
515-545, 594-664 and 756-788); the `.inc` files were re-extracted from 1.4.0 on
2026-11-06 and the rebuilt probe reproduced `expected.txt` byte for byte. `probe.c` supplies their globals and a fixed raw-word tape in
place of SFMT (`genrand_real2` is one word times 2^-32), and draws one marker
word after each non-rejected candidate, as the C samplers draw the Metropolis
word only when `rejectFlag` is false.

C draws a destination, then tests `icnt > Nsite*Nsite`, then tests occupancy:
the draw made after the cap is rejected even when it names an empty site. The
ordinary source-selection loop has no cap in C. `result` rows are
`mi ri rj s t reject candidate_words words_after_metropolis`.

The probe and its generator are maintained in mvmc-rs
(`c_toolbox/hopping_retry_559/`, AtelierArith/mvmc-rs#559). Generated on Linux
x86_64 with GCC 13.3.0 (`cc -std=c11 -O0`). To rebuild from an mVMC 1.4.0
checkout, extract the line ranges above into the `.inc` files and run
`cc -std=c11 -O0 probe.c -o probe && ./probe > expected.txt`.

SHA-256:

```text
a9518a622b4c170e0f0d6bca98270425a0042ccbe1c1deb1106c8bbec0a28fc9  mVMC-1.4.0/src/mVMC/vmcmake.c
f6322ed48eb186c097074ea9287906fcf58129a44d91b932d0a6750ea9daba66  mVMC-1.4.0/src/mVMC/vmcmake_fsz.c
417af45c4352588f60d2d62219783bd3ef918aaab0e08f7ff72c35402829ab2d  hopping_upstream.inc
f1b120b5111b91282595878400284066b0fcbe3c46f17e96364b0b7a801aff93  fsz_upstream.inc
8b175389b16ef651cf38c4782e4c118915505f7c4443b9741fc4f4637aa32446  probe.c
1b3354bafaf3567b4777994b28d117cd7e528d4f8f541a74e9f060dbd236330c  expected.txt
```

Historical 1.3.0 extraction (`d73d06bd`): `vmcmake.c` SHA-256
`8431b58eaec53e5325c801b69aa8a56306ec6a154b54e099266f4e21e5dd42ed`,
`vmcmake_fsz.c` `2b43a47a85ec13833d453a010962334147024a6d6aae01b99bb8b2a0b80c15db`,
`hopping_upstream.inc` `29d7b83e…`, `fsz_upstream.inc` `8d8c96ab…` (same
function bodies; only the line-number comment differed).
