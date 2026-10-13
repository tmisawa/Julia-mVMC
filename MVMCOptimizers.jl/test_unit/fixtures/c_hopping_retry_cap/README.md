# C candidate retry-cap fixture

`expected.txt` comes from a standalone probe of the authoritative C candidate
functions in mVMC 1.3.0 (`d73d06bd`): `makeCandidate_hopping`
(`src/mVMC/vmcmake.c` lines 515-545) and `makeCandidate_hopping_fsz`,
`makeCandidate_hopping_csz`, `makeCandidate_LocalSpinFlip_conduction`
(`src/mVMC/vmcmake_fsz.c` lines 594-664 and 756-788), copied verbatim into the
`.inc` files. `probe.c` supplies their globals and a fixed raw-word tape in
place of SFMT (`genrand_real2` is one word times 2^-32), and draws one marker
word after each non-rejected candidate, as the C samplers draw the Metropolis
word only when `rejectFlag` is false.

C draws a destination, then tests `icnt > Nsite*Nsite`, then tests occupancy:
the draw made after the cap is rejected even when it names an empty site. The
ordinary source-selection loop has no cap in C. `result` rows are
`mi ri rj s t reject candidate_words words_after_metropolis`.

The probe and its generator are maintained in mvmc-rs
(`c_toolbox/hopping_retry_559/`, AtelierArith/mvmc-rs#559). Generated on Linux
x86_64 with GCC 13.3.0 (`cc -std=c11 -O0`). To rebuild from an mVMC 1.3.0
checkout, extract the line ranges above into the `.inc` files and run
`cc -std=c11 -O0 probe.c -o probe && ./probe > expected.txt`.

SHA-256:

```text
8431b58eaec53e5325c801b69aa8a56306ec6a154b54e099266f4e21e5dd42ed  mVMC-1.3.0/src/mVMC/vmcmake.c
2b43a47a85ec13833d453a010962334147024a6d6aae01b99bb8b2a0b80c15db  mVMC-1.3.0/src/mVMC/vmcmake_fsz.c
29d7b83e7c5707aa071025ee4b7ae886bf2de8c7b95b2a2b6cc915a8ebc25b7e  hopping_upstream.inc
8d8c96ab0b0a22b8d569ff248d57ff258986df94cc95cc53ccc54573fbdec6db  fsz_upstream.inc
8b175389b16ef651cf38c4782e4c118915505f7c4443b9741fc4f4637aa32446  probe.c
1b3354bafaf3567b4777994b28d117cd7e528d4f8f541a74e9f060dbd236330c  expected.txt
```
