// Standalone C-authority probe for issue #559. See README.md; never built by Cargo.
// The upstream candidate functions are included verbatim; this driver only supplies
// their globals and a raw-word stream in place of SFMT.
#include <assert.h>
#include <stdint.h>
#include <stdio.h>

static int Nsite, Ne, Nsize;
static int *LocSpn;
static uint32_t tape[256];
static int size, cursor;
static uint32_t gen_rand32(void) { assert(cursor < size); return tape[cursor++]; }
static double genrand_real2(void) { return gen_rand32() * (1.0 / 4294967296.0); }

#include "hopping_upstream.inc"
#include "fsz_upstream.inc"

enum kind { HOPPING, HOPPING_FSZ, HOPPING_CSZ, CONDUCTION };
static const char *kind_name[] = {"hopping", "hopping_fsz", "hopping_csz", "conduction"};

static void print_ints(const char *label, const int *v, int n) {
  printf("%s", label);
  for (int i = 0; i < n; i++) printf(" %d", v[i]);
  puts("");
}

// One candidate, then the Metropolis word that the sampling engines draw only for a
// non-rejected candidate (`if(rejectFlag) continue;` in vmcmake*.c). The extra word
// is a marker; the decision itself is outside this probe.
static void run_case(const char *name, enum kind k, int nsite, int nelec, int *loc,
                     int *idx, int *cfg, int *spn, const uint32_t *words, int nwords) {
  Nsite = nsite; LocSpn = loc;
  int nidx;
  if (k == HOPPING) { Ne = nelec; Nsize = 2 * nelec; nidx = 2 * nelec; }
  else { Nsize = nelec; Ne = nelec / 2; nidx = nelec; }
  size = 0; cursor = 0;
  for (int i = 0; i < nwords; i++) tape[size++] = words[i];
  tape[size++] = UINT32_C(0x7fffffff); // Metropolis marker
  int num[2 * 8] = {0};
  for (int r = 0; r < 2 * nsite; r++) num[r] = cfg[r] != -1;
  int mi = -1, ri = -1, rj = -1, s = -1, t = -1, reject = -1;
  switch (k) {
  case HOPPING: makeCandidate_hopping(&mi, &ri, &rj, &s, &reject, idx, cfg); t = s; break;
  case HOPPING_FSZ: makeCandidate_hopping_fsz(&mi, &ri, &rj, &s, &t, &reject, idx, cfg, num, spn); break;
  case HOPPING_CSZ: makeCandidate_hopping_csz(&mi, &ri, &rj, &s, &t, &reject, idx, cfg, num, spn); break;
  case CONDUCTION: makeCandidate_LocalSpinFlip_conduction(&mi, &ri, &rj, &s, &t, &reject, idx, cfg, num, spn); break;
  }
  const int candidate_words = cursor;
  if (!reject) (void)genrand_real2();
  printf("case %s %s\n", name, kind_name[k]);
  printf("nsite %d\nnelec %d\n", nsite, nelec);
  print_ints("loc_spn", loc, nsite);
  print_ints("ele_idx", idx, nidx);
  print_ints("ele_cfg", cfg, 2 * nsite);
  print_ints("ele_spn", spn, k == HOPPING ? 0 : nidx);
  printf("words");
  for (int i = 0; i < nwords; i++) printf(" %u", words[i]);
  puts("");
  printf("result %d %d %d %d %d %d %d %d\nend\n", mi, ri, rj, s, t, reject, candidate_words, cursor);
}

#define N(a) (int)(sizeof(a) / sizeof((a)[0]))
#define UP 0u
#define DN UINT32_C(0x80000000)

int main(void) {
  // Ordinary hopping, issue #559 minimal state: Nsite=2, Ne=1, both electrons on
  // site 0, so only site 1 is empty in either spin sector (cap Nsite^2 = 4).
  {
    int loc[] = {0, 0}, idx[] = {0, 0}, cfg[] = {0, -1, 0, -1};
    uint32_t final_empty[] = {0, UP, 0, 0, 0, 0, 0, 1};
    run_case("final_empty_destination_rejected", HOPPING, 2, 1, loc, idx, cfg, NULL, final_empty, N(final_empty));
    uint32_t last_permitted[] = {0, UP, 0, 0, 0, 0, 1};
    run_case("last_permitted_destination", HOPPING, 2, 1, loc, idx, cfg, NULL, last_permitted, N(last_permitted));
    uint32_t first_empty[] = {0, DN, 1};
    run_case("first_destination_empty", HOPPING, 2, 1, loc, idx, cfg, NULL, first_empty, N(first_empty));
    uint32_t all_occupied[] = {0, UP, 0, 0, 0, 0, 0, 0};
    run_case("all_destinations_occupied", HOPPING, 2, 1, loc, idx, cfg, NULL, all_occupied, N(all_occupied));
  }
  // Same contract at Nsite=4 (cap 16): 17 occupied draws then an empty one.
  {
    int loc[] = {0, 0, 0, 0}, idx[] = {0, 0}, cfg[] = {0, -1, -1, -1, 0, -1, -1, -1};
    uint32_t w[2 + 18] = {0, DN};
    w[2 + 17] = 3;
    run_case("nsite4_final_empty_destination_rejected", HOPPING, 4, 1, loc, idx, cfg, NULL, w, N(w));
    uint32_t v[2 + 17] = {0, DN};
    v[2 + 16] = 3;
    run_case("nsite4_last_permitted_destination", HOPPING, 4, 1, loc, idx, cfg, NULL, v, N(v));
  }
  // Ordinary source selection has no retry cap in C: eleven local-spin draws
  // (more than Nsite^2 + 1 = 10) are followed by an itinerant electron.
  {
    int loc[] = {1, 0, 0}, idx[] = {0, 1}, cfg[] = {0, -1, -1, -1, 0, -1};
    uint32_t w[2 * 12 + 1];
    for (int i = 0; i < 11; i++) { w[2 * i] = 0; w[2 * i + 1] = UP; }
    w[22] = 0; w[23] = DN; w[24] = 2;
    run_case("source_selection_beyond_cap", HOPPING, 3, 1, loc, idx, cfg, NULL, w, N(w));
  }
  // FSZ siblings (audit): Nsize=2, doublon on site 0, site 1 empty (cap 4).
  {
    int loc[] = {0, 0}, idx[] = {0, 0}, cfg[] = {0, -1, 1, -1}, spn[] = {0, 1};
    uint32_t fsz_final[] = {0, 0, UP, 0, DN, 0, UP, 0, DN, 0, UP, 1, UP};
    run_case("fsz_final_empty_destination_rejected", HOPPING_FSZ, 2, 2, loc, idx, cfg, spn, fsz_final, N(fsz_final));
    uint32_t fsz_last[] = {0, 0, UP, 0, DN, 0, UP, 0, DN, 1, UP};
    run_case("fsz_last_permitted_destination", HOPPING_FSZ, 2, 2, loc, idx, cfg, spn, fsz_last, N(fsz_last));
    uint32_t csz_final[] = {1, 0, 0, 0, 0, 0, 1};
    run_case("csz_final_empty_destination_rejected", HOPPING_CSZ, 2, 2, loc, idx, cfg, spn, csz_final, N(csz_final));
    uint32_t csz_last[] = {1, 0, 0, 0, 0, 1};
    run_case("csz_last_permitted_destination", HOPPING_CSZ, 2, 2, loc, idx, cfg, spn, csz_last, N(csz_last));
  }
  // FSZ conduction spin flip (audit): electron 0 sits on a local-spin site,
  // electron 1 is itinerant with an empty opposite-spin slot (cap 4).
  {
    int loc[] = {1, 0}, idx[] = {0, 1}, cfg[] = {0, 1, -1, -1}, spn[] = {0, 0};
    uint32_t final_valid[] = {0, 0, 0, 0, 0, 1};
    run_case("conduction_final_valid_draw_rejected", CONDUCTION, 2, 2, loc, idx, cfg, spn, final_valid, N(final_valid));
    uint32_t last[] = {0, 0, 0, 0, 1};
    run_case("conduction_last_permitted_draw", CONDUCTION, 2, 2, loc, idx, cfg, spn, last, N(last));
  }
  return 0;
}
