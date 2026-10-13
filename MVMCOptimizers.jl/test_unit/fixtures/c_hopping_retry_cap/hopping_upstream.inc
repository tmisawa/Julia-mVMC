/* Verbatim: extern/mVMC-1.3.0/src/mVMC/vmcmake.c lines 515-545 (GPLv3, see that file). */
void makeCandidate_hopping(int *mi_, int *ri_, int *rj_, int *s_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg) {
  const int icnt_max = Nsite*Nsite;
  int icnt;
  int mi, ri, rj, s, flag;

  flag = 0; // FALSE
  do {
    mi = gen_rand32()%Ne;
    s = (genrand_real2()<0.5) ? 0 : 1;
    ri = eleIdx[mi+s*Ne];
  } while (LocSpn[ri] == 1);

  icnt = 0;
  do {
    rj = gen_rand32()%Nsite;
    if(icnt> icnt_max){
      flag = 1; // TRUE
      break;
    }
    icnt+=1;
  } while (eleCfg[rj+s*Nsite] != -1 || LocSpn[rj]==1);

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_ = s;
  *rejectFlag_ = flag;

  return;
}
