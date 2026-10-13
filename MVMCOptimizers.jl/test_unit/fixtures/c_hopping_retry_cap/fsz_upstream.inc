/* Verbatim: extern/mVMC-1.3.0/src/mVMC/vmcmake_fsz.c lines 594-664 (GPLv3, see that file). */
void makeCandidate_hopping_fsz(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn) {
  const int icnt_max = Nsite*Nsite;
  int icnt;
  int mi, ri, rj, s, flag;
  int t; //fsz

  flag = 0; // FALSE
  do {
    mi = gen_rand32()%Nsize;
    s  = eleSpn[mi] ; //fsz 
    //t  = (genrand_real2()<0.5) ? s : 1-s; //fsz
    ri = eleIdx[mi];  //fsz
  } while (LocSpn[ri] == 1);

  icnt = 0;
  do {
    rj = gen_rand32()%Nsite;
    t  = (genrand_real2()<0.5) ? 0 : 1; //fsz
    if(icnt> icnt_max){
      flag = 1; // TRUE
      break;
    }
    icnt+=1;
  } while (eleCfg[rj+t*Nsite] != -1 || LocSpn[rj]==1);

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_  = s;
  *t_  = t;
  *rejectFlag_ = flag;

  return;
}
// mi (ri,s) -> mi (rj,s)
void makeCandidate_hopping_csz(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn) {
  const int icnt_max = Nsite*Nsite;
  int icnt;
  int mi, ri, rj, s, flag;
  int t; //fsz

  flag = 0; // FALSE
  do {
    mi = gen_rand32()%Nsize;
    s  = eleSpn[mi] ; //fsz 
    //t  = (genrand_real2()<0.5) ? s : 1-s; //fsz
    t  = s;//csz
    ri = eleIdx[mi];  //fsz
  } while (LocSpn[ri] == 1);

  icnt = 0;
  do {
    rj = gen_rand32()%Nsite;
    if(icnt> icnt_max){
      flag = 1; // TRUE
      break;
    }
    icnt+=1;
  } while (eleCfg[rj+t*Nsite] != -1 || LocSpn[rj]==1);

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_  = s;
  *t_  = t;
  *rejectFlag_ = flag;

  return;
}
/* Verbatim: extern/mVMC-1.3.0/src/mVMC/vmcmake_fsz.c lines 756-788 (GPLv3, see that file). */
void makeCandidate_LocalSpinFlip_conduction(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn) {
  const int icnt_max = Nsite*Nsite;
  int icnt=0;
  int mi, ri, rj, s, flag;
  int t; //fsz

  flag = 0; // FALSE
  do {
    mi = gen_rand32()%Nsize;
    s  = eleSpn[mi] ; //fsz 
    t  = 1-s;
    //t  = (genrand_real2()<0.5) ? s : 1-s; //fsz
    ri = eleIdx[mi];  //fsz
    rj = ri;  //fsz // note ! we assume local spin
    if(icnt> icnt_max){ // all doublons can not be accepted
      flag = 1; // TRUE
      break;
    }
    icnt+=1;
  } while (LocSpn[ri] == 1 || eleCfg[ri+t*Nsite] != -1);

  //if(t==s) flag=1;

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_  = s;
  *t_  = t;
  *rejectFlag_ = flag;

  return;
}
