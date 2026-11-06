/* Unchanged mVMC 1.4.0 StochasticOptCG_Main definition (stcopt_cg_impl.c:484-591).
 * Copyright (C) 2016 University of Tokyo; GPL-3.0-or-later.
 * Origin https://github.com/issp-center-dev/mVMC ; see c_cg_small_spd.md.
 */
#include <stdio.h>
#include <math.h>
typedef int MPI_Comm;
#define USE_IMAG 0
int NSROptCGMaxIter=1, NVMCSample=1;
double DSROptCGTol=1e-30;
void MPI_Comm_rank(MPI_Comm c,int*r){(void)c;*r=0;}
void MPI_Comm_size(MPI_Comm c,int*s){(void)c;*s=1;}
double xdot(const int n,const double* p,const double* q){double z=0;for(int i=0;i<n;i++)z+=p[i]*q[i];return z;}
/* Adapter is the one-dimensional sampled operator: O*(O*d), W=1,
 * mean O=0, stabilization=0. No BLAS/MPI executable claim. */
void fn_operate_by_S(int n,const double*d,double*q,double*v,MPI_Comm c){(void)n;(void)c; q[0]=v[4]*(v[4]*d[0]);}
int fn_StochasticOptCG_Main(const int nSmat, double *VecCG, MPI_Comm comm) {
  int si;
  int rank, size;
  int iter;
  int max_iter = (NSROptCGMaxIter > 0 ? NSROptCGMaxIter : nSmat);
  double delta, beta;
  double alpha;
  double delta_new, dq;
  double cg_thresh = DSROptCGTol*DSROptCGTol * (double)nSmat * (double)nSmat;
  //double cg_thresh = DSROptRedCut;

  double *x, *g, *sdiag, *stcO, *stcOs_real;
  double *stcOs_imag, *y_real, *y_imag, *z_local;
  double *q, *d, *r;

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): Start stcOptCG_Main\n", __FILE__, __LINE__);
#endif
  
  x = VecCG;
  g = x + nSmat;
  sdiag = g + nSmat;
  stcO = sdiag + nSmat;
  stcOs_real = stcO + nSmat;
  stcOs_imag = stcOs_real + NVMCSample*nSmat;
  y_real = stcOs_imag + USE_IMAG*NVMCSample*nSmat;
  y_imag = y_real + NVMCSample;
  z_local = y_imag + USE_IMAG*NVMCSample;
  q = z_local + nSmat;
  d = q + nSmat;
  r = d + nSmat;
  
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    d[si] = r[si] = g[si];
  }

  delta = xdot(nSmat, r, r);
  if (!isfinite(delta)) return -(0+1);

  for(iter=0; iter < max_iter; iter++){
    //check convergence 
    //if(rank==0) printf("%d: %d %.3e\n",rank, iter,delta);
#ifdef _DEBUG_STCOPT_CG
    fprintf(stderr, "delta = %lg, cg_thresh = %lg\n", delta, cg_thresh);
#endif
    if (!isfinite(delta)) return -(iter+1);
    if (fabs(delta) < cg_thresh) break;

    // compute vector q=S*d
    fn_operate_by_S(nSmat, d, q, VecCG, comm);
    dq = xdot(nSmat,d,q);
    if (!isfinite(dq) || fabs(dq) < 1.0e-20) return -(iter+1);
    alpha = delta/dq;
    if (!isfinite(alpha)) return -(iter+1);
  
    // update solution vector x=x+alpha*d
    #pragma omp parallel for default(shared) private(si)
    #pragma loop noalias
    for(si=0;si<nSmat;++si) {
      x[si] = x[si] + alpha*d[si];
    }
    // update residual vector r=r-alpha*q, q=S*d
    if((iter+1) % 20 == 0){
      fn_operate_by_S(nSmat, x, r, VecCG, comm);
      #pragma omp parallel for default(shared) private(si)
      #pragma loop noalias
      for(si=0;si<nSmat;++si) {
        r[si] = g[si] - r[si];
      }
    }else{
      #pragma omp parallel for default(shared) private(si)
      #pragma loop noalias
      for(si=0;si<nSmat;++si) {
        r[si] = r[si] - alpha*q[si];
      }
    }
    delta_new = xdot(nSmat,r,r);
    if (!isfinite(delta_new)) return -(iter+1);
    if (fabs(delta) < 1.0e-20) return -(iter+1);
    beta = delta_new/delta;
    if (!isfinite(beta)) return -(iter+1);

    //update the norm of residual vector r
    /* Use delta=beta*delta (mathematically delta_new) to preserve the legacy
       FP order for SR-CG reference data; the DiagScale branch uses
       delta=delta_new directly. */
    delta = beta*delta;
    // update direction vector d
    #pragma omp parallel for default(shared) private(si)
    #pragma loop noalias
    for(si=0;si<nSmat;++si) {
      d[si] = r[si] + beta*d[si];
    }
  }
  if (iter == max_iter && fabs(delta) >= cg_thresh) return -(max_iter+1);

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): iter = %d\n", __FILE__, __LINE__, iter);
  fprintf(stderr, "DEBUG in %s (%d): End stcOptCG_Main\n", __FILE__, __LINE__);
#endif

  return iter;
}

int main(void){double v[10]={0}; v[1]=ldexp(1.0,-54); v[4]=ldexp(1.0,-54); int it=fn_StochasticOptCG_Main(1,v,0); printf("%d %.17g %.17g %.17g\n",it,v[0],v[9],v[8]);return 0;}
