/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

his program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. 

You should have received a copy of the GNU General Public License 
along with this program. If not, see http://www.gnu.org/licenses/. 
*/
/* calculate the parameter change r[nSmat] from SOpt.
   Solve S*x = g */
int fn_StochasticOptPreconCG_DiagScale_Main(const int nSmat, double *VecCG, MPI_Comm comm) {
  int si;
  int rank, size;
  int iter;
  int max_iter = (NSROptCGMaxIter > 0 ? NSROptCGMaxIter : nSmat);
  double delta, beta;
  double alpha;
  double delta_new, dq;
  double cg_thresh = DSROptCGTol*DSROptCGTol * (double)nSmat * (double)nSmat;
  const double precon_factor = 1.0 + DSROptStaDel;
  const double precon_diag_eps = 1.0e-12;
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

  fn_operate_by_S(nSmat, x, r, VecCG, comm);
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    r[si] = g[si]-r[si];
  }

  #pragma omp parallel for default(shared) private(si,dq)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    dq = precon_factor * sdiag[si];
    if (!isfinite(dq) || fabs(dq) < precon_diag_eps) {
      dq = (dq < 0.0 ? -precon_diag_eps : precon_diag_eps);
    }
    d[si] = r[si]/dq;
  }

  delta = xdot(nSmat, r, d);
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

    #pragma omp parallel for default(shared) private(si,dq)
    #pragma loop noalias
    for(si=0;si<nSmat;++si) {
      dq = precon_factor * sdiag[si];
      if (!isfinite(dq) || fabs(dq) < precon_diag_eps) {
        dq = (dq < 0.0 ? -precon_diag_eps : precon_diag_eps);
      }
      q[si] = r[si]/dq;
    }
    delta_new = xdot(nSmat,r,q);
    if (!isfinite(delta_new)) return -(iter+1);
    if (fabs(delta) < 1.0e-20) return -(iter+1);
    beta = delta_new/delta;
    if (!isfinite(beta)) return -(iter+1);

    //update the norm of residual vector r
    delta = delta_new;
    // update direction vector d
    #pragma omp parallel for default(shared) private(si)
    #pragma loop noalias
    for(si=0;si<nSmat;++si) {
      //d[si] = r[si] + beta*d[si];
      d[si] = q[si] + beta*d[si];
    }
  }
  if (iter == max_iter && fabs(delta) >= cg_thresh) return -(max_iter+1);

#ifdef _DEBUG_STCOPT_CG
  fprintf(stderr, "DEBUG in %s (%d): iter = %d\n", __FILE__, __LINE__, iter);
  fprintf(stderr, "DEBUG in %s (%d): End stcOptCG_Main\n", __FILE__, __LINE__);
#endif

  return iter;
}

/* calculate the parameter change r[nSmat] from SOpt.
   Solve S*x = g */
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



/* calculate  z = S*x */
/* S is the overlap matrix*/
/* S[i][j] = OO[i+1][j+1] - OO[i+1][0] * OO[0][j+1]; */
int fn_operate_by_S(int nSmat, double *x, double *z, double *VecCG, MPI_Comm comm) {
  int rank,size,info=0;
  int si;
  int incx = 1, incy = 1;
  double coef;
  double one = 1.0, zero = 0.0;
  double invW = 1.0/Wc;
  char transT='T';
  char transN='N';

  double *sdiag, *stcO, *stcOs_real, *stcOs_imag;
  double *y_real, *y_imag;
  double *z_local;

  sdiag = VecCG + 2*nSmat;
  stcO = sdiag + nSmat;
  stcOs_real = stcO + nSmat;
  stcOs_imag = stcOs_real + NVMCSample*nSmat;
  y_real = stcOs_imag + USE_IMAG*NVMCSample*nSmat;
  y_imag = y_real + NVMCSample;
  z_local = y_imag + USE_IMAG*NVMCSample;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  MPI_Bcast(x, nSmat, MPI_DOUBLE, 0, comm);
  
  StartTimer(53);
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    z[si] = 0.0;
    z_local[si] = 0.0;
  }

  // y_real[sample] = sum{si} x[si] * O_real[si][sample]
  M_DGEMV(&transT, &nSmat, &NVMCSample, &one, stcOs_real, &nSmat, x, &incx, &zero, y_real, &incy);
#ifndef MVMC_SRCG_REAL
  // y_imag[sample] = sum{si} x[si] * O_imag[si][sample]
  M_DGEMV(&transT, &nSmat, &NVMCSample, &one, stcOs_imag, &nSmat, x, &incx, &zero, y_imag, &incy);
#endif

  // z_local[si] = sum{sample} O_real[si][sample] * y_real[sample] + O_imag[si][sample] * y_imag[sample]
  M_DGEMV(&transN, &nSmat, &NVMCSample, &one, stcOs_real, &nSmat, y_real, &incx, &zero, z_local, &incy);

#ifndef MVMC_SRCG_REAL
  M_DGEMV(&transN, &nSmat, &NVMCSample, &one, stcOs_imag, &nSmat, y_imag, &incx, &one, z_local, &incy);
#endif

  MPI_Barrier(comm);
  /* compute <OO>*x */
  SafeMpiAllReduce(z_local, z, nSmat, comm);
 
  /* compute <O>*x */
  coef = xdot(nSmat, stcO, x);

  /* y = S*x */
  #pragma omp parallel for default(shared) private(si)
  #pragma loop noalias
  for(si=0;si<nSmat;++si) {
    z[si] = invW*z[si] - coef*stcO[si];
    /* modify diagonal elements */
    // z[si] += DSROptStaDel * x[si];

    z[si] += sdiag[si]*DSROptStaDel*x[si];
  }

  StopTimer(53);

  return info;
}
