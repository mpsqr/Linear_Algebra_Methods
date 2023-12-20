/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la) {
}

double eigmax_poisson1D(int *la) {
  double h = 1.0/((*la) + 1.0);
  double temp = sin((M_PI*h*(*la))*0.5);
  return 4.0*temp*temp;
}

double eigmin_poisson1D(int *la) {
  double h = 1.0/((*la) + 1.0);
  double temp = sin((M_PI*h)*0.5);
  return 4.0*temp*temp;
}

double richardson_alpha_opt(int *la) {
  return 2/(eigmax_poisson1D(la) + eigmin_poisson1D(la)); // = 1/2
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv) {
  for (int i = 0; i < *la; i++) {
    for (int j = 0; j < *lab; j++) {
      MB[((*lab)*i) + j] = 0.0;
    }
    // Fill the diagonal
    MB[((*lab)*i) + (*ku)] = AB[((*lab)*i) + (*ku)];
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

