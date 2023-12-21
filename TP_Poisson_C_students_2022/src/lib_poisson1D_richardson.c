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

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite) {
  
  double* rk = malloc(sizeof(double) * (*la));
  double norm = 0.0;
  double norm_B = cblas_dnrm2(*la, RHS, 1); // L2 norm
  double inv_norm = 1 / norm_B;

  for ((*nbite) = 0; (*nbite) < (*maxit); (*nbite)++) {
    
    for (int i = 0; i < (*la); i++) {
      rk[i] = RHS[i];
    }


    // Calcul de b = b-Ax
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, rk, 1);
    // Résidu
    norm = cblas_dnrm2(*la, rk, 1) * inv_norm;
    resvec[(*nbite)] = norm;

    cblas_daxpy(*la, *alpha_rich, rk, 1, X, 1);

    if (norm <= (*tol))
      break;
  }

  free(rk);
  
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv) {
  for (int i = 0; i < *la; i++) {
    for (int j = 0; j < *lab; j++) {
      int ind = (*lab)*i;
      MB[ind + j] = 0.0;
    }
    // Fill the diagonal
    int ind = ((*lab)*i) + (*ku);
    MB[ind] = AB[ind];
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv) {
  for (int i = 0; i < *la; i++) {
    for (int j = 0; j < *lab; j++) {
      int ind = (*lab)*i;
      MB[ind + j] = 0.0;
    }
    // Fill the diagonal
    int ind = ((*lab)*i) + (*ku);
    MB[ind] = AB[ind];
    MB[ind + 1] = AB[ind + 1];
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite) {
  double* rk = malloc(sizeof(double) * (*la));
  double norm_B = cblas_dnrm2(*la, RHS, 1); // L2 norm
  double inv_norm = 1 / norm_B;
  double norm = 0.0;
  int* ipiv = malloc(sizeof(int) * (*la)); // Pivots
  int info = 0;// 0->success
  int NHRS = 1; // Right-hand sides
  int kuu = (*ku)-1;



  dgbtrf_(la, la, kl, &kuu, MB, lab, ipiv, &info);
  for ((*nbite) = 0; (*nbite) < (*maxit); (*nbite)++) {

    for (int i = 0; i < (*la); i++) {
      rk[i] = RHS[i];
    }

    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, rk, 1);

    norm = cblas_dnrm2(*la, rk, 1) * inv_norm;
    resvec[(*nbite)] = norm;

    dgbtrs_("N", la, kl, &kuu, &NHRS, MB, lab, ipiv, rk, la, &info, 1);

    cblas_daxpy(*la, 1, rk, 1, X, 1);

    if (norm <= (*tol))
      break;
  }

  free(rk);
  free(ipiv);

}

