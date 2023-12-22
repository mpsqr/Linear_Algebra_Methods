/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *RHS, *MY_RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;

  //double temp, relres;

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  MY_RHS = (double *) malloc(sizeof(double)*la); // For validating our results
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  // TODO : add validation tests

  struct timespec start, end;
  double time = 0.0;

  // DGBMV
  printf("TEST WITH DGBMV\n");

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  //write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  cblas_dgbmv(CblasColMajor, CblasTrans, la, la, kl, ku, 1, AB+1, lab, EX_SOL, 1, 0, MY_RHS, 1);
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);
  time = ((double)end.tv_sec + (double)end.tv_nsec/1e9) - ((double) start.tv_sec + (double)start.tv_nsec/1e9);

  //write_vec(MY_RHS, &la, "MY_RHS.dat");
  // Validating
  for (int i = 0; i < la; i++) {
    printf("[RHS-MY_RHS]: %lf\n", RHS[i]-MY_RHS[i]);
  }

  printf("Time taken by DGBMV: %lfs \n", time);
  printf("Forward error is: %e\n", forward_error(&la, EX_SOL, MY_RHS));

  printf("\n\n");



  // DGBTRF
  printf("TEST WITH DGBTRF\n");

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  info=0;
  ipiv = (int *) calloc(la, sizeof(int));

  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);


  // Solution (Triangular) 
  if (info==0){
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info, 1);
    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
  }else{
    printf("\n INFO = %d\n",info);
  }

  clock_gettime(CLOCK_MONOTONIC_RAW, &end);
  time = ((double)end.tv_sec + (double)end.tv_nsec/1e9) - ((double) start.tv_sec + (double)start.tv_nsec/1e9);

  printf("Time taken by DGBTRF: %lfs \n", time);
  //printf("Forward error is: %e\n", forward_error(&la, EX_SOL, RHS));

  printf("\n\n");


  
  // DGBTRFTRIDIAG
  printf("TEST WITH DGBTRFTRIDIAG\n");


  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  info = 0;

  // LU for tridiagonal matrix  (can replace dgbtrf_) 
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

  if (info==0){
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info, 1);
    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
  } else{
    printf("\n INFO = %d\n",info);
  }

  clock_gettime(CLOCK_MONOTONIC_RAW, &end);
  time = ((double)end.tv_sec + (double)end.tv_nsec/1e9) - ((double) start.tv_sec + (double)start.tv_nsec/1e9);


  //write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");

  printf("Time taken by DGBTRFTRIDIAG: %lfs \n", time);
  //printf("Forward error is: %e\n", forward_error(&la, EX_SOL, RHS));

  printf("\n\n");
  

  // DGBSV
  printf("TEST WITH DGBSV\n");

  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  info = 0;
  
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);
  time = ((double)end.tv_sec + (double)end.tv_nsec/1e9) - ((double) start.tv_sec + (double)start.tv_nsec/1e9);

  write_xy(RHS, X, &la, "SOL.dat");
  
  printf("Time taken by DGBTRFTRIDIAG: %lfs \n", time);
  //printf("Forward error is: %e\n", forward_error(&la, EX_SOL, RHS));


  free(RHS);
  free(MY_RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
