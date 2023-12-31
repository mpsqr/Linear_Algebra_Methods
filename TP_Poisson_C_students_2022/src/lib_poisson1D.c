/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){

  int k = 0;

  for (int i = 0; i < *kv; i++) {
    for (int j = 0; j <  *la; j++) {
      k = (j*(*lab));
      AB[i + k] = 0.0;
    }
  }

  for (int i = 0; i < *la; i++) {
    k = (i*(*lab));
    AB[*kv + k] = -1.0;
    AB[*kv + k + 1] = 2.0;
    AB[*kv + k + 2] = -1.0;
  }

  AB[*kv] = 0.0;
  AB[((*lab) * (*la)) - 1] = 0.0;

}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){

  for (int i = 0; i < *la; i++) {

    int k = i*(*lab);

    if (*kv >= 0) {
      for (int j = 0; j < *kv; j++) {
        AB[j+k] = 0.0;
      }
    }

    

    AB[*kv + k] = 0.0;
    AB[*kv + k + 1] = 1.0;
    AB[*kv + k + 2] = 0.0;

  }

  AB[1] = 0.0;
  AB[((*la) * (*lab)) - 1] = 0.0;

}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  RHS[0] = *BC0;
  
  for (int i = 1; i < *la-1; i++) {
    RHS[i] = 0.0;
  }
  
  RHS[*la-1] = *BC1;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){

  double diff = *BC1 - *BC0; // T1 - T0

  for (int i = 0; i < *la; i++) {
    EX_SOL[i] = *BC0 + (X[i]*diff);
  }

}  

void set_grid_points_1D(double* x, int* la){

  double h = 1.0/(*la+1.0); // Step

  for (int i = 0; i < *la; i++) {
    x[i] = h*(double)(i+1);
  }

}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

int indexABCol(int i, int j, int *lab) {
  return i + j*(*lab);
}


int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){

  ipiv[0]=1;

  for(int i = 1; i < *la; i++){
    if(!AB[ *lab*i-2 ]){
      *info = 1;
      return *info;
    }
    AB[(*lab*i)-1] *= 1/AB[*lab*i-2];
    AB[*lab*(i+1)-2] -= AB[*lab*i-1] * AB[*lab*(i+1)-3];
    ipiv[i]=i+1;
  }
  return *info;
}

double forward_error(int* la, double* EX_SOL, double* SOL) {
  
  double* temp = malloc(sizeof(double)*(*la));
  cblas_dcopy(*la, SOL, 1, temp, 1);

  double exsolution_norm = 0.0; // ||x-x'||
  // 2-norm
  for (int i = 0; i < *la; i++) {
    exsolution_norm += EX_SOL[i]*EX_SOL[i];
    temp[i] *= -1.0;
  }
  exsolution_norm = sqrt(exsolution_norm);


  // Daxpy
  for (int i = 0; i < *la; i++) {
    temp[i] += EX_SOL[i];
  }


  double solution_norm = 0.0; // ||x||
  // 2-norm
  for (int i = 0; i < *la; i++) {
    solution_norm += temp[i]*temp[i];
  }
  solution_norm = sqrt(solution_norm);

  free(temp);
  return solution_norm / exsolution_norm;
}