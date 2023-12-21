/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "../include/lib_poisson1D.h"
#include "../include/atlas_headers.h"

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
  double *RHS, *EX_SOL, *X, *Y;
  double **AAB;
  double *AB;

  double temp, relres;

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  
  set_grid_points_1D(X, &la);

  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  //printf("RHS\n");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  //printf("EX_SOL\n");
  write_vec(X, &la, "X_grid.dat");
  //printf("X\n");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
  //printf("AB\n");
  
  //void cblas_dgbmv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const int KL, const int KU, const double alpha, const double *A, const int lda, const double *X, const int incX, const double beta, double *Y, const int incY);
  
  cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1.0, AB, lab, EX_SOL, NRHS, 1.0, Y, NRHS);
  write_vec(Y, &la, "Y.dat");
  //printf("Y\n");

  /* Verify validity */
  temp = 0.0;
  for (int i = 0; i < la; i ++)
  {
    Y[i] = pow(Y[i] - RHS[i], 2);
    temp += Y[i];
  }
  temp = sqrt(temp);

  printf("Solution with LAPACK\n");


  /* LU Factorization */
  info=0;
  ipiv = (int *) calloc(la, sizeof(int));
  if (ipiv == NULL)
  {
    printf("\nFailed to allocate memory for ipiv\n");
    exit(1);
  }
  //printf("ipiv\n");
  if (info == 0)
  {
    LAPACK_dgbtrf(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    if (info!=0)
    {
      printf("\n INFO DGBTRF = %d\n",info);
    }
  }
  else
  {
    printf("\n DGBTRF: INFO = %d\n",info);
  }
  

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  //printf("LU\n");

  /* Solution (Triangular) */
  if (info == 0)
  {
    //void LAPACK_dgbtrs(char* trans, lapack_int* n, lapack_int* kl, lapack_int* ku, lapack_int* nrhs, const double* ab, lapack_int* ldab, const lapack_int* ipiv, double* b, lapack_int* ldb, lapack_int *info );
    int ldb_dgbtrs = la;
    LAPACK_dgbtrs('N', &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &ldb_dgbtrs, &info);
    if (info!=0)
    {
      printf("\n INFO DGBTRS = %d\n",info);
    }
  }
  else
  {
    printf("\n DGBTRS: INFO = %d\n",info);
  }
  

  /* It can also be solved with dgbsv */
  if (info == 0)
  {
    //void LAPACK_dgbsv( lapack_int* n, lapack_int* kl, lapack_int* ku, lapack_int* nrhs, double* ab, lapack_int* ldab, lapack_int* ipiv, double* b, lapack_int* ldb, lapack_int *info );
    int ldb_dgbsv = la;
    LAPACK_dgbsv(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &ldb_dgbsv, &info);
    if (info != 0)
    {
      printf("\n INFO DGBSV = %d\n",info);
      }
  }
  else
  {
    printf("\n DGBSV: INFO = %d\n",info);
  }
  write_xy(RHS, X, &la, "SOL.dat"); // col1 x, col2 rhs

  /* Relative forward error */
  // TODO : Compute relative norm of the residual
  
  relres = temp / dnrm2_(&la, RHS, &NRHS);
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
