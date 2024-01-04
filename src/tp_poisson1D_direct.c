/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "../include/lib_poisson1D.h"

int main(int argc, char *argv[])
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

  set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
 
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
  // result here: Y = A * EX_SOL, will be used as RHS for the solver
  cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1.0, AB, lab, EX_SOL, NRHS, 1.0, Y, NRHS);
  write_vec(Y, &la, "Y.dat");
  //printf("Y\n");


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

  info = LAPACKE_dgbtrf(LAPACK_COL_MAJOR, la, la, kl, ku, AB, lab, ipiv);
  

  /* LU for tridiagonal matrix  (can replace dgbtrf_) 
  ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  */

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  //printf("LU\n");

  /* Solution (Triangular) */
  
  int ldb_dgbtrs = la;
  info = LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N', la, kl, ku, NRHS, AB, lab, ipiv, RHS, ldb_dgbtrs);
  // solution is stored in RHS
  

  /* It can also be solved with dgbsv */
  // can jump the LU factorization step
  /*
  int ldb_dgbsv = la;
  info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku, NRHS, AB, lab, ipiv, RHS, ldb_dgbsv);
  // solution is stored in RHS
  */
  write_xy(RHS, X, &la, "SOL.dat"); // col1 x, col2 rhs

  /* Relative forward error */
  // Compute relative norm of the residual
  
  relres = relative_forward_error(RHS, Y, &NRHS);
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
