/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "../include/lib_poisson1D.h"
#include <time.h>

#define TRF 0
#define TRI 1
#define SV 2

int main(int argc, char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info;
  int NRHS;

  int IMPLEM;

  double T0, T1;
  double *RHS, *EX_SOL, *X;

  double *AB;

  double temp, relres;

  clock_t start_timer, end_timer;
  double elapsed_time;

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  if (argc == 2) 
  {
    IMPLEM = atoi(argv[1]);
  } 
  else if (argc > 2) 
  {
    perror("Application takes at most one argument");
    exit(1);
  }

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1); 
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "./RHS.dat");
  write_vec(EX_SOL, &la, "./EX_SOL.dat");
  write_vec(X, &la, "./X_grid.dat");


  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "./AB.dat");


  printf("Solution with LAPACK\n");

  start_timer = clock();

  /* LU Factorization */
  
  ipiv = (int *) calloc(la, sizeof(int));
  if (ipiv == NULL)
  {
    printf("\nFailed to allocate memory for ipiv\n");
    exit(1);
  }

  if (IMPLEM == TRF)
  {
    info = LAPACKE_dgbtrf(LAPACK_COL_MAJOR, la, la, kl, ku, AB, lab, ipiv);
    printf("TRF\n");
  }
  
  /* LU for tridiagonal matrix  (can replace LAPACKE_dgbtrf) */
  if (IMPLEM == TRI)
  {
    ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    printf("TRI\n");
    printf("ierr = %d\n", ierr);
  }
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "./LU.dat");


  /* Solution (Triangular) */
  if (IMPLEM == TRF || IMPLEM == TRI)
  {
    int ldb_dgbtrs = la;
    info = LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N', la, kl, ku, NRHS, AB, lab, ipiv, RHS, ldb_dgbtrs);
    printf("with LU factorization\n");
    // solution is stored in RHS
  }
  
   

  /* It can also be solved with dgbsv */
  // can jump the LU factorization step
  if (IMPLEM == SV)
  {
    int ldb_dgbsv = la;
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku, NRHS, AB, lab, ipiv, RHS, ldb_dgbsv);
    printf("SV without LU facorization\n");
    // solution is stored in RHS
  }
  

  write_xy(RHS, X, &la, "./SOL.dat"); // col1 x_grid, col2 solution

  /* Relative forward error */ 
  relres = relative_forward_error(RHS, EX_SOL, &la);
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(ipiv);

  end_timer = clock();
  elapsed_time = (double)(end_timer - start_timer) / (double)CLOCKS_PER_SEC;
  printf("Elapsed time for direct method: %lf seconds\n", elapsed_time);

  printf("\n\n--------- End -----------\n");
  return 0;
}
