/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"


void eig_poisson1D(double* eigval, int *la)
{ 
  for (int i = 0; i < *la; i ++)
  {
    eigval[i] = 2.0 * (1 + cos(i * M_PI / (*la + 1)));
  }  
}

double eigmax_poisson1D(double* eigval, int *la)
{
  double max_eig = eigval[0];

  for (int i = 1; i < *la; i ++) 
  {
    if (eigval[i] >= max_eig) 
    {
      max_eig = eigval[i];
    }
  }
  return max_eig;
}

double eigmin_poisson1D(double* eigval, int *la)
{ 
  double min_eig = eigval[0];

  for (int i = 1; i < *la; i ++) 
  {
    if (eigval[i] <= min_eig) 
    {
      min_eig = eigval[i];
    }
  }
  return min_eig;
}

double richardson_alpha_opt(double* eigval, int *la)
{ 
  double alpha_opt = 1.0 / (eigmax_poisson1D(eigval, *la) + eigmin_poisson1D(eigval, *la));
  return alpha_opt;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) 
{
    // Allocate memory for temporary vectors
    double *temp1 = (double *)malloc(sizeof(double) * (*la));
    //double *temp2 = (double *)malloc(sizeof(double) * (*la));

    // Initialize variables
    double alpha = *alpha_rich;
    int iter;

    // Initial residual: r_0 = b - Ax_0
    //void cblas_dcopy(const int N, const double *X, const int incX, double *Y, const int incY);
    cblas_dcopy(*la, RHS, 1, temp1, 1);  // temp1 = RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp1, 1); //temp1 = temp1 - AB * X
    double norm_res = cblas_dnrm2(*la, temp1, 1);

    // Iterate until convergence or max iterations reached
    for (iter = 0; iter < *maxit && norm_res > *tol; ++ iter) 
    {
        // Richardson iteration: x_{k+1} = x_k + alpha * (b - Ax_k)
        cblas_dcopy(*la, RHS, 1, temp1, 1);  // temp1 = RHS
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp1, 1);
        cblas_daxpy(*la, alpha, temp1, 1, X, 1);

        // Update residual: r_{k+1} = b - Ax_{k+1}
        cblas_dcopy(*la, RHS, 1, temp1, 1);  // temp1 = RHS
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp1, 1);
        norm_res = cblas_dnrm2(*la, temp1, 1);

        // Save the residual norm for each iteration
        resvec[iter] = norm_res;
        nbite[0] = iter + 1;
    }

    // Free allocated memory
    free(temp1);
    //free(temp2);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, int *kv, double *tol, int *maxit, double *resvec, int *nbite)
{
  double *temp1 = (double *)malloc(sizeof(double) * (*la));
  double *temp2 = (double *)malloc(sizeof(double) * (*la));

  // inverse of diagonal elements of AB
  for (int i = 0; i < *la; i ++)
  {
    MB[i] = 1.0 / AB[i * (*lab) + *ku + *kl];
  }

  int iter;
  double norm_res;

  // Initial residual: r_0 = b - Ax_0
  cblas_dcopy(*la, RHS, 1, temp1, 1);
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp1, 1);
  norm_res = cblas_dnrm2(*la, temp1, 1);

  // Iterate until convergence or max iterations reached
  for (iter = 0; iter < *maxit && norm_res > *tol; ++ iter) 
  {
    cblas_dcopy(*la, RHS, 1, temp1, 1);
    cblas_dcopy(*la, X, 1, temp2, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp1, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, *la, *la, 1.0, MB, *la, temp1, 1, 0.0, temp2, 1); //temp2 = temp2 + MB * temp1

    cblas_dcopy(*la, RHS, 1, temp1, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, temp2, 1, 1.0, temp1, 1);
    norm_res = cblas_dnrm2(*la, temp1, 1);

    resvec[iter] = norm_res;
    nbite[0] = iter + 1;
  }
  free(temp1);
  free(temp2);
}

// without inversing (D - E)
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, int *kv, double *tol, int *maxit, double *resvec, int *nbite)
{
  double *temp1 = (double *)malloc(sizeof(double) * (*la));
  double *temp2 = (double *)malloc(sizeof(double) * (*la));

  // MB = D - E
  for (int i = 0; i < *la; i++)
  {
    MB[i] = AB[i * (*lab) + *ku + *kl]; // Diagonal elements
    if (i < *la - 1)
    {
      MB[i] += AB[i * (*lab) + *ku + *kl + 1]; // Sub-diagonal elements
    }
  }
  
  int iter;
  double norm_res;

  cblas_dcopy(*la, RHS, 1, temp1, 1);
  cblas_dcopy(*la, X, 1, temp2, 1);
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp1, 1);
  norm_res = cblas_dnrm2(*la, temp1, 1);

  for (iter = 0; iter < *maxit && norm_res > *tol; ++iter)
  {
    cblas_dcopy(*la, RHS, 1, temp1, 1);
    cblas_dcopy(*la, X, 1, temp2, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp1, 1);
        
    for (int i = 0; i < *la; i++)
    {
      X[i] = (RHS[i] - temp1[i]) / MB[i];
    }

    norm_res = cblas_dnrm2(*la, temp1, 1);

    resvec[iter] = norm_res;
    nbite[0] = iter + 1;
  }

  free(temp1);
  free(temp2);
}


void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite)
{
    // Allocate memory for temporary vectors
    double *temp1 = (double *)malloc(sizeof(double) * (*la));
    double *temp2 = (double *)malloc(sizeof(double) * (*la));

    // Initialize variables
    int iter;
    double norm_res;

    // Initial residual: r_0 = b - Ax_0
    cblas_dcopy(*la, RHS, 1, temp1, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp1, 1);
    norm_res = cblas_dnrm2(*la, temp1, 1);

    // Iterate until convergence or max iterations reached
    for (iter = 0; iter < *maxit && norm_res > *tol; ++iter)
    {
        // Richardson iteration: x_{k+1} = x_k + MB * (b - Ax_k)
        cblas_dcopy(*la, RHS, 1, temp1, 1);
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp1, 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, *la, *la, 1.0, MB, *la, temp1, 1, 0.0, temp2, 1); // temp2 = temp2 + MB * temp1
        cblas_daxpy(*la, 1.0, temp2, 1, X, 1);                                                // X = X + temp2

        // Update residual: r_{k+1} = b - Ax_{k+1}
        cblas_dcopy(*la, RHS, 1, temp1, 1);
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp1, 1);
        norm_res = cblas_dnrm2(*la, temp1, 1);

        // Save the residual norm for each iteration
        resvec[iter] = norm_res;
        nbite[0] = iter + 1;
    }

    // Free allocated memory
    free(temp1);
    free(temp2);
}