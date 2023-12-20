/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la)
{
  for (int i = 0; i < *la; i ++)
    eigval[i] = 2.0 * (1 +cos(i * M_PI / (*la + 1)));
}

double eigmax_poisson1D(int *la)
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

double eigmin_poisson1D(int *la)
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

double richardson_alpha_opt(int *la)
{ 
  double alpha_opt = 1.0 / (eigmax_poisson1D(*la) + eigmin_poisson1D(*la));
  return alpha_opt;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) 
{
    // Allocate memory for temporary vectors
    double *temp1 = (double *)malloc(sizeof(double) * (*la));
    double *temp2 = (double *)malloc(sizeof(double) * (*la));

    // Initialize variables
    double alpha = *alpha_rich;
    int iter;

    // Initial residual: r_0 = b - Ax_0
    //void cblas_dcopy(const int N, const double *X, const int incX, double *Y, const int incY);
    cblas_dcopy(*la, RHS, 1, temp1, 1);  // temp1 = RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp1, 1);
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
    free(temp2);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

