/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv) 
{
    int i, j;
    int ku = *kv;  // Nombre de diagonales au-dessus de la diagonale principale
    int kl = *kv;  // Nombre de diagonales en dessous de la diagonale principale

    // Initialisation de la matrice à zéro
    for (j = 0; j < *lab; j ++) 
    {
        for (i = 0; i < *la; i ++) 
        {
            AB[i * (*lab) + j] = 0.0;
        }
    }

    // Remplissage de la bande principale
    for (i = 0; i < *la; i ++) 
    {
        AB[i * (*lab) + ku + kl] = 2.0;  // Éléments sur la diagonale principale
    }

    // Remplissage des diagonales au-dessus et en dessous de la diagonale principale
    for (i = 0; i < *la - 1; i ++) 
    {
        AB[(i + 1) * (*lab) + ku + kl - 1] = -1.0;  // Diagonale au-dessus
        AB[i * (*lab) + ku + kl + 1] = -1.0;         // Diagonale en dessous
    }
}


void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv) 
{
    int i, j;
    int ku = *kv;  // Nombre de diagonales au-dessus de la diagonale principale
    int kl = *kv;  // Nombre de diagonales en dessous de la diagonale principale

    // Initialisation de la matrice à zéro
    for (j = 0; j < *lab; j ++) 
    {
        for (i = 0; i < *la; i ++) 
        {
            AB[i * (*lab) + j] = 0.0;
        }
    }

    // Remplissage de la diagonale principale
    for (i = 0; i < *la; i ++) 
    {
        AB[i * (*lab) + ku + kl] = 1.0;  // Éléments sur la diagonale principale
    }
}


void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1)
{
    int i;
    for (i = 1; i < *la; i++) 
    {
        RHS[i] = 0.0;
    }
    RHS[0] = *BC0;
    RHS[*la - 1] = *BC1;
}


void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1) 
{
    int i;
    double T0 = *BC0;
    double T1 = *BC1;
    for (i = 0; i < *la; i++) 
    {
        EX_SOL[i] = X[i] * (T1 - T0) + T0;
    }
    EX_SOL[0] = T0;
    EX_SOL[*la - 1] = T1;
} 

void set_grid_points_1D(double* x, int* la) 
{
    int i;
    
    // Initialisation du vecteur à zéro
    for (i = 0; i < *la; i ++) 
    {
        x[i] = 0.0;
    }

    // Remplissage du vecteur avec des valeurs équidistantes
    for (i = 0; i < *la; i ++) 
    {
        x[i] = (double) (i + 1) / (double) ((*la) + 1);
    }
}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename)
{
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL)
  {
    for (ii=0;ii<(*lab);ii++)
    {
      for (jj=0;jj<(*la);jj++)
      {
	      fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else
  {
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename)
{
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL)
  {
    for (ii=0;ii<(*la);ii++)
    {
      for (jj=0;jj<(*lab);jj++)
      {
	      fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else
  {
    perror(filename);
  }
}

// output: row_index column_index value
void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename)
{
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL)
  {
    for (jj=1;jj<(*la);jj++)
    {
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++)
    {
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++)
    {
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else
  {
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename)
{
  int jj = 0; 
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL)
  {
    for (jj=0;jj<(*la);jj++)
    {
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else
  {
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename)
{
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL)
  {
    for (jj=0;jj<(*la);jj++)
    {
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else
  {
    perror(filename);
  } 
}  

int indexABCol(int i, int j, int *lab)
{
  return j * (*lab) + i;
}


// Function to perform LU factorization of a tridiagonal matrix
// Input:
//   la: Number of diagonals in the matrix (should be 3 for tridiagonal)
//   n: Order of the matrix
//   kl: Number of subdiagonals (should be 1 for tridiagonal)
//   ku: Number of superdiagonals (should be 1 for tridiagonal)
//   AB: Input matrix (stored in column-major order)
//   lab: Leading dimension of the matrix (lab >= n)
// Output:
//   AB: LU factorization of the matrix (overwritten on the input matrix)
//   ipiv: Array of pivot indices
//   info: Status indicator
int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info) 
{
    int i, j, k;
    double temp;
    // Check for invalid input
    if (*la != 3 || *kl != 1 || *ku != 1 || *lab < *n) 
    {
        *info = -1;
        return *info;
    }
    // Initialize info
    *info = 0;
    // Perform LU factorization
    for (k = 0; k < *n - 1; k ++) 
    {
        if (AB[k + (*ku) * (*lab) + k] == 0.0) 
        {
            *info = k + 1;  // Matrix is singular
            return *info;
        }
        // Compute multiplier
        temp = -AB[k + (*ku) * (*lab) + k + 1] / AB[k + (*ku) * (*lab) + k];
        
        // Store multiplier in the subdiagonal
        AB[k + 1 + (*kl) * (*lab) + k] = -temp;
        
        // Apply the multiplier to the remaining submatrix
        for (i = k + 1; i < k + 3 && i < *n; i ++) 
        {
            for (j = k + 1; j < k + 3 && j < *n; j ++) 
            {
                AB[i + (*ku) * (*lab) + j] += temp * AB[k + (*ku) * (*lab) + j];
            }
        }
    }
    // Check the last pivot element for singularity
    if (AB[*n - 1 + (*ku) * (*lab) + *n - 1] == 0.0) 
    {
        *info = *n;
    }
    return *info;
}


double relative_forward_error(double* x, double* y, int* la)
{
  int i;
  double temp1 = 0.0;
  double temp2 = 0.0;
  for (i = 0; i < *la; i ++)
  {
    temp1 += pow(x[i] - y[i], 2);
    temp2 += pow(y[i], 2);
  }
  return sqrt(temp1 / temp2);
}




// ********* replace lib_poisson1D_richardson.c ********

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
  double alpha_opt = 1.0 / (eigmax_poisson1D(eigval, la) + eigmin_poisson1D(eigval, la));
  return alpha_opt;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) 
{
    // Allocate memory for temporary vectors
    double *temp = (double *)malloc(sizeof(double) * (*la));
    if (temp == NULL) 
    {
    // Handle memory allocation failure
    fprintf(stderr, "Memory allocation failed for temp\n");
    exit(EXIT_FAILURE);
    }

    // Initialize variables
    double alpha = *alpha_rich;
    int iter;
    double norm_res;
  
    // Initial residual: r_0 = b - Ax_0
    cblas_dcopy(*la, RHS, 1, temp, 1);  // temp = RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp, 1); //temp = temp - AB * X
    norm_res = cblas_dnrm2(*la, temp, 1);
    resvec[0] = norm_res;
    *nbite = 0;

    // Iterate until convergence or max iterations reached
    for (iter = 1; iter < *maxit && norm_res > *tol; ++ iter) 
    {
        cblas_daxpy(*la, alpha, temp, 1, X, 1); //X = X + alpha * temp
        cblas_dcopy(*la, RHS, 1, temp, 1);  // temp = RHS
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp, 1); //temp = temp - AB * X

        // Save the residual norm for each iteration
        norm_res = cblas_dnrm2(*la, temp, 1);
        resvec[iter] = norm_res;
        *nbite = iter;
    }

    // Free allocated memory
    free(temp);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, int *kv, double *tol, int *maxit, double *resvec, int *nbite)
{
  double *temp = (double *)malloc(sizeof(double) * (*la));
  if (temp == NULL) 
  {
  // Handle memory allocation failure
  fprintf(stderr, "Memory allocation failed for temp\n");
  exit(EXIT_FAILURE);
  }

  // inverse of diagonal elements of AB
  for (int i = 0; i < *la; i ++)
  {
    MB[i * (*lab) + *ku + *kl] = 1.0 / AB[i * (*lab) + *ku + *kl];
  }

  int iter;
  double norm_res;

  // Initial residual: r_0 = b - Ax_0
  cblas_dcopy(*la, RHS, 1, temp, 1);
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp, 1);
  norm_res = cblas_dnrm2(*la, temp, 1);
  resvec[0] = norm_res;
  *nbite = 0;

  // Iterate until convergence or max iterations reached
  for (iter = 1; iter < *maxit && norm_res > *tol; ++ iter) 
  {
    cblas_dgemv(CblasColMajor, CblasNoTrans, *la, *la, 1.0, MB, *la, temp, 1, 1.0, X, 1); //X = X + MB * temp
    cblas_dcopy(*la, RHS, 1, temp, 1);  // temp = RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp, 1); //temp = temp - AB * X
    
    norm_res = cblas_dnrm2(*la, temp, 1);
    resvec[iter] = norm_res;
    *nbite = iter;
  }
  free(temp);
}

// without inversing (D - E)
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, int *kv, double *tol, int *maxit, double *resvec, int *nbite)
{
  double *temp = (double *)malloc(sizeof(double) * (*la));

  // MB = D - E
  for (int i = 0; i < *la; i ++)
  {
    MB[i * (*lab) + *ku + *kl] = AB[i * (*lab) + *ku + *kl]; // Diagonal elements
    if (i < *la - 1)
    {
      MB[i * (*lab) + *ku + *kl + 1] += AB[i * (*lab) + *ku + *kl + 1]; // Sub-diagonal elements
    }
  }
  
  int iter;
  double norm_res;

  cblas_dcopy(*la, RHS, 1, temp, 1);
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp, 1);
  norm_res = cblas_dnrm2(*la, temp, 1);
  resvec[0] = norm_res;
  *nbite = 0;

  for (iter = 1; iter < *maxit && norm_res > *tol; ++ iter)
  {
    for (int i = 0; i < *la; i++)
    {
      X[i] = (RHS[i] - temp[i]) / MB[i];
    }
    cblas_dcopy(*la, RHS, 1, temp, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp, 1);

    norm_res = cblas_dnrm2(*la, temp, 1);
    resvec[iter] = norm_res;
    *nbite = iter;
  }
  free(temp);
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