/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv) 
{
    int i, j;

    // Initialisation
    for (j = 0; j < *lab; j ++) 
    {
        for (i = 0; i < *la; i ++) 
        {
            AB[i * (*lab) + j] = 0.0;
        }
    }

    // Main diagonal
    for (i = 0; i < *la; i ++) 
    {
        AB[i * (*lab) + (*kv) + 1] = 2.0; 
    }

    // Upper and lower diagonals
    for (i = 0; i < *la - 1; i ++) 
    {
        AB[(i + 1) * (*lab) + (*kv)] = -1.0;  // Upper diagonal
        AB[i * (*lab) + (*kv) + 2] = -1.0;    // Lower diagonal
    }
}


void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv) 
{
    int i, j;

    // Initialisation
    for (j = 0; j < *lab; j ++) 
    {
        for (i = 0; i < *la; i ++) 
        {
            AB[i * (*lab) + j] = 0.0;
        }
    }

    // Main diagonal
    for (i = 0; i < *la; i ++) 
    {
        AB[i * (*lab) + (*kv) + 1] = 1.0;  
    }
}


void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1)
{
    int i;
    RHS[0] = *BC0;
    for (i = 1; i < *la; i ++) 
    {
        RHS[i] = 0.0;
    }
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
} 

void set_grid_points_1D(double* x, int* la) 
{
    int i;
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

// i, j from 1 to la
int indexABCol(int i, int j, int *lab)
{
  return j * (*lab - 1) + i - (*lab) + 2;
}


// Function to perform LU factorization of a tridiagonal matrix
// Input:
//   la: Leading dimension of the matrix (lab >= n)
//   n: Order of the matrix
//   kl: Number of subdiagonals (should be 1 for tridiagonal)
//   ku: Number of superdiagonals (should be 1 for tridiagonal)
//   AB: Input matrix (stored in column-major order)
//   lab: Number of diagonals in the matrix (should be 3 for tridiagonal)
// Output:
//   AB: LU factorization of the matrix (overwritten on the input matrix)
//   ipiv: Array of pivot indices
//   info: Status indicator
int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info) 
{
    int i, j, k;
    double temp;
   
    if (*lab != 4 || *kl != 1 || *ku != 1 || *la < *n) 
    {
        *info = -1;
        return *info;
    }
    
    *info = 0;

    for (k = 0; k < *n - 1; k ++) 
    {
        int p_i = k;
        double p_v = AB[indexABCol(k, k, lab)];
        for (int i = k + 1; i < *n; i ++) 
        {
            if (fabs(AB[indexABCol(i, k, lab)]) > fabs(p_v)) 
            {
                p_i = i;
                p_v = AB[indexABCol(i, k, lab)];
            }
        }
        if (p_i != k) 
        {
            for (j = 0; j < *lab; j ++) 
            {
                temp = AB[indexABCol(p_i, j, lab)];
                AB[indexABCol(p_i, j, lab)] = AB[indexABCol(k, j, lab)];
                AB[indexABCol(k, j, lab)] = temp;
            }
            ipiv[k] = p_i;
        }
        else 
        {
            ipiv[k] = k + 1;
        }

        for (i = k + 1; i < *n; i ++) 
        {
            AB[indexABCol(i, k, lab)] /= AB[indexABCol(k, k, lab)];
            for (j = k + 1; j < *n; j ++) 
            {
                AB[indexABCol(i, j, lab)] -= AB[indexABCol(i, k, lab)] * AB[indexABCol(k, j, lab)];
            }
        }
    }
    ipiv[*n - 1] = *n;

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
    eigval[i] = 2.0 * (1 + cos((i + 1) * M_PI / (*la + 1)));
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
  double alpha_opt = 2.0 / (eigmax_poisson1D(eigval, la) + eigmin_poisson1D(eigval, la));
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
    double norm_res;
  
    // Initial residual: r_0 = b - Ax_0
    cblas_dcopy(*la, RHS, 1, temp, 1);  // temp = RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp, 1); //temp = temp - AB * X
    norm_res = cblas_dnrm2(*la, temp, 1) / cblas_dnrm2(*la, RHS, 1);
    resvec[0] = norm_res;

    // Iterate until convergence or max iterations reached
    for (*nbite = 1; *nbite < *maxit && norm_res > *tol; ++ *nbite) 
    {
        cblas_daxpy(*la, alpha, temp, 1, X, 1); //X = X + alpha * temp
        cblas_dcopy(*la, RHS, 1, temp, 1);  // temp = RHS
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp, 1); //temp = temp - AB * X

        // Save the residual norm for each iteration
        norm_res = cblas_dnrm2(*la, temp, 1) / cblas_dnrm2(*la, RHS, 1);
        resvec[*nbite] = norm_res;
        
    }

    // Free allocated memory
    free(temp);
}



// MB = D^{-1}
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int*kl, int *kv)
{
  for (int i = 0; i < *la; i ++)
  {
    MB[i * (*lab) + *kv + *ku] = 1.0 / AB[i * (*lab) + *kv + *ku];
  }
}



// MB = (D - E)^{-1}
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int*kl, int *kv)
{
  // MB = D - E
  for (int i = 0; i < *la; i ++)
  {
    MB[i * (*lab) + *kv + *ku] = AB[i * (*lab) + *kv + *ku];                        // Diagonal elements
    MB[i * (*lab) + *kv + *ku + 1] = AB[i * (*lab) + *kv + *ku + 1];              // Sub-diagonal elements    
  } 
  // Inverse of MB
  for (int j = 0; j < *la; j ++) 
  {
    MB[indexABCol(j, j, lab)] = 1.0 / MB[indexABCol(j, j, lab)];

    for (int i = j + 1; i < *la; i ++) 
    {
      double sum = 0.0;
      for (int k = i; k <= j; k ++) 
      {
        sum += MB[indexABCol(i, k, lab)] * MB[indexABCol(k, j, lab)];
      }
      MB[indexABCol(i, j, lab)] = -sum / MB[indexABCol(i, i, lab)];
    }
  }
}



void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  double *temp = (double *)malloc(sizeof(double) * (*la));
  if (temp == NULL) 
  {
  // Handle memory allocation failure
  fprintf(stderr, "Memory allocation failed for temp\n");
  exit(EXIT_FAILURE);
  }

  double norm_res;
  
  // Initial residual: r_0 = b - Ax_0
  cblas_dcopy(*la, RHS, 1, temp, 1);  // temp = RHS
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp, 1); //temp = temp - AB * X
  norm_res = cblas_dnrm2(*la, temp, 1) / cblas_dnrm2(*la, RHS, 1);
  resvec[0] = norm_res;

  // Iterate until convergence or max iterations reached
  for (*nbite = 1; *nbite < *maxit && norm_res > *tol; ++ *nbite) 
  {
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1.0, MB, *lab, temp, 1, 1.0, X, 1); //X = X + MB * temp

    cblas_dcopy(*la, RHS, 1, temp, 1);  // temp = RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp, 1); //temp = temp - AB * X

    // Save the residual norm for each iteration
    norm_res = cblas_dnrm2(*la, temp, 1) / cblas_dnrm2(*la, RHS, 1);
    resvec[*nbite] = norm_res;
  }
  
  free(temp);
}