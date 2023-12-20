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
};


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
    for (i = 0; i < *la; i++) 
    {
        x[i] = (double) i / (double) ((*la) - 1);
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
  int jj;
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
  return i * (*lab) + j;
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
    for (k = 0; k < *n - 1; k++) 
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
        for (i = k + 1; i < k + 3 && i < *n; i++) 
        {
            for (j = k + 1; j < k + 3 && j < *n; j++) 
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