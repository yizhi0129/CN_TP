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
    for (j = 0; j < *lab; j++) 
    {
        for (i = 0; i < *la; i++) 
        {
            AB[i * (*lab) + j] = 0.0;
        }
    }

    // Remplissage de la bande principale
    for (i = 0; i < *la; i++) 
    {
        AB[i * (*lab) + ku + kl] = 2.0;  // Éléments sur la diagonale principale
    }

    // Remplissage des diagonales au-dessus et en dessous de la diagonale principale
    for (i = 0; i < *la - 1; i++) 
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
    for (j = 0; j < *lab; j++) 
    {
        for (i = 0; i < *la; i++) 
        {
            AB[i * (*lab) + j] = 0.0;
        }
    }

    // Remplissage de la diagonale principale
    for (i = 0; i < *la; i++) 
    {
        AB[i * (*lab) + ku + kl] = 1.0;  // Éléments sur la diagonale principale
    }
}


void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1)
{
  int ii;
  for (ii = 0; ii < (* la); ii ++)
  {
    RHS[ii] = 0.0;
  }
  RHS[0] = (* BC0);
  RHS[(* la) - 1] = (* BC1);
};


void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1)
{
  int ii;
  for (ii = 0; ii < (* la); ii ++)
  {
    EX_SOL[ii] = (* BC0) + (X[ii] - X[0]) * ((* BC1) - (* BC0)) / (X[(* la) - 1] - X[0]);
  }
};  

void set_grid_points_1D(double* x, int* la)
{
  int ii;
  for (ii = 0; ii < (* la); ii ++)
  {
    x[ii] = (double) ii / (double) ((* la) - 1);
  }
};

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
  return 0;
}
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info)
{
  return *info;
}
