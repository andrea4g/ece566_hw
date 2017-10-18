#include <stdio.h>
#include <stdlib.h>

#define N 10


typedef int** matrix;

int** allocate_zero_matrix(int dimension);
void print_matrix(matrix* A_add, int dimension);
void LU_decomposition(matrix* A_add ,matrix* L_add,matrix* U_add);

int main(int argc,char** argv) {

  matrix A,L,U;
  int i,j;

  A = allocate_zero_matrix(N);

  for ( i = 0; i < N; i++ ) {
    for (j = 0; j < N; j++ ) {
      A[i][j] = rand() % 10;
    }
  }


  LU_decomposition( &A,&L, &U); 

  print_matrix(&A, N);
  printf("\n");
  print_matrix(&L, N);
  printf("\n");
  print_matrix(&U, N);


  return 0;

}


void print_matrix(matrix* A_add, int dimension) {

  matrix A = *A_add;
  int i,j;


  for (i = 0; i < dimension; i++ ) {
    for ( j = 0; j < dimension; j++ ) {
      printf("%d\t", A[i][j]);
    }
    printf("\n");
  }
}


int** allocate_zero_matrix(int dimension) {
  
  int** mat;
  int i,j;

  mat = (int **) malloc(dimension*sizeof(int*));

  for ( i = 0; i < dimension; i++ ) {
    mat[i] = (int*) malloc(dimension*sizeof(int));
    for ( j = 0; j < dimension; j++ ) {
      mat[i][j] = 0;
    }
  }
  return mat;
}


void LU_decomposition(matrix* A_add ,matrix* L_add,matrix* U_add) {

  matrix A,L,U;
  int i,j,k;
  int sum;

  A = *A_add;
  L = allocate_zero_matrix(N);
  U = allocate_zero_matrix(N);

  for ( i = 0; i < N; i++ ) {
    L[i][i] = 1;
  }

  for ( i = 0; i < N; i++ ) {
    for ( j = 0; j < N; j++ ) {
      sum = 0;
      if ( j < i ) {
        for ( k = 0; k < j-1; k++ ) {
          sum += L[i][k]*U[k][j];
        }
        L[i][j] = (A[i][j] - sum)/U[j][j];
      } else {
        for ( k = 0; k < i-1; k++ ) {
          sum += L[i][k]*U[k][j];
        }
        U[i][j] = A[i][j] - sum;
      }
    }
  }

  *U_add = U;
  *L_add = L;

  return;
}

