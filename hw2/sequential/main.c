#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef float** matrix;

matrix allocate_zero_matrix(int dimension);
void print_matrix(matrix A, int dimension);
void LU_decomposition_serial(matrix A, int n);
float compute_det_serial(matrix A, int n);

int main(int argc,char** argv) {

  matrix A;
  int i,j,n;

  n = atoi(argv[1]);

  A = allocate_zero_matrix(n);
  
  srand(time(0));
  for ( i = 0; i < n; i++ ) {
    for (j = 0; j < n; j++ ) {
      A[i][j] = rand() % 20 - 10;
    }
  }

  print_matrix(A,n);

  printf("Det: %f\n", compute_det(A,n));

  return 0;

}


void print_matrix(matrix A, int dimension) {

  int i,j;

  for (i = 0; i < dimension; i++ ) {
    for ( j = 0; j < dimension; j++ ) {
      printf("%.2f\t", A[i][j]);
    }
    printf("\n");
  }
}

matrix allocate_zero_matrix(int dimension) {

  matrix mat;
  int i,j;

  mat = (matrix) malloc(dimension*sizeof(float*));

  for ( i = 0; i < dimension; i++ ) {
    mat[i] = (float*) malloc(dimension*sizeof(float));
    for ( j = 0; j < dimension; j++ ) {
      mat[i][j] = 0;
    }
  }
  return mat;
}


void LU_decomposition_serial(matrix A, int n) {

  int i,j,k;
  float sum;

  for ( i = 1; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      sum = 0;
      if ( j < i ) {
        for ( k = 0; k < j; k++ ) {
          sum += A[k][j]*A[i][k];
        }
        A[i][j] = (A[i][j] - sum)/A[j][j];
      } else {
        for ( k = 0; k < i; k++ ) {
          sum += A[i][k]*A[k][j];
        }
        A[i][j] = A[i][j] - sum;
      }
    }
  }


  return;
}


float compute_det_serial(Matrix A, int n) {

  int i,j;
  matrix B;
  float det;
  B = allocate_zero_matrix(n);

  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      B[i][j] = A[i][j];
    }
  }

  LU_decomposition_serial(B, n);
  det = 1;
  for (i = 0; i < n; i++ ) {
    det = det*B[i][i]; 
  }

  free(B);
  return det;
}


