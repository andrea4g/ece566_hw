#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define N_ITERATIONS 20
#define DEBUG 0

/*----------------------------------------------------------------------------*/
/*-------------------------------TYPES DEFINITION-----------------------------*/
/*----------------------------------------------------------------------------*/

typedef float** Matrix;
typedef float*  Flat_matrix;

/*----------------------------------------------------------------------------*/
/*-------------------------------FUNCTION PROTOTYPES--------------------------*/
/*----------------------------------------------------------------------------*/
void parse_input(int argc, char** argv, Matrix X, int n );
Matrix allocate_zero_matrix(int rows, int cols);
Flat_matrix flattenize_matrix(Matrix A, int rows, int cols);
Matrix deflattenize_matrix(Flat_matrix fmat, int rows, int cols);
void matrix_multiply(Matrix a, Matrix b, Matrix dest, int r, int c);
Matrix copy_matrix(Matrix src, int rows, int cols);
void print_matrix(Matrix A, int rows, int cols);
float compute_det_serial(Matrix A, int n);
int square_root(int p);
void cpy_mat(Matrix dest, Matrix src, int rows, int cols);
Flat_matrix flat_block_matrix(int num_1D_blocks, int num_elements_1D_block, Matrix A);
void LU_decomposition_serial(Matrix A, int n);
float compute_det_serial(Matrix A, int n);
void cannon(Matrix src, Matrix dest, MPI_Comm comm, int exponent, int n, int p, int sr_p, int row_per_block, int cols_per_block, int* my_cord, int my_rank, int root_rank, int* sendcounts, int* displs);

/*----------------------------------------------------------------------------*/
/*------------------------------------MAIN------------------------------------*/
/*----------------------------------------------------------------------------*/

int main(int argc, char** argv) {

  // variable declaration
  MPI_Comm mesh_comm;
  int reorder = 1;
  int period[2];
  int root_rank,my_rank;
  int my_cord[2], root_cord[2];
  int n,p;
  int exp;
  int i,j,iteration;
  //float result,partial_det;
  double time_vector[N_ITERATIONS],deviation;
  double average_time, final_time, initial_time;
  // pointer declaration
  int* sendcounts;
  int* displs;
  int dims[2];
  Matrix A, C;
  int sr_p;
  int row_per_block;
  int cols_per_block;

  int p_cord[2];
  char** command_line;
  
  command_line = (char **) malloc(argc*sizeof(char*));
  for ( i = 0; i < argc; i++ ) {
    command_line[i] = strdup(argv[i]);
  }
  // save in n the dimension of the Matrix
  n = atoi(argv[1]);
  // save the exponent
  exp = atoi(argv[2]);
  // Initialize MPI environment
  MPI_Init(&argc, &argv);
  // save the number of processors
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  // Compute the square root.
  sr_p = square_root(p);
  // save the number of row_per_block and cols_per_block for which each processor is in charge
  row_per_block = n / sr_p;
  cols_per_block = n / sr_p;
  // create a virtual topology for the 2-D mesh of k-elements
  dims[0] = sr_p;
  dims[1] = sr_p;
  period[0] = 1;
  period[1] = 1;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, reorder, &mesh_comm);
  // save the rank of each processor in rank
  MPI_Comm_rank(mesh_comm, &my_rank);
  // save the coordinates of each processor given the rank
  MPI_Cart_coords(mesh_comm, my_rank, 2, my_cord);

  // Allocate vectors for SCATTERV primitive
  sendcounts = (int*) malloc(p*sizeof(int));
  displs = (int*) malloc(p*sizeof(int));
  // Assign to the root the reminder
  for (i = 0; i < p; i++) {
    sendcounts[i] = n*n/p;
    MPI_Cart_coords(mesh_comm, i, 2, p_cord);
    displs[i] = (p_cord[0]*sr_p + p_cord[1])*n*n/p;
  }

  // Cordinates of the root process
  root_cord[0] = 0;
  root_cord[1] = 0;
  // get the rank of root processor in the topology
  MPI_Cart_rank(mesh_comm, root_cord, &root_rank);

  // if it is the root processor
  if (my_rank == root_rank) {
    // allocate the data
    A = allocate_zero_matrix(n,n);
    C = allocate_zero_matrix(n,n);
    parse_input(argc, command_line, A, n );
#if DEBUG==1
    print_matrix(A, n, n);
#endif
   
  }

  average_time = 0;
  for ( iteration = 0; iteration < N_ITERATIONS; iteration++ ) {
    // save initial time of the task
    initial_time = MPI_Wtime();

    if (exp == 1) {
    #if DEBUG
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
          printf("%f\t", A[i][j]);
        }
        printf("\n");
      }
    #endif
    }
    else {
      cannon(A, C, mesh_comm, exp, n, p, sr_p, row_per_block, cols_per_block, my_cord, my_rank, root_rank, sendcounts, displs);
    }
    final_time = MPI_Wtime();
    time_vector[iteration] = final_time - initial_time;
    average_time += time_vector[iteration];
  }
  average_time = average_time/N_ITERATIONS;

  if ( my_rank == root_rank ) {
#if DEBUG == 1
    printf("Det A:%f\n",compute_det_serial(A,n));
    printf("Det C with k=%d : %f\n", k, compute_det_serial(C,n));
#endif
    deviation = 0;
    for ( i = 0; i < N_ITERATIONS; i++ ) {
      deviation += (time_vector[i] - average_time)*(time_vector[i] - average_time);
    }
    // compute and print the rank of the processor and the time it took to complete the task
    printf("%f, %f\n", average_time, deviation);
  }

  // close the MPI environment
  MPI_Finalize();

  return 0;
}

/*----------------------------------------------------------------------------*/
/*-------------------------------FUNCTIONS------------------------------------*/
/*----------------------------------------------------------------------------*/

void parse_input(int argc, char** argv, int** X, int n ) {

  int sel;
  int p1,p2; // p1 is the prob of -1 and p2 is prob of +1
  int seq[4];
  int i,j,k;
  int r,idx; 
  int* A;


  A = (int*) malloc(n*sizeof(int));
  sel = atoi(argv[3]);

  if ( sel == 0 ) {
    p1 = atoi(argv[4]);
    p2 = atoi(argv[5]);

    srand(time(NULL));
    for (i = 0; i < n; i++) {
      for ( j = 0; j < n; j++) {
        r = rand() % 100;
        if ( r < p1 ) {
          X[i][j] = -1;
        } else {
          if ( r < p1 + p2) {
            X[i][j] = 1; 
          } else {
            X[i][j] = 0;
          }
        }
      } 
    }
  } else {
    for ( i = 0; i < 4; i++ ) {
      seq[i] = atoi(argv[4 + i]);
    }
    for ( i = 0; i < n; i++ ) {
      A[i] = seq[i % 4];
    }
    for ( i = 0; i < n; i++ ) {
      k = (i+1)*(i+2)/2;
      for ( j = 0; j < n; j++ ) {
        idx = (j - k) % n;
        if ( idx < 0 )
          X[i][j] = A[n + idx];
        else
          X[i][j] = A[idx];
      }
    }
  }

}
void cannon(Matrix src, // source matrix
    Matrix dest,        // destination matrix
    MPI_Comm comm,      // communicator
    int exponent,       // exponent
    int n,              // matrix dimension
    int p,              // number of processors
    int sr_p,           // square root of the processors
    int row_per_block,  // number of rows in the block
    int cols_per_block, // number of columns in the block
    int* my_cord,      // coordinates of the topology
    int my_rank,        // rank
    int root_rank,      // rank of the root
    int* sendcounts,
    int* displs
    )
{

  int q, i, j, k, l, m;
  int uprank, downrank, leftrank, rightrank;
  int shiftsource, shiftdest;
  Matrix B, D, copy, T;
  Flat_matrix A_flat, B_flat, T_flat, C_serial;
  Flat_matrix flatB, flatD;
  MPI_Status status;

  if (my_rank == root_rank) {
    C_serial = (Flat_matrix) malloc(n*n*sizeof(float));
    A_flat = flat_block_matrix(sr_p, n/sr_p, src);
  }
  B_flat = (Flat_matrix) malloc(n*n/p*sizeof(float));
  T = allocate_zero_matrix(row_per_block, cols_per_block);
  D = allocate_zero_matrix(row_per_block, cols_per_block);
  // scatter the data from source to all the processors
  MPI_Scatterv(A_flat, sendcounts, displs, MPI_FLOAT, B_flat, sendcounts[my_rank], MPI_FLOAT, root_rank, comm);
  B = deflattenize_matrix(B_flat, row_per_block, cols_per_block);
  //D = deflattenize_matrix(D_flat, n/sr_p, n/sr_p);
  D = copy_matrix(B, row_per_block, cols_per_block);
  copy = copy_matrix(B, row_per_block, cols_per_block);

  for (q = 0; q < exponent - 1; q++) {
    // reset T matrix
    for (i = 0; i < row_per_block; i++) {
      for (j = 0; j < cols_per_block; j++) {
        T[i][j] = 0;
      }
    }
    // Shift each row_per_block by i on left
    MPI_Cart_shift(comm, 1, -my_cord[0], &shiftsource, &shiftdest);
    flatD = flattenize_matrix(D, row_per_block, cols_per_block);
    MPI_Sendrecv_replace(flatD, row_per_block*cols_per_block, MPI_FLOAT, shiftdest, 1, shiftsource, 1, comm, &status);
    //for (i = 0; i < row_per_block*cols_per_block; i++) {
    //  printf("D_flat --(%d,%d) - %f\n", my_cord[0], my_cord[1], flatD[i]);
    //}
    D = deflattenize_matrix(flatD, row_per_block, cols_per_block);
    // Shift each column up by j
    MPI_Cart_shift(comm, 0, -my_cord[1], &shiftsource, &shiftdest);
    flatB = flattenize_matrix(B, row_per_block, cols_per_block);
    MPI_Sendrecv_replace(flatB, row_per_block*cols_per_block, MPI_FLOAT, shiftdest, 1, shiftsource, 1, comm, &status);
    //for (i = 0; i < row_per_block*cols_per_block; i++) {
    //  printf("B_flat --(%d,%d) - %f\n", my_cord[0], my_cord[1], flatB[i]);
    //}
    B = deflattenize_matrix(flatB, row_per_block, cols_per_block);
    //printf("QUI zio\n");
    for (j = 0; j < sr_p; j++) {
      matrix_multiply(D, B, T, row_per_block, cols_per_block); // D += D*B
      //left circ by 1
      MPI_Cart_shift(comm, 1, -1, &rightrank, &leftrank);
      flatD = flattenize_matrix(D, row_per_block, cols_per_block);
      //printf("ciaone\n");
      MPI_Sendrecv_replace(flatD, row_per_block*cols_per_block, MPI_FLOAT, leftrank, 1, rightrank, 1, comm, &status);
      //for (i = 0; i < row_per_block*cols_per_block; i++) {
      //  printf("D_flat --(%d,%d) - %f\n", my_cord[0], my_cord[1], flatD[i]);
      //}
      D = deflattenize_matrix(flatD, row_per_block, cols_per_block);
      //up circ by 1
      MPI_Cart_shift(comm, 0, -1, &downrank, &uprank);
      flatB = flattenize_matrix(B, row_per_block, cols_per_block);
      MPI_Sendrecv_replace(flatB, row_per_block*cols_per_block, MPI_FLOAT, uprank, 1, downrank, 1, comm, &status);
      B = deflattenize_matrix(flatB, row_per_block, cols_per_block);
    }
    // copy T into D
    for (i = 0; i < row_per_block; i++) {
      for (j = 0; j < cols_per_block; j++) {
        B[i][j] = copy[i][j];
        D[i][j] = T[i][j];
      }
    }
  }
  T_flat = flattenize_matrix(T, row_per_block, cols_per_block);
  MPI_Gather(T_flat, row_per_block*cols_per_block, MPI_FLOAT, C_serial, row_per_block*cols_per_block, MPI_FLOAT, 0, comm);
  // reoder the matrc in the root rank
  if (my_rank == 0) {
    k = 0;
    for (i = 0; i < sr_p; i++ ) {
      for (j = 0; j < sr_p; j++ ) {
        for (l = 0; l < row_per_block; l++) {
          for (m = 0; m < cols_per_block; m++) {
            dest[i*row_per_block + l][j*cols_per_block + m] = C_serial[k];
            k++;
          }
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

void matrix_multiply(Matrix a, Matrix b, Matrix dest, int r, int c) {
  int i, j, k;

  for (i = 0; i < r; i++) {
    for (j = 0; j < c; j++) {
      for (k = 0; k < c; k++) {
        dest[i][j] += a[i][k] * b[k][j];
      }
    }
  }
}
/*----------------------------------------------------------------------------*/

void print_matrix(Matrix A, int rows, int cols) {

  int i,j;

  //printf("R: %d, C: %d\n", rows, cols);
  for (i = 0; i < rows; i++ ) {
    for ( j = 0; j < cols; j++ ) {
      fprintf(stderr,"%.2f\t", A[i][j]);
    }
    fprintf(stderr,"\n");
  }
  return;
}

/*----------------------------------------------------------------------------*/

Matrix allocate_zero_matrix(int rows, int cols) {

  Matrix mat;
  int i,j;

  mat = (Matrix) malloc(rows*sizeof(float*));

  for ( i = 0; i < rows; i++ ) {
    mat[i] = (float*) malloc(cols*sizeof(float));
    for ( j = 0; j < cols; j++ ) {
      mat[i][j] = 0;
    }
  }
  return mat;
}

/*----------------------------------------------------------------------------*/

Matrix copy_matrix(Matrix src, int rows, int cols) {

  Matrix dest;
  int i, j;

  dest = (Matrix) malloc(rows*sizeof(float*));

  for (i = 0; i < rows; i++) {
    dest[i] = (float*) malloc(cols*sizeof(float));
    for (j = 0; j < cols; j++) {
      dest[i][j] = src[i][j];
    }
  }
  return dest;
}

/*----------------------------------------------------------------------------*/

void cpy_mat(Matrix dest, Matrix src, int rows, int cols) {

  int i, j;

  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      dest[i][j] = src[i][j];
    }
  }
}
/*----------------------------------------------------------------------------*/

Matrix deflattenize_matrix(Flat_matrix fmat, int rows, int cols) {

  Matrix mat;
  int i,j;

  mat = allocate_zero_matrix(rows,cols);

  for ( i = 0; i < rows; i++ ) {
    for (j = 0; j < cols; j++) {
      mat[i][j] = fmat[i*cols + j];
    }
  }

  return mat;
}

/*----------------------------------------------------------------------------*/

Flat_matrix flattenize_matrix(Matrix A, int rows, int cols) {

  Flat_matrix fmat;
  int i,j;

  fmat = (Flat_matrix) malloc(rows*cols*sizeof(float));

  for ( i = 0; i < rows; i++ ) {
    for (j = 0; j < cols; j++) {
      fmat[i*cols + j] = A[i][j];
    }
  }

  return fmat;
}

/*----------------------------------------------------------------------------*/

int square_root(int p) {

  int flag,result;
  result = 1;
  flag = 1;

  while (flag) {
    if (result * result == p) {
      flag = 0;
    }
    else {
      result++;
    }
  }
  return result;
}

/*----------------------------------------------------------------------------*/
Flat_matrix flat_block_matrix(int num_1D_blocks, int num_elements_1D_block, Matrix A) {

  int i,j,k,l,m;

  Flat_matrix A_flat;


  A_flat = (Flat_matrix) malloc(
      num_1D_blocks*num_1D_blocks*num_elements_1D_block*num_elements_1D_block*sizeof(float));

  k = 0;
  for ( i = 0; i < num_1D_blocks; i++ ) {
    for ( j = 0; j < num_1D_blocks; j++ ) {
      for ( m = 0; m < num_elements_1D_block; m++ ) {
        for ( l = 0; l < num_elements_1D_block; l++ ) {
          A_flat[k] = A[i*num_elements_1D_block + m][j*num_elements_1D_block + l];
          k++;
        }
      }
    }
  }

  return A_flat;
}

float compute_det_serial(Matrix A, int n) {

  int i,j;
  Matrix B;
  float det;
  B = allocate_zero_matrix(n,n);

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


/* void LU_decomposition_serial(Matrix A, int n) { */

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
      }
      else {
        for ( k = 0; k < i; k++ ) {
          sum += A[i][k]*A[k][j];
        }
        A[i][j] = A[i][j] - sum;
      }
    }
  }
  return;
}
