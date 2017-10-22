#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define N_ITERATIONS 1
#define DEBUG 1

/*-------------------------------TYPES DEFINITION-----------------------------*/
typedef float** Matrix;
typedef float*  Flat_matrix;
/*----------------------------------------------------------------------------*/


/*-------------------------------FUNCTION PROTOTYPES--------------------------*/
Matrix allocate_zero_matrix(int rows, int cols);
Flat_matrix flattenize_matrix(Matrix A, int rows, int cols);
Matrix deflattenize_matrix(Flat_matrix fmat, int rows, int cols );
void print_matrix(Matrix A, int rows, int cols);
void compute_intern(Matrix A, int rows_A, int n, int head_offset_row, int my_cord);
void compute_extern( Matrix A, int rows_A, int n, Matrix B, int head_offset_row, int rows_B,int my_cord );
void LU_decomposition( int p, Matrix A, int my_cord, int n, int* rows_division, MPI_Comm comm);
void LU_decomposition_serial(Matrix A, int n);
float compute_det_serial(Matrix A, int n);
/*----------------------------------------------------------------------------*/



int main(int argc, char** argv) {

  // variable declaration
  MPI_Comm ring_comm;
  int reorder = 1;
  int period = 1;
  int root_rank,root_cord;
  int n,p;
  int my_rank, my_cord;
  int i,j,iteration;
  float result,partial_det;
  int rows_per_proc,reminder;
  double time_vector[N_ITERATIONS],deviation;
  double average_time, final_time, initial_time;
  // pointer declaration
  int* rows_division;
  int* sendcounts;
  int* displs;
  Matrix A,B;  
  Flat_matrix A_flat, B_flat;

  // save in n the dimension of theMatrix
  n = atoi(argv[1]);

  // Initialize MPI environment
  MPI_Init(&argc, &argv);
  // save the number of processors
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  // save the number of rows for which each processor is in charge
  rows_per_proc = n / p;
  reminder = n % p; //ignored for now

  // create a virtual topology for the ring of k-elements
  MPI_Cart_create(MPI_COMM_WORLD, 1, &p, &period, reorder, &ring_comm);
  // save the rank of each processor in rank
  MPI_Comm_rank(ring_comm, &my_rank);
  // save the coordinates of each given the rank of the processor
  MPI_Cart_coords(ring_comm, my_rank, 1, &my_cord);

  // Allocate vectors for SCATTERV primitive
  sendcounts = (int*) malloc(p*sizeof(int));
  displs = (int*) malloc(p*sizeof(int));
  rows_division = (int*) malloc(p*sizeof(int));
  // Assign to the root the reminder
  for ( i = 0; i < p; i++ ) {
    sendcounts[i] = rows_per_proc*n;
    rows_division[i] = rows_per_proc;
    if ( i < reminder ) {
      sendcounts[i] += n;
      rows_division[i] += 1;
    }
    if ( i == 0 ) {
      displs[i] = 0;
    } else {
      displs[i] = displs[i-1] + sendcounts[i-1];
    }
  }

  root_cord = 0;
  // get the rank of root processor in the topology
  MPI_Cart_rank(ring_comm, &root_cord, &root_rank);


  // if it is the root processor
  if ( my_cord == root_cord) {
    // allocate the data
    A = allocate_zero_matrix(n,n);
    // declare the seed
    srand(time(0));
    for (i = 0; i < n; i++) {
      for ( j = 0; j < n; j++ ) {
        A[i][j] = rand() % 6;
      }
    }
    A_flat = flattenize_matrix(A, n,n);
  #if DEBUG
    print_matrix(A,n,n);
  #endif
  }
  
  B_flat = (Flat_matrix) malloc(sendcounts[my_cord]*sizeof(float));

  average_time = 0;
  for ( iteration = 0; iteration < N_ITERATIONS; iteration++ ) {
    // save initial time of the task
    initial_time = MPI_Wtime();
    // scatter the data from source to all the processors
    MPI_Scatterv(A_flat, sendcounts, displs, MPI_FLOAT, B_flat, sendcounts[my_cord], MPI_FLOAT, root_rank, ring_comm);
    B = deflattenize_matrix(B_flat,rows_division[my_cord],n);
    LU_decomposition(p,B,my_cord,n, rows_division, ring_comm); 
    partial_det = 1;
    for ( i = 0; i < rows_division[my_cord]; i++ ) {
      partial_det = partial_det*B[i][i + displs[my_cord]/n];
    }
    // apply reduce operation (MPI_SUM) on the root processor
    MPI_Reduce(&partial_det, &result, 1, MPI_FLOAT, MPI_PROD, root_rank, ring_comm);
    // save final time of the task
    final_time = MPI_Wtime();
    if ( my_rank == root_rank) {
    #if DEBUG
      printf("DET: %f\n", result);
    #endif
      // and free the data array dinamically allocated
    }
    time_vector[iteration] = final_time - initial_time;
    average_time += time_vector[iteration];
  }
  average_time = average_time/N_ITERATIONS;


  if ( my_rank == root_rank ) {
    deviation = 0;
    for ( i = 0; i < N_ITERATIONS; i++ ) {
      deviation += (time_vector[i] - average_time)*(time_vector[i] - average_time);
    }
    // compute and print the rank of the processor and the time it took to complete the task
    printf("Av_time: %f ,dev: %f\n", average_time, deviation);
  }


  // close the MPI environment
  MPI_Finalize();

  return 0;

}


void print_matrix(Matrix A, int rows, int cols) {

  int i,j;

  printf("R: %d, C: %d\n", rows, cols);
  for (i = 0; i < rows; i++ ) {
    for ( j = 0; j < cols; j++ ) {
      printf("%.2f\t", A[i][j]);
    }
    printf("\n");
  }
  return;
}


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


Matrix deflattenize_matrix(Flat_matrix fmat, int rows, int cols ) {

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



void LU_decomposition(
    int p,                      // Number of processors.
    Matrix A,                   // A matrix.
    int my_cord,                // Cordinates of this processor.
    int n,                      // Number of cols of A.
    int* rows_division,         // How many rows for each processor.
    MPI_Comm comm) {            // Communicator.

  int k,i, rank_src, my_rank;
  int max_rows = rows_division[0];
  int head_offset_row = 0;
  Matrix B;
  Flat_matrix A_flat, B_flat;
  
  for ( i = 1; i < p; i++ ) {
    if ( rows_division[i] > max_rows ) {
      max_rows = rows_division[i];
    }
  }
  B_flat = (Flat_matrix) malloc(max_rows*n*sizeof(float));

  MPI_Cart_rank(comm, &my_cord, &my_rank);

  for ( k = 0; k < my_cord; k++ ) {
    if ( k != 0 ) {
      head_offset_row += rows_division[k-1];
    }
    MPI_Cart_rank(comm, &k, &rank_src);
    MPI_Bcast(B_flat, rows_division[k]*n, MPI_FLOAT, rank_src, comm);
    B = deflattenize_matrix(B_flat, rows_division[k], n);
    for ( i = 0; i < n; i++ )
      printf("%d %d %.2f\n", my_cord,k, B[0][i]);
    compute_extern(A, rows_division[my_cord], n, B, head_offset_row, rows_division[k],my_cord);
    free(B);
  }
  head_offset_row += rows_division[k-1];
  compute_intern(A,rows_division[my_cord],n,head_offset_row,my_cord);
  if ( my_cord != p - 1 ) {
    free(B_flat);
    B_flat = flattenize_matrix(A,rows_division[my_cord],n);
    MPI_Bcast(B_flat, rows_division[my_cord]*n, MPI_FLOAT, my_rank, comm);
    free(B_flat);
  }

  B_flat = (Flat_matrix) malloc(max_rows*n*sizeof(float));
  for ( k = my_cord + 1; k < p - 1; k++ ) {
    MPI_Cart_rank(comm, &k, &rank_src);
    MPI_Bcast(B_flat, rows_division[k]*n, MPI_FLOAT, rank_src, comm);
  }


  print_matrix(A,rows_division[my_cord],n);

  return;

}


void compute_extern(
    Matrix A,                   // A Matrix = internal LU combination. See report.
    int rows_A,                 // Number of rows of A, it depends on the processor coordinate.
    int n,                      // Number of cols of A and B. The partition is done on the rows.
    Matrix B,                   // B Matrix = received sub-matrix.
    int head_offset_row,        // Virtual index of the first row in B.
    int rows_B,
    int my_cord) {              // Number of rows of B.

  int i, j, k;
  int physical_row_index;
  int tail_offset_row = head_offset_row + rows_B;
  

  for ( k = head_offset_row; k < tail_offset_row; k++ ) {
    physical_row_index = k - head_offset_row;
    for ( i = 0; i < rows_A; i++ ) {
      printf("%d: A[%d][%d] = A[%d][%d]/B[%d][%d]\n",my_cord,i,k,i,k,physical_row_index,k);
      A[i][k] = A[i][k]/B[physical_row_index][k];
      for ( j = k+1; j < n; j++ ) {
        printf("%d: A[%d][%d] = A[%d][%d] - A[%d][%d]*B[%d][%d]\n",my_cord,i,j,i,j,i,k,physical_row_index,j);
        A[i][j] = A[i][j] - A[i][k]*B[physical_row_index][j];
      }
    }
  }

}

void compute_intern(
    Matrix A,                   // A matrix = internal LU combination. See report.
    int rows_A,                 // Number of rows of A. 
    int n,                      // Number of cols of A.
    int head_offset_row,
    int my_cord) {      // Virtual index of the first row of A

  int i,j,k;
  float sum;


  for ( i = 1 + head_offset_row; i < head_offset_row + rows_A; i++ ) {
    for ( j = head_offset_row; j < n; j++ ) {
      sum = 0;
      if ( j < i ) {
        for ( k = head_offset_row; k < j; k++ ) {
          printf("%d: sum += A[%d][%d]*A[%d][%d]\n", my_cord, i - head_offset_row,k,k-head_offset_row,j);
          sum += A[i - head_offset_row][k]*A[k-head_offset_row][j];
        }
        printf("%d: A[%d][%d] = (A[%d][%d] - sum)/A[%d][%d]\n",
                my_cord,i - head_offset_row,j, i - head_offset_row,j, j - head_offset_row,j);
        A[i - head_offset_row][j] 
          = (A[i - head_offset_row][j] - sum)/A[j - head_offset_row][j];
      } else {
        for ( k = head_offset_row; k < i; k++ ) {
          printf("%d: sum += A[%d][%d]*A[%d][%d]\n", my_cord, i - head_offset_row,k,k-head_offset_row,j);
          sum += A[i- head_offset_row][k]*A[k - head_offset_row][j];
        }
        printf("%d: A[%d][%d] = (A[%d][%d] - sum)\n",
                my_cord,i - head_offset_row,j, i - head_offset_row,j);
        A[i - head_offset_row][j] = A[i - head_offset_row][j] - sum;
      }
    }
  }
}



float compute_det_serial(Matrix A, int n) {

  int i,j;
  Matrix B;
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


void LU_decomposition_serial(Matrix A, int n) {

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

