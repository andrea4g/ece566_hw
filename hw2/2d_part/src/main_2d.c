#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define N_ITERATIONS 50
#define DEBUG 0

/*-------------------------------TYPES DEFINITION-----------------------------*/
typedef float** Matrix;
typedef float*  Flat_matrix;
/*----------------------------------------------------------------------------*/


/*-------------------------------FUNCTION PROTOTYPES--------------------------*/
Matrix allocate_zero_matrix(int rows, int cols);
Flat_matrix flattenize_matrix(Matrix A, int rows, int cols);
Matrix deflattenize_matrix(Flat_matrix fmat, int rows, int cols );
void print_matrix(Matrix A, int rows, int cols);
void LU_decomposition( int p, int sr_p, Matrix A, int* my_cord, int n, int* rows_division, MPI_Comm comm);

void compute_intern(Matrix A,int n);
void compute_only_left(Matrix A, Matrix B, int n);
void compute_up_left(Matrix A, Matrix B, Matrix C, int n);
void compute_only_up(Matrix A, Matrix B, int n);
void send_on_col(MPI_Comm mesh_comm, Matrix A, int n, int sr_p, int srt_row, int srt_col);
void receive(MPI_Comm mesh_comm, int col, int row, Matrix mailbox, int n);
void send_on_row(MPI_Comm mesh_comm, Matrix A, int n, int sr_p, int srt_row, int srt_col);

void LU_decomposition_serial(Matrix A, int n);
float compute_det_serial(Matrix A, int n);
int square_root(int p);
/*----------------------------------------------------------------------------*/




int main(int argc, char** argv) {

  // variable declaration
  MPI_Comm mesh_comm;
  int reorder = 1;
  int period[2];
  int root_rank,my_rank;
  int my_cord[2], root_cord[2];
  int n,p;
  int i,j,iteration,m,l,k;
  float result,partial_det;
  int rows_per_proc,reminder;
  double time_vector[N_ITERATIONS],deviation;
  double average_time, final_time, initial_time;
  // pointer declaration
  int* rows_division;
  int* sendcounts;
  int* displs;
  int dims[2];
  Matrix A,B;
  Flat_matrix A_flat, B_flat;

  int sr_p;
  
  int p_cord[2];

  // save in n the dimension of theMatrix
  n = atoi(argv[1]);

  // Initialize MPI environment
  MPI_Init(&argc, &argv);
  // save the number of processors
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  // Compute the square root.
  sr_p = square_root(p);
  // save the number of rows and cols for which each processor is in charge
  rows_per_proc = n / sr_p;

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
  for ( i = 0; i < p; i++ ) {
    sendcounts[i] = n*n/p;
    MPI_Cart_coords(mesh_comm, i, 2, p_cord);
    displs[i] = (p_cord[0]*3 + p_cord[1])*n*n/p;
  }

  // Cordinates of the root process
  root_cord[0] = 0;
  root_cord[1] = 0;
  // get the rank of root processor in the topology
  MPI_Cart_rank(mesh_comm, root_cord, &root_rank);

  // if it is the root processor
  if ( my_rank == root_rank) {
    // allocate the data
    A = allocate_zero_matrix(n,n);
    // declare the seed
    srand(time(NULL));
    for (i = 0; i < n; i++) {
      for ( j = 0; j < n; j++ ) {
        A[i][j] = rand() % 10 + 1;
      }
    }
    A_flat = (Flat_matrix) malloc(n*n*sizeof(float));
    k = 0;
    for ( i = 0; i < sr_p; i++ ) {
      for ( j = 0; j < sr_p; j++ ) {
        for ( m = 0; m < n/sr_p; m++ ) {
          for ( l = 0; l < n/sr_p; l++ ) {
            A_flat[k] = A[i*n/sr_p + m][j*n/sr_p + l];
            k++;
          }
        }
      }
    }

  #if DEBUG
    print_matrix(A,n,n);
  #endif
  }
  
  B_flat = (Flat_matrix) malloc(n*n/p*sizeof(float));

  average_time = 0;
  for ( iteration = 0; iteration < N_ITERATIONS; iteration++ ) {
    // save initial time of the task
    initial_time = MPI_Wtime();
    // scatter the data from source to all the processors
    MPI_Scatterv(A_flat, sendcounts, displs, MPI_FLOAT, B_flat, sendcounts[my_rank], MPI_FLOAT, root_rank, mesh_comm);
    B = deflattenize_matrix(B_flat,n/sr_p,n/sr_p);
    LU_decomposition(p,sr_p,B,my_cord,ni/sr_p, rows_division, mesh_comm); 
    partial_det = 1;

/*---------------------------------------------------------------------------------------------------------------------------*/
    // apply reduce operation (MPI_SUM) on the root processor
   // MPI_Reduce(&partial_det, &result, 1, MPI_FLOAT, MPI_PROD, root_rank, ring_comm);
    // save final time of the task
    final_time = MPI_Wtime();
    if ( my_rank == root_rank) {
    #if DEBUG
      printf("DET serial: %f\n", compute_det_serial(A,n));   
      printf("DET parallel: %f\n", result);
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
    printf("%f, %f\n", average_time, deviation);
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
    int sr_p,                   // Square root number of procs.
    Matrix A,                   // A matrix.
    int* my_cord,               // Cordinates of this processor.
    int n,                      // Number of cols of A.
    int* rows_division,         // How many rows for each processor.
    MPI_Comm comm) {            // Communicator.

  int my_row, my_col;
  int k;

  Matrix mailbox_up, mailbox_left;

  mailbox_up = allocate_zero_matrix(n,n);
  mailbox_left = allocate_zero_matrix(n,n);

  my_row = my_cord[0];
  my_col = my_cord[1];

  for ( k = 0; (k < my_row) && (k <= my_col); k++ ) {
    receive(comm, k, my_row, mailbox_up, n);
    if ( my_col == k ) {
      compute_only_up(A, mailbox_up, n);
      send_on_row(comm, A, n, sr_p, my_row, my_col);
    } else {
      receive(comm, my_row,k,mailbox_left,n);
      compute_up_left(A,mailbox_left,mailbox_up,n);
    }
  }
  if ( my_col == k ) {
    compute_intern(A,n);
    if ( (my_row < (p - 1)) && (my_col < p -1 )) {
      send_on_row(comm, A, n, sr_p, my_row, my_col);
      send_on_col(comm, A, n, sr_p, my_row, my_col);
    }
  }
  if ( my_col > k ) {
    receive(comm, k, my_row, mailbox_left, n);
    compute_only_left(A,mailbox_left,n);
    if ( my_row < p - 1 ) {
      send_on_col(comm, A, n, sr_p,my_row, my_col);
    }
  }

  return;
}

void send_on_row(MPI_Comm mesh_comm, Matrix A, int n, int sr_p, int srt_row, int srt_col) {

  Flat_matrix A_flat;
  int i, dest_rank;
  int dest_cord[2];

  A_flat = flattenize_matrix(A,n,n);

  dest_cord[0] = srt_row;
  for ( i = srt_col + 1; i < sr_p; i++ ) {
    dest_cord[1] = i;
    MPI_Cart_rank(mesh_comm, dest_cord, &dest_rank);
    MPI_Send(A_flat, n*n, MPI_FLOAT, dest_rank, 0, mesh_comm);
  }
  free(A_flat);

  return;
}


void receive(MPI_Comm mesh_comm, int col, int row, Matrix mailbox, int n) {

  Flat_matrix A_flat;
  int i, j, k, src_rank;
  int src_cord[2];
  MPI_Status status;

  A_flat = (Flat_matrix) malloc(n*n*sizeof(float));

  src_cord[0] = row;
  src_cord[1] = col;
  MPI_Cart_rank(mesh_comm, src_cord, &src_rank);
  MPI_Recv(A_flat, n*n, MPI_FLOAT, src_rank, 0, mesh_comm, &status);

  k = 0;
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      mailbox[i][j] = A_flat[k++];
    }
  }

  free(A_flat);

  return;
}








void send_on_col(MPI_Comm mesh_comm, Matrix A, int n, int sr_p, int srt_row, int srt_col) {

  Flat_matrix A_flat;
  int i,dest_rank;
  int dest_cord[2];

  A_flat = flattenize_matrix(A,n,n);

  dest_cord[1] = srt_col;
  for ( i = srt_row + 1; i < sr_p; i++ ) {
    dest_cord[0] = i;
    MPI_Cart_rank(mesh_comm, dest_cord, &dest_rank);
    MPI_Send(A_flat, n*n, MPI_FLOAT, dest_rank, 0, mesh_comm);
  }
  free(A_flat);
  return;
}







void compute_only_up(Matrix A, Matrix B, int n) {

  int i,j,k;
  float sum;

  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      sum = 0;
      for ( k = 0; k < j; k++ ) {
        sum += A[i][k]*B[k][j];
      }
      A[i][j] = A[i][j]/B[j][j];
    }
  }

  return;
}


void compute_up_left(Matrix A, Matrix B, Matrix C, int n) {

  int i,j,k;
  float sum;

  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      sum = 0; 
      for ( k = 0; k < n; k++) {
        sum += B[i][k]*C[k][j];
      }
      A[i][j] = A[i][j] - sum;
    }
  }
  return;
}


void compute_only_left(Matrix A, Matrix B, int n) {

  int i,j,k;
  float sum;

  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      sum = 0;
      for ( k = 0; k < i; k++ ) {
        sum += B[i][k]*A[k][j];
      }
        A[i][j] = A[i][j] - sum;
    }
  }
 

  return;

}


void compute_intern(Matrix A,int n) {

  int i,j,k;
  float sum;

  for ( i = 1; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      sum = 0;
      if ( i < j ) {
        for ( k = 0; k < j; k++ ) {
          sum += A[i][k]*A[k][j];
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




int square_root(int p) {

  int flag,result;
  result = 1;
  flag = 1;

  while ( flag ) {
    if ( result*result == p ) {
      flag = 0;
    } else {
      result++;
    }
  }

  return result;

}



