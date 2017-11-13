#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define N_ITERATIONS 1
#define DEBUG 1

#define DIM_i 0
#define DIM_j 1
#define DIM_k 2


/*-------------------------------TYPES DEFINITION-----------------------------*/
typedef float** Matrix;
typedef float*  Flat_matrix;

/*-------------------------------FUNCTION PROTOTYPES--------------------------*/
Matrix allocate_zero_matrix(int rows, int cols);
Flat_matrix flattenize_matrix(Matrix A, int rows, int cols);
Matrix deflattenize_matrix(Flat_matrix fmat, int rows, int cols );
Flat_matrix flat_block_matrix(int num_1D_blocks, int num_elements_1D_block, Matrix A);

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

/*------------------------------------MAIN------------------------------------*/
int main(int argc, char** argv) {

  // variable declaration
  MPI_Comm mesh_comm;
  int reorder = 1;
  int period[3];
  int root_rank,my_rank;
  int my_cord[3], root_cord[3];
  int n,p;
  int i,j,iteration,m,l,k;
  float result,partial_det,mail_box;
  int rows_per_proc,reminder;
  double time_vector[N_ITERATIONS],deviation;
  double average_time, final_time, initial_time;
  // pointer declaration
  int* rows_division;
  int* sendcounts;
  int* displs;
  int dims[3];
  Matrix A,B;
  Flat_matrix A_flat, B_flat;
  MPI_Status status;
  int sr_p;

  int p_cord[2];


  MPI_comm ring_i;
  MPI_comm ring_j;
  MPI_comm ring_k;
  MPI_comm mesh_ij;



  // save in n the dimension of theMatrix
  n = atoi(argv[1]);

  // Initialize MPI environment
  MPI_Init(&argc, &argv);
  // save the number of processors
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  // Compute the square root.
  cr_p = cube_root(p);
  // save the number of rows and cols of each subblock of the initial matrix
  rows_per_proc = n / cr_p;

  // create a virtual topology for the 3-D mesh of p^(1/3)-elements
  dims[0] = cr_p;
  dims[1] = cr_p;
  dims[2] = cr_p;
  period[0] = 1;
  period[1] = 1;
  period[2] = 1;
  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, period, reorder, &mesh_comm);
  // save the rank of each processor in rank
  MPI_Comm_rank(mesh_comm, &my_rank);
  // save the coordinates of each processor given the rank
  MPI_Cart_coords(mesh_comm, my_rank, 3, my_cord);

  dims[DIM_i] = 0;
  dims[DIM_j] = 1;
  dims[DIM_k] = 0;
  MPI_Cart_sub(mesch_comm, dims, &ring_j);
  dims[DIM_i] = 0;
  dims[DIM_j] = 0;
  dims[DIM_k] = 1;
  MPI_Cart_sub(mesch_comm, dims, &ring_k);
  dims[DIM_i] = 1;
  dims[DIM_j] = 0;
  dims[DIM_k] = 0;
  MPI_Cart_sub(mesch_comm, dims, &ring_i);
  dims[DIM_i] = 1;
  dims[DIM_j] = 1;
  dims[DIM_k] = 0;
  MPI_Cart_sub(mesch_comm, dims, &mesh_ij);

  // Cordinates of the root process
  root_cord[0] = 0;
  root_cord[1] = 0;
  root_cord[2] = 0;
  // get the rank of root processor in the topology
  MPI_Cart_rank(mesh_comm, root_cord, &root_rank);

  if ( my_cord[DIM_k] == 0 ) {
    my_flat_block_A = (Flat_matrix) malloc(num_elements_per_block*sizeof(float));
    my_flat_block_B = (Flat_matrix) malloc(num_elements_per_block*sizeof(float));
  }

  // if it is the root processor
  if ( my_rank == root_rank) {
    // allocate the data
    A = allocate_zero_matrix(n,n);
    B = allocate_zero_matrix(n,n);
    // declare the seed
    srand(time(NULL));
    for (i = 0; i < n; i++) {
      for ( j = 0; j < n; j++ ) {
        A[i][j] = 2*(rand() % 2) - 1;i
        B[i][j] = 2*(rand() % 2) - 1;i
      }
    }
    A_flat = flat_block_matrix(cr_p, n/cr_p, A);
    B_flat = flat_block_matrix(cr_p, n/cr_p, B);
    MPI_Scatter(A_flat, num_elements_per_block , MPI_FLOAT,
               my_flat_block_A, num_elements_per_block, MPI_FLOAT, root_rank,
               mesh_ij);
    MPI_Scatter(B_flat, num_elements_per_block, MPI_FLOAT,
               my_flat_block_B, num_elements_per_block, MPI_FLOAT, root_rank,
               mesh_ij);
  } else {
    if ( my_cord[DIM_k] == 0 ) {
      MPI_Scatter(A_flat, num_elements_per_block, MPI_FLOAT,
                my_flat_block_A, num_elements_per_block, MPI_FLOAT, root_rank,
                mesh_ij);
      MPI_Scatter(B_flat, num_elements_per_block, MPI_FLOAT,
                my_flat_block_B, num_elements_per_block, MPI_FLOAT, root_rank,
                mesh_ij);
    }
  }

  if ( my_cord[DIM_k] == 0) {
    dest_cord[DIM_i] = my_cord[DIM_i];
    dest_cord[DIM_j] = my_cord[DIM_j];
    dest_cord[DIM_k] = my_cord[DIM_j];
    MPI_Cart_rank(mesh_comm, dest_cord, &dest_rank);
    MPI_Send(my_flat_block_A,num_elements_per_block,MPI_FLOAT,dest_rank,0,mesh_comm);
    dest_cord[DIM_i] = my_cord[DIM_i];
    dest_cord[DIM_j] = my_cord[DIM_j];
    dest_cord[DIM_k] = my_cord[DIM_i];
    MPI_Cart_rank(mesh_comm, dest_cord, &dest_rank);
    MPI_Send(my_flat_block_B,num_elements_per_block,MPI_FLOAT,dest_rank,0,mesh_comm);
  } else {
    src_cord[DIM_i] = my_cord[DIM_i];
    src_cord[DIM_j] = my_cord[DIM_j];
    src_cord[DIM_k] = 0;
    MPI_Cart_rank(mesh_comm, src_cord, &src_rank);
    if ( my_cord[DIM_j] == my_cord[DIM_k] ) {
      MPI_Recv(my_flat_block_A, num_elements_per_block,MPI_FLOAT,src_rank,0,mesh_comm);
    }
    if ( my_cord[DIM_i] == my_cord[DIM_k]) {
      MPI_Recv(&my_flat_block_B,num_elements_per_block,MPI_FLOAT,src_rank,0,mesh_comm);
    }
  }

  MPI_Comm_group(ring_j,&group_ring_j);
  MPI_Comm_group(ring_i,&group_ring_i);
  MPI_Comm_group(mesh_comm,&group_main);

  src_cord[DIM_j] = my_cord[DIM_k];
  src_cord[DIM_i] = my_cord[DIM_i];
  src_cord[DIM_k] = my_cord[DIM_k];
  MPI_Cart_rank(mesh_comm, src_cord, &src_rank);
  MPI_Group_translate_ranks(group_main, 1, &src_rank,
                            group_ring_j, &src_rank_ring);

  MPI_Bcast(my_flat_block_A, num_elements_per_block, MPI_FLOAT, src_rank_ring,
               ring_j);

  src_cord[DIM_i] = my_cord[DIM_k];
  src_cord[DIM_j] = my_cord[DIM_j];
  src_cord[DIM_k] = my_cord[DIM_k];
  MPI_Cart_rank(mesh_comm, src_cord, &src_rank);
  MPI_Group_translate_ranks(group_main, 1, &src_rank,
                            group_ring_i, &src_rank_ring);

  MPI_Bcast(my_flat_block_B, num_elements_per_block, MPI_FLOAT, src_rank_ring,
               ring_i);


  partial_C = internal_mul(A_block,B_block,n/cr_p,n/cr_p);
  

  src_cord[DIM_i] = my_cord[DIM_i];
  src_cord[DIM_j] = my_cord[DIM_j];
  src_cord[DIM_k] = 0;
  MPI_Cart_rank(mesh_comm, src_cord, &src_rank);
  MPI_Group_translate_ranks(group_main, 1, &src_rank,
                            group_ring_k, &src_rank_ring);


  MPI_Reduce(
      my_flat_block_C,
      mail_box,
      num_elements_per_block,
      MPI_FLOAT,
      MPI_SUM,
      src_rank,
      group_ring_k);


  if ( my_cord[k] == 0 ) {
    
  }


  // close the MPI environment
  MPI_Finalize();

  return 0;
}

/*-------------------------------FUNCTIONS------------------------------------*/

Matrix internal_mul(Matrix A, Matrix B, int rows, int cols) {

  int i,j,k;
  Matrix C = allocate_zero_matrix(rows,cols);

  for ( i = 0; i < rows; i++ ) {
    for ( j = 0; j < rows; j++ ) {
      for ( k = 0; k < cols; k++ ) {
          C[i][j] += A[i][k]*B[k][j];
        }
    }
  }

  return C;

}







void print_matrix(Matrix A, int rows, int cols) {

  int i,j;

  printf("R: %d, C: %d\n", rows, cols);
  for (i = 0; i < rows; i++ ) {
    for ( j = 0; j < cols; j++ ) {
      fprintf(stderr,"%.2f\t", A[i][j]);
    }
    fprintf(stderr,"\n");
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

  mailbox_up = allocate_zero_matrix(n/sr_p,n/sr_p);
  mailbox_left = allocate_zero_matrix(n/sr_p,n/sr_p);

  my_row = my_cord[0];
  my_col = my_cord[1];

  for ( k = 0; (k < my_row) && (k <= my_col); k++ ) {
    receive(comm, my_col, k, mailbox_up, n/sr_p);
    if ( my_col == k ) {
      compute_only_up(A, mailbox_up, n/sr_p);
      send_on_row(comm, A, n/sr_p, sr_p, my_row, my_col);
    } else {
      receive(comm, k, my_row,mailbox_left,n/sr_p);
      compute_up_left(A,mailbox_left,mailbox_up,n/sr_p);
      //printf("P[%d][%d] here\n", my_row, my_col);
    }
  }
  if ( my_col == k ) {
  #if DEBUG
    printf("HERE\n");
  #endif
    compute_intern(A,n/sr_p);
  #if DEBUG
    print_matrix(A, n/sr_p, n/sr_p);
  #endif
    if ( (my_row < (sr_p - 1)) && (my_col < sr_p -1 )) {
      send_on_row(comm, A, n/sr_p, sr_p, my_row, my_col);
      send_on_col(comm, A, n/sr_p, sr_p, my_row, my_col);
    }
  }
  if ( my_col > k ) {
    receive(comm, k, my_row, mailbox_left, n/sr_p);
    compute_only_left(A,mailbox_left,n/sr_p);
    if ( my_row < sr_p - 1 ) {
      send_on_col(comm, A, n/sr_p, sr_p,my_row, my_col);
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
      A[i][j] = (A[i][j]-sum)/B[j][j];
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
#if DEBUG
  print_matrix(A,n,n);
#endif
  for ( i = 1; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      sum = 0;
      if ( j < i ) {
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




