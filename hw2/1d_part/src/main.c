#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N_ITERATIONS 10
#define DEBUG 0

typedef float** Matrix;
typedef float*  Flat_matrix;

Matrix allocate_zero_matrix(int dimension);
Flat_matrix flattenize_matrix(Matrix* A_add, int dimension);

void print_matrix(Matrix* A_add, int dimension);


int main(int argc, char** argv) {

  // variable declaration
  MPI_Comm ring_comm;
  int reorder = 1;
  int period = 1;
  int root_rank,root_cord;

  int n,p;
  int my_rank;
  int i,iteration;
  int result,partial_sum;
  int rows_per_proc,reminder;
  double time_vector[N_ITERATIONS],deviation;
  double average_time, final_time, initial_time;
  // pointer declaration
  Matrix A; 
  Flat_matrix A_flat;


  

  int my_cord;
  int* data;
  int* partial_data;
  int* sendcounts;
  int* displs;

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
  // Assign to the root the reminder
  for ( i = 0; i < p; i++ ) {
    sendcounts[i] = rows_per_proc*n;
    if ( i < reminder ) {
      sendcounts[i] += n;
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
    srand((unsigned int) 0);
    for (i = 0; i < n; i++) {
      for ( j = 0; j < n; j++ ) {
        A[i][j] = rand() % 100;
      }
    }
    A_flat = flattenize_matrix(&A, n);
  #if DEBUG
    print_matrix(&A,n);
  #endif
  } else {
    B_flat = (Flat_matrix) malloc(sendcounts[my_cord]*sizeof(float));
  }



  average_time = 0;
  for ( iteration = 0; iteration < N_ITERATIONS; iteration++ ) {
    // save initial time of the task
    initial_time = MPI_Wtime();
    // scatter the data from source to all the processors
    MPI_Scatterv(A_flat, sendcounts, displs, MPI_FLOAT, B_flat, sendcounts[my_cord], MPI_FLOAT, root_rank, ring_comm);
    B = deflattenize_matrix(B_flat,sendcounts[my_cord]/n,n);
    LU_decomposition();
    
    // apply reduce operation (MPI_SUM) on the root processor
    MPI_Reduce(&partial_sum, &result, 1, MPI_INT, MPI_SUM, root_rank, ring_comm);
    // save final time of the task
    final_time = MPI_Wtime();
    if ( my_rank == root_rank) {
    #if DEBUG
      printf("Sum: %d\n", result);
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
    free(data);
  }

  // free the dynamic memory allocated
  free(partial_data);

  // close the MPI environment
  MPI_Finalize();

  return 0;

}


void print_matrix(Matrix* A_add, int dimension) {

  Matrix A = *A_add;
  int i,j;

  for (i = 0; i < dimension; i++ ) {
    for ( j = 0; j < dimension; j++ ) {
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


Flat_matrix flattenize_matrix(Matrix* A_add, int dimension) {

  Matrix mat;
  Flat_matrix fmat;
  int i,j;

  mat = *A_add;
  fmat = (Flat_matrix) malloc(dimension*dimension*sizeof(float));

  for ( i = 0; i < dimension; i++ ) {
    for (j = 0; j < dimension; j++) {
      fmat[i*dimension + j] = mat[i][j];
    }
  }

  return fmat;

}



void LU_decomposition() {

  if ( my_cord == 0 ) {
    compute_adjustment();
    send_block();
  } else {
    for ( k = 0; k < my_cord - 1; k++ ) {
      receive_block();
      compute_adjustment();
    }
    send_block();
  }

  return;

}



void adjustment() {

  for ( i = 1 + first_row_index; i < number_of_rows[0]; i++ ) {
    for ( j = first_row_index; j < n; j++ ) {
      sum = 0;
      if ( j < i) {
        for ( k = first_row_index; k < j; k++ ) {
          sum += A[i - first_row_index][k]*A[k-first_row_index][j];
        }
        A[i - first_row_index][j] = (A[i - first_row_index][j] - sum)/A[j][j];
      } else {
        for ( k = first_row_index; k < i + first_row_index; k++ ) {
          sum += A[i-first_row_index][k]*A[k - first_row_index][j];
        }
        A[i][j] = A[i][j] - sum;
      }
    }
  }
}





void adjustment_0() {

  for ( i = 1; i < number_of_rows[0]; i++ ) {
    for ( j = 0; j < n; j++ ) {
      sum = 0;
      if ( j < i) {
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
}

void adjustemnt_extern() {

  

  for ( i = 1; i < number_of_rows[0]; i++ ) {
    for ( j = 0; j < n; j++ ) {
      sum = 0;
      if ( j < i) {
        for ( k = 0; k < j; k++ ) {
          sum += A[i][k]*A[k][j];
        }
      } else {
        for ( k = 0; k < i; k++ ) {
          sum += A[i][k]*A[k][j];
        }
        A[i][j] = A[i][j] - sum;
      }
    }
  }

}











