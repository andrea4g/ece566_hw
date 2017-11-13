#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define N_ITERATIONS 1
#define DEBUG 1

/*----------------------------------------------------------------------------*/
/*-------------------------------TYPES DEFINITION-----------------------------*/
/*----------------------------------------------------------------------------*/

typedef float** Matrix;
typedef float*  Flat_matrix;

/*----------------------------------------------------------------------------*/
/*-------------------------------FUNCTION PROTOTYPES--------------------------*/
/*----------------------------------------------------------------------------*/

Matrix allocate_zero_matrix(int rows, int cols);
Flat_matrix flattenize_matrix(Matrix A, int rows, int cols);
Matrix deflattenize_matrix(Flat_matrix fmat, int rows, int cols);
void cannon(Matrix dest, Matrix src, MPI_Comm comm, int rows, int cols, int* coordsi, int num_blocks);
void matrix_multiply(Matrix a, Matrix b, Matrix dest, int r, int c);
Matrix copy_matrix(Matrix src, int rows, int cols);
void print_matrix(Matrix A, int rows, int cols);
float compute_det_serial(Matrix A, int n);
int square_root(int p);

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
  int i,j,iteration,m,l,k;
  //float result,partial_det;
  int rows_per_proc;
  double time_vector[N_ITERATIONS],deviation;
  double average_time, final_time, initial_time;
  // pointer declaration
  int* sendcounts;
  int* displs;
  int dims[2];
  Matrix A, B, D;
  Flat_matrix A_flat, B_flat;
  int sr_p;
  MPI_Status status;

  int p_cord[2];

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
//#if DEBUG
//  printf("p = %d, sr_p = %d\n", p, sr_p);
//#endif
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
    displs[i] = (p_cord[0]*sr_p + p_cord[1])*n*n/p;
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
      #if DEBUG
        A[i][j] = 1;
      #else
        num = rand() % 100;
        if (num % 3 == 0) {
          A[i][j] = 1;
        }
        else if (num % 3 == 1){
          A[i][j] = 0;
        }
        else if (num % 3 == 2){
          A[i][j] = -1;
        }
      #endif
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

  //#if DEBUG
  //  printf("Original matrix\n");
  //  print_matrix(A,n,n);
  //#endif
  }

  B_flat = (Flat_matrix) malloc(n*n/p*sizeof(float));

  average_time = 0;
  for ( iteration = 0; iteration < N_ITERATIONS; iteration++ ) {
    // save initial time of the task
    initial_time = MPI_Wtime();
    // scatter the data from source to all the processors
    MPI_Scatterv(A_flat, sendcounts, displs, MPI_FLOAT, B_flat, sendcounts[my_rank], MPI_FLOAT, root_rank, mesh_comm);
    B = deflattenize_matrix(B_flat, n/sr_p, n/sr_p);
    D = allocate_zero_matrix(n/sr_p, n/sr_p);
    D = copy_matrix(B, n/sr_p, n/sr_p);
  //#if DEBUG
  //  printf("(%d,%d) rank %d: ", my_cord[0], my_cord[1], my_rank);
  //  for (i = 0; i < n/sr_p; i++) {
  //    for (j = 0; j < n/sr_p; j++) {
  //      printf("%.2f\t", D[i][j]+1);
  //    }
  //    printf("\n");
  //  }
  //#endif

  int uprank, downrank, leftrank, rightrank;

  int rows = n/sr_p;
  int cols = n/sr_p;
  int i, j, k;
  Flat_matrix flatB, flatD;

  int shiftsource, shiftdest;
  // Shift each rows by i on left
  MPI_Cart_shift(mesh_comm, 0, -my_cord[0], &shiftsource, &shiftdest);
  //printf("prima flat D\n");
  flatD = flattenize_matrix(D, rows, cols);
  //printf("quasi D\n");
  MPI_Sendrecv_replace(&flatD, rows*cols, MPI_FLOAT, shiftdest, 1, shiftsource, 1, mesh_comm, &status);
  //printf("dopo D\n");
  D = deflattenize_matrix(flatD, rows, cols);
  //printf("dopo deflat D\n");
  // Shift each column up by j
  MPI_Cart_shift(mesh_comm, 1, -my_cord[1], &shiftsource, &shiftdest);
  //printf("prima flat B\n");
  flatB = flattenize_matrix(B, rows, cols);
  //printf("quasi B\n");
  MPI_Sendrecv_replace(&flatB, rows*cols, MPI_FLOAT, shiftdest, 1, shiftsource, 1, mesh_comm, &status);
  //printf("dopo B\n");
  B = deflattenize_matrix(flatB, rows, cols);
  //printf("dopo deflat B\n");

  printf("QUI\n");

//  for (k = 0; k < sr_p; k++) {
//    for (i = 0; i < sr_p; i++) {
      for (j = 0; j < sr_p; j++) {
        matrix_multiply(D, B, D, rows, cols); // D += D*B
        //left circ by 1
        MPI_Cart_shift(mesh_comm, 0, -1, &rightrank, &leftrank);
        flatD = flattenize_matrix(D, rows, cols);
        printf("ciao\n");
        MPI_Sendrecv_replace(&flatD, rows*cols, MPI_FLOAT, leftrank, 1, rightrank, 1, mesh_comm, &status);
        D = deflattenize_matrix(flatD, rows, cols);
        printf("meta'\n");
        //up circ by 1
        MPI_Cart_shift(mesh_comm, 1, -1, &downrank, &uprank);
        flatB = flattenize_matrix(B, rows, cols);
        MPI_Sendrecv_replace(&flatB, rows*cols, MPI_FLOAT, uprank, 1, downrank, 1, mesh_comm, &status);
        B = deflattenize_matrix(flatB, rows, cols);
      }
//    }
//  }

  //src_flat = flattenize_matrix(src_cp, rows, cols);

  //MPI_Gather(&src_flat, rows*columns, MPI_FLOAT, dest, 1, MPI_FLOAT, 0, comm);



  //  printf("Cannon start\n");
  //  cannon(D, B, mesh_comm, n/sr_p, n/sr_p, my_cord, sr_p);
  //  printf("dopo cannon\n");
  //#if DEBUG
  //  for (i = 0; i < n/sr_p; i++) {
  //    printf("(%d,%d): ", my_cord[0], my_cord[1]);
  //    for (j = 0; j < n/sr_p; j++) {
  //      printf("%.2f\t", D[i][j]);
  //    }
  //    printf("\n");
  //  }
  //#endif



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
//    printf("%f, %f\n", average_time, deviation);
  }

  // close the MPI environment
  MPI_Finalize();

  return 0;
}

/*----------------------------------------------------------------------------*/
/*-------------------------------FUNCTIONS------------------------------------*/
/*----------------------------------------------------------------------------*/

void cannon(Matrix dest, Matrix src, MPI_Comm comm, int rows, int cols, int* coords, int num_blocks) {
}

/*----------------------------------------------------------------------------*/

void matrix_multiply(Matrix a, Matrix b, Matrix dest, int r, int c) {
  int i, j, k;

  for (i = 0; i < r; i++) {
    for (j = 0; j < c; j++) {
      for (k = 0; k < r; k++) {
        dest[i][j] += a[i][k] + b[k][j];
      }
    }
  }
}
/*----------------------------------------------------------------------------*/

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
