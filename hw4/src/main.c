#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>

#define N_ITERATIONS 1
#define DEBUG 1

/*----------------------------------------------------------------------------*/
/*-------------------------------TYPES DEFINITION-----------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*-------------------------------FUNCTION PROTOTYPES--------------------------*/
/*----------------------------------------------------------------------------*/

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
  double time_vector[N_ITERATIONS],deviation;
  double average_time, final_time, initial_time;
  // pointer declaration
  int dims[2];
  Matrix A, C;

  int p_cord[2];
  // save in n the dimension of the Matrix
  n = atoi(argv[1]);
  // Initialize MPI environment
  MPI_Init(&argc, &argv);
  // save the number of processors
  MPI_Comm_size(MPI_COMM_WORLD, &p);
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

  // Cordinates of the root process
  root_cord[0] = 0;
  root_cord[1] = 0;
  // get the rank of root processor in the topology
  MPI_Cart_rank(mesh_comm, root_cord, &root_rank);

  // if it is the root processor
  if (my_rank == root_rank) {

  }

  average_time = 0;
  for ( iteration = 0; iteration < N_ITERATIONS; iteration++ ) {
    // save initial time of the task
    initial_time = MPI_Wtime();

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

/*----------------------------------------------------------------------------*/
/*-------------------------------FUNCTIONS------------------------------------*/
/*----------------------------------------------------------------------------*/

