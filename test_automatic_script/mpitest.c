/* MPI program sample */
/* HPCC, Sept. 2015 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char **argv) {
  int input = 0;
  int world_size;
  int rank;
  int period = 1;
  int reorder = 1;
  int dim = 32;
  double initial_time,final_time;
  MPI_Comm ring_comm;

  input = atoi(argv[1]);


  // Initialize the MPI environment
  MPI_Init(&argc,&argv);
  initial_time = MPI_Wtime();


  // get the total number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // get the rank of current process
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  final_time = MPI_Wtime();
  printf("%d of %d processors. Time: %f\n", rank, world_size, final_time - initial_time);

  MPI_Finalize();

  return 0;
}
