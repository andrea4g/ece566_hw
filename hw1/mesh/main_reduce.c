#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

//#define N 4096000
#define N 32
#define P 16

int main(int argc, char** argv) {

  MPI_Comm mesh_comm;

  double side = P;
  side = sqrt(side);

  int dim[2] = {(int)side, (int)side};
  int period = 0;
  int reorder = 1;
  int cord;
  int world_size;
  int rank;
  int data[N];
  int partial_data[N/P];
  int i;
  int result;
  int partial_sum = 0;
  double initial_time, final_time;

  MPI_Init(&argc,&argv);

  initial_time = MPI_Wtime();

  MPI_Cart_create(MPI_COMM_WORLD, 1, &dim, &period, &reorder, &mesh_comm);

  MPI_Comm_size(mesh_comm, &world_size);

  MPI_Comm_rank(mesh_comm, &rank);
  MPI_Cart_coords(mesh_comm, rank, 1, &cord);

  if (cord[0] == 0 && coord[1] == 0) {
    //srand(NULL);
    for (i = 0; i < N; i++) {
      //data[i] = rand() % 20;
      data[i] = i;
    }
  }

  MPI_Scatter(data, N/P, MPI_INT, partial_data, N/P, MPI_INT, 0, mesh_comm);

  for ( i = 0 ; i < N/P; i++ )
    partial_sum += partial_data[i];

  MPI_Reduce(&partial_sum, &result, 1, MPI_INT, MPI_SUM, 0, mesh_comm);

  if ( rank == 0 )
    printf("Sum: %d\n", result);

  final_time = MPI_Wtime();
  printf("Proc: %d, %f\n", cord, final_time - initial_time);

  MPI_Finalize();

  return 0;

}
