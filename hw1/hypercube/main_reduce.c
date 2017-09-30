#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 4096000
#define D 4
#define P (1 << D)


int main(int argc, char** argv) {

  MPI_Comm hc_comm;
  int dim[D];
  int period[D];
  int reorder = 1;
  int cord[D];
  int world_size;
  int rank;
  int data[N];
  int partial_data[N/(P)];
  int i;
  int result;
  int partial_sum;
  double final_time, initial_time;
  char flag;

  MPI_Init(&argc,&argv);

  for ( i = 0; i < D; i++ ) {
    dim[i] = 2;
    period[i] = 1;
  }

  initial_time = MPI_Wtime();

  MPI_Cart_create(MPI_COMM_WORLD, D, dim, period, reorder, &hc_comm);

  MPI_Comm_size(hc_comm, &world_size);

  MPI_Comm_rank(hc_comm, &rank);
  MPI_Cart_coords(hc_comm, rank, D, cord);

  flag = 1;
  for ( i = 0; i < D; i++ ) {
    if ( cord[i] == 1 )
      flag = 0;
  }

  if (flag) {
    //srand(NULL);
    for (i = 0; i < N; i++) {
      //data[i] = rand() % 20;
      data[i] = i;
    }
  }

  MPI_Scatter(data, N/(P), MPI_INT, partial_data, N/(P), MPI_INT, 0, hc_comm);

  partial_sum = 0;
  for ( i = 0 ; i < N/(P); i++ )
    partial_sum += partial_data[i];

  MPI_Reduce(&partial_sum, &result, 1, MPI_INT, MPI_SUM, 0, hc_comm);

  if ( rank == 0 )
    printf("Sum: %d\n", result);

  final_time = MPI_Wtime();
  printf("Proc: %d, time: %f\n", rank, final_time-initial_time);

  MPI_Finalize();

  return 0;

}
