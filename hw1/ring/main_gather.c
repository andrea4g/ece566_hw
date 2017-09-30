#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define P 16
#define N 32

int main(int argc, char** argv) {

  MPI_Comm ring_comm;
  int dim[] = {4,4};
  int period[] = {1,1};
  int reorder = 1;
  int cord[2];
  int world_size;
  int rank;
  int data[N];
  int partial_data[N/P];
  int i;
  int *result = NULL;
  int partial_sum;
  int total;
  double initial_time, final_time;

  /* */
  MPI_Init(&argc,&argv);

  initial_time = MPI_Wtime();

  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &ring_comm);

  MPI_Comm_size(ring_comm, &world_size);

  MPI_Comm_rank(ring_comm, &rank);
  MPI_Cart_coords(ring_comm, rank, 2, cord);

  if (cord[0] == 0 && cord[1] == 0) {
    //srand(NULL);
    for (i = 0; i < N; i++) {
      //data[i] = rand() % 20;
      data[i] = i;
    }
  }

  /* Implement the broadcast (one-to-all communication) */
  MPI_Bcast(&data, N, MPI_INT, 0, ring_comm);

  printf("rank: %d, cord: [%d][%d]", rank, cord[0], cord[1]);
  for (i = 0; i < N; i++) {
    printf("%d\t", data[i]);
  }
  printf("\n");

  /* sum the data */
  partial_sum = data[rank*2] + data[rank*2+1];

  if (rank == 0) {
    result = malloc(sizeof(int) * P);
  }
  /* Reduce the data to source node */
  MPI_Gather(&partial_sum, 1, MPI_INT, result, 1, MPI_INT, 0, ring_comm);

  if (rank == 0){
    total = 0;
    for (i = 0; i < P; i++) {
      printf("res[%d] = %d\n", i, result[i]);
      total = total + result[i];
    }
    printf("Sum: %d\n", total);
    free(result);
  }

  final_time = MPI_Wtime();
  printf("Proc: %d, %f\n", cord, final_time - initial_time);

  MPI_Finalize();

  return 0;
}
