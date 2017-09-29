#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define P 16
#define N 32
#define D 3

int main(int argc, char** argv) {

  MPI_Comm hc_comm;
  int dim[] = {2, 2, 2};
  int period[] = {1, 1, 1};
  int reorder = 1;
  int cord[3];
  int world_size;
  int rank;
  int data[N];
  int partial_data[N/P];
  int i;
  int *result = NULL;
  int partial_sum;
  int total;

  /* */
  MPI_Init(&argc,&argv);

  MPI_Cart_create(MPI_COMM_WORLD, D, dim, period, reorder, &hc_comm);

  MPI_Comm_size(hc_comm, &world_size);

  MPI_Comm_rank(hc_comm, &rank);
  MPI_Cart_coords(hc_comm, rank, D, cord);

  if (cord[0] == 0 && cord[1] == 0 && cord[2] == 0) {
    //srand(NULL);
    for (i = 0; i < N; i++) {
      //data[i] = rand() % 20;
      data[i] = i;
    }
  }
  /* Implement the scatter (one-to-all personalized) */
  MPI_Scatter(data, N/P, MPI_INT, partial_data, N/P, MPI_INT, 0, hc_comm);

  printf("rank: %d, cord: [%d][%d][%d], %d, %d\n", rank, cord[0], cord[1], cord[2], partial_data[0], partial_data[1]);
  /* sum the data */
  partial_sum = partial_data[0] + partial_data[1];
  if (rank == 0) {
    result = malloc(sizeof(int) * P);
  }
  /* Reduce the data to source node */
  MPI_Gather(&partial_sum, 1, MPI_INT, result, 1, MPI_INT, 0, hc_comm);

  if (rank == 0){
    total = 0;
    for (i = 0; i < P; i++) {
      printf("res[%d] = %d\n", i, result[i]);
      total = total + result[i];
    }
    printf("Sum: %d\n", total);
    free(result);
  }

  MPI_Finalize();

  return 0;
}
