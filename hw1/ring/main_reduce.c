#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


int main(int argc, char** argv) {

  MPI_Comm ring_comm;
  int period = 1;
  int reorder = 1;
  int cord;
  int k;
  int rank;
  int n;
  int* data;
  int* partial_data;
  int card_partial_data;
  int i;
  int root_cord,root_rank;
  int result;
  int partial_sum = 0;
  double initial_time, final_time;

  n = atoi(argv[1]);                // Number of elements of the array A

  MPI_Init(&argc,&argv);

  initial_time = MPI_Wtime();

  MPI_Comm_size(MPI_COMM_WORLD, &k);

  MPI_Cart_create(MPI_COMM_WORLD, 1, &k, &period, reorder, &ring_comm);
  

  card_partial_data = n/k;
  data = (int *) malloc(n*sizeof(int));
  partial_data = (int *) malloc(card_partial_data*sizeof(int));
  

  MPI_Comm_rank(ring_comm, &rank);

  MPI_Cart_coords(ring_comm, rank, 1, &cord);

  if (cord == 0) {
    srand(NULL);
    for (i = 0; i < n; i++) {
      data[i] = rand();
    }
  }

  root_cord = 0;
  MPI_Cart_rank(ring_comm,&root_cord,&root_rank);
  MPI_Scatter(data, card_partial_data, MPI_INT, partial_data, card_partial_data, MPI_INT,root_rank,ring_comm);

  printf("After scatter\n");
  for ( i = 0 ; i < card_partial_data; i++ )
    partial_sum += partial_data[i];

  MPI_Reduce(&partial_sum, &result, 1, MPI_INT, MPI_SUM, root_rank, ring_comm);

  if ( cord == 0 )
    printf("Sum: %d\n", result);

  final_time = MPI_Wtime();
  printf("Proc: %d, %f\n", cord, final_time - initial_time);

  free(data);
  free(partial_data);
  MPI_Finalize();

  return 0;

}
