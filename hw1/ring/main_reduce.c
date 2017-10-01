#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {

  // Communicator declaration
  MPI_Comm ring_comm;
  // variable declarations
  int period = 1;
  int reorder = 1;
  int cord;
  int k;
  int rank;
  int n;
  int card_partial_data;
  int i;
  int root_cord,root_rank;
  int result;
  int partial_sum = 0;
  // poiner declarations
  int* data;
  int* partial_data;
  double initial_time, final_time;

  // save the number of elements of the array
  n = atoi(argv[1]);
  // inizialization of MPI environment
  MPI_Init(&argc,&argv);
  // save the start time of the task
  initial_time = MPI_Wtime();
  // save in k the size of the virtual topology (#processors)
  MPI_Comm_size(MPI_COMM_WORLD, &k);
  // create the topology (1 dimension, k processors in the dimension,
  // wraparound yes, reorder yes)
  MPI_Cart_create(MPI_COMM_WORLD, 1, &k, &period, reorder, &ring_comm);
  // save the number of elements to send to each processor
  card_partial_data = n / k;

  // allocate the data of the starting array (sender buffer)
  data = (int*) malloc(n * sizeof(int));
  // allocate the data of the receiver buffer
  partial_data = (int*) malloc(card_partial_data * sizeof(int));

  // save in rank the rank of each processor (assigned with MPI_Cart_create)
  MPI_Comm_rank(ring_comm, &rank);
  // save in cord the coordinates given a rank
  MPI_Cart_coords(ring_comm, rank, 1, &cord);

  // if it is the source processor
  if (cord == 0) {
    // declare the seed
    srand(NULL);
    // fill the array data with random numbers
    for (i = 0; i < n; i++) {
      data[i] = rand();
    }
  }
  // initialize root_cord
  root_cord = 0;
  // given the root_cord, save in root_rank the rank (save in root_rank the
  // rank of the root processor)
  MPI_Cart_rank(ring_comm, &root_cord, &root_rank);
  // send the data with a one-to-all personalized broadcast communication
  // the root send n/k elements to each i-th processor
  MPI_Scatter(data, card_partial_data, MPI_INT, partial_data, card_partial_data, MPI_INT,root_rank,ring_comm);

  // sum the data received with the scatter
  for (i = 0 ; i < card_partial_data; i++)
    partial_sum += partial_data[i];

  // apply reduce operation on the root node and save the sum in result
  MPI_Reduce(&partial_sum, &result, 1, MPI_INT, MPI_SUM, root_rank, ring_comm);

  // if it is the root processor
  if (cord == 0)
    printf("Sum: %d\n", result);

  // save final time
  final_time = MPI_Wtime();
  // print coordinates of each processor and the time it took to complete the task
  printf("Proc: %d, %f\n", cord, final_time - initial_time);

  // free the memory allocated with the malloc
  free(data);
  free(partial_data);

  // terminate the MPI environment execution
  MPI_Finalize();

  return 0;
}
