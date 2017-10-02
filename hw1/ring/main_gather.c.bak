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
  int total;
  int partial_sum = 0;
  double initial_time, final_time;
  // poiner declarations
  int* data;
  int* partial_data;
  int* result;

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

  // allocate the data of the receiver buffer of the scatter
  partial_data = (int*) malloc(card_partial_data * sizeof(int));
  // allocate the data of the receiver buffer of the gather
  result = (int*) malloc(k * sizeof(int));

  // save in rank the rank of each processor (assigned with MPI_Cart_create)
  MPI_Comm_rank(ring_comm, &rank);
  // save in cord the coordinates given a rank
  MPI_Cart_coords(ring_comm, rank, 1, &cord);

  // if it is the source processor
  if (cord == 0) {
    // allocate the data of the starting array (sender buffer)
    data = (int*) malloc(n * sizeof(int));
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

  // gather the data
  MPI_Gather(&partial_sum, 1, MPI_INT, result, 1, MPI_INT, 0, ring_comm);
  // if is is root processor
  if (rank == 0){
    // initialize total to 0
    total = 0;
    // sum all the data gathered in result by all processors
    for (i = 0; i < k; i++) {
      total = total + result[i];
    }
    printf("Sum: %d\n", total);
  }

  // save final time
  final_time = MPI_Wtime();
  // print coordinates of each processor and the time it took to complete the task
  printf("Proc: %d, %f\n", cord, final_time - initial_time);

  // free the memory allocated with the malloc
  if (cord == 0) {
    free(data);
  }
  free(partial_data);
  free(result);

  // terminate the MPI environment execution
  MPI_Finalize();

  return 0;
}
