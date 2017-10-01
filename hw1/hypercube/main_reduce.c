#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {

  // variable declaration
  MPI_Comm hc_comm;
  int reorder = 1;
  int root_rank;
  int world_size;
  int rank;
  int i;
  int result;
  int partial_sum;
  unsigned int n;
  unsigned int temp;
  unsigned int card_partial_data;
  unsigned int d;
  char flag;
  double final_time, initial_time;
  // pointer declaration
  int* root_cord;
  int* dim;
  int* period;
  int* cord;
  int* data;
  int* partial_data;

  // save in n the length of the array
  n = atoi(argv[1]);

  // Initialize MPI environment
  MPI_Init(&argc, &argv);

  // save initial time of the task
  initial_time = MPI_Wtime();

  // save the number of processors
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // save the number of elements that have to be delivered to each processor
  card_partial_data = n / world_size;
  // initialize d to 0
  d = 0;
  temp = world_size;
  // compute the log2 of word_size and save the value in d
  while ( !(temp & 0x01) ) {
    d++;
    temp = temp >> 1;
  }

  // allocate the memory
  dim          = (int*) malloc(d * sizeof(int));
  root_cord    = (int*) malloc(d * sizeof(int));
  cord         = (int*) malloc(d * sizeof(int));
  period       = (int*) malloc(d * sizeof(int));
  partial_data = (int*) malloc(card_partial_data * sizeof(int));

  // fill the dim array with 2 (2 processor per dimension)
  // each dimension has the wraparound
  for (i = 0; i < d; i++) {
    dim[i] = 2;
    period[i] = 1;
  }

  // create a virtual topoly for the hypercube of d-dimensions
  MPI_Cart_create(MPI_COMM_WORLD, d, dim, period, reorder, &hc_comm);
  // save the rank of each processor in rank
  MPI_Comm_rank(hc_comm, &rank);
  // save the coordinates of each given the rank of the processor
  MPI_Cart_coords(hc_comm, rank, d, cord);

  // initialize a flag for the root process (root processor - flag = 1)
  flag = 1;
  for (i = 0; i < d; i++) {
    root_cord[i] = 0;
    // if it not the root processor the flag will be = 0
    if (cord[i] == 1)
      flag = 0;
  }
  // if it is the root processor
  if (flag) {
    // allocate the data
    data = (int*) malloc(n * sizeof(int));
    // declare the seed
    srand(NULL);
    // fill data with random values
    for (i = 0; i < n; i++) {
      data[i] = rand();
    }
  }
  // get the rank of each processor in the topology
  MPI_Cart_rank(hc_comm, root_cord, &root_rank);
  // scatter the data from souce to all the processors
  MPI_Scatter(data, card_partial_data, MPI_INT, partial_data, card_partial_data, MPI_INT, root_rank, hc_comm);
  // initialize the partial sum of each processor
  partial_sum = 0;
  // compute the sum of all the values received with the scatter
  for ( i = 0 ; i < card_partial_data; i++ )
    partial_sum += partial_data[i];

  // apply reduce operation (MPI_SUM) on the root processor
  MPI_Reduce(&partial_sum, &result, 1, MPI_INT, MPI_SUM, root_rank, hc_comm);

  // if it is the root processor print the sum
  if (flag){
    printf("Sum: %d\n", result);
    // and free the data array dinamically allocated
    free(data);
  }
  // save final time of the task
  final_time = MPI_Wtime();
  // compute and print the rank of the processor and the time it took to complete the task
  printf("Proc: %d, time: %f\n", rank, final_time-initial_time);

  // free the dynamic memory allocated
  free(dim);
  free(root_cord);
  free(cord);
  free(period);
  free(partial_data);

  // close the MPI environment
  MPI_Finalize();

  return 0;
}
