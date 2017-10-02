#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {

  // variable declaration
  MPI_Comm hc_comm;
  int acc;
  int i, index;
  int mask, msg_dest, msg_source, my_virtual_id;
  int partial_sum;
  int rank,result, root_rank;
  int reorder = 1;
  int virtual_source, virtual_dest;
  int world_size;
  unsigned int card_partial_data;
  unsigned int d;
  unsigned int n;
  unsigned int temp;
  char flag;
  double initial_time, final_time;
  // pointer declaration
  int* cord;
  int* data;
  int* dim;
  int* partial_data;
  int* period;
  int* root_cord;

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
  data         = (int*) malloc(n * sizeof(int));
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
    // declare the seed
    srand(NULL);
    // fill data with random values
    for (i = 0; i < n; i++) {
      //data[i] = rand();
      data[i] = i;
      //printf("data[i] = %d\n", data[i]);
    }
  }

  // get the rank of each processor in the topology
  MPI_Cart_rank(hc_comm, root_cord, &root_rank);
  // 1-to-all personalized from source to all the processors
  // implemented using send and receive
  my_virtual_id = rank ^ root_rank;
  mask = world_size - 1;
  index = n / 2;
  for (i = d - 1; i >= 0; i--) {
    mask = mask ^ (1 << i);
    if ((my_virtual_id & mask) == 0) {
      if ((my_virtual_id & (1 << i)) == 0) {
        virtual_dest = my_virtual_id ^ (1 << i);
        MPI_Send(&data[index], index, MPI_INT, (virtual_dest ^ root_rank), 0 , hc_comm);
      }
      else {
        virtual_source = my_virtual_id ^ (1 << i);
        MPI_Recv(data, index, MPI_INT, (virtual_source ^ root_rank), 0, hc_comm, MPI_STATUS_IGNORE);
      }
    }
    index = index >> 1;
  }
  // save the final value of the scatter inside partial_data
  for (i = 0; i < card_partial_data; i++) {
    partial_data[i] = data[i];
  }
  // initialize the partial sum of each processor
  partial_sum = 0;
  // compute the sum of all the values received with the scatter
  for (i = 0 ; i < card_partial_data; i++)
    partial_sum += partial_data[i];

  // apply reduce operation (MPI_SUM) on the root processor
  // implemented using send and receive
  mask = 0;
  for (i = 0; i < d; i++) {
    if ((rank & mask) == 0) {
      if ((rank & (1 << i)) != 0) {
        msg_dest = rank ^ (1 << i);
        MPI_Send(&partial_sum, 1, MPI_INT, msg_dest, 0, hc_comm);
      }
      else {
        msg_source = rank ^ (1 << i);
        MPI_Recv(&acc, 1, MPI_INT, msg_source, 0, hc_comm, MPI_STATUS_IGNORE);
        partial_sum += acc;
      }
    }
    mask = mask ^ (1 << i);
  }
  // save the final result of the reduction inside result
  result = partial_sum;

  // if it is the root processor print the sum
  if (flag){
    printf("Sum: %d\n", result);
  }
  // save final time of the task
  final_time = MPI_Wtime();
  // compute and print the rank of the processor and the time it took to complete the task
  printf("Proc: %d, time: %f\n", rank, final_time-initial_time);

  // free the dynamic memory allocated
  free(data);
  free(dim);
  free(root_cord);
  free(cord);
  free(period);
  free(partial_data);

  // close the MPI environment
  MPI_Finalize();

  return 0;
}
