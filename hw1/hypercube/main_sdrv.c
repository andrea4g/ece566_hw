#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define DEBUG 0 
#define N_ITERATIONS 1000


int main(int argc, char** argv) {

  // variable declaration
  MPI_Comm hc_comm;
  int acc;
  int i, index;
  int temp1;
  int k;
  int mask, msg_dest, msg_source;
  int my_cord;
  int partial_sum;
  int rank,result, root_rank;
  int reorder = 1;
  int virtual_dest_cord, virtual_source_cord;
  int world_size;
  int actual_rank;
  int save;
  int iteration;
  unsigned int card_partial_data;
  unsigned int d;
  unsigned int bit;
  unsigned int n;
  unsigned int temp;
  char flag;
  double time_vector[N_ITERATIONS],deviation;
  double average_time, final_time, initial_time;
  // pointer declaration
  int* cord;
  int* c;
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
  c            = (int*) malloc(d * sizeof(int));
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
    srand((unsigned) 0);
    // fill data with random values
    for (i = 0; i < n; i++) {
      data[i] = rand();
    }
  }
  my_cord = 0;
  for (i = d - 1; i >= 0; i--) {
    save = 0;
    save = cord[i] << i;
    my_cord += save;
  }

  // get the rank of the root in the topology
  MPI_Cart_rank(hc_comm, root_cord, &root_rank);
  // 1-to-all personalized from source to all the processors
  // implemented using send and receive
 
  average_time = 0;
  for ( iteration = 0; iteration < N_ITERATIONS; iteration++ ) { 
    
    initial_time = MPI_Wtime();
    mask = world_size - 1;
    index = n / 2;
    for (i = d - 1; i >= 0; i--) {
      mask = mask ^ (1 << i);
      if ((my_cord & mask) == 0) {
        if ((my_cord & (1 << i)) == 0) {
          virtual_dest_cord = my_cord ^ (1 << i);
          temp1 = virtual_dest_cord;
          for (k = 0; k < d; k++) {
            bit = virtual_dest_cord & 0x01;
            c[k] = bit;
            virtual_dest_cord = virtual_dest_cord >> 1;
          }
          MPI_Cart_rank(hc_comm, c, &actual_rank);
          MPI_Send(&data[index], index, MPI_INT, actual_rank, 0 , hc_comm);
          }
        else {
          virtual_source_cord = my_cord ^ (1 << i);
          temp1 = virtual_source_cord;
          for (k = 0; k < d; k++) {
            bit = virtual_source_cord  & 0x01;
            c[k] = bit;
            virtual_source_cord = virtual_source_cord >> 1;
          }
          MPI_Cart_rank(hc_comm, c, &actual_rank);
          MPI_Recv(data, index, MPI_INT, actual_rank, 0, hc_comm, MPI_STATUS_IGNORE);
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
      if ((my_cord & mask) == 0) {
        if ((my_cord & (1 << i)) != 0) {
          msg_dest = my_cord ^ (1 << i);
          for (k = 0; k < d; k++) {
            bit =  msg_dest & 0x01;
            c[k] = bit;
            msg_dest = msg_dest >> 1;
          }
          MPI_Cart_rank(hc_comm, c, &actual_rank);
          MPI_Send(&partial_sum, 1, MPI_INT, actual_rank, 0, hc_comm);
        }
        else {
          msg_source = my_cord ^ (1 << i);
          for (k = 0; k < d; k++) {
            bit =  msg_source & 0x01;
            c[k] = bit;
            msg_source = msg_source >> 1;
          }
          MPI_Cart_rank(hc_comm, c, &actual_rank);
          MPI_Recv(&acc, 1, MPI_INT, actual_rank, 0, hc_comm, MPI_STATUS_IGNORE);
          partial_sum += acc;
        }
      }
      mask = mask ^ (1 << i);
    }
    // save the final result of the reduction inside result
    result = partial_sum;

    final_time = MPI_Wtime();
    // if it is the root processor print the sum
#if DEBUG
    if (flag){
      printf("Sum: %d\n", result);
    }
#endif
    
    // save final time of the task
    final_time = MPI_Wtime();
    time_vector[iteration] = final_time - initial_time;
    average_time += time_vector[iteration];
    // compute and print the rank of the processor and the time it took to complete the task 
  } 
  average_time = average_time/N_ITERATIONS;
  
  if ( flag ) {
    deviation = 0;
    for ( i = 0; i < N_ITERATIONS; i++ ) {
      deviation += (time_vector[i] - average_time)*(time_vector[i] - average_time);
    }
    // compute and print the rank of the processor and the time it took to complete the task
    printf("Av_time: %f ,dev: %f\n", average_time, deviation);
  }

  // free the dynamic memory allocated
  free(cord);
  free(c);
  free(data);
  free(dim);
  free(partial_data);
  free(period);
  free(root_cord);

  // close the MPI environment
  MPI_Finalize();

  return 0;
}
