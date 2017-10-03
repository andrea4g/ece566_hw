#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N_ITERATIONS 1000
#define DEBUG 0


int main(int argc, char** argv) {

  // variable declaration
  MPI_Comm hc_comm;
  int reorder = 1;
  int root_rank;
  int n,k;
  int my_rank;
  int i,iteration;
  int result,partial_sum;
  int temp;
  int card_partial_data,reminder;
  int d;
  char flag;
  double time_vector[N_ITERATIONS],deviation;
  double average_time, final_time, initial_time;
  // pointer declaration
  int* root_partial_sum;
  int* root_cord;
  int* dim;
  int* period;
  int* my_cord;
  int* data;
  int* partial_data;
  int* sendcounts;
  int* displs;

  // save in n the length of the array
  n = atoi(argv[1]);

  // Initialize MPI environment
  MPI_Init(&argc, &argv);
  // save the number of processors
  MPI_Comm_size(MPI_COMM_WORLD, &k);
  // save the number of elements that have to be delivered to each processor
  card_partial_data = n / k;
  reminder = n % k;
  // initialize d to 0
  d = 0;
  temp = k;
  // compute the log2 of word_size and save the value in d
  while ( !(temp & 0x01) ) {
    d++;
    temp = temp >> 1;
  }

  // allocate the memory
  dim          = (int*) malloc(d * sizeof(int));
  root_cord    = (int*) malloc(d * sizeof(int));
  my_cord      = (int*) malloc(d * sizeof(int));
  period       = (int*) malloc(d * sizeof(int));


  // fill the dim array with 2 (2 processor per dimension)
  // each dimension has the wraparound
  for (i = 0; i < d; i++) {
    dim[i] = 2;
    period[i] = 1;
  }

  // create a virtual topoly for the hypercube of d-dimensions
  MPI_Cart_create(MPI_COMM_WORLD, d, dim, period, reorder, &hc_comm);
  // save the rank of each processor in rank
  MPI_Comm_rank(hc_comm, &my_rank);
  // save the coordinates of each given the rank of the processor
  MPI_Cart_coords(hc_comm, my_rank, d, my_cord);

  // initialize a flag for the root process (root processor -> flag = 1)
  flag = 1;
  for (i = 0; i < d; i++) {
    root_cord[i] = 0;
    // if it not the root processor the flag will be = 0
    if (my_cord[i] == 1)
      flag = 0;
  }
  // if it is the root processor
  if (flag) {
    // allocate the data
    data = (int*) malloc(n * sizeof(int));
    // declare the seed
    srand((unsigned int) 0);
  #if DEBUG
    // fill data with consecutive values
    for (i = 0; i < n; i++) {
      data[i] = i;
    }
  #else
    // fill data with random values
    for (i = 0; i < n; i++) {
      data[i] = rand();
    }
  #endif
  }
  // get the rank of each processor in the topology
  MPI_Cart_rank(hc_comm, root_cord, &root_rank);
  
  // Allocate vectors for SCATTERV primitive
  sendcounts = (int*) malloc(k*sizeof(int));
  displs = (int*) malloc(k*sizeof(int));
  // Assign to the root the reminder
  sendcounts[root_rank] = reminder + card_partial_data;
  displs[root_rank] = 0;
  for ( i = root_rank + 1; i < k; i++ ) {
    sendcounts[i] = card_partial_data;
    displs[i] = (reminder) + (i-root_rank)*card_partial_data;
  }
  for ( i = 0; i < root_rank; i++ ) {
    sendcounts[i] = card_partial_data;
    displs[i] = (reminder) + (k - root_rank + i)*card_partial_data;
  }
  if ( flag ) {
    root_partial_sum = (int* ) malloc( k * sizeof(int));
    partial_data = (int*) malloc((card_partial_data+reminder) * sizeof(int));
  } else {
    partial_data = (int*) malloc(card_partial_data * sizeof(int));
  }

  average_time = 0;
  for ( iteration = 0; iteration < N_ITERATIONS; iteration++ ) {
    // save initial time of the task
    initial_time = MPI_Wtime();
    // scatter the data from source to all the processors
    MPI_Scatterv(data, sendcounts, displs, MPI_INT, partial_data, sendcounts[my_rank], MPI_INT, root_rank, hc_comm);
    // initialize the partial sum of each processor
    partial_sum = 0;
    // compute the sum of all the values received with the scatter
    for ( i = 0 ; i < sendcounts[my_rank]; i++ )
      partial_sum += partial_data[i];
    // apply reduce operation (MPI_SUM) on the root processor
    MPI_Gather(&partial_sum, 1, MPI_INT, root_partial_sum, 1, MPI_INT, root_rank, hc_comm); 

    if ( flag ) {
      result = 0;
      for ( i = 0 ; i < k; i++ ) {
        result += root_partial_sum[i];
      }
    }

    // save final time of the task
    final_time = MPI_Wtime();
    if (flag) {
    #if DEBUG
      printf("Sum: %d\n", result);
    #endif
      // and free the data array dinamically allocated
    }
    time_vector[iteration] = final_time - initial_time;
    average_time += time_vector[iteration];
  }
  average_time = average_time/N_ITERATIONS;




  if ( flag ) {
    deviation = 0;
    for ( i = 0; i < N_ITERATIONS; i++ ) {
      deviation += (time_vector[i] - average_time)*(time_vector[i] - average_time);
    }
    // compute and print the rank of the processor and the time it took to complete the task
    printf("Av_time: %f ,dev: %f\n", average_time, deviation);
    free(data);
  }

  // free the dynamic memory allocated
  free(dim);
  free(root_cord);
  free(my_cord);
  free(period);
  free(partial_data);

  // close the MPI environment
  MPI_Finalize();

  return 0;

}
