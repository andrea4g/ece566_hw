#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N_ITERATIONS 1000
#define DEBUG 1


int main(int argc, char** argv) {

  // variable declaration
  MPI_Comm ring_comm;
  int reorder = 1;
  int root_rank,root_cord;
  int n,k;
  int my_rank;
  int i,iteration;
  int result,partial_sum;
  int card_partial_data,reminder;
  double time_vector[N_ITERATIONS],deviation;
  double average_time, final_time, initial_time;
  // pointer declaration
  int period = 1;
  int my_cord;
  int* root_partial_sum;
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

  // create a virtual topology for the ring of k-elements
  MPI_Cart_create(MPI_COMM_WORLD, 1, &k, &period, reorder, &ring_comm);
  // save the rank of each processor in rank
  MPI_Comm_rank(ring_comm, &my_rank);
  // save the coordinates of each given the rank of the processor
  MPI_Cart_coords(ring_comm, my_rank, 1, &my_cord);

  root_cord = 0;
  // if it is the root processor
  if ( my_cord == root_cord) {
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
  MPI_Cart_rank(ring_comm, &root_cord, &root_rank);
  
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
  if ( my_rank == root_rank ) {
    root_partial_sum = (int *) malloc(k * sizeof(int));
    partial_data = (int*) malloc((card_partial_data+reminder) * sizeof(int));
  } else {
    partial_data = (int*) malloc(card_partial_data * sizeof(int));
  }

  average_time = 0;
  for ( iteration = 0; iteration < N_ITERATIONS; iteration++ ) {
    // save initial time of the task
    initial_time = MPI_Wtime();
    // scatter the data from source to all the processors
    MPI_Scatterv(data, sendcounts, displs, MPI_INT, partial_data, sendcounts[my_rank], MPI_INT, root_rank, ring_comm);
    // initialize the partial sum of each processor
    partial_sum = 0;
    // compute the sum of all the values received with the scatter
    for ( i = 0 ; i < sendcounts[my_rank]; i++ )
      partial_sum += partial_data[i];
    // apply reduce operation (MPI_SUM) on the root processor
    MPI_Gather(&partial_sum, 1, MPI_INT, root_partial_sum, 1, MPI_INT, 0, ring_comm);

    if ( my_rank == root_rank) {
      result = 0;
      for ( i = 0; i < k; i++ ) {
        result += root_partial_sum[i];
      }
    }

    // save final time of the task
    final_time = MPI_Wtime();
    if ( my_rank == root_rank) {
    #if DEBUG
      printf("Sum: %d\n", result);
    #endif
      // and free the data array dinamically allocated
    }
    time_vector[iteration] = final_time - initial_time;
    average_time += time_vector[iteration];
  }
  average_time = average_time/N_ITERATIONS;


  if ( my_rank == root_rank ) {
    deviation = 0;
    for ( i = 0; i < N_ITERATIONS; i++ ) {
      deviation += (time_vector[i] - average_time)*(time_vector[i] - average_time);
    }
    // compute and print the rank of the processor and the time it took to complete the task
    printf("Av_time: %f ,dev: %f\n", average_time, deviation);
    free(data);
  }

  // free the dynamic memory allocated
  free(partial_data);

  // close the MPI environment
  MPI_Finalize();

  return 0;

}
