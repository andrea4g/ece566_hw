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
  int my_cord,sender_mask,number_of_steps,new_edge,edge;
  int dest_cord,dest_rank,error_value,source_cord;
  int temp;
  MPI_Status status;


  n = atoi(argv[1]);                // Number of elements of the array A

  MPI_Init(&argc,&argv);

  initial_time = MPI_Wtime();

  MPI_Comm_size(MPI_COMM_WORLD, &k);

  MPI_Cart_create(MPI_COMM_WORLD, 1, &k, &period, reorder, &ring_comm);
  
  MPI_Comm_rank(ring_comm, &rank);
  MPI_Cart_coords(ring_comm, rank, 1, &my_cord);
 
  number_of_steps = 0;
  while ( (1 << number_of_steps) < k ) {
    number_of_steps++;
  }  
  printf("%d\n",number_of_steps);
  
  data = (int *)malloc(n*sizeof(int));
  if ( my_cord == 0 ) {
    for ( i = 0; i < n; i++ )
      data[i] = i;
  }
 
  edge = n;
  sender_mask = (1 << number_of_steps) - 1;        // mask = 11...11_ceil(log2(k))
  for ( i = 0; i < number_of_steps; i++ ) {
    new_edge = edge >> 1;
    if ( (my_cord & sender_mask) == 0 ) {
      dest_cord = (1 << (number_of_steps - i - 1))^my_cord; 
      MPI_Cart_rank(ring_comm,&dest_cord,&dest_rank);
      printf("Proc: %d send to %d from %d to %d\n",my_cord,dest_cord,new_edge,edge-1); 
      error_value = 
        MPI_Send(&data[new_edge], edge - new_edge , MPI_INT, dest_cord, 0,ring_comm);
    }
    if (((my_cord & (sender_mask >> 1)) == 0) && (( (my_cord >> (number_of_steps - 1 - i )) & 0x01) == 1 )) {
      source_cord = my_cord^(1 << (number_of_steps - 1 - i));
      printf("Proc: %d rcv from %d\n",my_cord,source_cord); 
      error_value = 
        MPI_Recv(data, edge - new_edge, MPI_INT, source_cord, 0,ring_comm, &status);
    }
    if ( i != number_of_steps - 1) {
      edge = new_edge;
    }
    sender_mask = sender_mask >> 1;
  } 
 

  MPI_Cart_rank(ring_comm,&root_cord,&root_rank);
  
  
  MPI_Finalize();

  return 0;

}
