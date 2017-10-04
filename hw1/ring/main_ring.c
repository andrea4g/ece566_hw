#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


int main(int argc, char** argv) {

  MPI_Comm ring_comm;
  int mailbox;
  int period = 1;
  int reorder = 1;
  int my_rank, my_cord;
  int k,n;
  int* data;
  int* partial_data;
  int card_partial_data;
  int i;
  int root_cord,root_rank;
  int result;
  int partial_sum = 0;
  double initial_time, final_time;
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
  
  data = (int *)malloc(n*sizeof(int));
  if ( my_cord == 0 ) {
    for ( i = 0; i < n; i++ )
      data[i] = i;
  }


  if ( my_coord == 0 ) {
    MPI_Cart_shift(ring_comm, 0, 1, &my_rank, &dest_rank);
    MPI_Send(&data[card_partial_data + reminder], 
        (k-1)*card_partial_data, MPI_INT, dest_rank, 0, ring_comm);
    partial_sum = 0;
    for ( i = 0 ; i < card_partial_data + reminder; i++ ) {
      partial_sum += data[i];
    }
  } else {
    MPI_Cart_shift(ring_comm, 0, -1, &my_rank, &src_rank);
    MPI_Recv(partial_data, 
        (k-my_coord)*card_partial_data, MPI_INT, src_rank, 0, ring_comm, &status);
    if ( my_coord < k ) {
      MPI_Send(&partial_data[card_partial_data], 
          (k - my_cord - 1)*card_partial_data, MPI_INT, dest_rank, 0, ring_comm);
    }
    partial_sum = 0;
    for ( i = 0 ; i < card_partial_data; i++ ) {
      partial_sum += partial_data[i];
    }
  }

  MPI_Cart_rank(ring_comm,&root_cord,&root_rank);

  if ( my_rank != 0 && my_rank < N/2 ) {
    
  }




  if ( my_cord == 0 ) {
    printf("%d\n", partial_sum);
  }


  MPI_Finalize();

  return 0;

}
