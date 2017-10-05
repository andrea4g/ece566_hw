#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


int main(int argc, char** argv) {

  MPI_Comm ring_comm;
  unsigned int mailbox,reminder;
  int period = 1;
  int reorder = 1;
  int my_rank, my_coord;
  int k,n;
  unsigned int* data;
  unsigned int* partial_data;
  unsigned int card_partial_data;
  int i;
  int root_cord,root_rank;
  long unsigned int partial_sum = 0;
  double initial_time, final_time;
  int dst_coord,src_coord;
  int src_rank,dst_rank;
  MPI_Status status;


  n = atoi(argv[1]);                // Number of elements of the array A

  MPI_Init(&argc,&argv);

  initial_time = MPI_Wtime();

  MPI_Comm_size(MPI_COMM_WORLD, &k);

  card_partial_data = n/k;
  reminder = n % k;

  MPI_Cart_create(MPI_COMM_WORLD, 1, &k, &period, reorder, &ring_comm);
  
  MPI_Comm_rank(ring_comm, &my_rank);
  MPI_Cart_coords(ring_comm, my_rank, 1, &my_coord);
  

  if ( my_coord == 0 ) {
    data = (unsigned int *)malloc(n*sizeof(int));
    for ( i = 0; i < n; i++ )
      data[i] = i;
  } else {
    partial_data = (unsigned int *) malloc ( (k-my_coord)*card_partial_data*sizeof(int) );
  }

  MPI_Cart_shift(ring_comm, 0, +1, &my_rank, &dst_rank);
  MPI_Cart_shift(ring_comm, 0, -1, &my_rank, &src_rank);

  if ( my_coord == 0 ) {
    dst_coord = 1;
    MPI_Cart_rank(ring_comm, &dst_coord, &dst_rank);
    MPI_Send(&data[card_partial_data + reminder], 
        (k-1)*card_partial_data, MPI_INT, dst_rank, 0, ring_comm);
    partial_sum = 0;
    for ( i = 0 ; i < card_partial_data + reminder; i++ ) {
      partial_sum += data[i];
    }
  } else {
    src_coord = my_coord - 1;
    MPI_Cart_rank(ring_comm, &src_coord, &src_rank);
    MPI_Recv(partial_data, 
        (k-my_coord)*card_partial_data, MPI_INT, src_rank, 0, ring_comm, &status);
    if ( my_coord < (k-1) ) {  
      dst_coord = (my_coord + 1) % k;
      MPI_Cart_rank(ring_comm, &dst_coord, &dst_rank);
      MPI_Send(&partial_data[card_partial_data], 
          (k - my_coord - 1)*card_partial_data, MPI_INT, dst_rank, 0, ring_comm);
    }
    partial_sum = 0;
    for ( i = 0 ; i < card_partial_data; i++ ) {
      partial_sum += partial_data[i];
    }
  }

  if ( my_coord == k/2  ) {
    dst_coord = my_coord - 1;
    MPI_Cart_rank(ring_comm, &dst_coord, &dst_rank);
    MPI_Send(&partial_sum, 1, MPI_INT, dst_rank, 0, ring_comm);
  } else if ( my_coord == k/2 + 1 ) {
    dst_coord = my_coord + 1;
    MPI_Cart_rank(ring_comm, &dst_coord, &dst_rank);
    MPI_Send(&partial_sum, 1, MPI_INT, dst_rank, 0, ring_comm);
  } else if ( my_coord == 0 ) {
    src_coord = my_coord + 1;
    MPI_Cart_rank(ring_comm, &src_coord, &src_rank);
    MPI_Recv(&mailbox, 1, MPI_INT, src_rank, 0, ring_comm, &status);
    partial_sum += mailbox;
    src_coord = (my_coord - 1) % k;
    MPI_Cart_rank(ring_comm, &src_coord, &src_rank);
    MPI_Recv(&mailbox, 1, MPI_INT, src_rank, 0, ring_comm, &status);
    partial_sum += mailbox;
  } else {
    src_coord = my_coord < k/2 ? (my_coord + 1) % k : (my_coord - 1) % k;
    dst_coord = my_coord < k/2 ? (my_coord - 1) % k : (my_coord + 1) % k;
    MPI_Cart_rank(ring_comm, &src_coord, &src_rank);
    MPI_Cart_rank(ring_comm, &dst_coord, &dst_rank);
    MPI_Recv(&mailbox, 1, MPI_INT, src_rank, 0, ring_comm, &status);
    partial_sum += mailbox;
    MPI_Send(&partial_sum, 1, MPI_INT, dst_rank, 0, ring_comm);
  }

 if ( my_coord == 0 ) {
    printf("%lu\n", partial_sum);
  }

  MPI_Finalize();

  return 0;

}
