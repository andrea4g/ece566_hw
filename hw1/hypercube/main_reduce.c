#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {

  MPI_Comm hc_comm;
  unsigned int d;
  int* dim;
  int* period;
  int reorder = 1;
  int* cord;
  int world_size;
  int rank;
  int* data;
  int* partial_data;
  unsigned int card_partial_data;
  int i;
  int result;
  int partial_sum;
  double final_time, initial_time;
  char flag;
  unsigned int n;
  unsigned int temp;
  int* root_cord;
  int root_rank;


  n = atoi(argv[1]);
  
  MPI_Init(&argc,&argv);
  initial_time = MPI_Wtime();

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  card_partial_data = n/world_size;
  printf("%d\n", card_partial_data);
  d = 0;
  temp = world_size;
  while ( !(temp & 0x01) ) {
    d++;
    temp = temp >> 1;
  } 
  printf("%d\n", n);

  dim = (int* ) malloc(d*sizeof(int));
  root_cord = (int* ) malloc(d*sizeof(int));
  cord = (int* ) malloc(d*sizeof(int));
  period = (int* ) malloc(d*sizeof(int));
  partial_data = (int* ) malloc(card_partial_data*sizeof(int));

  printf("%x\n", partial_data);

  for ( i = 0; i < d; i++ ) {
    dim[i] = 2;
    period[i] = 1;
  }
  
  MPI_Cart_create(MPI_COMM_WORLD, d, dim, period, reorder, &hc_comm);

  MPI_Comm_rank(hc_comm, &rank);
  MPI_Cart_coords(hc_comm, rank, d, cord);

  flag = 1;
  for ( i = 0; i < d; i++ ) {
    root_cord[i] = 0;
    if ( cord[i] == 1 )
      flag = 0;
  }


  if (flag) {
    data = (int* ) malloc(n*sizeof(int));
    srand(NULL);
    for (i = 0; i < n; i++) {
      data[i] = rand();
    }
  }

  MPI_Cart_rank(hc_comm,root_cord,&root_rank);

  printf("%d\n",root_rank);
  MPI_Scatter(data, card_partial_data, MPI_INT, partial_data, card_partial_data, MPI_INT, root_rank, hc_comm);

  printf("ciaoneeee\n");

  partial_sum = 0;
  for ( i = 0 ; i < card_partial_data; i++ )
    partial_sum += partial_data[i];

  MPI_Reduce(&partial_sum, &result, 1, MPI_INT, MPI_SUM, root_rank, hc_comm);

  if ( flag )
    printf("Sum: %d\n", result);

  final_time = MPI_Wtime();
  printf("Proc: %d, time: %f\n", rank, final_time-initial_time);

  free(dim);
  free(root_cord);
  free(cord);
  free(period);
  free(data);
  free(partial_data);

  MPI_Finalize();

  return 0;

}
