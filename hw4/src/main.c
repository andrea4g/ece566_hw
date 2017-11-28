#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include "header/path.h"
#include "header/stack.h"
#include "header/graph.h"


#define N_ITERATIONS      1
#define DEBUG             1

#define ROOT_RANK         0
#define HOMETOWN          0

#define TERMINATION       0
#define REQUEST_WORK      1
#define REQUEST_REJECTED  2
#define REQUEST_ACCEPTED  3
#define PBSC              4

#define TOK_COLOR_WHITE   0
#define TOK_COLOR_BLACK   1
#define TOK_COLOR_GREEN   2

#define EXP_THR 3


/*----------------------------------------------------------------------------*/
/*-------------------------------TYPES DEFINITION-----------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*-------------------------------FUNCTION PROTOTYPES--------------------------*/
/*----------------------------------------------------------------------------*/
int cost( Graph g, Path p,int cost_edge, int current_node, 
          int act_best_sol_cost, int* new_act_cost, int* new_est_cost);
int** read_matrix_from_file(char* filename,int* n);
int tsp_best_solution(Graph g, Stack s, int p, int my_rank, int n);
int verify_request(Stack s, int n);
void send_request_work(int my_rank, int p);
void rcv_pbsc(int* act_best_sol_cost_ptr);
void serve_pendant_requests(Stack s, int* proc_color, int n, int my_rank);
int terminate(int n,int p, int my_rank, int* proc_color_ptr,
              int* act_best_sol_cost_ptr, Stack s);
int check_termination(int p, int my_rank, int* proc_color_ptr);
int work(Graph g,int n, Stack s, int act_best_sol_cost, Path* best_tour_ptr);
void broadcast_act_best_sol_cost(int my_rank, int p, int act_best_sol_cost);
/*----------------------------------------------------------------------------*/
/*------------------------------------MAIN------------------------------------*/
/*----------------------------------------------------------------------------*/

int main(int argc, char** argv) {

  // variable declaration
  int my_rank;
  int procs_number;
  int i,j,iteration;
  double time_vector[N_ITERATIONS],deviation;
  double average_time, final_time, initial_time;
  int n;
  int** adj_matrix;
  int* buffer;
  int est_tour_cost;
  Stack s;
  Graph g;
  Path p;


  

  // Initialize MPI environment
  MPI_Init(&argc, &argv);
  // save the number of processors
  MPI_Comm_size(MPI_COMM_WORLD, &procs_number);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


  // if it is the root processor
  if (my_rank == ROOT_RANK) {
    adj_matrix = read_matrix_from_file(argv[1],&n);
    buffer = malloc(n*n*sizeof(int));
    for ( i = 0; i < n; i++ ) {
      for (j = 0; j < n; j++ ) {
        buffer[i*n + j] = adj_matrix[i][j];
      }
    }
  }

  MPI_Bcast(&n,1,MPI_INT, ROOT_RANK,MPI_COMM_WORLD);
  if ( my_rank != ROOT_RANK) {
    buffer = malloc(n*n*sizeof(int));
    adj_matrix = malloc(n*sizeof(int *));
    for ( i = 0; i < n; i++ ) {
      adj_matrix[i] = malloc(n*sizeof(int));
    }
  }

  MPI_Bcast(buffer, n*n, MPI_INT, ROOT_RANK,MPI_COMM_WORLD);

  for ( i = 0; i < n; i++ ) {
    for (j = 0; j < n; j++ ) {
      adj_matrix[i][j] = buffer[i*n + j];
    }
  }

  g = init_graph(n,adj_matrix);
  s = init_stack();

  average_time = 0;
  for ( iteration = 0; iteration < N_ITERATIONS; iteration++ ) {
    // save initial time of the task
    initial_time = MPI_Wtime();
    if ( my_rank == ROOT_RANK ) {
      est_tour_cost = 0;
      for ( i = 0; i < n; i++ ) {
        est_tour_cost += (min_edge(g,i) + sec_min_edge(g,i))/2;
      }
      p = init_path(n,est_tour_cost);
      add_node_path(p, HOMETOWN);
      push(s,p);
    }
    tsp_best_solution(g,s,procs_number,my_rank,n);
  }
  average_time = average_time/N_ITERATIONS;

  if ( my_rank == ROOT_RANK ) {

    deviation = 0;
    for ( i = 0; i < N_ITERATIONS; i++ ) {

      deviation += (time_vector[i] - average_time)*(time_vector[i] - average_time);
    }
    // compute and print the rank of the processor and the time it took to complete the task
    printf("%f, %f\n", average_time, deviation);
  }

  // close the MPI environment
  MPI_Finalize();

  return 0;
}

/*----------------------------------------------------------------------------*/
/*-------------------------------FUNCTIONS------------------------------------*/
/*----------------------------------------------------------------------------*/
int tsp_best_solution(Graph g, Stack s, int p, int my_rank, int n) {

  int new_act_best_sol_cost, act_best_sol_cost;
  int proc_color;
  int flag;
  int first_time, first_color, first_rank_dst;
  Path best_tour;
  MPI_Request req;

  best_tour = NULL;
  first_time = 1;
  
  
  proc_color  = TOK_COLOR_WHITE;
  first_color = TOK_COLOR_WHITE; 
  first_rank_dst = (ROOT_RANK + 1) % p;

  flag = 1;
  act_best_sol_cost = -1;

  while( flag ) {
    if ( stack_empty(s) ) {
      if ( my_rank == ROOT_RANK && first_time == 1 ) { 
        MPI_Isend(&first_color,
                  1,
                  MPI_INT,
                  first_rank_dst,
                  TERMINATION,
                  MPI_COMM_WORLD,
                  &req);
        first_time = 0;
      }
      flag = !terminate(n,p,my_rank,&proc_color,&act_best_sol_cost,s);
    } else {
      new_act_best_sol_cost = work(g,n,s,act_best_sol_cost,&best_tour);
      if ( new_act_best_sol_cost != act_best_sol_cost ) {
        act_best_sol_cost = new_act_best_sol_cost;
        broadcast_act_best_sol_cost(my_rank,p,act_best_sol_cost);
      }
      serve_pendant_requests(s, &proc_color, n, my_rank);
    }
  }

  //cleanup_messages();

  return act_best_sol_cost;

}


void serve_pendant_requests(Stack s, int* proc_color, int n, int my_rank) {

  MPI_Status status;
  int flag;
  int trash;
  int work_amount;
  int dim_buffer;
  int requester_rank;
  Stack work_to_send;
  char* buffer;


  MPI_Iprobe(MPI_ANY_SOURCE,
            REQUEST_WORK,
            MPI_COMM_WORLD,
            &flag,
            &status);

  while ( flag ) {
    requester_rank = status.MPI_SOURCE;
    MPI_Recv(&trash, 0, MPI_INT, requester_rank, REQUEST_WORK,
             MPI_COMM_WORLD, &status);
    work_amount = dimension_stack(s);
    if ( work_amount >= 3 ) {
      work_to_send = split_stack(s);
      buffer = serialize_stack(s, n, &dim_buffer);
      MPI_Send( &buffer, dim_buffer, MPI_BYTE, requester_rank,
                REQUEST_ACCEPTED, MPI_COMM_WORLD);
      if ( requester_rank < my_rank ) {
        *proc_color = TOK_COLOR_BLACK;
      }
    } else {
      MPI_Send(&trash, 0, MPI_BYTE, requester_rank,
                REQUEST_REJECTED, MPI_COMM_WORLD);
    }
    MPI_Iprobe(MPI_ANY_SOURCE,
              REQUEST_WORK,
              MPI_COMM_WORLD,
              &flag,
              &status);
  }

  return;

}


void broadcast_act_best_sol_cost(int my_rank, int p, int act_best_sol_cost) {

  int i;
  MPI_Request req;

  for ( i = 0; i < p; i++ ) {
    if ( i != my_rank ) {
      MPI_Isend(&act_best_sol_cost, 1, MPI_INT, i, 
      PBSC, MPI_COMM_WORLD, &req);
    }
  }

  return;

}


int** read_matrix_from_file(char* filename,int* n) {

  FILE* fid;
  int num_nodes;
  int** mat;
  int i,j;

  fid = fopen(filename,"r");

  fscanf(fid,"%d",&num_nodes);
  mat = malloc(num_nodes*sizeof(int*));
  for ( i = 0; i < num_nodes; i++ ) {
    mat[i] = malloc(num_nodes*sizeof(int));
    for ( j = 0; j < num_nodes; j++ ) {
      fscanf(fid,"%d", &mat[i][j]);
    }
  }

  *n = num_nodes;
  return mat;

}


int work(Graph g,int n,Stack s, int act_best_sol_cost, Path* best_tour_ptr) {

  Path p, new_p;
  int new_act_cost, new_est_cost;
  int current_node,first_node;
  int i, count;
  int edge_cost;

  count = 0;
  while ( !stack_empty(s) && count < EXP_THR ) {
    pop(s,&p);
    current_node = extract_last_node(p);
    if ( get_dimension_path(p) == n ) {
      first_node = extract_first_node(p);
      edge_cost = get_edge_cost(g,current_node,first_node);
      if ( edge_cost != -1 ) {
        new_act_cost = get_act_tour_cost(p) + edge_cost;
        if (  act_best_sol_cost == -1 || 
              new_act_cost < act_best_sol_cost ) {
          finalize_path(*best_tour_ptr);
          *best_tour_ptr = copy_path(p);
          act_best_sol_cost = new_act_cost;
        }
      }
    } else {
      for ( i = 0; i < n; i++ ) {
        edge_cost = get_edge_cost(g,current_node,i);
        if ( ( edge_cost > -1 ) &&
              !visited(p,i) ) {
          if ( cost(g,p,edge_cost,i,
                    act_best_sol_cost,&new_act_cost,&new_est_cost) ) {
            new_p = copy_path(p);
            add_node_path(new_p,i);
            set_act_tour_cost(new_p,new_act_cost);
            set_est_tour_cost(new_p,new_est_cost);
            push(s,new_p);
          }
        }
      }
    }
    count++;
    finalize_path(p);
  }

  return act_best_sol_cost;

}

int cost(
    Graph g,
    Path p,
    int cost_edge,
    int current_node,
    int act_best_sol_cost,
    int* new_act_cost,
    int* new_est_cost){

  int act_cost,est_cost;
  int last_node;
  int total_cost;

  act_cost = get_act_tour_cost(p);
  est_cost = get_est_tour_cost(p);
  last_node = extract_last_node(p);

  act_cost = act_cost + cost_edge;

  if ( get_dimension_path(p) == 1 ) {
    est_cost =
      est_cost - (min_edge(g,last_node) + min_edge(g,current_node))/2;
  } else {
    est_cost =
      est_cost - (sec_min_edge(g,last_node) + min_edge(g,current_node))/2;
  }

  total_cost = est_cost + act_cost;

  if (  (act_best_sol_cost > -1) &&
        (total_cost >= act_best_sol_cost ) )
    return 0;

  *new_est_cost = est_cost;
  *new_act_cost = act_cost;

  return 1;
}

int check_termination(int p, int my_rank, int* proc_color_ptr){

  int flag, value;
  int proc_color;
  MPI_Request req;

  int rank_src = (my_rank - 1) % p;
  int rank_dst = (my_rank + 1) % p;
  proc_color = *proc_color_ptr;


  MPI_Iprobe(rank_src,
            TERMINATION,
            MPI_COMM_WORLD,
            &flag,
            MPI_STATUS_IGNORE);

  if ( flag ) {
    MPI_Recv( &value,
              1,
              MPI_INT,
              rank_src,
              TERMINATION,
              MPI_COMM_WORLD,
              MPI_STATUS_IGNORE);
    if( value == TOK_COLOR_GREEN ) {
      return 1;
    } else {
      if ( proc_color == TOK_COLOR_WHITE ) {
        if ( my_rank == ROOT_RANK ) {
          value = TOK_COLOR_GREEN;
        }
        MPI_Isend(&value,
                  1,
                  MPI_INT,
                  rank_dst,
                  TERMINATION,
                  MPI_COMM_WORLD,
                  &req);
      } else {
        value = TOK_COLOR_BLACK;
        MPI_Isend(&value,
                  1,
                  MPI_INT,
                  rank_dst,
                  TERMINATION,
                  MPI_COMM_WORLD,
                  NULL);
        *proc_color_ptr = TOK_COLOR_WHITE;
      }
    }
  }

  return 0;
}


int terminate(int n, int p, int my_rank, int* proc_color_ptr,
              int* act_best_sol_cost_ptr, Stack s) {

  int request_on_fly = 0;
  int response;

  while ( 1 ) {
    if ( check_termination(p,my_rank,proc_color_ptr) )
      return 1;
    rcv_pbsc(act_best_sol_cost_ptr);
    if ( !request_on_fly ) {
      send_request_work(my_rank,p);
      request_on_fly = 1;
    } else {
      response = verify_request(s,n); // Returns a value < 0 in case until now
                                    // no response received.
      if ( response >= 0   ) {  // Request Rejected or Accepted
        request_on_fly = 0;
        if ( response == 0 )    // Request Accepted
          return 0;
      }
    }
  }

  return 255; // Never here
}


void rcv_pbsc(int* act_best_sol_cost_ptr) {

  MPI_Status status;
  int flag;
  int poss_best_sol_cost;
  int act_best_sol_cost;

  act_best_sol_cost = *act_best_sol_cost_ptr;

  MPI_Iprobe(MPI_ANY_SOURCE,
            PBSC,
            MPI_COMM_WORLD,
            &flag,
            &status);

  while ( flag ) {
    MPI_Recv(&poss_best_sol_cost, 1, MPI_INT, status.MPI_SOURCE, PBSC,
             MPI_COMM_WORLD, &status);
    if ( poss_best_sol_cost < act_best_sol_cost ) {
      act_best_sol_cost = poss_best_sol_cost;
    }
    MPI_Iprobe(MPI_ANY_SOURCE,
              PBSC,
              MPI_COMM_WORLD,
              &flag,
              &status);
  }

  *act_best_sol_cost_ptr = act_best_sol_cost;
  return;

}


void send_request_work(int my_rank, int p) {

  int dest_rank;
  MPI_Request req;

  srand(time(NULL));
  dest_rank = (my_rank + (rand() % (p-1)) + 1) % p;

  MPI_Isend(&dest_rank,
            1,
            MPI_INT,
            dest_rank,
            REQUEST_WORK,
            MPI_COMM_WORLD,
            &req);

}
/*  < 0 in case no request
 *  = 0 in case accepted
 *  > 0 in case rejected
*/
int verify_request(Stack s, int n){

  int trash;
  int flag;
  char* buffer;
  int count;
  MPI_Status status;
  Stack new_s;

  MPI_Iprobe(MPI_ANY_SOURCE,
            REQUEST_REJECTED,
            MPI_COMM_WORLD,
            &flag,
            &status);

  if ( flag ) {
    MPI_Recv(&trash, 1, MPI_INT, status.MPI_SOURCE, REQUEST_REJECTED,
             MPI_COMM_WORLD, &status);
    return 1;
  }

  MPI_Iprobe(MPI_ANY_SOURCE,
            REQUEST_ACCEPTED,
            MPI_COMM_WORLD,
            &flag,
            &status);

  if ( flag ) {
    MPI_Get_count(&status, MPI_BYTE, &count );
    buffer = malloc(count*sizeof(char));
    MPI_Recv(buffer, count, MPI_BYTE, status.MPI_SOURCE, REQUEST_ACCEPTED,
             MPI_COMM_WORLD, &status);
    new_s = deserialize_stack(buffer, n);
    insert_stack(new_s,s);
    free(buffer);
    return 1;
  }


  return -1;

}
/*
void cleanup_messages(){
 
  int flag;
  MPI_Status status;

  MPI_Iprobe(MPI_ANY_SOURCE,
            MPI_ANY_TAG,
            MPI_COMM_WORLD,
            &flag,
            &status);

  while ( flag ) {

    MPI_Recv(&poss_best_sol_cost, 1, MPI_INT, MPI_ANY, PBSC,
             MPI_COMM_WORLD, &status);
    if ( poss_best_sol_cost < act_best_sol_cost ) {
      act_best_sol_cost = poss_best_sol_cost;
    }
    MPI_Iprobe(MPI_ANY_SOURCE,
              PBSC,
              MPI_COMM_WORLD,
              &flag,
              &status);
  }
}*/













