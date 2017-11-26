#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>

#define N_ITERATIONS      1
#define DEBUG             1

#define ROOT_RANK         0

#define TERMINATION       0
#define REQUEST_WORK      1
#define REQUEST_REJECTED  2
#define REQUEST_ACCEPTED  3
#define PBSC              4

#define TOK_COLOR_WHITE   0
#define TOK_COLOR_RED     1
#define TOK_COLOR_GREEN   2




/*----------------------------------------------------------------------------*/
/*-------------------------------TYPES DEFINITION-----------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*-------------------------------FUNCTION PROTOTYPES--------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*------------------------------------MAIN------------------------------------*/
/*----------------------------------------------------------------------------*/

int main(int argc, char** argv) {

  // variable declaration
  int root_rank,my_rank;
  int p;
  int i,j,iteration;
  double time_vector[N_ITERATIONS],deviation;
  double average_time, final_time, initial_time;
  // pointer declaration
  int dims[2];
  Matrix A, C;

  // Initialize MPI environment
  MPI_Init(&argc, &argv);
  // save the number of processors
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  root_rank = 0;

  // if it is the root processor
  if (my_rank == root_rank) {
    
  }

  average_time = 0;
  for ( iteration = 0; iteration < N_ITERATIONS; iteration++ ) {
    // save initial time of the task
    initial_time = MPI_Wtime();
    
  }
  average_time = average_time/N_ITERATIONS;

  if ( my_rank == root_rank ) {

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
int Work(Graph g, Stack s, int act_best_sol_cost, Path* best_tour) {

  Path p;
  int new_act_cost, new_est_cost;
  int first_node;

  count = 0;
  while ( !empty_stack(s) && count < EXP_THR ) {
    p = pop(s);
    if ( get_dimension_path(p) == n ) {
      first_node = extract_first_node(p);
      edge_cost = get_edge_cost(g,current_node,first_node);
      new_act_cost = get_act_tour_cost(p) + edge_cost;
      if ( new_act_cost < act_best_sol_cost ) {
        finalize_path(*best_tour);
        *best_tour = p;
        act_best_sol_cost = new_act_cost;
      }
    } else {
      current_node = extract_last_node(p);
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
      count++;
    }
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

  if ( get_dim_path(p) == 1 ) {
    est_cost =
      est_cost - (min_edge(G,last_node) + min_edge(G,current_node))/2;
  } else {
    est_cost =
      est_cost - (sec_min_edge(G,last_node) + min_edge(G,current_node))/2;
  }

  total_cost = est_cost + act_cost;

  if (  (act_best_sol_cost > -1) &&
        (total_cost >= act_best_sol_cost ) )
    return 0;

  *new_est_cost = est_cost;
  *new_act_cost = act_cost;

  return 1;
}

int check_termination(){

  int flag, value;

  MPI_Probe(pb->rank - 1,
            TERMINATION,
            MPI_COMM_WORLD,
            &flag,
            MPI_STATUS_IGNORE);

  if ( flag ) {
    MPI_Recv( &value,
              1,
              MPI_INT,
              pb->rank - 1,
              TERMINATION,
              MPI_COMM_WORLD,
              MPI_STATUS_IGNORE);
    if( value == TOK_COLOR_GREEN ) {
      return 1;
    } else {
      if ( pb->col == TOK_COLOR_WHITE ) {
        if ( pb->rank == ROOT_RANK ) {
          value = TOK_COLOR_GREEN;
        }
        MPI_Isend(&value,
                  1,
                  MPI_INT,
                  pb->rank + 1,
                  TERMINATION,
                  MPI_COMM_WORLD,
                  NULL);
      } else {
        value = TOK_COLOR_BLACK;
        MPI_Isend(&value,
                  1,
                  MPI_INT,
                  pb->rank + 1,
                  TERMINATION,
                  MPI_COMM_WORLD,
                  NULL);
        pb->col = TOK_COLOR_WHITE;
      }
    }
  }

  return 0;
}


int terminate() {

  int request_on_fly = 0;

  while ( 1 ) {
    flag = check_termination();
    if ( flag )
      return 1;
    rcv_pbsc();
    if ( !request_on_fly ) {
      send_request_work();
      request_on_fly = 1;
    } else {
      verify_request();
      if ( response == 1  ) {
        insert_stack();
        return 0;
      } else {
        request_on_fly = 0;
      }
    }
  }

  return 255; // Never here
}























