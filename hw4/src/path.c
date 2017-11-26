#include "header/path.h"


struct path {

  int max_dim;
  int* nodes;
  int* visited;
  int est_tour_cost;
  int act_tour_cost;
  int dim;

};

int finalize_path(Path p) {

  int error;
  error = 1;

  if ( p == NULL ) 
    return error;

  error += 2*free(p->nodes);
  error += 4*free(p->visited);
  error += 8*free(p);

  return error;

}


Path init_path(int num_nodes, int est_tour_cost) {

  Path p;

  p = malloc(sizeof(struct path));

  p->max_dim = num_nodes;
  p->dim = 0;
  p->nodes   = malloc(p->max_dim*sizeof(int));
  p->visited = malloc(p->max_dim*sizeof(int));
  for ( i = 0; i < num_nodes; i++ ) {
    p->visited[i] = 0;
  }
  p->est_tour_cost = est_tour_cost;
  p->act_tour_cost = 0;

  return p;

}

int get_est_tour_cost(Path p) {
  return p->est_tour_cost;
}

int get_act_tour_cost(Path p) {
  return p->act_tour_cost;
}


int set_est_tour_cost(Path p, int est_tour_cost) {
  
  if ( p == NULL )
    return 1;
  p->est_tour_cost = est_tour_cost;
  return 0;

}


int set_act_tour_cost(Path p, int act_tour_cost) {

  if ( p == NULL )
    return 1;
  p->act_tour_cost = act_tour_cost;
  return 0;

}


int add_node(Path p, int node) {

  if ( p->dim == p->max_dim )
    return 1;

  p->nodes[p->dim] = node;
  p->visited[node] = 1;
  p->dim++;

  return 1;

}

void print_path(Path p) {

  int i = 0;

  if ( p == NULL )
    return;

  for ( i = 0; i < p->dim; i++ ) {
    printf("%d, ", p->nodes[i]);
  }
  printf("\n");

  return;

}

int visited(Path p, int node) {

  return p->visited[node];

}

int extract_last_node(Path p) {

  return p->dim?p->nodes[p->dim - 1]:-1;

}

int extract_fist_node(Path p) {

  return p->dim?p->nodes[0]:-1;

}

