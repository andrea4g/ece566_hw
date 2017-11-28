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

  if ( p == NULL ) 
    return 1;

  free(p->nodes);
  free(p->visited);
  free(p);

  return 0;

}


Path init_path(int num_nodes, int est_tour_cost) {

  Path p;
  int i;

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


int add_node_path(Path p, int node) {

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

int extract_first_node(Path p) {

  return p->dim?p->nodes[0]:-1;

}


char* serialize_path(Path p) {

  char* buffer;
  int index;
  int n;

  buffer = malloc((4 + 2*(p->dim))*sizeof(int));
  index = 0;
  n = p->max_dim;

  memcpy(&buffer[index], &(p->max_dim),       sizeof(int));
  index = index + sizeof(int);
  memcpy(&buffer[index], &(p->est_tour_cost), sizeof(int));
  index = index + sizeof(int);
  memcpy(&buffer[index], &(p->act_tour_cost), sizeof(int));
  index = index + sizeof(int);
  memcpy(&buffer[index], &(p->dim),           sizeof(int));
  index = index + sizeof(int);
  memcpy(&buffer[index], p->nodes,            n*sizeof(int));
  index = index + n*sizeof(int);
  memcpy(&buffer[index], p->visited,          n*sizeof(int));

  return buffer;
}



Path deserialize_path(char* buffer) {

  Path p;
  int index;
  int n;

  p = malloc(sizeof(struct path));
  
  index = 0;
  memcpy( &(p->max_dim),&buffer[index],       sizeof(int));
  n = p->max_dim;
  p->nodes    = malloc(n*sizeof(int));
  p->visited  = malloc(n*sizeof(int));
  index = index + sizeof(int);
  memcpy( &(p->est_tour_cost),&buffer[index], sizeof(int));
  index = index + sizeof(int);
  memcpy( &(p->act_tour_cost),&buffer[index], sizeof(int));
  index = index + sizeof(int);
  memcpy( &(p->dim), &buffer[index],          sizeof(int));
  index = index + sizeof(int);
  memcpy( p->nodes,&buffer[index],            n*sizeof(int));
  index = index + n*sizeof(int);
  memcpy( p->visited, &buffer[index],         n*sizeof(int));

  return p;


}

int get_dimension_path(Path p) {

  if ( p == NULL ) {
    return -1;
  }

  return p->dim;

}


Path copy_path(Path p) {

  Path new_p;
  int i;

  new_p = init_path(p->max_dim,0); 
  
  new_p->est_tour_cost = p->est_tour_cost;
  new_p->act_tour_cost = p->act_tour_cost;
  new_p->dim = p->dim;

  for ( i = 0; i < p->dim; i++ ) {
    new_p->nodes[i]   = p->nodes[i];
  }
  for ( i = 0; i < p->max_dim; i++ ) {
    new_p->visited[i] = p->visited[i];
  }


  return new_p;


}




