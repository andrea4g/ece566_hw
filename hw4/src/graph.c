#include "header/graph.h"


struct graph {

  int num_nodes;
  int** adj_matrix;
  int* min_edge;
  int* sec_min_edge;

};

int finalize_graph(Graph g) {

  int i;

  if ( g == NULL )
    return 1;

  free(g->min_edge);
  free(g->sec_min_edge);

  for ( i = 0; i < g->num_nodes; i++ ) {
    free(g->adj_matrix[i]);
  }
  free(g->adj_matrix);
  free(g);

  return 0;

}


Graph init_graph(int num_nodes, int** adj_matrix) {

  Graph g;
  int i,min;

  g = malloc(sizeof(struct graph));

  g->num_nodes = num_nodes;
  g->adj_matrix = malloc(num_nodes*sizeof(int*));
  g->min_edge = malloc(num_nodes*sizeof(int));
  g->sec_min_edge = malloc(num_nodes*sizeof(int));

  for ( i = 0; i < num_nodes; i++ ) {
    g->adj_matrix[i] = malloc(num_nodes*sizeof(int));
    for ( j = 0; j < num_nodes; j++ ) {
      g->adj_matrix[i][j] = adj_matrix[i][j];
    }
  }

  for ( i = 0; i < num_nodes; i++ ) {
    min = adj_matrix[i][0];
    for (j = 0; j < num_nodes; j++ ) {
      if ( adj_matrix[i][j] < min )
        min = adj_matrix[i][j];
    }
    g->min_edge[i] = min;
  }

  for ( i = 0; i < num_nodes; i++ ) {
    sec_min = -1;
    for (j = 0; j < num_nodes; j++ ) {
      if ((adj_matrix[i][j] > g->min_edge[i]) &&
          (sec_min == -1 || adj_matrix[i][j] < sec_min)) {
        sec_min = adj_matrix[i][j];
      }
    }
    g->sec_min_edge[i] = sec_min;
  }

  return g;

}

int get_edge_cost(Graph g, int node_src, int node_dst) {

  return g->adj_matrix[node_src][node_dst];

}


int min_edge(Graph g, int node) {

  return g->min_edge[node];
}

int sec_min_edge(Graph g,int node) {
  return g->sec_min_edge[node];
}
