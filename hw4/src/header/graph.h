#ifndef GRAPH_H
#define GRAPH_H
#include <stdlib.h>
#include <stdio.h>

typedef struct graph* Graph;

Graph init_graph(int num_nodes, int** adj_matrix);
int finalize_graph(Graph g);
int get_edge_cost(Graph g, int node_src, int node_dst);




#endif
