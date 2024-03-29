#ifndef PATH_H
#define PATH_H
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct path* Path;

Path init_path(int num_nodes, int est_tour_cost);
int finalize_path(Path p);
int get_est_tour_cost(Path p);
int get_act_tour_cost(Path p);
int set_est_tour_cost(Path p, int est_tour_cost);
int set_act_tour_cost(Path p, int act_tour_cost);
int add_node_path(Path p, int node);
int extract_last_node(Path p);
int extract_first_node(Path p);
int visited(Path p, int node);
void print_path(Path p);
int get_dimension_path(Path p);
char* serialize_path(Path p);
Path deserialize_path(char* buffer);
Path copy_path(Path p);
#endif
