#ifndef PATH_H
#define PATH_H
#include <stdlib.h>
#include <stdio.h>

typedef struct path* Path;

Path init_path(int num_nodes, int est_tour_cost);
int get_est_tour_cost(Path p);
int get_act_tour_cost(Path p);
int set_est_tour_cost(Path p, int est_tour_cost);
int set_act_tour_cost(Path p, int act_tour_cost);
int add_node(Path p, int node);
void print_path(Path p);


#endif
