#include "header/stack.h"

typedef struct node* link;

struct stack{
  link tos;   // point to the top of the stack
  int count;  // Number of elements of the stack
};

struct node {
  Path p;
  int cost;
  link next_node;
};


link create_new_node(Path p, int c);

  
Stack init_stack() {
  Stack s;
  s = malloc(sizeof(struct stack));
  if ( s == NULL )
    return NULL;

  s->count = 0;
  s->tos = NULL;
  return s;
}

link create_new_node(Path p, int c) {
  link l;

  l = malloc(sizeof(struct node));
  
  if ( l == NULL )
    return NULL;

  l->p = p;
  l->cost = c;
  return l;
}


int push(Stack s, Path p,int c) {

  link l;

  if ( s == NULL ) 
    return 1;

  l = create_new_node(p,c);
  if ( l == NULL )
    return 1;

  l->next_node = s->tos;
  s->tos = l;
  s->count = s->count + 1;

  return 0;

}

int stack_empty(Stack s) {
  return (s->count == 0);
}

int pop(Stack s, Path* p, int* c) {

  link l;

  if ( stack_empty(s) )
    return 1;

  l = s->tos;
  s->tos = l->next_node;
  s->count--;

  *p = l->p;
  *c = l->cost;

  free(l);

  return 0;
}


int finalize_stack(Stack s) {

  int i;
  link l,prev_l;

  if ( s == NULL )
    return 0;

  l = s->tos;
  for ( i = 0; i < s->count; i++ ) {
    finalize_path(l->p);
    prev_l = l;
    l = l->next_node;
    free(prev_l);
  }

  free(s);

  return 0;
}


int insert_stack(Stack src, Stack dest) {

  link l;

  if ( src == NULL )
    return 1;
  if ( src->count == 0 )
    return 0;

  dest->count = src->count + dest->count;

  l = src->tos;

  while( l->next_node != NULL )
    l = l->next_node;

  l->next_node = dest->tos;
  dest->tos = src->tos;

  free(src);
  return 0;
}

Stack split_stack(Stack org_s) {

  Stack new_s;
  int n_skip_nodes,i;
  link l;

  new_s = init_stack();
  
  n_skip_nodes = org_s->count*2/3;
  l = org_s->tos;

  for ( i = 0; i < n_skip_nodes - 1; i++ ) {
    l = l->next_node;
  }

  new_s->tos = l->next_node;
  l->next_node = NULL;

  new_s->count = org_s->count - n_skip_nodes;
  org_s->count = n_skip_nodes;

  return new_s;

}




