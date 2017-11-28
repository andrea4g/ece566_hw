#include <stdio.h>
#include <stdlib.h>
#include "header/stack.h"
#include "header/path.h"


int main(int argc, char** argv) {

  Stack s, new;
  Path p;
  int dimension;
  char* buffer;

  s = init_stack();
  p = init_path(5,0);
  add_node_path(p,0);
  add_node_path(p,3);
  add_node_path(p,2);
  push(s,p);
  p = copy_path(p);
  add_node_path(p,4);
  push(s,p);
  p = copy_path(p);
  add_node_path(p,1);
  push(s,p);
  p = copy_path(p);
  push(s,p);

  buffer = serialize_stack(s,5,&dimension);
  new = deserialize_stack(buffer,5);

  new = split_stack(s);



  return 0;

}
