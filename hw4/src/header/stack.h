#ifndef STACK_H
#define STACK_H
#include <stdlib.h>
#include "path.h"

typedef struct stack* Stack;

Stack init_stack();
int push(Stack s, Path p);
int stack_empty(Stack s);
int pop(Stack s, Path* p);
int finalize_stack(Stack s);
int insert_stack(Stack src, Stack dest);
Stack split_stack(Stack org_s);
int dimension_stack(Stack s);
char* serialize_stack(Stack s,int n);
Stack deserialize_stack(char* buffer,int n);

#endif
