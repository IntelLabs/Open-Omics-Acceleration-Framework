#ifndef _STACK_
#define _STACK_

typedef struct {
    int *base;
    int top;
    int size;
    FILE *trace;
} integer_stack_t;

typedef integer_stack_t * stack;

stack stack_create(const int size);
int   stack_pop(stack s);
void  stack_push(stack s, int i);
int   stack_size(const stack& s);
void  stack_trace(const stack& s, FILE *f);
int   stack_test(void);
int   stack_depth(const stack& s);

#endif
