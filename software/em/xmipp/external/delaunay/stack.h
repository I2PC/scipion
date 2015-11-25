#ifndef STACK_H
#define STACK_H


//*****************************************************************************
//                           DATA TYPES DEFINITION
//*****************************************************************************
struct Convex_T
{
	int		vertex_Index;
	int		edge_ID;
};

struct stack_T
{
	int					size;
	int					nElements;
	struct Convex_T	*data;
};

int 	initialize_Stack(struct stack_T *stack, int size);
struct Convex_T		peak(struct stack_T *stack);
struct Convex_T		peak_Prev(struct stack_T *stack);
void	update_Peak(struct stack_T *stack, int edge);
void	push(struct stack_T *stack, struct Convex_T edge);
void	pop(struct stack_T *stack);
int		get_Stack_nElements(struct stack_T *stack);
void	reset_Stack(struct stack_T *stack);
void	print_Stack(struct stack_T *stack);
void	finalize_Stack(struct stack_T *stack);

#endif
