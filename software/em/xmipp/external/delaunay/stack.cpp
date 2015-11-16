#include "defines.h"
#ifdef LOGGING
#include "log.h"
#endif
#include "stack.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int initialize_Stack(struct stack_T *stack, int size)
{
    int     ret=SUCCESS;        // Return value.

	// Initialize stack attributes.
	stack->size = size;
	stack->nElements = 0;
	stack->data = (struct Convex_T *) malloc(sizeof(struct Convex_T)*size);

	// Check error allocating memory.
	if (stack->data == NULL)
	{
	    ret = FAILURE;
	    printf("Error allocating memory in initialize_Stack");
#ifdef LOGGING
		sprintf( log_Text, "Error allocating memory in initialize_Stack\n");
		write_Log( log_Text);
#endif
	}

	return(ret);
}

struct Convex_T		peak(struct stack_T *stack)
{
	// Get peak element of the stack.
	return(stack->data[stack->nElements-1]);
}

struct Convex_T		peak_Prev(struct stack_T *stack)
{
	// Get the top 2 elements of the stack.
	return(stack->data[stack->nElements-2]);
}

void	update_Peak(struct stack_T *stack, int edge)
{
	// Update edge of peak of the stack.
	stack->data[stack->nElements-1].edge_ID	= edge;
}

void	push(struct stack_T *stack, struct Convex_T edge)
{
	// Insert a new element in the top of the stack.
	stack->data[stack->nElements] = edge;
	stack->nElements++;
}


void	pop(struct stack_T *stack)
{
	// Check if stack is not empty.
	if (stack->nElements > 0)
	{
		// Decrease number of elements.
		stack->nElements--;
	}
}

int		get_Stack_nElements(struct stack_T *stack)
{
	// Return number of elements.
	return(stack->nElements);
}

void	reset_Stack(struct stack_T *stack)
{
	// Reset data of the stack (without resizing).
	stack->nElements = 0;
	memset( stack->data, 0, sizeof(struct Convex_T)*stack->size);
}

void	print_Stack(struct stack_T *stack)
{
	int	i=0;					// Loop counter.

	// Print all elements in stack.
	for (i=0; i<stack->nElements ;i++)
	{
		printf("Position: %d. Edge id: %d with origin vertex %d.\n",
											i,
											stack->data[i].edge_ID,
											stack->data[i].vertex_Index);
	}
}

void	finalize_Stack(struct stack_T *stack)
{
	// Deallocate memory.
	stack->size = 0;
	stack->nElements = 0;
	if (stack->data != NULL)
	{
		free(stack->data);
		stack->data = NULL;
	}
}
