#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include "dcel.h"
#include "stack.h"

/**************************************************************************
* Data types definition
**************************************************************************/
struct Triangulation_T
{
	struct stack_T 	stack;
	struct DCEL_T	*dcel;
};


/**************************************************************************
* Public functions declaration
**************************************************************************/
int		initialize_Triangulation(struct Triangulation_T *triang, struct DCEL_T *dcel);
void	finalize_Triangulation(struct Triangulation_T *triang);

void 	build_Triangulation( struct Triangulation_T *triang);


#endif
