#include "sorting.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//#define DEBUG_SORTING

/****************************************************************************
                           PRIVATE FUNCTION HEADERS
*****************************************************************************/
int		move_Pivot( struct Dcel_Vertex_T *origin, struct Dcel_Vertex_T *list, int init_Index, int final_Index);
void	quicksort( struct Dcel_Vertex_T *origin, struct Dcel_Vertex_T *list, int init_Index, int final_Index);


/****************************************************************************
                           PUBLIC FUNCTION BODIES
*****************************************************************************/
void	sort(struct DCEL_T *dcel)
{
	// Call quicksort with DCEL structure.
	quicksort( &dcel->vertex[0], &dcel->vertex[1], 0, dcel->nVertex-2);
}

/***************************************************************************
* Name: clutter
* IN:		N/A
* OUT:		N/A
* IN/OUT:	dcel			dcel data
* RETURN:	N/A
* Description: clutters the set of point of the DCEL.
***************************************************************************/
void	clutter(struct DCEL_T *dcel)
{
	int	i=0;					// Loop counter.
	int	index1=0, index2=0;		// Array indexes.

	// Set seed.
	srand(time(NULL));

	// Loop set of points.
	for (i=0; i<dcel->nVertex ;i++)
	{
		// Generate random indexes to swap.
		index1 = rand() %  dcel->nVertex;
		index2 = rand() %  dcel->nVertex;

		// Swap elements.
		swap_Vertex( dcel, index1, index2);
	}
}

void    set_Highest_First(struct DCEL_T *dcel, int  (*f)(struct Point_T *, struct Point_T *))
{
    int     i=0;                // Loop counter.
    int     length=0;           // Number of vertexes in DCEL.
    int     highest_Index=0;    // Index of highest point.

    // Get number of vertexes.
    length = get_Number_Vertex( dcel);
    highest_Index=0;

    // Find index of highest point.
    for (i=1; i<length ;i++)
    {
        if ((*f) ( &dcel->vertex[i].vertex, &dcel->vertex[highest_Index].vertex))
        {
            highest_Index = i;
        }
    }

    // Swap highest vertex and vertex at position 0.
    swap_Vertex( dcel, 0, highest_Index);

#ifdef DEBUG_SORTING
    printf("Highest point at position %d.", highest_Index);
    print_Point( &dcel->vertex[0].vertex);
    fflush(NULL);
#endif
}

/****************************************************************************
                           PRIVATE FUNCTION BODIES
*****************************************************************************/
int		move_Pivot( struct Dcel_Vertex_T *origin, struct Dcel_Vertex_T *list, int init_Index, int final_Index)
{
	int i=0, first_Index=0;			// Loop variables.
    int pivot_Index=0;				// Pivot index.
	struct Dcel_Vertex_T pivot_Point;	// Pivot point.
    struct Dcel_Vertex_T temp;			// Register to exchange elements in list.

	// Get position and value of pivot point.
    pivot_Index = init_Index;
    pivot_Point = list[pivot_Index];

	// Set first index.
	first_Index = init_Index + 1;

	// Loop.
    for (i=first_Index; i<=final_Index; i++)
	{
		// If turn is right then angle is lower in list[i] than in pivot element.
		if ((check_Turn( &origin->vertex, &pivot_Point.vertex, &list[i].vertex) == RIGHT_TURN) ||
			((check_Turn( &origin->vertex, &pivot_Point.vertex, &list[i].vertex) == COLINEAR) &&
			(pivot_Point.vertex.x > list[i].vertex.x)))
		{
			// Increase pivot index.
			pivot_Index++;

			// Move current element to lower position in set.
            temp = list[i];
            list[i] = list[pivot_Index];
            list[pivot_Index] = temp;
		}
	}

	// Move pivot point.
    temp = list[init_Index];
    list[init_Index] = list[pivot_Index];
    list[pivot_Index] = temp;

	return pivot_Index;
}

void	quicksort( struct Dcel_Vertex_T *origin, struct Dcel_Vertex_T *list, int init_Index, int final_Index)
{
	int pivot_Index=0;			// Index of the element used as pivot.

	// Check if list is one element long.
    if(init_Index < final_Index)
	{
		// Move first element to its position.
		pivot_Index = move_Pivot( origin, list, init_Index, final_Index);

		// Order from initial element to previous element to pivot.
		quicksort( origin, list, init_Index, pivot_Index-1);

		// Order from next element to pivot until end of set.
		quicksort( origin, list, pivot_Index+1, final_Index);
	}
}



