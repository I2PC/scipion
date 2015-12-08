#include "defines.h"
#include "float.h"
#include <math.h>
#include "polygon.h"
#include <stdio.h>
#include <stdlib.h>

#define MIN_SLOPE	0.0005

/**************************************************************************
* Private function headers
**************************************************************************/
int    compute_Segment( struct Point_T *p, struct Point_T *q, struct Segment_T *out);
enum Segment_Direction compute_Direction( struct Segment_T *segment, struct Point_T *p, struct Point_T *q);

/**************************************************************************
* Public function bodies
**************************************************************************/
void        middle_Point(struct Point_T *p, struct Point_T *q, struct Point_T *middle)
{
    // Compute middle point of edge.
    middle->x = (p->x + q->x) / (TYPE) 2.0;
    middle->y = (p->y + q->y) / (TYPE) 2.0;
}


//#define DEBUG_EXTEND_SEGMENT
/********************************************************************
FUNCTION:   extend_Segment
INPUT:      p   Point_T
            q   Point_T
OUTPUT:     r   Point_T
PRE:        p != NULL, q != NULL and r != NULL
POST:
DESCRIPTION:Computes r using p and q to compute the m and n of a line
            (y=mx + n). Once it has computed m and n, the uses them to
            compute r that is far away in the direction of pq line.
COMPLEXITY: O(1)
*********************************************************************/
void extend_Segment( struct Point_T *p, struct Point_T *q, struct Point_T *r)
{
	int		error=SUCCESS;
    struct Segment_T segment;         // Segment to store m and n values.
    enum Segment_Direction direction; // pq line direction.

    // Compute m and of pq line.
#ifdef DEBUG_EXTEND_SEGMENT
    printf("Px %f Py %f Qx %f Qy %f\n", p->x, p->y, q->x, q->y);
#endif
    error = compute_Segment( p, q, &segment);

    if (error == SUCCESS)
    {
		// Compute direction of pq line.
		direction = compute_Direction( &segment, p, q);

		// Line towards positive x coordinates values.
		if ((direction == FROM_0_TO_90) || (direction == FROM_270_TO_360))
		{
			r->x = MAX_X_COORD;

			// Set y value y=mx + n
			r->y = segment.m*r->x + segment.n;
		}
		// Line towards negative x coordinates values.
		else if ((direction == FROM_90_TO_180) || (direction == FROM_180_TO_270))
		{
			r->x = -MAX_X_COORD;

			// Set y value y=mx + n
			r->y = segment.m*r->x + segment.n;
		}
		// Line is parallel to X coordinates axis towards right.
		else if (direction == HORIZONTAL_0)
		{
			r->x = MAX_X_COORD;
			r->y = p->y;
		}
		// Line is parallel to X coordinates axis towards left.
		else if (direction == HORIZONTAL_180)
		{
			r->x = -MAX_X_COORD;
			r->y = p->y;
		}
		// Line is parallel to Y coordinates axis upwards.
		else if (direction == VERTICAL_90)
		{
			r->x = p->x;
			r->y = MAX_Y_COORD;
		}
		// Line is parallel to Y coordinates axis downwards.
		else
		{
			r->x = p->x;
			r->y = -MAX_Y_COORD;
		}
    }
    else
    {
    	if (error == -1)
    	{
    		r->x = p->x;
    		r->y = -MAX_Y_COORD;
    	}
    	else
    	{
    		r->x = p->x;
    		r->y = MAX_Y_COORD;
    	}
    }
}



//#define DEBUG_GET_CENTRE
/***************************************************************************
* Name: get_Centre
* IN:	triang			triangle whose circumcentre must be computed
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	triangle circumcentre
* Description: 	Computes the circumcentre of the triangle formed by the
* 				three points of the griangle
***************************************************************************/
struct Point_T 	get_Centre( struct Triangle_T *triang)
{
	struct 	Point_T point;					// Return value.
	TYPE 	x1=0.0, y1=0.0;
	TYPE  x2=0.0, y2=0.0;
	TYPE 	slope1=0.0, slope2=0.0;
	TYPE  n1=0.0, n2=0.0;
	int		valid1=0, valid2=0;

#ifdef DEBUG_GET_CENTRE
	//print_Triangle( triang);
#endif

	// Get point between two triangle vertex.
	x1 = (triang->vertex[0].x + triang->vertex[1].x) / 2.0;
	y1 = (triang->vertex[0].y + triang->vertex[1].y) / 2.0;

	// Compute slope.
	if ((triang->vertex[0].x - triang->vertex[1].x) != 0)
	{
		slope1 = (triang->vertex[0].y - triang->vertex[1].y) /
				(triang->vertex[0].x - triang->vertex[1].x);

#ifdef DEBUG_GET_CENTRE
		printf("Slope between points 0 and 1 is %lf.\n", slope1);
#endif
		if (slope1 == 0)
		{
			valid1 = 0;
		}
		else
		{
			valid1 = 1;
			slope1 = -(1/slope1);

			// Compute n1.
			n1 = y1 - slope1*x1;
		}
	}
	else
	{
		valid1 = 1;
		slope1 = 0;
		n1 = y1;
	}

	// Get point between two triangle vertex.
	x2 = (triang->vertex[1].x + triang->vertex[2].x) / 2.0;
	y2 = (triang->vertex[1].y + triang->vertex[2].y) / 2.0;

	// Compute slope.
	if ((triang->vertex[1].x - triang->vertex[2].x) != 0)
	{
		slope2 = (triang->vertex[1].y - triang->vertex[2].y) /
				(triang->vertex[1].x - triang->vertex[2].x);

#ifdef DEBUG_GET_CENTRE
		printf("Slope between points 1 and 2 is %lf.\n", slope2);
#endif
		if (slope2 == 0)
		{
			valid2 = 0;
		}
		else
		{
			valid2 = 1;
			slope2 = -(1/slope2);

			// Compute n2.
			n2 = y2 - slope2*x2;
		}
	}
	else
	{
		valid2 = 1;
		slope2 = 0;
		n2 = y2;
	}

	if (valid1)
	{
		if (valid2)
		{
			// Compute x.
			point.x = (n2 - n1) / (slope1 - slope2);

			// Compute y.
			point.y = point.x*slope1 + n1;
		}
		else
		{
			point.x = x2;
			point.y = point.x*slope1 + n1;
		}
	}
	else
	{
		point.x = x1;
		point.y = point.x*slope2 + n2;
	}

#ifdef DEBUG_GET_CENTRE
	printf("x %f y %f.\n", point.x, point.y);
#endif

	return(point);
}


//#define DEBUG_COMPUTE_SEGMENT
/********************************************************************
FUNCTION:   compute_Segment
INPUT:      p   Point_T
            q   Point_T
OUTPUT:     out Segment_T
PRE:        p != NULL, q != NULL
POST:       m = (q->y - p->y) / (q->x - p->x)
            n = p->y - (out->m*p->x)
DESCRIPTION:Computes the values m and n of the pq line (y=mx + n).
COMPLEXITY: O(1)
*********************************************************************/
int    compute_Segment( struct Point_T *p, struct Point_T *q, struct Segment_T *out)
{
	int		ret=SUCCESS;		// Return value.

    // Avoid division by 0.
    if ((q->x - p->x) != 0)
    {
        // Compute slope.
        out->m = (q->y - p->y) / (q->x - p->x);

        // Compute n.
        out->n = p->y - (out->m*p->x);
    }
    // m is +/-(infinite) but must be converted into a valid value.
    else
    {
    	if (q->y < p->y)
    	{
    		out->m = -FLT_MAX;
    		ret = -1;
    	}
    	else
    	{
    		out->m = -FLT_MAX;
    		ret = -2;
    	}
#ifdef DEBUG_COMPUTE_SEGMENT
    	printf("Slope is +/- inifinite and set to %f\n", out->m);
#endif
    }

    return(ret);
}


enum Segment_Direction compute_Direction( struct Segment_T *segment, struct Point_T *p, struct Point_T *q)
{
    enum Segment_Direction direction=FROM_0_TO_90;       // Return value.

    // Only first or third cuadrant.
    if (segment->m > 0)
    {
        // First cuadrant.
        if (p->y < q->y)
        {
            direction = FROM_0_TO_90;
        }
        // Third cuadrant.
        else
        {
            direction = FROM_180_TO_270;
        }
    }
    // Second or fourht cuadrant.
    else if (segment->m < 0)
    {
        // First cuadrant.
        if (p->y > q->y)
        {
            direction = FROM_270_TO_360;
        }
        // Third cuadrant.
        else
        {
            direction = FROM_90_TO_180;
        }
    }
    // Slope is zero -> line parallel to X axis.
    else if (segment->m == 0)
    {
        // Paralle to X axis to the right.
        if (p->x < q->x)
        {
            direction = HORIZONTAL_0;
        }
        // Paralle to X axis to the left.
        else
        {
            direction = HORIZONTAL_180;
        }
    }
    // Slope is infinite -> line parallel to Y axis.
    else
    {
        // Paralle to Y axis upwards.
        if (p->y < q->y)
        {
            direction = VERTICAL_90;
        }
        // Paralle to Y axis downwards.
        else
        {
            direction = VERTICAL_270;
        }
    }

    return(direction);
}



