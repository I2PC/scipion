#include <math.h>
#include "point.h"
#include <stdio.h>


/*****************************************************************************
 * Defines section
 ****************************************************************************/
#define 	COLLINEAR_THRESHOLD		0.0

/*****************************************************************************
* Private functions declaration
*****************************************************************************/
TYPE 	Euclidean( struct Point_T *p, struct Point_T *q);


 /*****************************************************************************
 * Public functions declaration
 *****************************************************************************/
/*****************************************************************************
* Name: 		distance
* Input: 		point p point q func
* Description: distance from p to q.
* Output: 		computes the distance from p to q using a distance
* 				function passed as parameter.
* Complexity: 	O(1)
*****************************************************************************/
TYPE distance( struct Point_T *p, struct Point_T *q)
{
	TYPE dist=0.0;

	// Compute distance.
	dist = Euclidean( p, q);

	return(dist);
}


int	equal_Point(struct Point_T *p1, struct Point_T *p2)
{
    int     	equal=FALSE;        // Return value.
    POINT_T   	diff=0.0;           // Coordinates difference value.

    // Check if points are too close.
    diff = p1->x - p2->x;
    if ((diff < DELTA_DIFF) && (diff > -DELTA_DIFF))
    {
        diff = p1->y - p2->y;
        if ((diff < DELTA_DIFF) && (diff > -DELTA_DIFF))
        {
            equal = TRUE;
        }
    }
    else
    {
        equal = FALSE;
    }

	return(equal);
}

int	higher_Point(struct Point_T *p1, struct Point_T *p2, int (*f)(struct Point_T *, struct Point_T *))
{
    return((*f)( p1, p2));
}

void	print_Point(struct Point_T *point)
{
	printf("Point (%f, %f)\n", point->x, point->y);
}

//#define DEBUG_CHECK_TURN
enum Turn_T		check_Turn(struct Point_T *p1, struct Point_T *p2, struct Point_T *p3)
{
	TYPE          area=0.0;       // Signed area.
	enum Turn_T     turn;           // Return value.

	// Compute signed area of the triangle formed by p1, p2 and p3.
	area = signed_Area( p1, p2, p3);

#ifdef DEBUG_CHECK_TURN
	printf("Area %f. P (%lf,%lf). Q (%lf,%lf). R (%lf,%lf)\n", area, p1->x, p1->y, p2->x, p2->y, p3->x, p3->y);
#endif

	// Higher than zero -> turn left.
	if (area > COLLINEAR_THRESHOLD)
	{
#ifdef DEBUG_CHECK_TURN
		printf("LEFT\n");
#endif
		turn = LEFT_TURN;
	}
	// Lower than zero -> turn right.
	else if (area < -COLLINEAR_THRESHOLD)
	{
#ifdef DEBUG_CHECK_TURN
		printf("RIGHT\n");
#endif
		turn = RIGHT_TURN;
	}
	// If area is close to zero then points are collinear.
	else
	{
#ifdef DEBUG_CHECK_TURN
		printf("COLLINEAR\n");
#endif
		turn = COLINEAR;
	}

	return(turn);
}


int	in_Circle(struct Point_T *p1, struct Point_T *p2, struct Point_T *p3, struct Point_T *q)
{
	int			ret=0;			// Return value.
	POINT_T 	value=0.0;		// Determinant value.
	POINT_T   	temp[9];		// Intermediate values.

	// Compute Ax - Dx, Ay - Dy and (Ax-Dx)² + (Ay-Dy)²
	temp[0] = (p1->x - q->x);
	temp[1] = (p1->y - q->y);
	temp[2] = (POINT_T) (pow((p1->x - q->x), 2) + pow((p1->y - q->y), 2));

	// Compute Bx - Dx, By - Dy and (Bx-Dx)² + (By-Dy)²
	temp[3] = (p2->x - q->x);
	temp[4] = (p2->y - q->y);
	temp[5] = (POINT_T) (pow((p2->x - q->x), 2) + pow((p2->y - q->y), 2));

	// Compute Cx - Dx, Cy - Dy and (Cx-Dx)² + (Cy-Dy)²
	temp[6] = (p3->x - q->x);
	temp[7] = (p3->y - q->y);
	temp[8] = (POINT_T) (pow((p3->x - q->x), 2) + pow((p3->y - q->y), 2));

	// Compute determinant.
	value = (temp[0]*temp[4]*temp[8]) +
			(temp[1]*temp[5]*temp[6]) +
			(temp[2]*temp[3]*temp[7]) -
			(temp[2]*temp[4]*temp[6]) -
			(temp[5]*temp[7]*temp[0]) -
			(temp[8]*temp[1]*temp[3]);

	// If positive then point q belongs to p1-p2-p3 circumference.
	if (value > 0.0)
	{
		ret = 1;
	}

	return(ret);
}

bool	interior_Triangle(struct Point_T *p1, struct Point_T *p2, struct Point_T *p3, struct Point_T *q)
{
    bool is_In_Triangle=false;   // Return value.

	// Check if new point is interior to triangle formed by p1, p2 and p3.
	if ((signed_Area( p1, p2, q) > 0) &&
        (signed_Area( p2, p3, q) > 0) &&
        (signed_Area( p3, p1, q) > 0))
	{
	    is_In_Triangle = true;
	}
	else
	{
	    is_In_Triangle = false;
	}

	return(is_In_Triangle);
}


int	has_Extreme_Coordinates(struct Point_T *p)
{
	int has_Extreme=FALSE;		// Return value.

	if ((p->x >= MAX_X_COORD) || (p->y >= MAX_Y_COORD) ||
		(p->x <= -MAX_X_COORD) || (p->y <= -MAX_Y_COORD))
	{
		has_Extreme = TRUE;
	}

	return(has_Extreme);
}

int	lower_X(struct Point_T *p1, struct Point_T *p2)
{
	// Return true if y coordinate is higher.
	return(p1->x < p2->x);
}

int	higher_X(struct Point_T *p1, struct Point_T *p2)
{
	// Return true if y coordinate is higher.
	return(p1->x > p2->x);
}

int	higher_Y(struct Point_T *p1, struct Point_T *p2)
{
	// Return true if y coordinate is higher.
	return(p1->y > p2->y);
}

/********************************************************************/
/* Name: 		lexicographic_Higher
 * Input:
 * Description:
 * Output:
 * Complexity:	O(1)
 * */
/********************************************************************/
int lexicographic_Higher(struct Point_T *p1, struct Point_T *p2)
{
    int     higher=FALSE;       // Return value.

    // Check if Y coordinate is higher.
    if (p1->y > p2->y)
    {
        higher = TRUE;
    }
    else
    {
        // If Y coordinate is equal then check X coordinate.
        if ((p1->y == p2->y) &&
            (p1->x > p2->x))
        {
            higher = TRUE;
        }
    }

    return(higher);
}

/********************************************************************
 * Name: 		Euclidean
 * Input: 		*p1     reference to Point_T
 * Output:      *p2     reference to Point_T
 * Description: copies p1 into p2.
 * Output: 		p2.x = p1.x and p2.y = p1.y
 * Complexity:	O(1)
********************************************************************/
void copy_Point(struct Point_T *p1, struct Point_T *p2)
{
    // Copy p1 x and y values into p2.
    p2->x = p1->x;
    p2->y = p1->y;
}

/********************************************************************/
/* Name: 		Euclidean
 * Input: 		point p point q
 * Description: computes the euclidean distance from p to q.
 * Output: 		euclidean distance from point p to point q
 * Complexity:	O(1)
 * */
/********************************************************************/
TYPE Euclidean( struct Point_T *p, struct Point_T *q)
{
	TYPE dist=0.0;

	// Compute euclidean distance.
	dist = sqrt(pow( p->x - q->x, 2) + pow( p->y - q->y, 2));

	return(dist);
}


TYPE signed_Area(struct Point_T *p1, struct Point_T *p2, struct Point_T *p3)
{
	double	area=0.0;			// Return value.

	area = (- (double) p2->x*p1->y + (double) p3->x*p1->y + (double) p1->x*p2->y -
			(double) p3->x*p2->y - (double) p1->x*p3->y + (double) p2->x*p3->y);

	return(area);
}

