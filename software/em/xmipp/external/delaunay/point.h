#ifndef POINT_H
#define POINT_H

#include "defines.h"

/**************************************************************************
* Data types definition
**************************************************************************/
// Type of turn of three points line.
enum Turn_T {COLINEAR=1, LEFT_TURN=2, RIGHT_TURN=3};

struct Point_T
{
	POINT_T		x;
	POINT_T		y;
};


/*****************************************************************************
* Public functions declaration
*****************************************************************************/
// Comparison interface.
int		    equal_Point(struct Point_T *p1, struct Point_T *p2);
int	        higher_Point(struct Point_T *p1, struct Point_T *p2, int (*f)(struct Point_T *, struct Point_T *));
int			lower_X(struct Point_T *p1, struct Point_T *p2);
int			higher_X(struct Point_T *p1, struct Point_T *p2);
int			higher_Y(struct Point_T *p1, struct Point_T *p2);
int         lexicographic_Higher(struct Point_T *p1, struct Point_T *p2);
void        copy_Point(struct Point_T *p1, struct Point_T *p2);

// Distance.
TYPE 	    distance( struct Point_T *p, struct Point_T *q);
TYPE		signed_Area(struct Point_T *p1, struct Point_T *p2, struct Point_T *p3);

// Location interface.
enum Turn_T	check_Turn(struct Point_T *p1, struct Point_T *p2, struct Point_T *p3);
int			in_Circle(struct Point_T *p1, struct Point_T *p2, struct Point_T *p3, struct Point_T *q);
bool		interior_Triangle(struct Point_T *p1, struct Point_T *p2, struct Point_T *p3, struct Point_T *q);
int			has_Extreme_Coordinates(struct Point_T *p);

// I/O interface
void	    print_Point(struct Point_T *point);

#endif
