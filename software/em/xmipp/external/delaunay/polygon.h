#ifndef POLYGON_H
#define POLYGON_H

#include "point.h"


/**************************************************************************
* Defines section
**************************************************************************/
#define N_POINTS_TRIANGLE	3


/**************************************************************************
* Data types definition
**************************************************************************/
// Triangle definition.
struct Triangle_T
{
	struct Point_T vertex[N_POINTS_TRIANGLE];       // Array of three points.
};


/*****************************************************************************
* Public functions declaration
*****************************************************************************/
// Lines and segments API.
void middle_Point(struct Point_T *p, struct Point_T *q, struct Point_T *middle);
void extend_Segment( struct Point_T *p, struct Point_T *q, struct Point_T *r);

// Triangles API.
struct Point_T get_Centre( struct Triangle_T *triang);

#endif
