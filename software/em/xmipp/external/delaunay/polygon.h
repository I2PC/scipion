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
// Directions of a line or segment (in degrees).
enum    Segment_Direction { HORIZONTAL_0,
                            VERTICAL_90,
                            HORIZONTAL_180,
                            VERTICAL_270,
                            FROM_0_TO_90,
                            FROM_90_TO_180,
                            FROM_180_TO_270,
                            FROM_270_TO_360};


// Segment definition Y= m*x + n.
struct Segment_T
{
    TYPE   	m;                              // Segment slope (m).
    TYPE   	n;                              // n value.
    enum    Segment_Direction direction;    // Segment direction.
};


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
