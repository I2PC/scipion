#ifndef VORONOI_H
#define VORONOI_H

#include "dcel.h"
#include "delaunay.h"

/**************************************************************************
* Data types definition
**************************************************************************/
struct Voronoi_T
{
    int				computed;
    struct DCEL_T  	dcel;
};


/*****************************************************************************
* Public functions declaration
*****************************************************************************/
int  initialize_Voronoi( struct Voronoi_T *voronoi_Set, int n);
int	 get_Voronoi_Centre( struct Voronoi_T *voronoi, int face,
													struct Point_T *centre);
void build_Voronoi( struct Voronoi_T *voronoi_Set, struct DCEL_T *dcel);
int  select_Closest_Point( struct Delaunay_T *delaunay, struct Voronoi_T *voronoi,
					struct Point_T *p, struct Point_T *q, double *lowest_Distance);
void finalize_Voronoi( struct Voronoi_T *voronoi_Set);


#endif
