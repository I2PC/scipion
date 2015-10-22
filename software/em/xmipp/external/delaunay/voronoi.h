#ifndef VORONOI_H
#define VORONOI_H

#include "dcel.h"

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
int  	initialize_Voronoi( struct Voronoi_T *voronoi_Set, int n);
int	 	get_Voronoi_Centre( struct Voronoi_T *voronoi, int face,
													struct Point_T *centre);
int		inner_To_Voronoi_Area( struct Voronoi_T *voronoi, int index,
													struct Point_T *q);
void 	build_Voronoi( struct Voronoi_T *voronoi_Set, struct DCEL_T *dcel);
void 	finalize_Voronoi( struct Voronoi_T *voronoi_Set);

#endif
