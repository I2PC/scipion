#ifndef DELAUNAY_H
#define DELAUNAY_H

#include "defines.h"
#include "dcel.h"
#include "graph.h"
#include "point.h"
#include "voronoi.h"

/*****************************************************************************
* Data types definition
*****************************************************************************/
struct Delaunay_T
{
	bool				voronoiComputed;	// Voronoi diagram computed.
	struct Graph_T 		graph;				// Incremental graph data.
	struct DCEL_T		*dcel;				// DCEL data.
	struct Voronoi_T 	voronoi;			// Voronoi data.

	int 				faceIndex;			// DCEL face iterator.
};

/**************************************************************************
* Public functions declaration
**************************************************************************/
int		init_Delaunay( struct Delaunay_T *delaunay, int nPoints);
void	delete_Delaunay( struct Delaunay_T *delaunay);
int		insert_Point( struct Delaunay_T *delaunay, double x, double y);
int		create_Delaunay_Triangulation( struct Delaunay_T *delaunay,
														int createVoronoi);
int 	get_Face_Points( struct Delaunay_T *delaunay, int face_ID,
												struct Point_T *p,
												struct Point_T *q,
												struct Point_T *r);
bool 	next_Face_Iterator(struct Delaunay_T *delaunay, struct Point_T *p1,
														struct Point_T *p2,
														struct Point_T *p3);

int 	initialize_Delaunay(struct Delaunay_T *delaunay, struct DCEL_T *dcel);
void 	finalize_Delaunay(struct Delaunay_T *delaunay);

void	incremental_Delaunay(struct Delaunay_T *delaunay);
void 	build_Delaunay_From_Triangulation(struct DCEL_T *dcel);
void    purge_Delaunay( struct Delaunay_T *delaunay);

int  	select_Closest( struct Delaunay_T *delaunay, int index);
void 	select_Two_Closest( struct Delaunay_T *delaunay, int *first, int *second);
bool 	select_Closest_Point( struct Delaunay_T *delaunay, struct Point_T *p,
														struct Point_T *q,
														double *lowest_Distance);
bool 	select_Closest_Point_DCEL( struct DCEL_T *dcel, int nAnchors,
														struct Point_T *p,
												   	   	struct Point_T *q,
														double *lowest_Distance);

#endif

