#ifndef DCEL_H
#define DCEL_H

#include "point.h"
#include "stack.h"
#include <stdio.h>


/**************************************************************************
* Defines
**************************************************************************/
#define     EXTERNAL_FACE       	0

#define 	POINTS_TYPE		    	1
#define 	DCEL_TYPE		    	2

#define     NO_UPDATE           	-1

#define 	CONVEX_HULL_LEN_FACTOR 	1

/**************************************************************************
* Data types definition
**************************************************************************/
enum Resize_DCEL_T {RESIZE_POINTS,
			  	  	RESIZE_EDGES,
					RESIZE_FACES};


struct	Dcel_Vertex_T
{
	int					origin_Edge;
	struct Point_T	    vertex;
};

struct Dcel_Edge_T
{
	int					origin_Vertex;
	int					twin_Edge;
	int					previous_Edge;
	int					next_Edge;
	int					face;
};

struct	Dcel_Face_T
{
	int					edge;
	int					imaginaryFace;
};

struct DCEL_T
{
	bool					incremental;

	// Vertex data.
	int 					nVertex;
	int						sizeVertex;
	struct	Dcel_Vertex_T	*vertex;

	// Edges data.
	int 					nEdges;
	int						sizeEdges;
	int						*edgeChecked;
	struct	Dcel_Edge_T	    *edges;

	// Vertex data.
	int 					nFaces;
	int						sizeFaces;
	struct	Dcel_Face_T	    *faces;
};


/**************************************************************************
* Public function headers.
**************************************************************************/
// General DCEL interface.
int 	initialize_DCEL( struct DCEL_T *dcel, int nPoints, int nEdges, int nFaces);
void	reset_DCEL( struct DCEL_T *dcel);
void    copy_Dcel( struct DCEL_T *in_dcel, struct DCEL_T *out_dcel);
int 	resize_DCEL( struct DCEL_T *dcel, enum Resize_DCEL_T resize_Type);
void    check_DCEL_Data( struct DCEL_T *dcel);
void	shake_Points_DCEL( struct DCEL_T *dcel);
void	finalize_DCEL( struct DCEL_T *dcel);

// Vertexes interface.
int 	get_Number_Vertex( struct DCEL_T *dcel);
int	 	insertPoint( struct DCEL_T *dcel, struct Point_T *point);
int 	insertVertex( struct DCEL_T *dcel, struct Dcel_Vertex_T vertex);
int	 	insert_Vertex_At( struct DCEL_T *dcel, struct Dcel_Vertex_T vertex, int index);
int 	update_Vertex_Edge_At( struct DCEL_T *dcel, int edge_ID, int index);
int		swap_Vertex( struct DCEL_T *dcel, int index1, int index2);
struct Dcel_Vertex_T *get_Vertex( struct DCEL_T *dcel, int index);
struct Point_T       *get_Vertex_Point( struct DCEL_T *dcel, int index);

// Edges interface.
int 	get_Number_Edges( struct DCEL_T *dcel);
int 	get_Number_Real_Edges( struct DCEL_T *dcel);
int 	insertEdge( struct DCEL_T *dcel, int origin, int twin, int prev, int next, int face);
int	    update_Edge( struct DCEL_T *dcel, int origin, int twin, int prev, int next, int face, int index);
int		are_Twins( struct DCEL_T *dcel, struct Dcel_Edge_T *edge1, struct Dcel_Edge_T *edge2);
int		is_External_Edge( struct DCEL_T *dcel, int index);
int		set_Edge_Not_Checked( struct DCEL_T *dcel, int index, int	*n);
int     get_Edge_Origin_Vertex( struct DCEL_T *dcel, int edge_Index);
int 	get_Edge_In_Convex_Hull( struct DCEL_T *dcel, int faceID);
struct  Dcel_Edge_T *get_Edge( struct DCEL_T *dcel, int index);
int 	copy_Edge( struct DCEL_T *dcel, int index, struct Dcel_Edge_T *edge);

// Faces interface.
int     get_Number_Faces(struct DCEL_T *dcel);
int     get_Number_Real_Faces(struct DCEL_T *dcel);
int 	insertFace( struct DCEL_T *dcel, int edge_ID);
int		update_Face(struct DCEL_T *dcel, int edge_ID, int index);
struct Dcel_Face_T *get_Face( struct DCEL_T *dcel, int index);
int     get_Face_Vertex( struct DCEL_T *dcel, int face_ID, struct Dcel_Vertex_T *v1,
									struct Dcel_Vertex_T *v2, struct Dcel_Vertex_T *v3);
int     is_Negative_Any_Vertex( struct DCEL_T *dcel, int edge_Index);
void	printFace( struct DCEL_T *dcel, int faceID);

// Points interface.
enum Turn_T return_Turn( struct DCEL_T *dcel, struct Point_T *p,
												int source_ID,
												int dest_ID);
void    shake_Dcel( struct DCEL_T *dcel);
void    get_Extreme_Point( struct DCEL_T *dcel, int (*f)(struct Point_T *, struct Point_T *),
																				struct Point_T *p);
void	get_Vertex_Of_Face( struct DCEL_T *dcel, int face, int *index1, int *index2, int *index3);
bool	is_Interior_To_Face( struct DCEL_T *dcel, struct Point_T *p, int face);

// Figures interface.
bool	get_Convex_Hull( struct DCEL_T *dcel, int *length, int *edges);
bool 	is_Interior_To_Convex_Hull( struct DCEL_T *dcel, struct Point_T *p, bool *error);

// I/O DCEL interface.
int 	read_DCEL( struct DCEL_T *dcel, char *fileName);
int     read_Points_Flat_File( struct DCEL_T *dcel, const char *fileName);
int		write_DCEL( struct DCEL_T *dcel, int type, const char *fileName);
void	print_DCEL( struct DCEL_T *dcel);
int 	read_Points_DCEL( FILE *fd, int nPoints, struct DCEL_T *dcel);
void	generate_Random_Points_DCEL( int nPoints, struct DCEL_T *dcel, TYPE maxX, TYPE maxY);
void	generate_Cluster_Points_DCEL( int nPoints, struct DCEL_T *dcel,
							int nClusters, int radius, TYPE maxX, TYPE maxY);

#define DCEL_STATISTICS
#ifdef DCEL_STATISTICS
void 	print_Dcel_Statistics( char *fileName, struct DCEL_T *dcel);
#endif

#endif
