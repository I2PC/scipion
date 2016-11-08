#include "dcel.h"
#include "defines.h"
#include "delaunay.h"
#include <float.h>
#include "graph.h"
#ifdef LOGGING
#include "log.h"
#endif
#ifdef LOCATION_STATISTICS
#include "statistics.h"
#endif
#include "point.h"
#include "polygon.h"
#include "sorting.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#ifdef LOCATION_STATISTICS
#include "/home/jnavas/workspace/Delaunay_Scipion/performance.h"
#endif

#ifdef DEBUG
#include <GL/glut.h>
#endif

/*****************************************************************************
* Defines
*****************************************************************************/
#define MAX_NEIGHBORS		50
#define EXTERNAL_FACE       0

/*****************************************************************************
* Variables declaration
*****************************************************************************/
struct Node_T 					root_Node;
struct Graph_T  				graph;

/*****************************************************************************
* Private functions declaration
*****************************************************************************/
void 		insert_First_Node( struct Delaunay_T *delaunay);
void		get_Vertex_Of_Node( struct Graph_T *graph, int node_Index, int *index1, int *index2, int *index3);
int			is_Interior_To_Node( struct DCEL_T *dcel, struct Graph_T *graph,
														struct Point_T *p,
														int node_Index);
int         is_Strictly_Interior_To_Node(struct DCEL_T *dcel, struct Graph_T *graph, struct Point_T *p, int node_Index);
bool		analyze_Face( struct Delaunay_T *delaunay, int point_Index);
void    	split_Triangle( struct DCEL_T *dcel, struct Graph_T *graph, int point_Index, int node_Index, int nTriangles);
double  	modified_Signed_Area( struct DCEL_T *dcel, struct  Node_T *node);
int     	select_Colinear_Edge( struct DCEL_T *dcel, int point_Index, int edge_ID);
void		check_Edge( struct DCEL_T *dcel, struct Graph_T *graph, int edge_ID);
void		flip_Edges_Dcel( struct DCEL_T *dcel, struct Graph_T *graph, int edge_ID);


/***************************************************************************
* Public functions declaration
***************************************************************************/
/***************************************************************************
* Name: init_Delaunay
* IN:		nPoints			# points in triangulation
* OUT:		delaunay		delaunay data
* IN/OUT:	N/A
* RETURN:	SUCCESS if delaunay allocated. FAILURE i.o.c
* Description: 	allocates structures to store "nPoints" points.
***************************************************************************/
int		init_Delaunay( struct Delaunay_T *delaunay, int nPoints)
{
	int		ret=SUCCESS; 		// Return value.

	// Allocate DCEL.
	delaunay->dcel = (struct DCEL_T *)malloc(sizeof(struct DCEL_T));
	if (delaunay->dcel == NULL)
	{
#ifdef LOGGING
		sprintf( log_Text, "Error allocating DCEL in init_Delaunay\n");
		write_Log( log_Text);
#endif
		printf("Error allocating DCEL in init_Delaunay\n");
		ret = FAILURE;
	}
	else
	{
		// Create DCEL data to store nPoints.
		if (initialize_DCEL( delaunay->dcel, nPoints, INVALID, INVALID) == FAILURE)
		{
#ifdef LOGGING
			sprintf( log_Text, "Error calling initialize_DCEL in init_Delaunay\n");
			write_Log( log_Text);
#endif
			printf("Error allocating DCEL in init_Delaunay\n");
			ret = FAILURE;
		}
		else
		{
			// Initialize graph.
			if (initialize_Graph( &delaunay->graph, delaunay->dcel->sizeVertex*10) == FAILURE)
			{
#ifdef LOGGING
				sprintf( log_Text, "Error allocating graph in init_Delaunay\n");
				write_Log( log_Text);
#endif
				printf("Error allocating graph in init_Delaunay\n");
				ret = FAILURE;

				delaunay->voronoiComputed = true;
			}
		}
	}

	return(ret);
}


/***************************************************************************
* Name: delete_Delaunay
* IN:		N/A
* OUT:		N/A
* IN/OUT:	delaunay		delaunay data
* RETURN:	N/A
* Description: 	Free delaunay data.
***************************************************************************/
void	delete_Delaunay( struct Delaunay_T *delaunay)
{
	// Deallocate DCEL.
	finalize_DCEL(delaunay->dcel);

	// Deallocate Delaunay.
	finalize_Delaunay(delaunay);
}


/***************************************************************************
* Name: insert_Point
* IN:		x				x coordinate of new point
* 			y				y coordinate of new point
* OUT:		N/A
* IN/OUT:	delaunay		delaunay data
* RETURN:	N/A
* Description: 	Inserts a new point in the Delaunay set of points.
***************************************************************************/
int		insert_Point( struct Delaunay_T *delaunay, double x, double y)
{
	int		ret=SUCCESS;			// Return value.

	// Check if dcel is full.
	if (delaunay->dcel->nVertex == delaunay->dcel->sizeVertex)
	{
		printf("Dcel is full. Resizing DCEL\n");

		// Resize vertex array.
		if (resize_DCEL( delaunay->dcel, RESIZE_POINTS) == FAILURE)
		{
#ifdef LOGGING
			sprintf("Error resizing DCEL in insert_Point.\n");
			write_Log( log_Text);
#endif
			printf("Error resizing DCEL in insert_Point.\n");
			ret = FAILURE;
		}
	}

	if (ret == SUCCESS)
	{
		// Update coordinates.
		delaunay->dcel->vertex[delaunay->dcel->nVertex].vertex.x = x;
		delaunay->dcel->vertex[delaunay->dcel->nVertex].vertex.y = y;

		// Increase # points.
		delaunay->dcel->nVertex++;
	}

	return(ret);
}


/***************************************************************************
* Name: create_Delaunay_Triangulation
* IN:		createVoronoi	determines if Voronoi areas must be computed.
* OUT:		N/A
* IN/OUT:	delaunay		delaunay data
* RETURN:	SUCCESS if Voronoi computed. FAILURE i.o.c.
* Description: 	Computes the Delaunay triangulation using points stored in
* 				the DCEL field of the "delaunay" data. If "createVoronoi"
* 				TRUE, then it also computes the Voronoi area of the
* 				triangulation
***************************************************************************/
int create_Delaunay_Triangulation( struct Delaunay_T *delaunay, int createVoronoi)
{
	int		ret=SUCCESS;		// Return value.

	// Build Delaunay triangulation from scratch.
	incremental_Delaunay( delaunay);

	// Set as imaginary faces with vertex points P_MINUS_2 or P_MINUS_1.
	purge_Delaunay( delaunay);

	delaunay->voronoiComputed = false;

	// Check if Voronoi areas must be computed.
	if (createVoronoi)
	{
		// Allocate Voronoi structures.
		if (initialize_Voronoi( &delaunay->voronoi, get_Number_Faces( delaunay->dcel)*2) == FAILURE)
		{
#ifdef LOGGING
			sprintf("Error allocating Voronoi data in create_Delaunay_Triangulation.\n");
			write_Log( log_Text);
#endif
			printf("Error allocating Voronoi data in create_Delaunay_Triangulation.\n");
			ret = FAILURE;
		}
		else
		{
			// Build Voronoi diagram.
			build_Voronoi( &delaunay->voronoi, delaunay->dcel);

			// Set Voronoi flag as computed.
			delaunay->voronoiComputed = true;
		}
	}

	return(ret);
}



/***************************************************************************
* Name: initialize_Delaunay
* IN:		dcel			dcel data
* OUT:		N/A
* IN/OUT:	delaunay		delaunay data
* RETURN:	N/A
* Description: 	Allocates graph and gets reference to input dcel
***************************************************************************/
int 	initialize_Delaunay(struct Delaunay_T *delaunay, struct DCEL_T *dcel)
{
	int		ret=SUCCESS;		// Return value.

	// Get reference to DCEL.
	delaunay->dcel = dcel;
	delaunay->voronoiComputed = false;

	// Initialize graph.
	if (initialize_Graph( &delaunay->graph, delaunay->dcel->nVertex*10) == FAILURE)
	{
#ifdef LOGGING
		sprintf("Error allocating graph in DCEL in insert_Point.\n");
		write_Log( log_Text);
#endif
		printf("Error allocating graph in DCEL in insert_Point.\n");
		ret = FAILURE;
	}

	return(ret);
}


/***************************************************************************
* Name: finalize_Delaunay
* IN:		N/A
* OUT:		N/A
* IN/OUT:	delaunay		delaunay data
* RETURN:	N/A
* Description: 	Free graph, Voronoi areas and dereferences dcel.
***************************************************************************/
void 	finalize_Delaunay(struct Delaunay_T *delaunay)
{
	// Delete graph.
	finalize_Graph( &delaunay->graph);

	// Delete Voronoi.
	if (delaunay->voronoiComputed)
	{
		finalize_Voronoi( &delaunay->voronoi);
	}

	// Dereference DCEL.
	delaunay->dcel = NULL;
}


//#define DEBUG_INCREMENTAL_DELAUNAY
/***************************************************************************
* Name: incremental_Delaunay
* IN:		N/A
* OUT:		N/A
* IN/OUT:	delaunay		delaunay data
* RETURN:	N/A
* Description:  Builds a Delaunay triangulation using the incremental
* 				algorithm.
***************************************************************************/
void	incremental_Delaunay(struct Delaunay_T *delaunay)
{
    int     point_Index=0;                  // Points loop counter.

#ifdef DELAUNAY_STATISTICS
    // Initialize graph statistics (first 2 points are located in first triangle).
    delaunay_Stat.trianglesFound[0] = 2;
    delaunay_Stat.trianglesFound[1] = 0;
    delaunay_Stat.trianglesFound[2] = 0;
    delaunay_Stat.nFlipped = 0;
#endif
    delaunay->dcel->incremental = true;

    // Set highest point at first position of the DCEL vertex array.
    set_Highest_First( delaunay->dcel, &lexicographic_Higher);

    // Insert initial nodes.
    insert_First_Node( delaunay);

    // Update edge from new point.
    update_Vertex_Edge_At( delaunay->dcel, 1, 0);

    // Insert first 6 edges due to first point.
    insertEdge( delaunay->dcel, 1, 4, 3, 2, 1);
    insertEdge( delaunay->dcel, P_MINUS_2, 6, 1, 3, 1);
    insertEdge( delaunay->dcel, P_MINUS_1, 5, 2, 1, 1);
    insertEdge( delaunay->dcel, P_MINUS_2, 1, 6, 5, 0);
    insertEdge( delaunay->dcel, 1, 3, 4, 6, 0);
    insertEdge( delaunay->dcel, P_MINUS_1, 2, 5, 4, 0);

    // Insert first internal face and external face.
    insertFace( delaunay->dcel, 4);
    insertFace( delaunay->dcel, 1);

    // Loop all other points.
    for (point_Index=1; point_Index<delaunay->dcel->nVertex ; point_Index++)
    {
        // Insert new point into triangle where it is located.
        analyze_Face( delaunay, point_Index);
    }
}


/***************************************************************************
* Name: build_Delaunay_From_Triangulation
* IN:		N/A
* OUT:		N/A
* IN/OUT:	delaunay		delaunay data
* RETURN:	N/A
* Description:  Builds a Delaunay triangulation using the incremental
* 				algorithm.
***************************************************************************/
void 	build_Delaunay_From_Triangulation(struct DCEL_T *dcel)
{
	int		i=0;							// Loop counters.
	int		nPending=0;						// Edges pending.
	int		temp=0;							// Intermediate variable.
	struct Point_T triangle[3];
	struct Dcel_Vertex_T *vertex=NULL;	// Current triangle.
	struct Dcel_Edge_T *prev=NULL;		// Previous of current edge.
	struct Dcel_Edge_T *twin=NULL;		// Twin of current edge.
	struct Dcel_Edge_T *neighbor=NULL;	// Twin of current edge.
#ifdef DEBUG_DCEL
	char	prueba;
	char	filename[20];
	int		step=0;
#endif
#ifdef DELAUNAY_STATISTICS
	delaunay_Stat.nFlipped = 0;
#endif

	// Initialize variables.
	i=0;
	nPending = dcel->nEdges;
	dcel->incremental = false;

	// Analyze all edges.
	while (nPending > 0)
	{
		// Check edge is still pending.
		if (!dcel->edgeChecked[i])
		{
#ifdef DEBUG_DCEL
			printf("Edge %d in face %d and its twin is %d in face %d\n", i+1, dcel->edges[i].face, dcel->edges[i].twin_Edge, dcel->edges[dcel->edges[i].twin_Edge-1].face);
#endif
			// Check edge is not in convex hull.
			if ((dcel->edges[i].face != 0) && (dcel->edges[dcel->edges[i].twin_Edge-1].face) != 0)
			{
				// Get origin of first edge of current face.
				vertex = get_Vertex( dcel, dcel->edges[i].origin_Vertex-1);
				triangle[0] = vertex->vertex;

				// Get origin of previous edge.
				prev = get_Edge( dcel, dcel->edges[i].previous_Edge-1);
				vertex = get_Vertex( dcel, prev->origin_Vertex-1);
				triangle[2] = vertex->vertex;

				// Get origin of twin of current edge of current face.
				twin = get_Edge( dcel, dcel->edges[i].twin_Edge-1);
				vertex = get_Vertex( dcel, twin->origin_Vertex-1);
				triangle[1] = vertex->vertex;

				// Get origin of previous edge of twin.
				neighbor = get_Edge( dcel, twin->previous_Edge-1);
				vertex = get_Vertex( dcel, neighbor->origin_Vertex-1);

				// Check if neighbor point is in circle.
				if (in_Circle( &triangle[0], &triangle[1], &triangle[2], &vertex->vertex))
				{
#ifdef DELAUNAY_STATISTICS
					delaunay_Stat.nFlipped++;
#endif
#ifdef DEBUG_DCEL
					printf("Edge %d wrong. Pending %d. \n", i+1, nPending);
					scanf( "%c", &prueba);
#endif
					// Update vertex.
					dcel->vertex[dcel->edges[i].origin_Vertex-1].origin_Edge = twin->next_Edge;
					dcel->vertex[twin->origin_Vertex-1].origin_Edge = dcel->edges[i].next_Edge;

					// Update origin of current and twin edges.
					temp = dcel->edges[dcel->edges[i].previous_Edge-1].origin_Vertex;
					dcel->edges[twin->twin_Edge-1].origin_Vertex = dcel->edges[twin->previous_Edge-1].origin_Vertex;
					dcel->edges[dcel->edges[i].twin_Edge-1].origin_Vertex = temp;

					// Update next edges.
					dcel->edges[dcel->edges[i].next_Edge-1].next_Edge = dcel->edges[i].twin_Edge;
					dcel->edges[twin->next_Edge-1].next_Edge = twin->twin_Edge;
					dcel->edges[dcel->edges[i].previous_Edge-1].next_Edge = twin->next_Edge;
					dcel->edges[twin->previous_Edge-1].next_Edge = dcel->edges[i].next_Edge;
					dcel->edges[twin->twin_Edge-1].next_Edge = dcel->edges[twin->twin_Edge-1].previous_Edge;
					dcel->edges[dcel->edges[i].twin_Edge-1].next_Edge = dcel->edges[dcel->edges[i].twin_Edge-1].previous_Edge;

					// Update previous edges.
					dcel->edges[twin->twin_Edge-1].previous_Edge = dcel->edges[dcel->edges[i].next_Edge-1].next_Edge;
					dcel->edges[dcel->edges[i].twin_Edge-1].previous_Edge = dcel->edges[twin->next_Edge-1].next_Edge;
					dcel->edges[dcel->edges[i].next_Edge-1].previous_Edge = dcel->edges[dcel->edges[i].previous_Edge-1].next_Edge;
					dcel->edges[twin->next_Edge-1].previous_Edge = dcel->edges[twin->previous_Edge-1].next_Edge;
					dcel->edges[dcel->edges[i].previous_Edge-1].previous_Edge = dcel->edges[twin->twin_Edge-1].next_Edge;
					dcel->edges[twin->previous_Edge-1].previous_Edge = dcel->edges[dcel->edges[i].twin_Edge-1].next_Edge;

					// Update faces in edges.
					dcel->edges[dcel->edges[i].previous_Edge-1].face = dcel->edges[i].face;
					dcel->edges[twin->previous_Edge-1].face = twin->face;

					// Update faces.
					dcel->faces[dcel->edges[i].face].edge = twin->twin_Edge;
					dcel->faces[twin->face].edge = dcel->edges[i].twin_Edge;

					// Update pending edges.
					dcel->edgeChecked[i] = 1;
					nPending--;

					// If twin not checked set as checked.
					if (!dcel->edgeChecked[dcel->edges[i].twin_Edge-1])
					{
						nPending--;
						dcel->edgeChecked[dcel->edges[i].twin_Edge-1] = 1;
					}

					// Update pending edges.
					set_Edge_Not_Checked( dcel, dcel->edges[i].previous_Edge-1, &nPending);
					set_Edge_Not_Checked( dcel, dcel->edges[i].next_Edge-1, &nPending);
					set_Edge_Not_Checked( dcel, dcel->edges[dcel->edges[i].twin_Edge-1].previous_Edge-1, &nPending);
					set_Edge_Not_Checked( dcel, dcel->edges[dcel->edges[i].twin_Edge-1].next_Edge-1, &nPending);
#ifdef DEBUG_DCEL
					printf("Edge %d flipped. Pending %d. \n", i+1, nPending);
					step++;
					sprintf( filename, "setp%d.txt", step);
					write_DCEL( dcel, filename);
					draw_DCEL_Points( &delaunay_Dcel, 0.1, WHITE);
					draw_DCEL_Triangulation( &delaunay_Dcel);
					glutSwapBuffers();
#endif
				}
				// Edge OK -> Not to be flipped.
				else
				{
					dcel->edgeChecked[i] = 1;
					nPending--;
#ifdef DEBUG_DCEL
					printf("Edge %d is OK. Pending %d. \n", i+1, nPending);
					scanf( "%c", &prueba);
#endif
				}
			}
			// Edge in convex hull -> Not to be flipped.
			else
			{
				dcel->edgeChecked[i] = 1;
				nPending--;
#ifdef DEBUG_DCEL
				printf("Edge %d in convex hull. Peding %d\n", i+1, nPending);
				scanf( "%c", &prueba);
#endif
			}
		}
		else
		{
#ifdef DEBUG_DCEL
			printf("Edge %d already checked. Pending %d\n", i+1, nPending);
			scanf( "%c", &prueba);
#endif
		}

		// Next edge.
		i++;
		i = i % dcel->nEdges;
	}
}

void    purge_Delaunay( struct Delaunay_T *delaunay)
{
    int edge_Index=0;           // Next edge to update.

    // Initialize edge and point indexes.
    edge_Index = delaunay->dcel->edges[edge_Index].previous_Edge-1;

    do
    {
        // Invalid current face.
        update_Face( delaunay->dcel, INVALID, delaunay->dcel->edges[edge_Index].face);

        // Get index of edge departing from next point in convex hull.
        edge_Index = delaunay->dcel->edges[edge_Index].previous_Edge-1;
        edge_Index = delaunay->dcel->edges[edge_Index].twin_Edge-1;

        // Get index of next edge to update.
        edge_Index = delaunay->dcel->edges[edge_Index].previous_Edge-1;

        // Check if new unbounded point is P_MINUS_1.
        if (delaunay->dcel->edges[edge_Index].origin_Vertex == P_MINUS_1)
        {
            // Invalid current face.
            update_Face( delaunay->dcel, INVALID, delaunay->dcel->edges[edge_Index].face);

            // Skip edges in current face.
            edge_Index = delaunay->dcel->edges[edge_Index].twin_Edge-1;
            edge_Index = delaunay->dcel->edges[edge_Index].previous_Edge-1;
        }

    } while(delaunay->dcel->edges[edge_Index].origin_Vertex != 1);
    // Until top most point found again.

    // Invalid current face.
    update_Face( delaunay->dcel, INVALID, delaunay->dcel->edges[edge_Index].face);

    // Invalid zero face.
    update_Face( delaunay->dcel, INVALID, 0);
}


int    select_Closest( struct Delaunay_T *delaunay, int index)
{
    int     finished=FALSE;                 // Loop control flag.
	int     point_Index=0;                  // Index of a point.
	int     closest_Index=0;                // Index of closest point.
	int     edge_Index=0;                   // Index of current edge.
	struct Dcel_Vertex_T *origin=NULL;      // Vertex data.
	struct Dcel_Vertex_T *new_Point=NULL;   // Vertex data.
	double  closest_Dist=0.0, new_Dist=0.0; // Distance to current closest point.

	// Get edge departing from point.
	origin = get_Vertex( delaunay->dcel, index);

	// Get first edge to analyze.
	edge_Index = origin->origin_Edge-1;

#ifdef DEBUG_TWO_CLOSEST
	printf("Searching closest point to %d. Starting edge %d\n", index, edge_Index);
#endif

    // Get index of current point.
    closest_Index = delaunay->dcel->edges[delaunay->dcel->edges[edge_Index].next_Edge-1].origin_Vertex-1;
    new_Point = get_Vertex( delaunay->dcel, closest_Index);
	closest_Dist = distance( &origin->vertex, &new_Point->vertex);

#ifdef DEBUG_TWO_CLOSEST
	printf("Initial point is %d whose distance is %lf\n", closest_Index, closest_Dist);
#endif

    // Initialize loop.
    finished = FALSE;
	while (!finished)
	{
	    // Select new edge.
	    edge_Index = delaunay->dcel->edges[delaunay->dcel->edges[edge_Index].previous_Edge-1].twin_Edge-1;
#ifdef DEBUG_TWO_CLOSEST
	    printf("Checking edge %d\n", edge_Index);
#endif

	    if (edge_Index == origin->origin_Edge-1)
	    {
            finished = TRUE;

#ifdef DEBUG_TWO_CLOSEST
            printf("Search finished.\n");
#endif
	    }
	    else
	    {
	        point_Index = delaunay->dcel->edges[delaunay->dcel->edges[edge_Index].next_Edge-1].origin_Vertex-1;

#ifdef DEBUG_TWO_CLOSEST
	        printf("Checking point %d\n", point_Index);
#endif

	        if (point_Index >= 0)
	        {
                new_Point = get_Vertex( delaunay->dcel, point_Index);
                new_Dist = distance( &origin->vertex, &new_Point->vertex);

#ifdef DEBUG_TWO_CLOSEST
                printf("New distance is %lf and previous is %lf.", new_Dist, closest_Dist);
#endif
                if (new_Dist < closest_Dist)
                {
                    closest_Dist = new_Dist;
                    closest_Index = point_Index;
#ifdef DEBUG_TWO_CLOSEST
                	printf("New point is closer. \n");
#endif
                }
	        }
#ifdef DEBUG_TWO_CLOSEST
	        printf("Not checked because it is an imaginary point\n");
#endif
	    }
	}

	return(closest_Index);
}


//#define #ifdef DEBUG_GET_FACE_POINTS
/***************************************************************************
* Name: get_Face_Points
* IN:	dcel			dcel data
* 		face_ID			face to recover its points
* OUT:	p				first point of the face
* 		q				second point of the face
* 		r				third point of the face
* Description: 	gets the three points that form a triangle associated to the
* 				face_ID face.
***************************************************************************/
int    get_Face_Points( struct Delaunay_T *delaunay, int face_ID, struct Point_T *p,
																struct Point_T *q,
																struct Point_T *r)
{
	int		i=1;
	int		index=0;
	int		found=FALSE;
	struct Dcel_Edge_T  *edge=NULL;	// Current edge.
	struct Dcel_Face_T  *face=NULL;	// Current face.
	int		faceIndex=0;

	i=1;
	while ((!found) && (i<get_Number_Faces( delaunay->dcel)))
	{
		// Get i-face.
		face = get_Face( delaunay->dcel, i);

		if (face->imaginaryFace != INVALID)
		{
			index = face->edge-1;
			edge = get_Edge( delaunay->dcel, index);

			faceIndex++;
#ifdef DEBUG_GET_FACE_POINTS
			printf("Face index %d. Target face %d\n", faceIndex, face_ID);
#endif
			if (faceIndex == face_ID)
			{
#ifdef DEBUG_GET_FACE_POINTS
				printf("%d \n ", dcel->edges[index].origin_Vertex-1);
				printf("%d \n", dcel->edges[edge->next_Edge-1].origin_Vertex-1);
				printf("%d \n", dcel->edges[edge->previous_Edge-1].origin_Vertex-1);
#endif
				if (delaunay->dcel->edges[index].origin_Vertex-1 >= 0)
				{
					p->x = delaunay->dcel->vertex[delaunay->dcel->edges[index].origin_Vertex-1].vertex.x;
					p->y = delaunay->dcel->vertex[delaunay->dcel->edges[index].origin_Vertex-1].vertex.y;
				}
				else
				{
					p->x = INVALID;
					p->y = INVALID;
				}
				if (delaunay->dcel->edges[edge->next_Edge-1].origin_Vertex-1 >= 0)
				{
					q->x = delaunay->dcel->vertex[delaunay->dcel->edges[edge->next_Edge-1].origin_Vertex-1].vertex.x;
					q->y = delaunay->dcel->vertex[delaunay->dcel->edges[edge->next_Edge-1].origin_Vertex-1].vertex.y;
				}
				else
				{
					q->x = INVALID;
					q->y = INVALID;
				}
				if (delaunay->dcel->edges[edge->previous_Edge-1].origin_Vertex-1 >= 0)
				{
					r->x = delaunay->dcel->vertex[delaunay->dcel->edges[edge->previous_Edge-1].origin_Vertex-1].vertex.x;
					r->y = delaunay->dcel->vertex[delaunay->dcel->edges[edge->previous_Edge-1].origin_Vertex-1].vertex.y;
				}
				else
				{
					r->x = INVALID;
					r->y = INVALID;
				}

				found=TRUE;
			}
		}
		i++;
	}

	return(found);
}


//#define DEBUG_NEXT_FACE_ITERATOR
/***************************************************************************
* Name: next_Face_Iterator
* IN:	delaunay		delaunay data
* OUT:	p1				first point of the face
* 		p2				second point of the face
* 		p3				third point of the face
* Description: 	gets the three points that form a triangle associated to a
* 				triangle. First call returns the first triangle and successive
* 				calls return the next triangle (the function remembers last
* 				triangle returned). The function returns FALSE if the end
* 				of the list of triangles was reached and a real triangle
* 				was not found.
***************************************************************************/
bool 	next_Face_Iterator(struct Delaunay_T *delaunay, struct Point_T *p1,
														struct Point_T *p2,
														struct Point_T *p3)
{
	bool 	validFace;				// Return value.
	int		edgeIndex;				// Edge index.
	bool	finished;				// Loop control flag.
	struct Dcel_Face_T  *face=NULL;	// Current face.
	struct Dcel_Edge_T  *edge=NULL;	// Current edge.
	int		nFaces;					// # faces in DCEL.

	// Get # faces in DCEL.
	nFaces = get_Number_Faces( delaunay->dcel);

	// Check if face index is out of bounds.
	if (delaunay->faceIndex < nFaces)
	{
		// Search next face.
		finished = false;
		while (!finished)
		{
			// Get i-face.
			face = get_Face( delaunay->dcel, delaunay->faceIndex);

#ifdef DEBUG_NEXT_FACE_ITERATOR
			printf("\nNew face %d of %d\n", delaunay->faceIndex, nFaces);
#endif

			delaunay->faceIndex++;

			// Check if current face is imaginary.
			if (face->imaginaryFace != INVALID)
			{
				// Update return value and loop control flag to exit.
				finished = true;
				validFace = true;

				// Get edge index and edge data.
				edgeIndex = face->edge - 1;
				edge = get_Edge( delaunay->dcel, edgeIndex);

#ifdef DEBUG_NEXT_FACE_ITERATOR
				printf("\t1st point index %d \n", delaunay->dcel->edges[edgeIndex].origin_Vertex-1);
				printf("\t2nd point index %d \n", delaunay->dcel->edges[edge->next_Edge-1].origin_Vertex-1);
				printf("\t3rd point index %d \n", delaunay->dcel->edges[edge->previous_Edge-1].origin_Vertex-1);
#endif

				// Update output points.
				p1->x = delaunay->dcel->vertex[delaunay->dcel->edges[edgeIndex].origin_Vertex-1].vertex.x;
				p1->y = delaunay->dcel->vertex[delaunay->dcel->edges[edgeIndex].origin_Vertex-1].vertex.y;
				p2->x = delaunay->dcel->vertex[delaunay->dcel->edges[edge->next_Edge-1].origin_Vertex-1].vertex.x;
				p2->y = delaunay->dcel->vertex[delaunay->dcel->edges[edge->next_Edge-1].origin_Vertex-1].vertex.y;
				p3->x = delaunay->dcel->vertex[delaunay->dcel->edges[edge->previous_Edge-1].origin_Vertex-1].vertex.x;
				p3->y = delaunay->dcel->vertex[delaunay->dcel->edges[edge->previous_Edge-1].origin_Vertex-1].vertex.y;
			}
			else
			{
#ifdef DEBUG_NEXT_FACE_ITERATOR
				// Get edge index and edge data.
				edgeIndex = face->edge - 1;
				edge = get_Edge( delaunay->dcel, edgeIndex);
				printf("Skipping imaginary face %d\n", delaunay->faceIndex-1);
				printf("\t1st point index %d \n", delaunay->dcel->edges[edgeIndex].origin_Vertex-1);
				printf("\t2nd point index %d \n", delaunay->dcel->edges[edge->next_Edge-1].origin_Vertex-1);
				printf("\t3rd point index %d \n", delaunay->dcel->edges[edge->previous_Edge-1].origin_Vertex-1);
#endif
				// Check if all faces searched.
				if (delaunay->faceIndex >= nFaces)
				{
					finished = true;
					validFace = false;
#ifdef DEBUG_NEXT_FACE_ITERATOR
					printf("\t\t\tLast face %d found. Resetting index.\n", delaunay->faceIndex);
#endif
					delaunay->faceIndex = 1;
				}
			}
		}
	}
	// Reset face index.
	else
	{
		delaunay->faceIndex = 1;
		validFace = false;
	}

	return(validFace);
}




void    select_Two_Closest( struct Delaunay_T *delaunay, int *first, int *second)
{
    int point_Index1=0, point_Index2=0;
    int edge_Index=0;                // Loop counter.
    double lowest_Dist=DBL_MAX;      // Distance between closest points.
    double new_distance=0.0;          // New distance.
	struct Dcel_Vertex_T *p1=NULL;   // Vertex data.
	struct Dcel_Vertex_T *p2=NULL;   // Vertex data.

    // Check all edges.
    for (edge_Index=0; edge_Index<delaunay->dcel->nEdges; edge_Index++)
    {
        // If twin edge higher than current edge then points have not yet been checked.
        if (delaunay->dcel->edges[edge_Index].twin_Edge-1 > edge_Index)
        {
            // Get two points indexes.
            point_Index1 = delaunay->dcel->edges[edge_Index].origin_Vertex-1;
            point_Index2 = delaunay->dcel->edges[delaunay->dcel->edges[edge_Index].twin_Edge-1].origin_Vertex-1;

            // Check if both points are valid (not P_MINUS_1 nor P_MINUS_2).
            if ((point_Index1 >= 0) && (point_Index2 >= 0))
            {
                // Get two points.
                p1 = get_Vertex( delaunay->dcel, point_Index1);
                p2 = get_Vertex( delaunay->dcel, point_Index2);

                // Compute distance.
                new_distance = distance( &p1->vertex, &p2->vertex);

                // Check if new distance is lower than current closest distance.
                if (new_distance < lowest_Dist)
                {
                    // Update lowest distance.
                    lowest_Dist = new_distance;

                    (*first) = point_Index1;
                    (*second) = point_Index2;
                }
            }
        }
    }
}


//#define DEBUG_SELECT_CLOSEST_POINT
/***************************************************************************
* Name: select_Closest_Point
* IN:	delaunay				delaunay triangulation
* 		voronoi					voronoi areas associated to delaunay triangulation
* 		p						input point
* OUT:		q					closest point to input "p" point
* 			lowest_Distance		distance from input "p" to output "q"
* RETURN:	True				TRUE if closest point found. False i.o.c.
* Description: finds the closest point in the DCEl to an input "p" point.
***************************************************************************/
bool select_Closest_Point( struct Delaunay_T *delaunay, struct Point_T *p,
														struct Point_T *q,
														double *lowest_Distance)
{
    int     i=0, j=0;			// Loop counters.
    int     node_Index=0;  		// Index of current node analyzed.
    bool	found=false;		// Loop control flag.
    bool	finished=false;		// Loop control flag.
	int		child_Index=0;		// Children node ID.
    int     nChildren=0;        // Number of children of current node.
    int		edge_Id=0;			// Edge identifier.

    int		current_Point=0;	// Current point checked.
    int		first_Point=0;		// First point checked.

    // Circular buffer to store point candidates.
    int		in=0;
    int		out=0;
    int		nItems=0;
    int		candidatesSize=MAX_NEIGHBORS;
    int		*candidates=NULL, *refCandidates=NULL;
    int		*inserted=NULL;
    int		nElems=0;					// # elements to copy.
    int		nDiscarded=0;
    bool	allocationError=false;		// Error allocating memory.

#ifdef LOCATION_STATISTICS
    // # triangles searched before a new point is located.
    location_Stat.current++;
	location_Stat.begin[location_Stat.current] = getTime();
    location_Stat.nTrianglesSearched[location_Stat.current] = 0;
#endif

#ifdef DEBUG_SELECT_CLOSEST_POINT
	printf("Point to search: ");
	print_Point( p);
#endif

    // Loop while not leaf node found and not end of set of points and not.
    finished = false;
    while (!finished)
    {
    	// If node is a leaf then triangle that surrounds point has been found.
    	if (is_Leaf_Node( &delaunay->graph, node_Index))
    	{
#ifdef DEBUG_SELECT_CLOSEST_POINT
    		printf("Leaf node found: %d\n", node_Index);
#endif
    		// Allocate neighbors array.
    		candidates = (int *) calloc( candidatesSize, sizeof(int));
    		inserted = (int *) calloc( delaunay->voronoi.dcel.nVertex, sizeof(int));

    		// Insert points from leaf node.
    		for (i=0; i<N_POINTS ;i++)
			{
				if (delaunay->graph.nodes[node_Index].points_Index[i] >= 0)
				{
					candidates[in] = delaunay->graph.nodes[node_Index].points_Index[i];
					inserted[candidates[in]] = TRUE;
#ifdef DEBUG_SELECT_CLOSEST_POINT
					printf("Inserted point %d. # %d. In %d. Out %d\n", candidates[in], nItems+1, in, out);
#endif
					nItems++;
					in++;
				}
			}

			found = false;
			allocationError = false;
			while ((!found) && (!allocationError))
			{
#ifdef DEBUG_SELECT_CLOSEST_POINT
				printf("Checking Point %d.\n", candidates[out]);
#endif
				if (candidates[out] >= 0)
				{
#ifdef LOCATION_STATISTICS
					// Update # triangles searched.
					location_Stat.nTrianglesSearched[location_Stat.current]++;
#endif
					if (inner_To_Voronoi_Area( &delaunay->voronoi, candidates[out]-1, p))
					{
						// End inner and outer loops.
						found = true;
						finished = true;

						// Update point and distance.
						q->x = delaunay->dcel->vertex[candidates[out]-1].vertex.x;
						q->y = delaunay->dcel->vertex[candidates[out]-1].vertex.y;
						(*lowest_Distance) = distance( p, q);

#ifdef DEBUG_SELECT_CLOSEST_POINT
						printf("Point %d found.\n", candidates[out]);
			    		getchar();
#endif
					}
					else
					{
#ifdef DEBUG_SELECT_CLOSEST_POINT
						printf("Discard point %d.\n", candidates[out]);
#endif
						// Get first neighbor.
						edge_Id = delaunay->dcel->vertex[candidates[out]-1].origin_Edge;
						edge_Id = delaunay->dcel->edges[edge_Id - 1].twin_Edge;
						current_Point = delaunay->dcel->edges[edge_Id - 1].origin_Vertex;
						first_Point = current_Point;

						do
						{
							// Check if circular buffer is full.
							if (nItems == candidatesSize)
							{
#ifdef DEBUG_SELECT_CLOSEST_POINT
								printf("Circular buffer is full. Size %d and current # elements %d\n",
																					nItems,
																					candidatesSize);
#endif
								// Compute # elements to copy.
								nElems = candidatesSize - out;

								// Double size of the candidates array.
								refCandidates = candidates;
								candidatesSize = candidatesSize*2;
								candidates = (int *) calloc( candidatesSize, sizeof(int));
								if (candidates != NULL)
								{
									// Copy current candidates.
									memcpy( candidates, &refCandidates[out], nElems*sizeof(int));
									memcpy( &candidates[nElems], refCandidates, (nItems-nElems)*sizeof(int));
									free( refCandidates);
								}
								else
								{
#ifdef LOGGING
									sprintf( log_Text, "Error allocating neighbors array in select_Closest_Point.\n");
									write_Log( log_Text);
#endif
									printf("Error allocating neighbors array in select_Closest_Point.\n");
									allocationError = true;
									break;
								}

								in = nItems;
								out = 0;
							}

							if (current_Point >= 0)
							{
								if (!inserted[current_Point])
								{
									inserted[current_Point] = TRUE;

									// Insert point.
									candidates[in] = current_Point;
									nItems++;
									in++;
									in = in % candidatesSize;
#ifdef DEBUG_SELECT_CLOSEST_POINT
									printf("Inserted point %d. # %d. In %d. Out %d\n", current_Point, nItems, in, out);
#endif
								}
#ifdef DEBUG_SELECT_CLOSEST_POINT
								else
								{
									printf("Point %d already checked.\n", current_Point);

								}
#endif
							}

							// Get next point.
							edge_Id = delaunay->dcel->edges[edge_Id - 1].next_Edge;
							edge_Id = delaunay->dcel->edges[edge_Id - 1].twin_Edge;
							current_Point = delaunay->dcel->edges[edge_Id - 1].origin_Vertex;
						}
						while (current_Point != first_Point);

						// Check if all vertices searched.
						nDiscarded++;
						if (nDiscarded == delaunay->dcel->nVertex)
						{
							allocationError = true;
							finished = true;
							found = false;
						}

#ifdef DEBUG_SELECT_CLOSEST_POINT
						printf("All %d neighbors inserted.\n", candidates[out]);
						getchar();
#endif
					}
				}
#ifdef DEBUG_SELECT_CLOSEST_POINT
				else
				{
					printf("Imaginary point %d.\n", candidates[out]);
					getchar();
				}
#endif
				// Check next point.
#ifdef DEBUG_SELECT_CLOSEST_POINT
				printf("Removed point %d. # %d. In %d. Out %d\n", candidates[out], nItems-1, in, out+1);
				getchar();
#endif
				out++;
				nItems--;
				out = out % candidatesSize;
			}

			// Free candidates circular buffer.
			free( candidates);
			free( inserted);
    	}
    	else
    	{
			// Search triangle in children nodes.
			i=0;
			found = false;
			nChildren = get_nChildren_Node( &delaunay->graph, node_Index);

#ifdef DEBUG_SELECT_CLOSEST_POINT
			printf("Start search in node %d\n", node_Index);
#endif
			while ((!found) && (i < nChildren))
			{
				// Get i-children.
				child_Index = get_iChildren_Node( &delaunay->graph, node_Index, i);

#ifdef LOCATION_STATISTICS
				// Update # triangles searched.
				location_Stat.nTrianglesSearched[location_Stat.current]++;
#endif

#ifdef DEBUG_SELECT_CLOSEST_POINT
				printf("Checking %d-child in node %d. Node %d\n", i, child_Index, node_Index);
#endif
				// Check if point is interior to i-child node.
				if (is_Interior_To_Node( delaunay->dcel, &delaunay->graph, p, child_Index))
				{
					// Search in next children node.
					node_Index = child_Index;

					// End loop.
					found = true;
#ifdef DEBUG_SELECT_CLOSEST_POINT
					printf("Point inside node %d\n", node_Index);
#endif
				}
				// Next child.
				else
				{
					i++;
					if (i == nChildren)
					{
						// Point must be equal to an existing point. Retry all points in current node.
						found = false;
						i=0;
						while ((!found) && (i < nChildren))
						{
							// Get i-children.
							child_Index = get_iChildren_Node( &delaunay->graph, node_Index, i);

							// Get closest point between the three points of the triangle.
							for (j=0; j<N_POINTS ;j++)
							{
								//printf("Point %d is (%lf,%lf)\n", delaunay->graph.nodes[child_Index].points_Index[j],
								//		delaunay->dcel->vertex[delaunay->graph.nodes[child_Index].points_Index[j]-1].vertex.x,
								//		delaunay->dcel->vertex[delaunay->graph.nodes[child_Index].points_Index[j]-1].vertex.y);
								if (equal_Point( p, &delaunay->dcel->vertex[delaunay->graph.nodes[child_Index].points_Index[j]-1].vertex))
								{
									(*lowest_Distance) = 0.0;
									q->x = delaunay->dcel->vertex[delaunay->graph.nodes[child_Index].points_Index[j]-1].vertex.x;
									q->y = delaunay->dcel->vertex[delaunay->graph.nodes[child_Index].points_Index[j]-1].vertex.y;
									found = true;
									//printf("Found point (%lf,%lf) equal to (%lf,%lf)\n", q->x, q->y, p->x, p->y);
#ifdef DEBUG_SELECT_CLOSEST_POINT
									printf("Index point %d. Distance %lf\n", delaunay->graph.nodes[child_Index].points_Index[i], (*lowest_Distance));
#endif
								}
							}

							i++;
						}

						// End main loop and force inner loop to finish.
						finished = true;
						i = nChildren;

						if (!found)
						{
							//exit(0);
							printf( "ERROR: No nodes surround new point.\n");
							print_Point( p);
#ifdef LOGGING
							sprintf( log_Text, "ERROR: No nodes surround new point.\n");
							write_Log( log_Text);
#endif
						}
					}
				}
			}
    	}
    }

#ifdef LOCATION_STATISTICS
    location_Stat.end[location_Stat.current] = getTime();
#endif

    return(found);
}

//#define DEBUG_SELECT_CLOSEST_DCEL
/***************************************************************************
* Name: select_Closest_Point_DCEL
* IN:	dcel			triangulation DCEL
* 		nAnchors		# anchors to be used as initial points
* 		p				input point
* OUT:	q				closest point to input "p" point
* 		lowest_Distance	distance from input "p" to output "q"
* RETURN:	True		TRUE if closest point found. False i.o.c.
* Description: finds the closest point in the DCEl to an input "p" point.
***************************************************************************/
bool 	select_Closest_Point_DCEL( struct DCEL_T *dcel, int nAnchors, struct Point_T *p,
												   	   	   	   	   	   struct Point_T *q,
																	   double *lowest_Distance)
{
	bool	error=false;			// Convex hull error flag.
	bool 	foundClosest=false;		// Return value.
	bool 	foundEdge=false;		// Found edge flag.
	int		i=0;					// Loop index.
	int		index=0;				// Array index.
	int		edgeIndex=0;			// DCEL edge index.
	int		selectedVertex=0;		// Vertex selected as initial point.
	int		faceID=0;				// Face identifier.
	double	newDistance=0.0;		// Distance to between points.
	struct Point_T *r1, *r2;		// Points references.
	int		nPoints=0;				// # points in convex hull.
	int		*points;

	// Initialize output value.
	(*lowest_Distance) = DBL_MAX;

#ifdef LOCATION_STATISTICS
    // # triangles searched before a new point is located.
    location_Stat.current++;
    location_Stat.nTrianglesSearched[location_Stat.current] = 0;
#endif

    // Check if point is exterior to convex hull.
    if (!is_Interior_To_Convex_Hull( dcel, p, &error))
    {
#ifdef DEBUG_SELECT_CLOSEST_DCEL
    	printf("Point (%f,%f) is exterior to convex hull\n", p->x, p->y);
#endif
    	// Get points in convex hull.
    	get_Convex_Hull( dcel, &nPoints, points);

    	// Check all points in convex hull.
    	for (i=0; i<nPoints ;i++)
    	{
			// Check if current anchor is closest point.
			newDistance = distance( p, &dcel->vertex[index].vertex);
#ifdef DEBUG_SELECT_CLOSEST_DCEL
			printf("Point %d is (%f,%f). Distance %f\n", index+1,
										dcel->vertex[index].vertex.x,
										dcel->vertex[index].vertex.y,
										newDistance);
#endif
			if (newDistance < (*lowest_Distance))
			{
				// Update output data.
				(*lowest_Distance) = newDistance;
				selectedVertex = index;
#ifdef DEBUG_SELECT_CLOSEST_DCEL
				printf("New point is closer\n");
#endif
			}
    	}
    }
    else
    {
#ifdef DEBUG_SELECT_CLOSEST_DCEL
    	printf("Point (%f,%f) is INTERIOR to convex hull\n", p->x, p->y);
    	printf("Selecting %d random anchors\n", nAnchors);
#endif
    	srand((int) time(NULL));

    	// Select random anchor points.
		for (i=0; i<nAnchors ;i++)
		{
	        // Select a new random point.
			index = rand() % dcel->nVertex;

			// Check if current anchor is closest point.
			newDistance = distance( p, &dcel->vertex[index].vertex);
#ifdef DEBUG_SELECT_CLOSEST_DCEL
			printf("Point %d is (%f,%f). Distance %f\n", index+1,
										dcel->vertex[index].vertex.x,
										dcel->vertex[index].vertex.y,
										newDistance);
#endif
			if (newDistance < (*lowest_Distance))
			{
				// Update output data.
				(*lowest_Distance) = newDistance;
				selectedVertex = index;
#ifdef DEBUG_SELECT_CLOSEST_DCEL
				printf("New point is closer\n");
#endif
			}
		}

		// Get face where selected point is the origin vertex.
		faceID = dcel->edges[dcel->vertex[selectedVertex].origin_Edge-1].face;

		// Main loop to start location from selected anchor.
		foundClosest = false;
		while (!foundClosest)
		{
#ifdef DEBUG_SELECT_CLOSEST_DCEL
			printf("Searching in face %d\n", faceID);
#endif
			// Check if point is internal to face.
			if (is_Interior_To_Face( dcel, p, faceID))
			{
#ifdef DEBUG_SELECT_CLOSEST_DCEL
				printf("Point (%f,%f) is interior to face %d\n", p->x, p->y, faceID);
#endif
				// Initialize output value.
				(*lowest_Distance) = DBL_MAX;

				// Get edge in current face.
				edgeIndex = dcel->faces[faceID].edge - 1;

				// Check all triangle vertex.
				for (i=0; i<N_POINTS_TRIANGLE ;i++)
				{
					// Get vertex.
					index = dcel->edges[edgeIndex].origin_Vertex-1;

					// Check if current vertex is closer.
					newDistance = distance( p, &dcel->vertex[index].vertex);
#ifdef DEBUG_SELECT_CLOSEST_DCEL
					printf("Point %d is (%f,%f). Distance %f\n", index+1,
												dcel->vertex[index].vertex.x,
												dcel->vertex[index].vertex.y,
												newDistance);
#endif
					if (newDistance < (*lowest_Distance))
					{
						// Update output data.
						(*lowest_Distance) = newDistance;
						selectedVertex = index;
#ifdef DEBUG_SELECT_CLOSEST_DCEL
						printf("New point %d is closer\n", index+1);
#endif
					}

					// Get next edge.
					edgeIndex = dcel->edges[edgeIndex].next_Edge - 1;
				}

				// Get reference to closest point.
				q->x = dcel->vertex[selectedVertex].vertex.x;
				q->y = dcel->vertex[selectedVertex].vertex.y;
				foundClosest = true;

#ifdef DEBUG_SELECT_CLOSEST_DCEL
				printf("Closest point index %d is (%f,%f)\n", selectedVertex+1, q->x, q->y);
#endif
			}
			else
			{
				// Get first edge.
				edgeIndex = dcel->faces[faceID].edge - 1;

				// Check edge to continue searching face.
				foundEdge = false;
				while(!foundEdge)
				{
					// Get edge extreme points.
					r1 = &dcel->vertex[dcel->edges[edgeIndex].origin_Vertex-1].vertex;
					r2 = &dcel->vertex[dcel->edges[dcel->edges[edgeIndex].twin_Edge-1].origin_Vertex-1].vertex;

#ifdef DEBUG_SELECT_CLOSEST_DCEL
					printf("Checking edge %d. Points %d %d\n", edgeIndex+1,
							dcel->edges[edgeIndex].origin_Vertex,
							dcel->edges[dcel->edges[edgeIndex].twin_Edge-1].origin_Vertex);
#endif
					// If right turn then this is the edge.
					if (check_Turn( r1, r2, p) == RIGHT_TURN)
					{
						// End inner loop.
						foundEdge = true;

						// Get next edge.
						faceID = dcel->edges[dcel->edges[edgeIndex].twin_Edge-1].face;
#ifdef DEBUG_SELECT_CLOSEST_DCEL
						printf("Right turn. New face is %d\n", faceID);
#endif
					}
					// Get next edge in triangle.
					else
					{
						edgeIndex = dcel->edges[edgeIndex].next_Edge - 1;
#ifdef DEBUG_SELECT_CLOSEST_DCEL
						printf("Left turn. New edge is %d\n", edgeIndex+1);
#endif
					}
				}
			}
		}
    }

	return(foundClosest);
}

//*****************************************************************************
//                      PRIVATE FUNCTION BODIES
//*****************************************************************************
void       insert_First_Node( struct Delaunay_T *delaunay)
{
    struct  Node_T node;            // Temp variable.

    // Insert root triangle.
	node.points_Index[0] = 1;
	node.points_Index[1] = P_MINUS_2;
	node.points_Index[2] = P_MINUS_1;
    node.face_ID = 1;
    insert_Node( &delaunay->graph, &node);
}




//#define DEBUG_IS_INTERIOR_TO_NODE
/***************************************************************************
* Name: is_Interior_To_Node
* IN:	dcel			triangulation DCEL
* 		graph			incremental Delaunay graph
* 		p				input point
* 		node_Index		graph node index
* OUT:	N/A
* RETURN:	True		TRUE if point is interior to node
* Description: check if input point p is interior to "node_Index" node.
***************************************************************************/
int     is_Interior_To_Node(struct DCEL_T *dcel, struct Graph_T *graph, struct Point_T *p, int node_Index)
{
    int     is_Interior=FALSE;				// Return value.
	int		id1=0, id2=0, id3=0;			// IDs of vertex points.

	//Get an edge of the face associated to the node.
	get_Vertex_Of_Node( graph, node_Index, &id1, &id2, &id3);

#ifdef DEBUG_IS_INTERIOR_TO_NODE
	printf(" Vertices %d %d %d\n", id1, id2, id3);
#endif

    // Check if there is not a right turn.
    if (!(return_Turn( dcel, p, id1, id2) == RIGHT_TURN) &&
        !(return_Turn( dcel, p, id2, id3) == RIGHT_TURN) &&
        !(return_Turn( dcel, p, id3, id1) == RIGHT_TURN))
    {
        // Point is interior.
        is_Interior = TRUE;
    }

    return(is_Interior);
}

int is_Strictly_Interior_To_Node(struct DCEL_T *dcel, struct Graph_T *graph, struct Point_T *p, int node_Index)
{
    int     is_Interior=FALSE;				// Return value.
	int		id1=0, id2=0, id3=0;			// IDs of vertex points.

	//Get an edge of the face associated to the node.
	get_Vertex_Of_Node( graph, node_Index, &id1, &id2, &id3);

	// Check if there is always a turn left.
    if ((return_Turn( dcel, p, id1, id2) == LEFT_TURN) &&
        (return_Turn( dcel, p, id2, id3) == LEFT_TURN) &&
        (return_Turn( dcel, p, id3, id1) == LEFT_TURN))
    {
        // Point is interior.
        is_Interior = TRUE;
    }

    return(is_Interior);
}



//#define DEBUG_ANALYZE_FACE
bool analyze_Face( struct Delaunay_T *delaunay, int point_Index)
{
    int     i=0;                            // Loop counter.
    int     node_Index=0;                   // Index of current node analyzed.
    bool	found=false;					// Loop control flag and return value.
	bool	finished=false;					// Loop control flag.
    struct  Dcel_Vertex_T   *point=NULL;   	// Pointer to points in DCEL.

    // Get new point to insert.
    point = get_Vertex( delaunay->dcel, point_Index);

#ifdef DEBUG_ANALYZE_FACE
    printf("*******************\nSearching point %d coordinates (%lf,%lf)\n", point_Index+1,
																				point->vertex.x,
																				point->vertex.y);
#endif

    // Initialize loop control flag.
    finished = FALSE;
    while (!finished)
    {
        // If node is a leaf then triangle found.
    	if (delaunay->graph.nodes[node_Index].nChildren == 0)
        {
#ifdef DEBUG_ANALYZE_FACE
           	printf("Found leaf node %d.\n", node_Index);
#endif

            // Check if point is strictly interior (not over an edge).
            if (is_Strictly_Interior_To_Node(delaunay->dcel, &delaunay->graph, &point->vertex, node_Index))
            {
#ifdef DEBUG_ANALYZE_FACE
            	printf("It is strictly interior\n");
#endif
                // Split current triangle creating 3 triangles.
                split_Triangle( delaunay->dcel, &delaunay->graph, point_Index, node_Index, 3);
            }
            // Point over an edge.
            else
            {
#ifdef DEBUG_ANALYZE_FACE
            	printf("It is over an edge\n");
#endif

                // Split current triangle creating 2 triangles.
                split_Triangle( delaunay->dcel, &delaunay->graph, point_Index, node_Index, 2);
            }

            // Finish main loop.
            finished = TRUE;
        }
        else
        {
            // Search triangle in children nodes.
            i=0;
            found = FALSE;
            while ((!found) && (i < delaunay->graph.nodes[node_Index].nChildren))
            {
#ifdef DEBUG_ANALYZE_FACE
            	printf("Trying %d-child from node %d: node id %d.\n", i, node_Index, delaunay->graph.nodes[node_Index].children_Index[i]);
#endif

                // Check if point is interior to i-child node.
                if (is_Interior_To_Node( delaunay->dcel, &delaunay->graph, &point->vertex,
                						delaunay->graph.nodes[node_Index].children_Index[i]))
                {
#ifdef DELAUNAY_STATISTICS
                	// Update triangle where point is located.
                	delaunay_Stat.trianglesFound[i]++;
#endif

                    // Search in next children node.
                    node_Index = delaunay->graph.nodes[node_Index].children_Index[i];

                    // End loop.
                    found = TRUE;

#ifdef DEBUG_ANALYZE_FACE
                    printf("Interior to node %d. Fetch node data.\n", node_Index);
                    printf("Node vertices are %d %d %d.\n", delaunay->graph.nodes[node_Index].points_Index[0],
															delaunay->graph.nodes[node_Index].points_Index[1],
															delaunay->graph.nodes[node_Index].points_Index[2]);
#endif
                }
                // Next child.
                else
                {
                    i++;

                    // Check if all children checked.
                    if (i == delaunay->graph.nodes[node_Index].nChildren)
                    {
                    	// Force exit function.
                    	finished = TRUE;
                    	i = delaunay->graph.nodes[node_Index].nChildren;
                        printf("CRITICAL ERROR: point (%lf,%lf) is not interior to any node.\n", point->vertex.x, point->vertex.y);
                        printf("Current node %d has %d children.\n", node_Index, delaunay->graph.nodes[node_Index].nChildren);
                        printf("Node vertices are %d %d %d.\n", delaunay->graph.nodes[node_Index].points_Index[0],
																delaunay->graph.nodes[node_Index].points_Index[1],
																delaunay->graph.nodes[node_Index].points_Index[2]);
                        exit(0);
                    }
                }
            }
        }
    }

	return(found);
}

//#define DEBUG_SPLIT_TRIANGLE
void split_Triangle( struct DCEL_T *dcel, struct Graph_T *graph, int point_Index, int node_Index, int nTriangles)
{
    int     new_Edge_ID=0;							// Edge identifier.
    int     new_Face_ID=0;							// Face identifier.
    int     prev_Edge_ID=0;							// Stores previous edge id.
    int     next_Edge_ID=0;							// Stores next edge id.
    int     collinear_Edge_ID=0, collinear_Index=0; // Edge identifier of collinear edge.
	int		flip_Candidates[2];						// Edges that must be cheked due to split operation.
    struct  Dcel_Face_T  *face=NULL;				// Pointer to current face.
	int		old_Node_ID1=0, old_Node_ID2=0;			// Old nodes id.
	struct  Node_T *node=NULL, new_Node[3];		// Nodes.

	double   area[3];

    // Get identifiers of next edge and face to be created.
    new_Edge_ID = get_Number_Edges( dcel) + 1;
    new_Face_ID = get_Number_Faces( dcel);

    // Update edge departing from new point.
    update_Vertex_Edge_At( dcel, new_Edge_ID, point_Index);

    // Get information of node where the new triangles are created.
    node = get_Node( graph, node_Index);

    // Get data of the face of the triangle to be splitted.
    face = get_Face( dcel, node->face_ID);

    // Check number of new triangles to create.
    if (nTriangles == 3)
    {
        // Save previous and next edges ID.
        prev_Edge_ID = dcel->edges[face->edge-1].previous_Edge;
        next_Edge_ID = dcel->edges[face->edge-1].next_Edge;

        // Insert two new edges: new_Edge_ID and new_Edge_ID+1.
        insertEdge( dcel, point_Index+1, new_Edge_ID+5, new_Edge_ID+1, face->edge, node->face_ID);
        insertEdge( dcel, dcel->edges[next_Edge_ID-1].origin_Vertex, new_Edge_ID+2, face->edge,
																	new_Edge_ID, node->face_ID);

        // Insert two new edges: new_Edge_ID+2 and new_Edge_ID+3.
        insertEdge( dcel, point_Index+1, new_Edge_ID+1, new_Edge_ID+3, next_Edge_ID, new_Face_ID);
        insertEdge( dcel, dcel->edges[prev_Edge_ID-1].origin_Vertex, new_Edge_ID+4, next_Edge_ID,
                                                                    new_Edge_ID+2, new_Face_ID);

        // Insert two new edges: new_Edge_ID+4 and new_Edge_ID+5.
        insertEdge( dcel, point_Index+1, new_Edge_ID+3, new_Edge_ID+5, prev_Edge_ID, new_Face_ID+1);
        insertEdge( dcel, dcel->edges[face->edge-1].origin_Vertex, new_Edge_ID, prev_Edge_ID,
                                                                    new_Edge_ID+4, new_Face_ID+1);

        // Update existing edges.
        update_Edge( dcel, NO_UPDATE, NO_UPDATE, new_Edge_ID, new_Edge_ID+1, NO_UPDATE, face->edge-1);
        update_Edge( dcel, NO_UPDATE, NO_UPDATE, new_Edge_ID+2, new_Edge_ID+3, new_Face_ID, next_Edge_ID-1);
        update_Edge( dcel, NO_UPDATE, NO_UPDATE, new_Edge_ID+4, new_Edge_ID+5, new_Face_ID+1, prev_Edge_ID-1);

        // Insert two new faces.
        insertFace( dcel, new_Edge_ID + 2);
        insertFace( dcel, new_Edge_ID + 4);

        // Update leaf node.
        node->children_Index[0] = get_Graph_Length( graph);
        node->children_Index[1] = get_Graph_Length( graph) + 1;
        node->children_Index[2] = get_Graph_Length( graph) + 2;
        update_Node( graph, node_Index, 3, node);

		// Insert three new nodes.
		new_Node[0].points_Index[0] = get_Edge_Origin_Vertex( dcel, new_Edge_ID);
		new_Node[0].points_Index[1] = get_Edge_Origin_Vertex( dcel, new_Edge_ID - 1);
		new_Node[0].points_Index[2] = get_Edge_Origin_Vertex( dcel, face->edge - 1);
        new_Node[0].face_ID = node->face_ID;
        area[0] = modified_Signed_Area( dcel, &new_Node[0]);
		new_Node[1].points_Index[0] = get_Edge_Origin_Vertex( dcel, new_Edge_ID + 2);
		new_Node[1].points_Index[1] = get_Edge_Origin_Vertex( dcel, new_Edge_ID + 1);
		new_Node[1].points_Index[2] = get_Edge_Origin_Vertex( dcel, next_Edge_ID - 1);
        new_Node[1].face_ID = new_Face_ID;
        area[1] = modified_Signed_Area( dcel, &new_Node[1]);
		new_Node[2].points_Index[0] = get_Edge_Origin_Vertex( dcel, new_Edge_ID + 4);
		new_Node[2].points_Index[1] = get_Edge_Origin_Vertex( dcel, new_Edge_ID + 3);
		new_Node[2].points_Index[2] = get_Edge_Origin_Vertex( dcel, prev_Edge_ID - 1);
        new_Node[2].face_ID = new_Face_ID + 1;
        area[2] = modified_Signed_Area( dcel, &new_Node[2]);
        if (area[0] > area[1])
        {
            // 2nd, 0th, 1st.
            if (area[2] > area[0])
            {
                insert_Node( graph, &new_Node[2]);
                insert_Node( graph, &new_Node[0]);
                insert_Node( graph, &new_Node[1]);
            }
            // 0th, 2nd, 1st.
            else
            {
                insert_Node( graph, &new_Node[0]);
                if (area[2] > area[1])
                {

                    insert_Node( graph, &new_Node[2]);
                    insert_Node( graph, &new_Node[1]);
                }
                // 0th, 1st, 2nd.
                else
                {
                    insert_Node( graph, &new_Node[1]);
                    insert_Node( graph, &new_Node[2]);
                }
            }
        }
        else
        {
            // 2nd, 1st, 0th.
            if (area[2] > area[1])
            {
                insert_Node( graph, &new_Node[2]);
                insert_Node( graph, &new_Node[1]);
                insert_Node( graph, &new_Node[0]);
            }
            else
            {
                insert_Node( graph, &new_Node[1]);
                if (area[2] > area[0])
                {
                    insert_Node( graph, &new_Node[2]);
                    insert_Node( graph, &new_Node[0]);
                }
                else
                {
                    insert_Node( graph, &new_Node[0]);
                    insert_Node( graph, &new_Node[2]);
                }
            }
        }

#ifdef DEBUG
		print_Graph(graph);
        print_DCEL( dcel);
#endif

		// Check if edges must be flipped.
		check_Edge( dcel, graph, face->edge);
		check_Edge( dcel, graph, prev_Edge_ID);
		check_Edge( dcel, graph, next_Edge_ID);
    }
    else
    {
        // Get edge identifier where new point is collinear.
        collinear_Edge_ID = select_Colinear_Edge( dcel, point_Index, face->edge);
		if (collinear_Edge_ID != -1)
		{
			collinear_Index = collinear_Edge_ID - 1;

			// Save previous and next edges ID.
			prev_Edge_ID = dcel->edges[collinear_Index].previous_Edge;
			next_Edge_ID = dcel->edges[collinear_Index].next_Edge;

#ifdef DEBUG_SPLIT_TRIANGLE
			printf("Collinear edges are %d and %d. Origins are %d and %d\n",
					collinear_Edge_ID,
					dcel->edges[collinear_Index].twin_Edge,
					dcel->edges[collinear_Index].origin_Vertex,
					dcel->edges[dcel->edges[collinear_Index].twin_Edge-1].origin_Vertex);
			printf("Faces that share edge are:\n");
			printFace( dcel, dcel->edges[collinear_Index].face);
			printFace( dcel, dcel->edges[dcel->edges[collinear_Index].twin_Edge-1].face);
#endif

			// Store edges that must be checked after split of first triangle operation.
			flip_Candidates[0] = next_Edge_ID;
			flip_Candidates[1] = prev_Edge_ID;

			// Store nodes ID that are going to be updated.
			old_Node_ID1 = get_Node_Assigned( graph, dcel->edges[collinear_Index].face);
			old_Node_ID2 = get_Node_Assigned( graph, dcel->edges[dcel->edges[collinear_Index].twin_Edge-1].face);

			// Update current face with new edge: new_Edge_ID.
			insertEdge( dcel, dcel->edges[prev_Edge_ID-1].origin_Vertex, new_Edge_ID+1, next_Edge_ID,
																				collinear_Edge_ID, node->face_ID);

			// Insert a new face with two new edges: new_Edge_ID+1 and new_Edge_ID+2.
			insertEdge( dcel, point_Index+1, new_Edge_ID, new_Edge_ID+2, prev_Edge_ID, new_Face_ID);
			insertEdge( dcel, dcel->edges[collinear_Index].origin_Vertex, new_Edge_ID+3, prev_Edge_ID,
																				new_Edge_ID+1, new_Face_ID);
			update_Vertex_Edge_At( dcel, new_Edge_ID+1, point_Index);
			update_Face( dcel, new_Edge_ID, node->face_ID);

			// Update existing edges.
			update_Edge( dcel, point_Index+1, NO_UPDATE, new_Edge_ID, NO_UPDATE, NO_UPDATE, collinear_Index);
			update_Edge( dcel, NO_UPDATE, NO_UPDATE, NO_UPDATE, new_Edge_ID, NO_UPDATE, next_Edge_ID-1);
			update_Edge( dcel, NO_UPDATE, NO_UPDATE, new_Edge_ID+1, new_Edge_ID+2, new_Face_ID, prev_Edge_ID-1);

			// Get node of current edge and update it.
			node = get_Node( graph, old_Node_ID1);
			node->children_Index[0] = get_Graph_Length( graph);
			node->children_Index[1] = get_Graph_Length( graph) + 1;
			node->children_Index[2] = INVALID;
			update_Node( graph, old_Node_ID1, 2, node);

			// Insert two new nodes in first node splitted.
			new_Node[0].points_Index[0] = get_Edge_Origin_Vertex( dcel, dcel->edges[collinear_Index].previous_Edge - 1);
			new_Node[0].points_Index[1] = get_Edge_Origin_Vertex( dcel, collinear_Index);
			new_Node[0].points_Index[2] = get_Edge_Origin_Vertex( dcel, dcel->edges[collinear_Index].next_Edge - 1);
			new_Node[0].face_ID = dcel->edges[collinear_Index].face;
			insert_Node( graph, &new_Node[0]);
#ifdef DEBUG_SPLIT_TRIANGLE
			printf("Splitting first triangle\n");
			printFace( dcel, new_Node[0].face_ID);
#endif
			new_Node[1].points_Index[0] = get_Edge_Origin_Vertex( dcel, dcel->edges[collinear_Index].previous_Edge - 1);
			new_Node[1].points_Index[1] = get_Edge_Origin_Vertex( dcel, dcel->edges[dcel->edges[dcel->edges[collinear_Index].previous_Edge-1].twin_Edge-1].previous_Edge-1);
			new_Node[1].points_Index[2] = get_Edge_Origin_Vertex( dcel, dcel->edges[dcel->edges[collinear_Index].previous_Edge - 1].twin_Edge-1);
			new_Node[1].face_ID = new_Face_ID;
			insert_Node( graph, &new_Node[1]);

			// Insert new face.
			insertFace( dcel, new_Edge_ID + 2);
#ifdef DEBUG_SPLIT_TRIANGLE
			printf("Splitting second triangle\n");
			printFace( dcel, new_Node[1].face_ID);
#endif

			// Update twin face.
			collinear_Edge_ID = dcel->edges[collinear_Index].twin_Edge;
			collinear_Index	 = collinear_Edge_ID-1;
			prev_Edge_ID	 = dcel->edges[collinear_Index].previous_Edge;
			next_Edge_ID	 = dcel->edges[collinear_Index].next_Edge;

			// Insert a new face with two new edges: new_Edge_ID+3 and new_Edge_ID+4.
			insertEdge( dcel, point_Index+1, new_Edge_ID+2, new_Edge_ID+4, next_Edge_ID, new_Face_ID+1);
			insertEdge( dcel, dcel->edges[prev_Edge_ID-1].origin_Vertex, new_Edge_ID+5, next_Edge_ID,
																				new_Edge_ID+3, new_Face_ID+1);

			// Update current face with new edge: new_Edge_ID+5.
			insertEdge( dcel, point_Index+1, new_Edge_ID+4, collinear_Edge_ID, prev_Edge_ID, dcel->edges[collinear_Index].face);

			// Update existing edges.
			update_Edge( dcel, NO_UPDATE, NO_UPDATE, new_Edge_ID+3, new_Edge_ID+4, new_Face_ID+1, next_Edge_ID-1);
			update_Edge( dcel, NO_UPDATE, NO_UPDATE, NO_UPDATE, new_Edge_ID+5, NO_UPDATE, collinear_Index);
			update_Edge( dcel, NO_UPDATE, NO_UPDATE, new_Edge_ID+5, NO_UPDATE, NO_UPDATE, prev_Edge_ID-1);

			// Get node of twin edge and update it.
			node = get_Node( graph, old_Node_ID2);
			node->children_Index[0] = get_Graph_Length( graph);
			node->children_Index[1] = get_Graph_Length( graph) + 1;
			node->children_Index[2] = INVALID;
			update_Node( graph, old_Node_ID2, 2, node);

			// Insert two new nodes in first node splitted.
			new_Node[0].points_Index[0] = get_Edge_Origin_Vertex( dcel, dcel->edges[collinear_Index].previous_Edge - 1);
			new_Node[0].points_Index[1] = get_Edge_Origin_Vertex( dcel, collinear_Index);
			new_Node[0].points_Index[2] = get_Edge_Origin_Vertex( dcel, dcel->edges[collinear_Index].next_Edge - 1);
			new_Node[0].face_ID = dcel->edges[collinear_Index].face;
			insert_Node( graph, &new_Node[0]);
			update_Face( dcel, collinear_Edge_ID, new_Node[0].face_ID);
#ifdef DEBUG_SPLIT_TRIANGLE
			printf("Splitting third triangle\n");
			printFace( dcel, new_Node[0].face_ID);
#endif
			new_Node[1].points_Index[0] = get_Edge_Origin_Vertex( dcel, dcel->edges[next_Edge_ID - 1].previous_Edge - 1);
			new_Node[1].points_Index[1] = get_Edge_Origin_Vertex( dcel, next_Edge_ID - 1);
			new_Node[1].points_Index[2] = get_Edge_Origin_Vertex( dcel, dcel->edges[next_Edge_ID - 1].next_Edge - 1);
			new_Node[1].face_ID = dcel->edges[next_Edge_ID - 1].face;
			insert_Node( graph, &new_Node[1]);

			// Insert new face.
			insertFace( dcel, new_Edge_ID + 4);

#ifdef DEBUG_SPLIT_TRIANGLE
			printf("Splitting forth triangle\n");
			printFace( dcel, new_Node[1].face_ID);
#endif

			// Check candidates from first triangle.
			check_Edge( dcel, graph, flip_Candidates[0]);
			check_Edge( dcel, graph, flip_Candidates[1]);

			// Check candidates from second triangle.
			check_Edge( dcel, graph, prev_Edge_ID);
			check_Edge( dcel, graph, next_Edge_ID);
		}
#ifdef DEBUG_SPLIT_TRIANGLE
		else
		{
			printf("Collinear edge is negative. Point index %d and face edge %d\n", point_Index, face->edge);
		}
#endif
    }
}


double   modified_Signed_Area( struct DCEL_T *dcel, struct  Node_T *node)
{
    double   area=0.0;           // Return value.

    // Check if any of the vertex is not real: P_MINUS_! or P_MINUS_2.
    if ((node->points_Index[0] < 0) ||
        (node->points_Index[1] < 0) ||
        (node->points_Index[2] < 0))
    {
        // Set area zero.
        area = 0.0;
    }
    else
    {
        // Compute area.
        area = signed_Area( &dcel->vertex[node->points_Index[0]-1].vertex,
                            &dcel->vertex[node->points_Index[1]-1].vertex,
                            &dcel->vertex[node->points_Index[2]-1].vertex);
    }

    return(area);
}

//#define DEBUG_SELECT_COLLINEAR_EDGE
int     select_Colinear_Edge( struct DCEL_T *dcel, int point_Index, int edge_ID)
{
	int		edge_Index=0;					// Index of edge of current face.
    int     collinear_Edge=0;				// Return value.
	int		id1=0, id2=0, id3=0;			// IDs of vertex points.
    struct  Point_T *p=NULL;				// Reference to new point.

	// Set edge index.
	edge_Index = edge_ID - 1;

	// Get index of the vertex of current face.
	id1 = get_Edge_Origin_Vertex( dcel, dcel->edges[edge_Index].previous_Edge-1);
	id2 = get_Edge_Origin_Vertex( dcel, edge_Index);
	id3 = get_Edge_Origin_Vertex( dcel, dcel->edges[edge_Index].next_Edge-1);

#ifdef DEBUG_SELECT_COLLINEAR_EDGE
	printf("Triangles points are %d, %d and %d\n", id1, id2, id3);
#endif

    // Get point coordinates.
    p  = get_Vertex_Point( dcel, point_Index);

	// Check if point is collinear to previous edge.
	if (return_Turn( dcel, p, id1, id2) == COLINEAR)
	{
#ifdef DEBUG_SELECT_COLLINEAR_EDGE
		printf("Collinear to points %d and %d\n", id1, id2);
#endif
		collinear_Edge = dcel->edges[edge_Index].previous_Edge;
	}
	// Check if point is collinear to input edge.
	else if (return_Turn( dcel, p, id2, id3) == COLINEAR)
	{
#ifdef DEBUG_SELECT_COLLINEAR_EDGE
		printf("Collinear to points %d and %d\n", id2, id3);
#endif
		collinear_Edge = edge_ID;
	}
	// Check if point is collinear to next edge.
	else if (return_Turn( dcel, p, id3, id1) == COLINEAR)
	{
#ifdef DEBUG_SELECT_COLLINEAR_EDGE
		printf("Collinear to points %d and %d\n", id3, id1);
#endif
		collinear_Edge = dcel->edges[edge_Index].next_Edge;
	}
	else
	{
		printf("Checking when point is collinear but it is not.\n");
#ifdef DEBUG
		sprintf( log_Text, "Checking when point is collinear but it is not.\n");
		write_Log( log_Text);
#endif
		collinear_Edge = -1;
	}

#ifdef DEBUG_SELECT_COLLINEAR_EDGE
	printf("Collinear edge is %d\n", collinear_Edge);
#endif

	return(collinear_Edge);
}



void	flip_Edges_Dcel( struct DCEL_T *dcel, struct Graph_T *graph, int edge_ID)
{
	int		temp=0;										// Temp variable.
	int		edge_Index=0;								// Edge index.
	struct  Node_T   *node=NULL, new_Node;
	int		old_Node_ID1=0, old_Node_ID2=0;				// Old nodes id.
	struct	Dcel_Edge_T *twin=NULL, *edge=NULL;

	// Get edge index.
	edge_Index = edge_ID-1;

	// Get edge and twin edge information.
	edge = get_Edge( dcel, edge_Index);
	twin = get_Edge( dcel, dcel->edges[edge_Index].twin_Edge-1);

#ifdef DEBUG
	printf("FLIPPING EDGES %d and %d\n", edge_Index+1, dcel->edges[edge_Index].twin_Edge);
	if ((edge->face > get_Number_Faces( dcel)) || (twin->face > get_Number_Faces( dcel)))
	{
		printf("Higher than number of faces\n");
		exit(0);
	}
#endif

	// Store nodes ID that are going to be updated to internal nodes.
	old_Node_ID1 = get_Node_Assigned( graph, edge->face);
	old_Node_ID2 = get_Node_Assigned( graph, twin->face);

#ifdef DEBUG
	if ((old_Node_ID1 <= 0) || (old_Node_ID1 > get_Graph_Length( graph)) ||
		(old_Node_ID2 <= 0) || (old_Node_ID2 > get_Graph_Length( graph)))
	{
		exit(0);
	}
#endif

	// Update vertex of flipped edge.
	if (edge->origin_Vertex > 0)
	{
		update_Vertex_Edge_At( dcel, twin->next_Edge, edge->origin_Vertex-1);
	}
	if (dcel->edges[edge->next_Edge-1].origin_Vertex > 0)
	{
		update_Vertex_Edge_At( dcel, edge->next_Edge, dcel->edges[edge->next_Edge-1].origin_Vertex-1);
	}

	// Update origin of current and twin edges.
	temp = dcel->edges[edge->previous_Edge-1].origin_Vertex;
	dcel->edges[edge_Index].origin_Vertex = dcel->edges[twin->previous_Edge-1].origin_Vertex;
	dcel->edges[edge->twin_Edge-1].origin_Vertex = temp;

	// Update next edges.
	dcel->edges[edge->next_Edge-1].next_Edge     = edge->twin_Edge;
	dcel->edges[twin->next_Edge-1].next_Edge     = edge_ID;
	dcel->edges[edge->previous_Edge-1].next_Edge = twin->next_Edge;
	dcel->edges[twin->previous_Edge-1].next_Edge = edge->next_Edge;
	dcel->edges[edge_Index].next_Edge            = edge->previous_Edge;
	dcel->edges[edge->twin_Edge-1].next_Edge     = twin->previous_Edge;

	// Update previous edges.
	dcel->edges[edge_Index].previous_Edge            = dcel->edges[dcel->edges[edge_Index].next_Edge-1].next_Edge;
	dcel->edges[edge->twin_Edge-1].previous_Edge     = dcel->edges[twin->next_Edge-1].next_Edge;
	dcel->edges[edge->next_Edge-1].previous_Edge     = edge_ID;
	dcel->edges[twin->next_Edge-1].previous_Edge     = edge->twin_Edge;
	dcel->edges[edge->previous_Edge-1].previous_Edge = edge->next_Edge;
	dcel->edges[twin->previous_Edge-1].previous_Edge = twin->next_Edge;

	// Update faces of edges that have moved to another face.
	dcel->edges[edge->previous_Edge-1].face = edge->face;
	dcel->edges[twin->previous_Edge-1].face = twin->face;

	// Update faces.
	dcel->faces[edge->face].edge = edge_ID;
	dcel->faces[twin->face].edge = edge->twin_Edge;

	// Get node of current edge and update it.
	node = get_Node( graph, old_Node_ID1);
    node->children_Index[0] = get_Graph_Length( graph);
    node->children_Index[1] = get_Graph_Length( graph) + 1;
    node->children_Index[2] = INVALID;
    update_Node( graph, old_Node_ID1, 2, node);

	// Get node of twin edge and update it.
	node = get_Node( graph, old_Node_ID2);
    node->children_Index[0] = get_Graph_Length( graph);
    node->children_Index[1] = get_Graph_Length( graph) + 1;
    node->children_Index[2] = INVALID;
    update_Node( graph, old_Node_ID2, 2, node);

    // Insert two new nodes.
	new_Node.points_Index[0] = get_Edge_Origin_Vertex( dcel, edge->previous_Edge - 1);
	new_Node.points_Index[1] = get_Edge_Origin_Vertex( dcel, edge_Index);
	new_Node.points_Index[2] = get_Edge_Origin_Vertex( dcel, edge->next_Edge - 1);
    new_Node.face_ID = edge->face;
    insert_Node( graph, &new_Node);
	new_Node.points_Index[0] = get_Edge_Origin_Vertex( dcel, twin->previous_Edge - 1);
	new_Node.points_Index[1] = get_Edge_Origin_Vertex( dcel, edge->twin_Edge - 1);
	new_Node.points_Index[2] = get_Edge_Origin_Vertex( dcel, twin->next_Edge - 1);
    new_Node.face_ID = twin->face;
    insert_Node( graph, &new_Node);

	// Check recursively edges that could be illegal.
	check_Edge( dcel, graph, edge->previous_Edge);
	check_Edge( dcel, graph, twin->next_Edge);
}

//#define DEBUG_CHECK_EDGES
void	check_Edge( struct DCEL_T *dcel, struct Graph_T *graph, int edge_ID)
{
	int		flip_Edges=FALSE;							// Flip needed flag.
	int		edge_Index=0;								// Edge index.
	struct	Point_T *common1=NULL, *common2=NULL, *p=NULL, *q=NULL;

	// Get edge inde.
	edge_Index = edge_ID-1;

	// Check if the edge is NOT in the external face.
	flip_Edges = FALSE;
	if (!is_External_Edge( dcel, edge_Index))
	{
		// Check if any of the vertex of current edge is P-2 or P-1.
		if (is_Negative_Any_Vertex( dcel, edge_ID))
		{
			// Check if any of the vertex of incident faces is P-2 or P-1.
			if ((dcel->edges[dcel->edges[edge_Index].previous_Edge-1].origin_Vertex > 0) &&
				(dcel->edges[dcel->edges[dcel->edges[edge_Index].twin_Edge-1].previous_Edge-1].origin_Vertex > 0))
			{
				// Get points of incident faces to current edge.
				p = get_Vertex_Point( dcel, dcel->edges[dcel->edges[edge_Index].previous_Edge-1].origin_Vertex-1);
				q = get_Vertex_Point( dcel, dcel->edges[dcel->edges[dcel->edges[edge_Index].twin_Edge-1].previous_Edge-1].origin_Vertex-1);

				// Set p as the vertex with highest y-coordinate.
				if (p->y < q->y)
				{
					p = get_Vertex_Point( dcel, dcel->edges[dcel->edges[dcel->edges[edge_Index].twin_Edge-1].previous_Edge-1].origin_Vertex-1);
					q = get_Vertex_Point( dcel, dcel->edges[dcel->edges[edge_Index].previous_Edge-1].origin_Vertex-1);
				}

				// Origin vertex is negative.
				if (dcel->edges[edge_Index].origin_Vertex < 0)
				{
					// Get destination point of current edge.
					common1 = get_Vertex_Point( dcel, dcel->edges[dcel->edges[edge_Index].twin_Edge-1].origin_Vertex-1);

					// Check if negative vertex is P-2.
					if (dcel->edges[edge_Index].origin_Vertex == P_MINUS_2)
					{
						// If turn LEFT_TURN then flip edge.
						if (check_Turn( p, q, common1) == LEFT_TURN)
						{
							flip_Edges = TRUE;
						}
					}
					else
					{
						// If turn RIGHT then flip edge.
						if (check_Turn( p, q, common1) == RIGHT_TURN)
						{
							flip_Edges = TRUE;
						}
					}
				}
				// Destination vertex is negative.
				else
				{
					// Get origin point of current edge.
					common1 = get_Vertex_Point( dcel, dcel->edges[edge_Index].origin_Vertex-1);

					// Check if negative vertex is P-2.
					if (dcel->edges[dcel->edges[edge_Index].twin_Edge-1].origin_Vertex == P_MINUS_2)
					{
						// If turn LEFT_TURN then flip edge.
						if (check_Turn( p, q, common1) == LEFT_TURN)
						{
							flip_Edges = TRUE;
						}
					}
					else
					{
						// If turn RIGHT then flip edge.
						if (check_Turn( p, q, common1) == RIGHT_TURN)
						{
							flip_Edges = TRUE;
						}
					}
				}
			}
		}
		// Vertex of candidate edge are positive.
		else
		{
			// Check if any of the other points of the triangles are P-2 or P-1.
			if ((dcel->edges[dcel->edges[edge_Index].previous_Edge-1].origin_Vertex > 0) &&
				(dcel->edges[dcel->edges[dcel->edges[edge_Index].twin_Edge-1].previous_Edge-1].origin_Vertex > 0))
			{
				// All points are positive -> normal in circle check.
				// Get points of edge to flip.
				common1 = get_Vertex_Point( dcel, dcel->edges[edge_Index].origin_Vertex-1);
				common2 = get_Vertex_Point( dcel, dcel->edges[dcel->edges[edge_Index].next_Edge-1].origin_Vertex-1);

				// Get points of faces.
				p = get_Vertex_Point( dcel, dcel->edges[dcel->edges[edge_Index].previous_Edge-1].origin_Vertex-1);
				q = get_Vertex_Point( dcel, dcel->edges[dcel->edges[dcel->edges[edge_Index].twin_Edge-1].previous_Edge-1].origin_Vertex-1);

				// Check if edge must be flipped.
				if (in_Circle( common1, common2, p, q))
				{
					flip_Edges = TRUE;
				}
			}
		}
	}

	// Check if edges must be flipped.
	if (flip_Edges)
	{
#ifdef DELAUNAY_STATISTICS
		// Update triangle where point is located.
		delaunay_Stat.nFlipped++;
#endif
#ifdef DEBUG_CHECK_EDGES
		printf("Flipping edge %d\n", edge_ID);
#endif
		// Flip edges.
		flip_Edges_Dcel( dcel, graph, edge_ID);
	}
}
