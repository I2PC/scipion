#include "defines.h"
#include "delaunay.h"
#include <float.h>
#include "graph.h"
#ifdef LOGGING
#include "log.h"
#endif
#include "sorting.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#ifdef DEBUG
#include <GL/glut.h>
#endif

#define     EXTERNAL_FACE       0


/*****************************************************************************
* Variables declaration
*****************************************************************************/
struct Node_T 	root_Node;
struct Graph_T  graph;

/*****************************************************************************
* Private functions declaration
*****************************************************************************/
void 		insert_First_Node( struct Delaunay_T *delaunay);
int			is_Interior_To_Node(struct DCEL_T *dcel, struct Graph_T *graph, struct Point_T *p, int node_Index);
int         is_Strictly_Interior_To_Node(struct DCEL_T *dcel, struct Graph_T *graph, struct Point_T *p, int node_Index);
void		analyze_Face( struct Delaunay_T *delaunay, int point_Index);
void    	split_Triangle( struct DCEL_T *dcel, struct Graph_T *graph, int point_Index, int node_Index, int nTriangles);
double  	modified_Signed_Area( struct DCEL_T *dcel, struct  Node_T *node);
int     	select_Colinear_Edge( struct DCEL_T *dcel, int point_Index, int edge_ID);
void		check_Edge( struct DCEL_T *dcel, struct Graph_T *graph, int edge_ID);
void		flip_Edges_Dcel( struct DCEL_T *dcel, struct Graph_T *graph, int edge_ID);


/*****************************************************************************
* Public functions declaration
*****************************************************************************/
int		init_Delaunay( struct Delaunay_T *delaunay, int nPoints)
{
	int		ret=SUCCESS; 		// Return value.

	// Allocate DCEL.
	delaunay->dcel = (struct DCEL_T *)malloc(sizeof(struct DCEL_T));

	// Create DCEL data to store nPoints.
	if (initialize_DCEL( delaunay->dcel, nPoints, nPoints*8, (nPoints+2)*2) == FAILURE)
	{
		printf("Error allocating memory when calling initialize_DCEL");
		ret = FAILURE;
	}
	else
	{
		// Initialize graph.
		ret = initialize_Graph( &delaunay->graph, delaunay->dcel->sizeVertex*10);
	}

	return(ret);
}


void	delete_Delaunay( struct Delaunay_T *delaunay)
{
	// Deallocate DCEL.
	finalize_DCEL(delaunay->dcel);
	free(delaunay->dcel);

	// Deallocate Delaunay.
	finalize_Delaunay(delaunay);
}

void	insert_Point( struct Delaunay_T *delaunay, double x, double y)
{
	// Update coordinates.
	delaunay->dcel->vertex[delaunay->dcel->nVertex].vertex.x = x;
	delaunay->dcel->vertex[delaunay->dcel->nVertex].vertex.y = y;

	// Increase # points.
	delaunay->dcel->nVertex++;
}


void	create_Delaunay_Triangulation( struct Delaunay_T *delaunay)
{
	// Build Delaunay triangulation from scratch.
	incremental_Delaunay( delaunay);

	// Remove P_MINUS_2 and P_MINUS_1 data from dcel.
	purge_Delaunay( delaunay);
}


int 	initialize_Delaunay(struct Delaunay_T *delaunay, struct DCEL_T *dcel)
{
	int		ret=SUCCESS;		// Return value.

	// Get reference to DCEL.
	delaunay->dcel = dcel;

	// Initialize graph.
	ret = initialize_Graph( &delaunay->graph, delaunay->dcel->nVertex*10);

	return(ret);
}


void 	finalize_Delaunay(struct Delaunay_T *delaunay)
{
	// Dereference DCEL.
	delaunay->dcel = NULL;

	// Delete graph.
	finalize_Graph( &delaunay->graph);
}


void	incremental_Delaunay(struct Delaunay_T *delaunay)
{
    int     point_Index=0;                  // Points loop counter.
#ifdef DEBUG
//    char    filename[20];                   // Debug filename.
#endif

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

#ifdef DEBUG
    print_Graph(graph);
    print_DCEL( dcel);
#endif

    // Loop all other points.
    for (point_Index=1; point_Index<delaunay->dcel->nVertex ; point_Index++)
    {
        // Insert new point into triangle where it is located.
        analyze_Face( delaunay, point_Index);
    }
}


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

	// Initialize variables.
	i=0;
	nPending = dcel->nEdges;

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
    //int saved_Index=0;          // Next edge to update.
    //int first_Index=0;          // Index of first edge updated.
    //int point_Index=0;          // Index of next point in convex hull.

    // Initialize edge and point indexes.
    //point_Index = 0;
    edge_Index = delaunay->dcel->edges[edge_Index].previous_Edge-1;
    //first_Index = edge_Index;

    do
    {
        // Update edge departing from point.
    	//delaunay->dcel->vertex[point_Index].origin_Edge = delaunay->dcel->edges[edge_Index].twin_Edge;

        // Invalid current face.
        update_Face( delaunay->dcel, INVALID, delaunay->dcel->edges[edge_Index].face);

        // Update edge so it belongs to convex hull.
        //delaunay->dcel->edges[edge_Index].face = EXTERNAL_FACE;

        // Get index of edge departing from next point in convex hull.
        //saved_Index = edge_Index;
        edge_Index = delaunay->dcel->edges[edge_Index].previous_Edge-1;
        edge_Index = delaunay->dcel->edges[edge_Index].twin_Edge-1;

        // Get point index.
        //point_Index = delaunay->dcel->edges[edge_Index].origin_Vertex-1;

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

        // Update next and previous links.
        //delaunay->dcel->edges[edge_Index].next_Edge = saved_Index + 1;
        //delaunay->dcel->edges[saved_Index].previous_Edge = edge_Index + 1;
    } while(delaunay->dcel->edges[edge_Index].origin_Vertex != 1);
    // Until top most point found again.

    // Invalid current face.
    update_Face( delaunay->dcel, INVALID, delaunay->dcel->edges[edge_Index].face);

    // Update edge so it belongs to convex hull.
    //delaunay->dcel->edges[edge_Index].face = EXTERNAL_FACE;

    // Set edge in external face.
    //delaunay->dcel->faces[EXTERNAL_FACE].edge = edge_Index+1;

    // Update next and previous links to close convex hull.
    //delaunay->dcel->edges[first_Index].next_Edge = edge_Index + 1;
    //delaunay->dcel->edges[edge_Index].previous_Edge = first_Index + 1;
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


//#define DEBUG_SELECT_CLOSEST_POINT
/***************************************************************************
* Name: select_Closest_Point
* IN:	delaunay				delaunay triangulation
* 		p						input point
* OUT:		q					closest point to input "p" point
* 			lowest_Distance		distance from input "p" to output "q"
* RETURN:	True				TRUE if closest point found. False i.o.c.
* Description: finds the closest point in the DCEl to an input "p" point.
***************************************************************************/
int select_Closest_Point( struct Delaunay_T *delaunay, struct Point_T *p, struct Point_T *q,
																	double *lowest_Distance)
{
    int     i=0;				// Loop counter.
    int     node_Index=0;  		// Index of current node analyzed.
    int		found=FALSE;		// Loop control flag.
    int		finished=FALSE;		// Loop control flag.
	int		child_Index=0;		// Children node ID.
    int     nChildren=0;        // Number of children of current node.
    double	new_Distance=0.0;	// Temporary distance.
    int		face_Id=0;			// Face identifier.
    int		edge_Id=0;			// Edge identifier.
    int 	nImaginary=0;
    struct 	Dcel_Face_T	*face;	// Face data.
	int 	currentPoint=0, firstPoint=0;

#ifdef DEBUG_SELECT_CLOSEST_POINT
	printf("Point to search: ");
	print_Point( p);
#endif

    // Loop while not leaf node found and not end of set of points and not.
    finished = FALSE;
    while (!finished)
    {
    	// If node is a leaf then triangle that surrounds point has been found.
    	if (is_Leaf_Node( &delaunay->graph, node_Index))
    	{
#ifdef DEBUG_SELECT_CLOSEST_POINT
    		printf("leaf node found: %d\n", node_Index);
#endif

			// Gets face assigned to node.
			face_Id = get_Face_Of_Node(&delaunay->graph, node_Index);
			face = get_Face( delaunay->dcel, face_Id);
			edge_Id = face->edge;

			// End main loop and update return value.
			finished = TRUE;
			found = TRUE;

			// Check if any point is imaginary.
			nImaginary = number_Imaginary_Points( &delaunay->graph, node_Index);
			if (nImaginary > 0)
			{
#ifdef DEBUG_SELECT_CLOSEST_POINT
				printf("Imaginary triangle with %d imaginary points\n", nImaginary);
#endif

				// Gets face assigned to node.
				face_Id = get_Face_Of_Node(&delaunay->graph, node_Index);
				face = get_Face( delaunay->dcel, face_Id);
				edge_Id = face->edge;

				// Get an edge with a positive origin point.
				while (delaunay->dcel->edges[edge_Id-1].origin_Vertex < 0)
				{
					edge_Id = delaunay->dcel->edges[edge_Id-1].next_Edge;
				}

#ifdef DEBUG_SELECT_CLOSEST_POINT
				printf("Found edge %d with non negative origin %d ", edge_Id, delaunay->dcel->edges[edge_Id-1].origin_Vertex);
				printf("and destination %d\n", delaunay->dcel->edges[delaunay->dcel->edges[edge_Id-1].twin_Edge-1].origin_Vertex);
#endif
				currentPoint = delaunay->dcel->edges[edge_Id-1].origin_Vertex;

				// Compute first point distance.
				(*lowest_Distance) = distance( p, &delaunay->dcel->vertex[currentPoint-1].vertex);
				q->x = delaunay->dcel->vertex[currentPoint-1].vertex.x;
				q->y = delaunay->dcel->vertex[currentPoint-1].vertex.y;

				// Check all neighbor points. Check first to set lowest value.
				edge_Id = delaunay->dcel->edges[edge_Id-1].twin_Edge;
				currentPoint = delaunay->dcel->edges[edge_Id-1].origin_Vertex;
				firstPoint = currentPoint;
				new_Distance = distance( p, &delaunay->dcel->vertex[currentPoint-1].vertex);
				if (new_Distance < (*lowest_Distance))
				{
					// Update new closest point and its distance.
					(*lowest_Distance) = new_Distance;
					q->x = delaunay->dcel->vertex[currentPoint-1].vertex.x;
					q->y = delaunay->dcel->vertex[currentPoint-1].vertex.y;
#ifdef DEBUG_SELECT_CLOSEST_POINT
					printf("Closest point %d with distance %lf\n", currentPoint, (*lowest_Distance));
#endif
				}

				edge_Id = delaunay->dcel->edges[edge_Id-1].next_Edge;
				edge_Id = delaunay->dcel->edges[edge_Id-1].twin_Edge;
				currentPoint = delaunay->dcel->edges[edge_Id-1].origin_Vertex;
#ifdef DEBUG_SELECT_CLOSEST_POINT
				printf("First point is %d.\nCurrent edge %d has point %d\n", firstPoint, edge_Id, currentPoint);
#endif
				// Check all points departing from first point.
				while (currentPoint != firstPoint)
				{
					// Skip imaginary points.
					if (currentPoint > 0)
					{
						new_Distance = distance( p, &delaunay->dcel->vertex[currentPoint-1].vertex);
						if (new_Distance < (*lowest_Distance))
						{
							// Update new closest point and its distance.
							(*lowest_Distance) = new_Distance;
							q->x = delaunay->dcel->vertex[currentPoint-1].vertex.x;
							q->y = delaunay->dcel->vertex[currentPoint-1].vertex.y;
#ifdef DEBUG_SELECT_CLOSEST_POINT
							printf("Closest point %d with distance %lf\n", currentPoint, (*lowest_Distance));
#endif
						}
					}

					// Get next edge before loop.
					edge_Id = delaunay->dcel->edges[edge_Id-1].next_Edge;
					edge_Id = delaunay->dcel->edges[edge_Id-1].twin_Edge;
					currentPoint = delaunay->dcel->edges[edge_Id-1].origin_Vertex;
#ifdef DEBUG_SELECT_CLOSEST_POINT
					printf("Next edge %d has point %d\n", edge_Id, currentPoint);
#endif
				}
			}
			else
			{
#ifdef DEBUG_SELECT_CLOSEST_POINT
				printf("Point is in an internal triangle.\n");
#endif
				// Get closest point between the three points of the triangle.
				(*lowest_Distance) = DBL_MAX;
				for (i=0; i<N_POINTS ;i++)
				{
					new_Distance = distance( p, &delaunay->dcel->vertex[delaunay->graph.nodes[node_Index].points_Index[i]-1].vertex);
#ifdef DEBUG_SELECT_CLOSEST_POINT
					printf("Index point %d. Distance %lf\n", delaunay->graph.nodes[node_Index].points_Index[i], new_Distance);
#endif
					if (new_Distance < (*lowest_Distance))
					{
						q->x = delaunay->dcel->vertex[delaunay->graph.nodes[node_Index].points_Index[i]-1].vertex.x;
						q->y = delaunay->dcel->vertex[delaunay->graph.nodes[node_Index].points_Index[i]-1].vertex.y;
						(*lowest_Distance) = new_Distance;
					}
				}
			}
    	}
    	else
    	{
			// Search triangle in children nodes.
			i=0;
			found = FALSE;
			nChildren = get_nChildren_Node( &delaunay->graph, node_Index);

#ifdef DEBUG_SELECT_CLOSEST_POINT
			printf("Start search in node %d\n", node_Index);
#endif
			while ((!found) && (i < nChildren))
			{
				// Get i-children.
				child_Index = get_iChildren_Node( &delaunay->graph, node_Index, i);

#ifdef DEBUG_SELECT_CLOSEST_POINT
				printf("Checking %d-child in node %d. Node %d\n", i, child_Index, node_Index);
#endif
				// Check if point is interior to i-child node.
				if (is_Interior_To_Node( delaunay->dcel, &delaunay->graph, p, child_Index))
				{
					// Search in next children node.
					node_Index = child_Index;

					// End loop.
					found = TRUE;
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
						// End loop.
						found = FALSE;
#ifdef DEBUG_SELECT_CLOSEST_POINT
						print_Point( p);
						printf("ERROR: No nodes surround new point.\n");
						finished = TRUE;
						finished = FALSE;
#endif
					}
				}
			}
    	}
    }

    return(found);
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
		index = face->edge-1;
		edge = get_Edge( delaunay->dcel, index);

		if (index >= 0)
		{
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
    if (!(return_Turn( dcel, p, id1, id2) == RIGHT) &&
        !(return_Turn( dcel, p, id2, id3) == RIGHT) &&
        !(return_Turn( dcel, p, id3, id1) == RIGHT))
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

//#define DEBUG_RETURN_TURN

enum Turn_T return_Turn( struct DCEL_T *dcel, struct Point_T *p, int source_ID, int dest_ID)
{
    enum Turn_T turn=LEFT_TURN;         // Return value.

    // Normal source point.
    if (source_ID > 0)
    {
        // Normal destination point.
        if (dest_ID > 0)
        {
            // If turn right then point is not in triangle.
			turn = check_Turn( &dcel->vertex[source_ID-1].vertex, &dcel->vertex[dest_ID-1].vertex, p);
        }
        // Destination point is P-2.
        else if (dest_ID == P_MINUS_2)
        {
        	if (source_ID == 1)
        	{
        		turn = LEFT_TURN;
        	}
            // Check if point is over line from source_Index point to P-2.
        	else if (higher_Point( p, &dcel->vertex[source_ID-1].vertex, &lexicographic_Higher))
            {
                turn = RIGHT;
            }
			else
			{
				turn = LEFT_TURN;
			}
        }
        // Destination point is P-1.
        else
        {
            // Check if point is over line from source_Index point to P-1.
            if (higher_Point( p, &dcel->vertex[source_ID-1].vertex, &lexicographic_Higher))
            {
				turn = LEFT_TURN;
            }
			else
			{
				turn = RIGHT;
			}
        }
    }
    else
    {
        // Source point is P-1 and destination cannot be p-2.
        if (source_ID == P_MINUS_1)
        {
        	if (dest_ID == 1)
        	{
        		turn = LEFT_TURN;
        	}
            // Check if point is over line from P-1 point to dest_Index point.
        	else if (higher_Point( p, &dcel->vertex[dest_ID-1].vertex, &lexicographic_Higher))
            {
				turn = RIGHT;
            }
			else
			{
				turn = LEFT_TURN;
			}
        }
        // Source point is P-2.
        else
        {
            // Check destination point.
            if (dest_ID != P_MINUS_1)
            {
				if (higher_Point( p, &dcel->vertex[dest_ID-1].vertex, &lexicographic_Higher))
				{
					turn = LEFT_TURN;
				}
				else
				{
					turn = RIGHT;
				}
            }
			else
			{
				// Points can only do a left turn.
				turn = LEFT_TURN;
			}
        }
    }

#ifdef DEBUG_RETURN_TURN
    printf("Turn between segment %d %d and point \n", source_ID, dest_ID);
    print_Point( p);
    if (turn == LEFT_TURN)
    {
    	printf(" is LEFT_TURN\n");
    }
    else if (turn == RIGHT)
    {
    	printf(" is RIGHT\n");
    }
    else
    {
    	printf(" is COLLINEAR\n");
    }
#endif

    return(turn);
}


void analyze_Face( struct Delaunay_T *delaunay, int point_Index)
{
    int     i=0;                            // Loop counter.
    int     node_Index=0;                   // Index of current node analyzed.
    int		found=FALSE;					// Loop control flag.
	int		finished=FALSE;					// Loop control flag.
	int		child_Index=0;					// Children node ID.
    int     nChildren=0;                    // Number of children of current node.
    struct  Dcel_Vertex_T   *point=NULL;   // Pointer to points in DCEL.

    // Get new point to insert.
    point = get_Vertex( delaunay->dcel, point_Index);

    // Initialize loop control flag.
    finished = FALSE;
    while (!finished)
    {
        // If node is a leaf then triangle found.
        if (is_Leaf_Node( &delaunay->graph, node_Index))
        {
            // Check if point is strictly interior (not over an edge).
            if (is_Strictly_Interior_To_Node(delaunay->dcel, &delaunay->graph, &point->vertex, node_Index))
            {
                // Split current triangle creating 3 triangles.
                split_Triangle( delaunay->dcel, &delaunay->graph, point_Index, node_Index, 3);
            }
            // Point over an edge.
            else
            {
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
            nChildren = get_nChildren_Node( &delaunay->graph, node_Index);
            while ((!found) && (i < nChildren))
            {
                // Get i-children.
                child_Index = get_iChildren_Node( &delaunay->graph, node_Index, i);

                // Check if point is interior to i-child node.
                if (is_Interior_To_Node( delaunay->dcel, &delaunay->graph, &point->vertex, child_Index))
                {
                    // Search in next children node.
                    node_Index = child_Index;

                    // End loop.
                    found = TRUE;
                }
                // Next child.
                else
                {
                    i++;
                    if (i == nChildren)
                    {
                        printf("POINTS ARE EQUAL\n");
                    }
                }
            }
        }
    }
}


void    split_Triangle( struct DCEL_T *dcel, struct Graph_T *graph, int point_Index, int node_Index, int nTriangles)
{
    int     new_Edge_ID=0;							// Edge identifier.
    int     new_Face_ID=0;							// Face identifier.
    int     prev_Edge_ID=0;							// Stores previous edge id.
    int     next_Edge_ID=0;							// Stores next edge id.
    int     colinear_Edge_ID=0, colinear_Index=0;   // Edge identifier of colinear edge.
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
        // Get edge identifier where new point is colinear.
        colinear_Edge_ID = select_Colinear_Edge( dcel, point_Index, face->edge);
		if (colinear_Edge_ID != -1)
		{
			colinear_Index = colinear_Edge_ID - 1;

			// Save previous and next edges ID.
			prev_Edge_ID = dcel->edges[colinear_Index].previous_Edge;
			next_Edge_ID = dcel->edges[colinear_Index].next_Edge;

			// Store edges that must be checked after split of first triangle operation.
			flip_Candidates[0] = next_Edge_ID;
			flip_Candidates[1] = prev_Edge_ID;

			// Store nodes ID that are going to be updated.
			old_Node_ID1 = get_Node_Assigned( graph, dcel->edges[colinear_Index].face);
			old_Node_ID2 = get_Node_Assigned( graph, dcel->edges[dcel->edges[colinear_Index].twin_Edge-1].face);

			// Update current face with new edge: new_Edge_ID.
			insertEdge( dcel, dcel->edges[prev_Edge_ID-1].origin_Vertex, new_Edge_ID+1, next_Edge_ID,
																				colinear_Edge_ID, node->face_ID);

			// Insert a new face with two new edges: new_Edge_ID+1 and new_Edge_ID+2.
			insertEdge( dcel, point_Index+1, new_Edge_ID, new_Edge_ID+2, prev_Edge_ID, new_Face_ID);
			insertEdge( dcel, dcel->edges[colinear_Index].origin_Vertex, new_Edge_ID+3, prev_Edge_ID,
																				new_Edge_ID+1, new_Face_ID);

			// Update existing edges.
			update_Edge( dcel, point_Index+1, NO_UPDATE, new_Edge_ID, NO_UPDATE, NO_UPDATE, colinear_Index);
			update_Edge( dcel, NO_UPDATE, NO_UPDATE, NO_UPDATE, new_Edge_ID, NO_UPDATE, next_Edge_ID-1);
			update_Edge( dcel, NO_UPDATE, NO_UPDATE, new_Edge_ID+1, new_Edge_ID+2, new_Face_ID, prev_Edge_ID-1);

			// Get node of current edge and update it.
			node = get_Node( graph, old_Node_ID1);
			node->children_Index[0] = get_Graph_Length( graph);
			node->children_Index[1] = get_Graph_Length( graph) + 1;
			node->children_Index[2] = INVALID;
			update_Node( graph, old_Node_ID1, 2, node);

			// Insert two new nodes in first node splitted.
			new_Node[0].points_Index[0] = get_Edge_Origin_Vertex( dcel, dcel->edges[colinear_Index].previous_Edge - 1);
			new_Node[0].points_Index[1] = get_Edge_Origin_Vertex( dcel, colinear_Index);
			new_Node[0].points_Index[2] = get_Edge_Origin_Vertex( dcel, dcel->edges[colinear_Index].next_Edge - 1);
			new_Node[0].face_ID = dcel->edges[colinear_Index].face;
			insert_Node( graph, &new_Node[0]);
			new_Node[1].points_Index[0] = get_Edge_Origin_Vertex( dcel, dcel->edges[dcel->edges[colinear_Index].twin_Edge - 1].previous_Edge - 1);
			new_Node[1].points_Index[1] = get_Edge_Origin_Vertex( dcel, dcel->edges[colinear_Index].twin_Edge - 1);
			new_Node[1].points_Index[2] = get_Edge_Origin_Vertex( dcel, dcel->edges[dcel->edges[colinear_Index].twin_Edge - 1].next_Edge - 1);
			new_Node[1].face_ID = new_Face_ID;
			insert_Node( graph, &new_Node[1]);

			// Update twin face.
			colinear_Edge_ID = dcel->edges[colinear_Index].twin_Edge;
			colinear_Index	 = colinear_Edge_ID-1;
			prev_Edge_ID	 = dcel->edges[colinear_Index].previous_Edge;
			next_Edge_ID	 = dcel->edges[colinear_Index].next_Edge;

			// Insert a new face with two new edges: new_Edge_ID+3 and new_Edge_ID+4.
			insertEdge( dcel, point_Index+1, new_Edge_ID+2, new_Edge_ID+4, next_Edge_ID, new_Face_ID+1);
			insertEdge( dcel, dcel->edges[prev_Edge_ID-1].origin_Vertex, new_Edge_ID+5, next_Edge_ID,
																				new_Edge_ID+3, new_Face_ID+1);

			// Update current face with new edge: new_Edge_ID+5.
			insertEdge( dcel, point_Index+1, new_Edge_ID+4, colinear_Edge_ID, prev_Edge_ID, dcel->edges[colinear_Index].face);

			// Update existing edges.
			update_Edge( dcel, NO_UPDATE, NO_UPDATE, new_Edge_ID+3, new_Edge_ID+4, new_Face_ID+1, next_Edge_ID-1);
			update_Edge( dcel, NO_UPDATE, NO_UPDATE, NO_UPDATE, new_Edge_ID+5, NO_UPDATE, colinear_Index);
			update_Edge( dcel, NO_UPDATE, NO_UPDATE, new_Edge_ID+5, NO_UPDATE, NO_UPDATE, prev_Edge_ID-1);

			// Insert two new faces.
			insertFace( dcel, new_Edge_ID + 2);
			insertFace( dcel, new_Edge_ID + 4);

			// Get node of twin edge and update it.
			node = get_Node( graph, old_Node_ID2);
			node->children_Index[0] = get_Graph_Length( graph);
			node->children_Index[1] = get_Graph_Length( graph) + 1;
			node->children_Index[2] = INVALID;
			update_Node( graph, old_Node_ID2, 2, node);

			// Insert two new nodes in first node splitted.
			new_Node[0].points_Index[0] = get_Edge_Origin_Vertex( dcel, dcel->edges[colinear_Index].previous_Edge - 1);
			new_Node[0].points_Index[1] = get_Edge_Origin_Vertex( dcel, colinear_Index);
			new_Node[0].points_Index[2] = get_Edge_Origin_Vertex( dcel, dcel->edges[colinear_Index].next_Edge - 1);
			new_Node[0].face_ID = dcel->edges[colinear_Index].face;
			insert_Node( graph, &new_Node[0]);
			new_Node[1].points_Index[0] = get_Edge_Origin_Vertex( dcel, dcel->edges[next_Edge_ID - 1].previous_Edge - 1);
			new_Node[1].points_Index[1] = get_Edge_Origin_Vertex( dcel, next_Edge_ID - 1);
			new_Node[1].points_Index[2] = get_Edge_Origin_Vertex( dcel, dcel->edges[next_Edge_ID - 1].next_Edge - 1);
			new_Node[1].face_ID = dcel->edges[next_Edge_ID - 1].face;
			insert_Node( graph, &new_Node[1]);

			// Check candidates from first triangle.
			check_Edge( dcel, graph, flip_Candidates[0]);
			check_Edge( dcel, graph, flip_Candidates[1]);

			// Check candidates from second triangle.
			check_Edge( dcel, graph, prev_Edge_ID);
			check_Edge( dcel, graph, next_Edge_ID);
		}
    }
}


double   modified_Signed_Area( struct DCEL_T *dcel, struct  Node_T *node)
{
    double   area=0.0;           // Return value.

    // Chec if any of the vertex is not real: P_MINUS_! or P_MINUS_2.
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
        area = signed_Area( &dcel->vertex[node->points_Index[0]].vertex,
                            &dcel->vertex[node->points_Index[1]].vertex,
                            &dcel->vertex[node->points_Index[2]].vertex);
    }

    return(area);
}


int     select_Colinear_Edge( struct DCEL_T *dcel, int point_Index, int edge_ID)
{
	int		edge_Index=0;					// Index of edge of current face.
    int     colinear_Edge=0;				// Return value.
	int		id1=0, id2=0, id3=0;			// IDs of vertex points.
    struct  Point_T *p=NULL;				// Reference to new point.

	// Set edge index.
	edge_Index = edge_ID - 1;

	// Get index of the vertex of current face.
	id1 = get_Edge_Origin_Vertex( dcel, dcel->edges[edge_Index].previous_Edge-1);
	id2 = get_Edge_Origin_Vertex( dcel, edge_Index);
	id3 = get_Edge_Origin_Vertex( dcel, dcel->edges[edge_Index].next_Edge-1);

    // Get point coordinates.
    p  = get_Vertex_Point( dcel, point_Index);

	// Check if point is colinear to previous edge.
	if (return_Turn( dcel, p, id1, id2) == COLINEAR)
	{
		colinear_Edge = dcel->edges[edge_Index].previous_Edge;
	}
	// Check if point is colinear to input edge.
	else if (return_Turn( dcel, p, id2, id3) == COLINEAR)
	{
		colinear_Edge = edge_ID;
	}
	// Check if point is colinear to next edge.
	else if (return_Turn( dcel, p, id3, id1) == COLINEAR)
	{
		colinear_Edge = dcel->edges[edge_Index].next_Edge;
	}
	else
	{
#ifdef DEBUG
		printf("Checking when point is colinear but it is not.\n");
		sprintf( log_Text, "Checking when point is colinear but it is not.\n");
		write_Log( log_Text);
#endif
		colinear_Edge = -1;
	}

	return(colinear_Edge);
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
						if (check_Turn( p, q, common1) == RIGHT)
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
						if (check_Turn( p, q, common1) == RIGHT)
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
		// Flip edges.
		flip_Edges_Dcel( dcel, graph, edge_ID);
	}
}
