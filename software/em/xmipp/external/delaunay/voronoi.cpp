#include "defines.h"
#include "graph.h"
#include <GL/glut.h>
#ifdef LOGGING
#include "log.h"
#endif
#include "polygon.h"
#include <stdlib.h>
#include <string.h>
#include "voronoi.h"


void insert_Voronoi( struct Voronoi_T *voronoi_Set, struct Point_T *centre);
int	 inner_To_Voronoi_Area( struct Voronoi_T *voronoi, int index,
												struct Point_T *q);


/***************************************************************************
* Name: initialize_Voronoi
* IN:		nPoints			# points in the Voronoi set
* OUT:		N/A
* IN/OUT:	voronoi			voronoi data
* GLOBAL:	N/A
* Description: Allocates voronoi data.
***************************************************************************/
int initialize_Voronoi( struct Voronoi_T *voronoi, int nPoints)
{
	int		ret=SUCCESS;		// Return value.

    // Initialize data.
	voronoi->computed = FALSE;

    // Initialize Voronoi DCEL.
	if (initialize_DCEL( &voronoi->dcel, nPoints, nPoints*8, (nPoints+2)*2) == FAILURE)
	{
#ifdef LOGGING
		sprintf( log_Text, "Error allocating memory in initialize_Voronoi\n");
		write_Log( log_Text);
#endif
		printf("Error allocating memory in initialize_Voronoi\n");
		ret = FAILURE;
    }

	return(ret);
}


/***************************************************************************
* Name: get_Voronoi_Centre
* IN:		voronoi			voronoi data
* 			face			face index
* OUT:		centre			face circumcentre
* IN/OUT:	N/A
* RETURN	FAILURE if face is out of bounds. SUCCESS i.o.c.
* GLOBAL:	N/A
* Description: gets the circumcentre of triangl stored at "face" position
***************************************************************************/
int    get_Voronoi_Centre( struct Voronoi_T *voronoi, int face, struct Point_T *centre)
{
	int		ret=SUCCESS;			// Loop counter.

	// Check if input face index is out of bounds.
    if (face < get_Number_Vertex( &voronoi->dcel))
    {
    	// Update output value.
        (*centre) = voronoi->dcel.vertex[face].vertex;
    }
    else
    {
    	// Print error message.
#ifdef LOGGING
		sprintf( log_Text, "Error in get_Voronoi_Centre. Face %d out of bounds.", face);
		write_Log( log_Text);
#endif
		printf( "Error in get_Voronoi_Centre. Face %d out of bounds.", face);
        ret = FAILURE;
    }

    return(ret);
}


/***************************************************************************
* Name: insert_Voronoi
* IN:		centre			new point
* OUT:		N/A
* IN/OUT:	voronoi			voronoi data
* GLOBAL:	N/A
* Description: Inserts a new point in the Voronoi set.
***************************************************************************/
void    insert_Voronoi( struct Voronoi_T *voronoi, struct Point_T *centre)
{
    // Insert centre.
    insertPoint( &voronoi->dcel, centre);
}


//#define DEBUG_BUILD_VORONOI_AREA
/***************************************************************************
* Name: build_Voronoi_Area
* IN:		dcel			dcel data
* 			pointIndex		index whose voronoi area is computed
* OUT:		N/A
* IN/OUT:	voronoi			voronoi data
* GLOBAL:	N/A
* Description: Computes the Voronoi area of the "pointIndex" point of the
* 				dcel data and store it in the "voronoi" data.
***************************************************************************/
void build_Voronoi_Area( struct Voronoi_T *voronoi, struct DCEL_T *dcel, int pointIndex)
{
	int		edgeID=0;					// Edge id.
	int		originPoint=0;				// origin point in a segment.
	int		firstFace=0;				// First face checked.
	int		currentFace=0;				// Adjacent face id.
	int		firstVoronoi_Edge=0;		// First Voronoi edge added to current Voronoi area.
	int		last_Face_Imaginary=FALSE;	// Check if last inserted point was in an imaginary face.
	int		voronoi_Edge=0;				// Current Voronoi edge to insert.
	int		voronoi_Face=0;				// Current Voronoi face to insert.
	struct  Dcel_Vertex_T *vertex;		// Temporary vertex.
	struct  Dcel_Vertex_T newVertex;	// Temporary vertex.
	struct  Dcel_Face_T *face;			// Temporary face.
	enum Turn_T turn;					// Turn type.

	// Temporary points.
	struct Point_T p1, p2, middle, centre, extreme;

	// Get vertex departing from input point index.
	vertex = get_Vertex( dcel, pointIndex);
	edgeID = vertex->origin_Edge;

	// Get adjacent face data.
	face = get_Face( dcel, dcel->edges[edgeID-1].face);

#ifdef DEBUG_BUILD_VORONOI_AREA
	printf("Starting with edge %d from face %d. Edge origin %d Destination %d\n",
									edgeID,
									dcel->edges[edgeID-1].face,
									dcel->edges[edgeID-1].origin_Vertex,
									dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex);
#endif

	// Skip external face.
	if (dcel->edges[edgeID-1].face != EXTERNAL_FACE)
	{
		// Skip imaginary faces to select as initial edge a real edge.
		while (face->imaginaryFace == INVALID)
		{
#ifdef DEBUG_BUILD_VORONOI_AREA
			printf("Skipping edge %d selected from face %d.\n", edgeID, dcel->edges[edgeID-1].face);
#endif
			// Get edge in adjacent face.
			edgeID = dcel->edges[dcel->edges[edgeID-1].previous_Edge-1].twin_Edge;

			// Get current face.
			face = get_Face( dcel, dcel->edges[edgeID-1].face);
		}

#ifdef DEBUG_BUILD_VORONOI_AREA
		printf("First edge %d selected from face %d.\n", edgeID, dcel->edges[edgeID-1].face);
		printf("Edge origin %d Destination %d\n", dcel->edges[edgeID-1].origin_Vertex,
										dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex);
#endif

		// Save current face as stop condition in current Voronoi area.
		firstFace = dcel->edges[edgeID-1].face;
		currentFace = firstFace;

		// Insert first edge of the Voronoi area.
		originPoint = firstFace-1;
		voronoi_Edge = get_Number_Edges( &voronoi->dcel) + 1;
		voronoi_Face = get_Number_Faces( &voronoi->dcel);
		insertEdge( &voronoi->dcel, originPoint, INVALID, INVALID, voronoi_Edge+1, voronoi_Face);

		// Save current edge to close Voronoi area.
		firstVoronoi_Edge = voronoi_Edge;

		// Save current circumentre.
		centre.x = voronoi->dcel.vertex[originPoint].vertex.x;
		centre.y = voronoi->dcel.vertex[originPoint].vertex.y;

		// Get next adjacent face.
		edgeID = dcel->edges[dcel->edges[edgeID-1].previous_Edge-1].twin_Edge;
		currentFace = dcel->edges[edgeID-1].face;

		// Loop to visit all face neighbors.
		while (currentFace != firstFace)
		{
			// Get current face data.
			face = get_Face( dcel, currentFace);

#ifdef DEBUG_BUILD_VORONOI_AREA
			printf("Analyzing face %d. Edge %d. Edge origin %d Destination %d\n", currentFace, edgeID,
					dcel->edges[edgeID-1].origin_Vertex,
					dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex);
#endif
			// Check if adjacent is NOT imaginary face.
			if ((face->imaginaryFace != INVALID) && (currentFace != 0))
			{
				// Get new circumcentre.
				originPoint = currentFace-1;

				// Save current circumentre.
				centre.x = voronoi->dcel.vertex[originPoint].vertex.x;
				centre.y = voronoi->dcel.vertex[originPoint].vertex.y;

				// Check if previous face was imaginary.
				if (last_Face_Imaginary)
				{
#ifdef DEBUG_BUILD_VORONOI_AREA
					printf("Last face was imaginary.\n");
#endif
					// Reset imaginary face flag.
					last_Face_Imaginary = FALSE;

					// Compute edge middle point.
					middle_Point( &dcel->vertex[dcel->edges[edgeID-1].origin_Vertex-1].vertex,
								  &dcel->vertex[dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex-1].vertex,
								  &middle);

#ifdef DEBUG_BUILD_VORONOI_AREA
					printf("Edge origin is %d. Edge destination %d\n", dcel->edges[edgeID-1].origin_Vertex,
												dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex);
					printf("Middle point (%lf, %lf) between points (%lf,%lf) and (%lf,%lf). Centre (%lf,%lf)\n",
							middle.x, middle.y,
							dcel->vertex[dcel->edges[edgeID-1].origin_Vertex-1].vertex.x,
							dcel->vertex[dcel->edges[edgeID-1].origin_Vertex-1].vertex.y,
							dcel->vertex[dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex-1].vertex.x,
							dcel->vertex[dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex-1].vertex.y,
							centre.x, centre.y);
#endif
					// Check turn between edge in convex hull and previous circumcentre.
					turn = check_Turn( &dcel->vertex[dcel->edges[edgeID-1].origin_Vertex-1].vertex,
									   &dcel->vertex[dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex-1].vertex,
									   &centre);
					if (turn == RIGHT)
					{
#ifdef DEBUG_BUILD_VORONOI_AREA
						printf("Turn RIGHT\n");
#endif
						copy_Point( &middle, &p1);
						copy_Point( &centre, &p2);
					}
					// Check if voronoi centre is into the convex hull.
					else if (turn == LEFT_TURN)
					{
#ifdef DEBUG_BUILD_VORONOI_AREA
						printf("Turn LEFT\n");
#endif
						copy_Point( &centre, &p1);
						copy_Point( &middle, &p2);
					}
					else
					{
#ifdef LOGGING
						sprintf( log_Text, "Collinear points in edge index %d\n", edgeID-1);
						write_Log( log_Text);
						// PENDING
						exit(0);
#endif
					}

					// Compute extreme of voronoi edge.
					extend_Segment( &p1, &p2, &extreme);

					// Add point to Voronoi set.
					originPoint = get_Number_Vertex( &voronoi->dcel);
					newVertex.vertex.x = extreme.x;
					newVertex.vertex.y = extreme.y;
					voronoi_Edge = get_Number_Edges( &voronoi->dcel) + 1;
					newVertex.origin_Edge = voronoi_Edge;
					insertVertex( &voronoi->dcel, newVertex);

					// Insert edge.
					insertEdge( &voronoi->dcel, originPoint, INVALID, voronoi_Edge-1, voronoi_Edge+1, voronoi_Face);
				}

#ifdef DEBUG_BUILD_VORONOI_AREA
				printf("Real face %d with edge %d. Edge origin %d Destination %d\n", currentFace, edgeID,
											dcel->edges[edgeID-1].origin_Vertex,
											dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex);
#endif
				// Restore circumcentre.
				originPoint = currentFace-1;

				// Insert new edge.
				voronoi_Edge = get_Number_Edges( &voronoi->dcel) + 1;
				insertEdge( &voronoi->dcel, originPoint, INVALID, voronoi_Edge-1, voronoi_Edge+1, voronoi_Face);

				// Save current circumentre.
				centre.x = voronoi->dcel.vertex[originPoint].vertex.x;
				centre.y = voronoi->dcel.vertex[originPoint].vertex.y;
			}
			// Current face is an imaginary face.
			else if (currentFace != 0)
			{
				if (!last_Face_Imaginary)
				{
#ifdef DEBUG_BUILD_VORONOI_AREA
					printf("Imaginary face %d found. Search start in edge %d\n", currentFace, edgeID);
#endif
					last_Face_Imaginary = TRUE;

					// Get real edge in current face.
					edgeID = get_Edge_In_Convex_Hull( dcel, face->edge-1);

					// Compute edge middle point.
					middle_Point( &dcel->vertex[dcel->edges[edgeID-1].origin_Vertex-1].vertex,
								  &dcel->vertex[dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex-1].vertex,
								  &middle);

#ifdef DEBUG_BUILD_VORONOI_AREA
					printf("Edge origin is %d. Edge destination %d\n", dcel->edges[edgeID-1].origin_Vertex,
												dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex);
					printf("Middle point (%lf, %lf) between points (%lf,%lf) and (%lf,%lf). Centre (%lf,%lf).\n",
							middle.x, middle.y,
							dcel->vertex[dcel->edges[edgeID-1].origin_Vertex-1].vertex.x,
							dcel->vertex[dcel->edges[edgeID-1].origin_Vertex-1].vertex.y,
							dcel->vertex[dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex-1].vertex.x,
							dcel->vertex[dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex-1].vertex.y,
							centre.x, centre.y);
#endif
					// Check turn between edge in convex hull and previous circumcentre.
					turn = check_Turn( &dcel->vertex[dcel->edges[edgeID-1].origin_Vertex-1].vertex,
									   &dcel->vertex[dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex-1].vertex,
									   &centre);
					if (turn == RIGHT)
					{
#ifdef DEBUG_BUILD_VORONOI_AREA
						printf("Turn RIGHT 2\n");
#endif
						copy_Point( &centre, &p1);
						copy_Point( &middle, &p2);
					}
					// Check if voronoi centre is into the convex hull.
					else if (turn == LEFT_TURN)
					{
#ifdef DEBUG_BUILD_VORONOI_AREA
						printf("Turn LEFT 2\n");
#endif
						copy_Point( &middle, &p1);
						copy_Point( &centre, &p2);
					}
					else
					{
						printf("Collinear points in edge index %d\n", edgeID-1);
#ifdef LOGGING
						sprintf( log_Text, "Collinear points in edge index %d\n", edgeID-1);
						write_Log( log_Text);
#endif
						// PENDING
						exit(0);
					}

					// Compute extreme of voronoi edge.
					extend_Segment( &p1, &p2, &extreme);

					// Add point to Voronoi set.
					originPoint = get_Number_Vertex( &voronoi->dcel);
					newVertex.vertex.x = extreme.x;
					newVertex.vertex.y = extreme.y;
					voronoi_Edge = get_Number_Edges( &voronoi->dcel) + 1;
					newVertex.origin_Edge = voronoi_Edge;
					insertVertex( &voronoi->dcel, newVertex);

					// Insert new edge.
					insertEdge( &voronoi->dcel, originPoint, INVALID, voronoi_Edge-1, voronoi_Edge+1, voronoi_Face);
				}
			}
#ifdef DEBUG_BUILD_VORONOI_AREA
			else
			{
				printf("Skipping 0 face.\n");
			}
#endif
			// Get next adjacent face and one edge in it.
			edgeID = dcel->edges[dcel->edges[edgeID-1].previous_Edge-1].twin_Edge;
			currentFace = dcel->edges[edgeID-1].face;
		}

		if (last_Face_Imaginary)
		{
			// Restore circumcentre.
			originPoint = currentFace-1;

			// Compute edge middle point.
			middle_Point( &dcel->vertex[dcel->edges[edgeID-1].origin_Vertex-1].vertex,
						  &dcel->vertex[dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex-1].vertex,
						  &middle);

			// Save current circumentre.
			centre.x = voronoi->dcel.vertex[originPoint].vertex.x;
			centre.y = voronoi->dcel.vertex[originPoint].vertex.y;

#ifdef DEBUG_BUILD_VORONOI_AREA
			printf("Edge origin is %d. Edge destination %d\n", dcel->edges[edgeID-1].origin_Vertex,
										dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex);
			printf("Middle point (%lf, %lf) between points (%lf,%lf) and (%lf,%lf). Centre (%lf,%lf)\n",
					middle.x, middle.y,
					dcel->vertex[dcel->edges[edgeID-1].origin_Vertex-1].vertex.x,
					dcel->vertex[dcel->edges[edgeID-1].origin_Vertex-1].vertex.y,
					dcel->vertex[dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex-1].vertex.x,
					dcel->vertex[dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex-1].vertex.y,
					centre.x, centre.y);
#endif
			// Check turn between edge in convex hull and previous circumcentre.
			turn = check_Turn( &dcel->vertex[dcel->edges[edgeID-1].origin_Vertex-1].vertex,
							   &dcel->vertex[dcel->edges[dcel->edges[edgeID-1].twin_Edge-1].origin_Vertex-1].vertex,
							   &centre);
			if (turn == RIGHT)
			{
#ifdef DEBUG_BUILD_VORONOI_AREA
				printf("Turn RIGHT\n");
#endif
				copy_Point( &middle, &p1);
				copy_Point( &centre, &p2);
			}
			// Check if voronoi centre is into the convex hull.
			else if (turn == LEFT_TURN)
			{
#ifdef DEBUG_BUILD_VORONOI_AREA
				printf("Turn LEFT\n");
#endif
				copy_Point( &centre, &p1);
				copy_Point( &middle, &p2);
			}
			else
			{
				printf("Collinear points in edge index %d\n", edgeID-1);
#ifdef LOGGING
				sprintf( log_Text, "Collinear points in edge index %d\n", edgeID-1);
				write_Log( log_Text);
#endif
				// PENDING
				exit(0);
			}

			// Compute extreme of voronoi edge.
			extend_Segment( &p1, &p2, &extreme);

			// Add point to Voronoi set.
			originPoint = get_Number_Vertex( &voronoi->dcel);
			newVertex.vertex.x = extreme.x;
			newVertex.vertex.y = extreme.y;
			voronoi_Edge = get_Number_Edges( &voronoi->dcel) + 1;
			newVertex.origin_Edge = voronoi_Edge;
			insertVertex( &voronoi->dcel, newVertex);

			// Insert edge.
			insertEdge( &voronoi->dcel, originPoint, INVALID, voronoi_Edge-1, voronoi_Edge+1, voronoi_Face);
		}

		// Update previous edge in first edge of this Voronoi area.
		update_Edge( &voronoi->dcel, INVALID, INVALID, voronoi_Edge-1, INVALID, INVALID, firstVoronoi_Edge-1);

#ifdef DEBUG_BUILD_VORONOI_AREA
		printf("Update edge %d with Prev %d.\n", firstVoronoi_Edge, voronoi_Edge);
#endif

		// Update next edge in last edge of this Voronoi area.
		voronoi_Edge = get_Number_Edges( &voronoi->dcel);
		update_Edge( &voronoi->dcel, INVALID, INVALID, INVALID, firstVoronoi_Edge, INVALID, voronoi_Edge-1);

#ifdef DEBUG_BUILD_VORONOI_AREA
		printf("Update edge %d with Next %d.\n", voronoi_Edge, firstVoronoi_Edge);
#endif

		// Insert face in voronoi dcel data.
		insertFace( &voronoi->dcel, firstVoronoi_Edge);
#ifdef DEBUG_BUILD_VORONOI_AREA
		printf("Insert face %d with edge %d.\n", get_Number_Faces( &voronoi->dcel), firstVoronoi_Edge);
		printf("------------------------------------------------------------------\n");
		//getchar();
		printf("Edge is %d\n", voronoi->dcel.faces[0].edge);
#endif
	}
}

//#define DEBUG_COMPUTE_VORONOI_CIRCUMCENTRES
/***************************************************************************
* Name: compute_Voronoi_Circumcentres
* IN:		dcel			dcel data
* OUT:		N/A
* IN/OUT:	voronoi			voronoi data
* GLOBAL:	N/A
* Description: For every face of the "dcel" input data computes its
* 				circumcentre.
***************************************************************************/
void compute_Voronoi_Circumcentres( struct Voronoi_T *voronoi, struct DCEL_T *dcel)
{
    int     i=0;                    // Loop counter.
    int		index=0;                // Array index.
    struct Triangle_T triang;	    // Triangle.
	struct Dcel_Face_T *face=NULL;	// Current face.
	struct	Point_T     point;
#ifdef DEBUG_COMPUTE_VORONOI_CIRCUMCENTRES
	int nReal=1, nImaginary=1;		// Statistics values.
#endif

    // Compute centres of triangulation but face 0 (external face).
    for (i=1; i<dcel->nFaces ;i++)
    {
		// Get i-face.
		face = get_Face( dcel, i);

		// Check it is a valid face.
		if (face->imaginaryFace != INVALID)
		{
#ifdef DEBUG_COMPUTE_VORONOI_CIRCUMCENTRES
			printf("Face %d is real.\n", i);
			nReal++;
#endif

			// Get edge index.
        	index = face->edge-1;

            // Build triangle.
            triang.vertex[0] = dcel->vertex[dcel->edges[index].origin_Vertex-1].vertex;
            triang.vertex[1] = dcel->vertex[dcel->edges[dcel->edges[index].next_Edge-1].origin_Vertex-1].vertex;
            triang.vertex[2] = dcel->vertex[dcel->edges[dcel->edges[index].previous_Edge-1].origin_Vertex-1].vertex;

            // Compute circumcentre of each internal face.
            point = get_Centre( &triang);
        }
        else
        {
#ifdef DEBUG_COMPUTE_VORONOI_CIRCUMCENTRES
			printf("Face %d is imaginary.\n", i);
			nImaginary++;
#endif
            // Face not valid. Insert dummy point to keep list of centers consistent with
            // incremental Delaunay structure.
            point.x = INVALID;
            point.y = INVALID;
        }

        // Add point to Voronoi set.
        insert_Voronoi( voronoi, &point);
    }

#ifdef DEBUG_COMPUTE_VORONOI_CIRCUMCENTRES
	printf("Total faces %d. Real %d. Imaginary %d.\n", dcel->nFaces-1, nReal, nImaginary);
	printf("--------------------------------------------------------\n");
#endif
}

//#define DEBUG_BUILD_VORONOI
void    build_Voronoi( struct Voronoi_T *voronoi, struct DCEL_T *dcel)
{
	int		pointIndex=0;		// Loop counter.
	int		nAreas=0;			// # Voronoi areas.

	// Compute Voronoi circumcentres for every non imaginary face.
	compute_Voronoi_Circumcentres( voronoi, dcel);

    // Compute Voronoi area for every point.
	nAreas = get_Number_Vertex( dcel);
    for (pointIndex=0; pointIndex<nAreas ;pointIndex++)
    {
#ifdef DEBUG_BUILD_VORONOI
    	printf("Building %d of %d Voronoi areas.\n", pointIndex, nAreas);
#endif

       	// Build Voronoi area for ith-point in the DCEL.
    	build_Voronoi_Area( voronoi, dcel, pointIndex);
    }

    // Update build flag.
    voronoi->computed = TRUE;

#ifdef DEBUG_BUILD_VORONOI
    printf("Voronoi diagram completed.\n");
#endif
}

#define MAX_NEIGHBORS	50

#define DEBUG_SELECT_CLOSEST_POINT
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
int select_Closest_Point( struct Delaunay_T *delaunay, struct Voronoi_T *voronoi,
							struct Point_T *p, struct Point_T *q, double *lowest_Distance)
{
    int     i=0, j=0;			// Loop counters.
    int     node_Index=0;  		// Index of current node analyzed.
    int		found=FALSE;		// Loop control flag.
    int		finished=FALSE;		// Loop control flag.
	int		child_Index=0;		// Children node ID.
    int     nChildren=0;        // Number of children of current node.
    int		edge_Id=0;			// Edge identifier.

    int		current_Point=0;	// Current point checked.
    int		first_Point=0;		// First point checked.

    // Circular buffer to store point candidates.
    int		in=0;
    int		out=0;
    int		nItems=0;
    int		candidates[MAX_NEIGHBORS];

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
    		printf("Leaf node found: %d\n", node_Index);
#endif
    		// Insert points from leaf node.
    		for (i=0; i<N_POINTS ;i++)
			{
				if (delaunay->graph.nodes[node_Index].points_Index[i] >= 0)
				{
					candidates[in] = delaunay->graph.nodes[node_Index].points_Index[i];
#ifdef DEBUG_SELECT_CLOSEST_POINT
					printf("Inserted point %d.\n", candidates[in]);
#endif
					nItems++;
					in++;
				}
			}

			found = FALSE;
			while (!found)
			{
				if (candidates[out] >= 0)
				{
					if (inner_To_Voronoi_Area( voronoi, candidates[out]-1, p))
					{
						// End inner and outer loops.
						found = TRUE;
						finished = TRUE;

						// Update point and distance.
						q->x = delaunay->dcel->vertex[candidates[out]-1].vertex.x;
						q->y = delaunay->dcel->vertex[candidates[out]-1].vertex.y;
						(*lowest_Distance) = distance( p, q);

#ifdef DEBUG_SELECT_CLOSEST_POINT
						printf("Point %d found.\n", candidates[out]);
			    		//getchar();
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
#ifdef DEBUG_SELECT_CLOSEST_POINT
							if (nItems ==  MAX_NEIGHBORS)
							{
#ifdef LOGGING
								sprintf( log_Text, "Error allocating neighbors array in select_Closest_Point.\n");
								write_Log( log_Text);
#endif
								// PENDING
								printf("Error allocating neighbors array in select_Closest_Point.\n");
								exit(0);
							}
							if (current_Point >= 0)
							{
								printf("Inserted point %d.\n", current_Point);
							}
#endif
							if (current_Point >= 0)
							{
								// Insert point.
								candidates[in] = current_Point;
								nItems++;
								in++;
								in = in % MAX_NEIGHBORS;
							}

							// Get next point.
							edge_Id = delaunay->dcel->edges[edge_Id - 1].next_Edge;
							edge_Id = delaunay->dcel->edges[edge_Id - 1].twin_Edge;
							current_Point = delaunay->dcel->edges[edge_Id - 1].origin_Vertex;
						}
						while (current_Point != first_Point);
					}
				}
#ifdef DEBUG_SELECT_CLOSEST_POINT
				else
				{
					printf("Imaginary point %d.\n", candidates[out]);
				}
#endif
				// Check next point.
				out++;
				nItems--;
				out = out % MAX_NEIGHBORS;
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
						// Point must be equal to an existing point. Retry all points in current node.
						found = FALSE;
						i=0;
						while ((!found) && (i < nChildren))
						{
							// Get i-children.
							child_Index = get_iChildren_Node( &delaunay->graph, node_Index, i);

							// Get closest point between the three points of the triangle.
							for (j=0; j<N_POINTS ;j++)
							{
								if (equal_Point( p, &delaunay->dcel->vertex[delaunay->graph.nodes[child_Index].points_Index[j]-1].vertex))
								{
									(*lowest_Distance) = 0.0;
									q->x = delaunay->dcel->vertex[delaunay->graph.nodes[child_Index].points_Index[j]-1].vertex.x;
									q->y = delaunay->dcel->vertex[delaunay->graph.nodes[child_Index].points_Index[j]-1].vertex.y;
									found = TRUE;
#ifdef DEBUG_SELECT_CLOSEST_POINT
									printf("Index point %d. Distance %lf\n", delaunay->graph.nodes[child_Index].points_Index[i], (*lowest_Distance));
#endif
								}
							}

							i++;
						}

						// End main loop and force inner loop to finish.
						finished = TRUE;
						i = nChildren;

						if (!found)
						{
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

    return(found);
}



#define DEBUG_INNER_TO_VORONOI_AREA
/***************************************************************************
* Name: inner_To_Voronoi_Area
* IN:	voronoi					voronoi areas associated to delaunay triangulation
* 		index					face index to be checked
* 		p						input point
* OUT:		q					closest point to input "p" point
* 			lowest_Distance		distance from input "p" to output "q"
* RETURN:	True				TRUE if closest point found. False i.o.c.
* Description: finds the closest point in the DCEl to an input "p" point.
***************************************************************************/
int	 inner_To_Voronoi_Area( struct Voronoi_T *voronoi, int index,
												struct Point_T *q)
{
	int	inner=TRUE;				// Return value.
	int	first_Edge_Index=0;		// First edge index of the area.
	int	edgeIndex=0;			// Current edge index.
	int	origin=0, dest=0;		// Origin and destination points.

#ifdef DEBUG_INNER_TO_VORONOI_AREA
	printf("Check if point (%lf,%lf) is inner to area %d.\n", q->x, q->y, index);
	//print_DCEL( &voronoi->dcel);
	//getchar();
#endif

	// Get first edge to control loop condition.
	edgeIndex = voronoi->dcel.faces[index].edge - 1;
	first_Edge_Index = edgeIndex;

	// Check all edges in current face.
	do
	{
		// get origin and destination points.
		origin 	= voronoi->dcel.edges[edgeIndex].origin_Vertex;
		dest 	= voronoi->dcel.edges[voronoi->dcel.edges[edgeIndex].next_Edge-1].origin_Vertex;

#ifdef DEBUG_INNER_TO_VORONOI_AREA
		printf("Checking turn between segment points %d %d.\n", origin, dest);
		printf("Coordinates (%lf,%lf) and (%lf,%lf)\n", voronoi->dcel.vertex[origin].vertex.x,
				voronoi->dcel.vertex[origin].vertex.y,
				voronoi->dcel.vertex[dest].vertex.x,
				voronoi->dcel.vertex[dest].vertex.y);
#endif

		// Skip segment between infinite points.
		if (!(has_Extreme_Coordinates( &voronoi->dcel.vertex[origin].vertex) &&
		      has_Extreme_Coordinates( &voronoi->dcel.vertex[dest].vertex)))
		{
			// Check edge turn between current edge and input point.
			if (check_Turn( &voronoi->dcel.vertex[origin].vertex,
							&voronoi->dcel.vertex[dest].vertex,
							q) == RIGHT)
			{
				inner = FALSE;
			}
		}
#ifdef DEBUG_INNER_TO_VORONOI_AREA
		else
		{
			printf("Skipping edge between infinite points.\n");
		}
#endif
		// Next edge.
		edgeIndex = voronoi->dcel.edges[edgeIndex].next_Edge - 1;

	} while ((inner) && (edgeIndex != first_Edge_Index));

#ifdef DEBUG_INNER_TO_VORONOI_AREA
	if (inner)
	{
		printf("Point is interior to Voronoi area.\n");
	}
	else
	{
		printf("RIGHT turn. Point is not interior to Voronoi area.\n");
	}
#endif

	return(inner);
}


/***************************************************************************
* Name: finalize_Voronoi
* IN:		N/A
* OUT:		N/A
* IN/OUT:	voronoi			voronoi data
* GLOBAL:	N/A
* Description: Deallocates voronoi data.
***************************************************************************/
void finalize_Voronoi(struct Voronoi_T *voronoi)
{
    // Reset number of elements.
	voronoi->computed = FALSE;

    // Free DCEL.
    finalize_DCEL( &voronoi->dcel);
}

