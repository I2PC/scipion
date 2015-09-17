#include "defines.h"
#include "delaunay.h"
#include <math.h>
#include "point.h"
#include "sorting.h"
#include "triangulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/**************************************************************************
* Public functions bodies
**************************************************************************/
int		initialize_Triangulation(struct Triangulation_T *triang, struct DCEL_T *dcel)
{
	int		ret=SUCCESS;			// Return value.

	// Initialize attributes.
	triang->dcel = dcel;

	// Allocate stack.
	ret = initialize_Stack( &triang->stack, triang->dcel->nVertex);

	return(ret);
}

void	finalize_Triangulation(struct Triangulation_T *triang)
{
	// Dereference dcel.
	triang->dcel = NULL;

	// Deallocate stack.
	finalize_Stack( &triang->stack);
}

void 	build_Triangulation( struct Triangulation_T *triang)
{
	int		i=0;							// Loop counter.
	int		len=0;							// Loop length.
	int		finished=0;						// Loop control flag.
#ifdef DEBUG_DELAUNAY
	int		expected_Edges=0;				// Number of edges in triangulation.
	int		expected_Triangles=0;			// Number of triangles in triangulation.
#endif

	int		edgeId=0, faceId=0;				// Edge and face Id counters.
	int		savedEdge=0;
	int		nextConvexEdge=0;				// Next edge to add to convex hull.
	int		originVertex1=0, originVertex2=0;
	int		originVertexIndex1=0, originVertexIndex2=0;
	int		convex_Peak[2];					// Last two points in convex hull.

	struct Dcel_Vertex_T *newVertex=NULL;	// New vertex to add to DCEL.
	struct Convex_T		convex_Edge;		//
	struct Convex_T		peak_Edge;			// Last edge in convex hull.

	// Order in angular order from bottom most point.
	sort( triang->dcel);

	// Insert invalid face. Updated when convex hull added.
	insertFace( triang->dcel, -1);

	// First vertex belongs to convex hull.
	convex_Edge.vertex_Index = 0;
	convex_Edge.edge_ID = 1;
	push( &triang->stack, convex_Edge);

	// Insert first edge.
	triang->dcel->vertex[0].origin_Edge = 1;
	insertEdge( triang->dcel, 1, -1, 3, 2, 1);

	// Second vertex belongs to convex hull.
	convex_Edge.vertex_Index = 1;
	convex_Edge.edge_ID = 2;
	push( &triang->stack, convex_Edge);

	// Insert second edge.
	triang->dcel->vertex[1].origin_Edge = 2;
	insertEdge( triang->dcel, 2, -1, 1, 3, 1);

	// Third element is inserted in convex hull at starting.
	convex_Edge.vertex_Index = 2;
	convex_Edge.edge_ID = 3;
	push( &triang->stack, convex_Edge);

	// Insert third edge.
	triang->dcel->vertex[2].origin_Edge = 3;
	insertEdge( triang->dcel, 3, -1, 2, 1, 1);

	// Insert first face.
	insertFace( triang->dcel, 1);

	// Initianlize variables before loop.
	faceId = 2;
	edgeId = 4;

	// Parse rest of points.
	for (i=3; i<triang->dcel->nVertex ;i++)
	{
		// Get current point.
		newVertex = get_Vertex( triang->dcel, i);

		// Initialize indexes and values.
		originVertexIndex1 = 0;
		originVertex1 = 1;
		originVertexIndex2 = i-1;
		originVertex2 = i;

		finished = 0;
		savedEdge = -1;
		nextConvexEdge = edgeId+2;
		peak_Edge = peak( &triang->stack);
		while (!finished)
		{
			// Get last two vertex from convex hull.
			convex_Edge = peak( &triang->stack);
			convex_Peak[0] = convex_Edge.vertex_Index;
			convex_Edge = peak_Prev( &triang->stack);
			convex_Peak[1] = convex_Edge.vertex_Index;

			// Insert new edge and update its twin.
			triang->dcel->vertex[originVertexIndex1].origin_Edge = edgeId;
			triang->dcel->edges[peak_Edge.edge_ID-1].twin_Edge = edgeId;
			insertEdge( triang->dcel, originVertex1, peak_Edge.edge_ID, edgeId+2, edgeId+1, faceId);
#ifdef DEBUG_DELAUNAY
			printf("A %d LE ASIGNO %d y A %d LE ASIGNO %d\n", edgeId, peak_Edge.edge_ID, peak_Edge.edge_ID, edgeId);
#endif

			// Insert new edge.
			triang->dcel->vertex[originVertexIndex2].origin_Edge = edgeId+1;
			insertEdge( triang->dcel, originVertex2, -1, edgeId, edgeId+2, faceId);
#ifdef DEBUG_DELAUNAY
			printf("A %d LE PONGO -1\n", edgeId+1);
#endif

			// Insert new edge.
			triang->dcel->vertex[i].origin_Edge = edgeId+2;
			insertEdge( triang->dcel, i+1, savedEdge, edgeId+1, edgeId, faceId);
			if (savedEdge != -1)
			{
				triang->dcel->edges[savedEdge-1].twin_Edge = edgeId+2;
#ifdef DEBUG_DELAUNAY
				printf("A %d PONGLE ASIGNO %d y a %d LEASIGNO %d\n", edgeId+2, savedEdge, savedEdge, edgeId+2);
			}
			else
			{
				printf("A %d LE ASIGNO -1\n", edgeId+2);
#endif
			}

			savedEdge = edgeId+1;

			// Check type of turn with new point.
			if ((check_Turn( &triang->dcel->vertex[convex_Peak[1]].vertex,
							&triang->dcel->vertex[convex_Peak[0]].vertex,
							&newVertex->vertex) == LEFT_TURN) ||
				(check_Turn( &triang->dcel->vertex[convex_Peak[1]].vertex,
							&triang->dcel->vertex[convex_Peak[0]].vertex,
							&newVertex->vertex) == COLINEAR))
			{
				// Update peak of convex hull.
				update_Peak( &triang->stack, edgeId+1);
#ifdef DEBUG_DELAUNAY
				printf("Actualizo PEAK con Edge %d\n", edgeId+1);
#endif

				// Add edge to convex hull.
				convex_Edge.vertex_Index = i;
				convex_Edge.edge_ID = nextConvexEdge;
				push( &triang->stack, convex_Edge);
				finished = 1;
#ifdef DEBUG_DELAUNAY
				printf("METO EN PILA PUNTO %d Y EDGE %d\n", convex_Edge.vertex_Index, convex_Edge.edge_ID);
#endif
			}
			else
			{
				// Update source points.
				convex_Edge = peak( &triang->stack);
				originVertexIndex1 = convex_Edge.vertex_Index;
				originVertex1 = convex_Edge.vertex_Index+1;

				// Pop convex hull.
#ifdef DEBUG_DELAUNAY
				printf("SACO DE LA PILA\n");
#endif
				pop( &triang->stack);

				peak_Edge = peak( &triang->stack);
				originVertexIndex2 = peak_Edge.vertex_Index;
				originVertex2 = peak_Edge.vertex_Index+1;
			}

			// Insert face.
			insertFace( triang->dcel, edgeId);

			// Update counters.
			faceId++;
			edgeId = edgeId + 3;
		}
	}

	// Save index of first edge from convex hull.
	convex_Edge = peak( &triang->stack);
	savedEdge = convex_Edge.edge_ID;

	pop( &triang->stack);
	peak_Edge = peak( &triang->stack);

	// Insert first edge from convex hull.
	triang->dcel->edges[peak_Edge.edge_ID-1].twin_Edge = edgeId;
	insertEdge( triang->dcel, convex_Edge.vertex_Index+1, peak_Edge.edge_ID,
								edgeId+get_Stack_nElements(&triang->stack),
								edgeId+1,
								0);

	// Next edge.
	edgeId = edgeId+1;

	// Update convex hull.
	len = get_Stack_nElements( &triang->stack)-1;
	for (i=0; i<len ;i++)
	{
		convex_Edge = peak( &triang->stack);
		pop( &triang->stack);
		peak_Edge = peak( &triang->stack);

		// Insert new edge fom convex hull
		triang->dcel->edges[peak_Edge.edge_ID-1].twin_Edge = edgeId;
		insertEdge( triang->dcel, convex_Edge.vertex_Index+1, peak_Edge.edge_ID, edgeId-1, edgeId+1, 0);

		// Next edge.
		edgeId++;
	}

	// Insert first edge from convx hull.
	triang->dcel->edges[savedEdge-1].twin_Edge = edgeId;
	insertEdge( triang->dcel, 1, savedEdge, edgeId-1, edgeId-len-1, 0);

	// Update convex hull face departing edge.
	triang->dcel->faces[0].edge = edgeId;

#ifdef DEBUG_DELAUNAY
	// Check error in number of triangles.
	expected_Edges = 3*triang->dcel->nVertex - triang->stack.nElements - 3;
	expected_Triangles = 2*triang->dcel->nVertex - triang->stack.nElements - 2;
	printf("Faces %d. Expected %d\nEdges %d. Expected %d\n",
											triang->dcel->nFaces,
											expected_Triangles,
											triang->dcel->nEdges,
											expected_Edges);
#endif
}
