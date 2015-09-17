#include "dcel.h"
#include "defines.h"
#include <errno.h>
#include <float.h>
#ifdef LOGGING
#include "log.h"
#endif
#include "point.h"
#include "sorting.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

/**************************************************************************
* Private functions declaration
**************************************************************************/
void	print_Vertex(struct Dcel_Vertex_T *vertex);
void	print_Edge(struct Dcel_Edge_T *edge);
void	print_Face(struct Dcel_Face_T *face);
int		get_Vertex_Index( struct DCEL_T *dcel, struct Point_T *point);
int		in_Convex_Hull( struct DCEL_T *dcel, struct Dcel_Edge_T *edge);
void	copy_Dcel(struct DCEL_T *in_Dcel, struct DCEL_T *out_Dcel);

/**************************************************************************
* Public functions bodies
**************************************************************************/
int 	initialize_DCEL( struct DCEL_T *dcel, int nPoints, int nEdges, int nFaces)
{
    int     ret=SUCCESS;        // Return value.

	// Initialize vertex attributes.
	dcel->nVertex = 0;
	dcel->sizeVertex = nPoints;

	// Allocate array of vertex.
	dcel->vertex = (struct Dcel_Vertex_T *) malloc(sizeof(struct Dcel_Vertex_T)*nPoints);

	// Initialize edges attributes.
	dcel->nEdges = 0;
	dcel->sizeEdges = nEdges;

	// Allocate array of edges.
	dcel->edgeChecked = (int*) calloc( nEdges, sizeof(int));
	dcel->edges = (struct Dcel_Edge_T *) malloc(sizeof(struct Dcel_Edge_T)*nEdges);

	// Initialize faces attributes.
	dcel->nFaces = 0;
	dcel->sizeFaces = nFaces;

	// Allocate array of vertex.
	dcel->faces = (struct Dcel_Face_T *) malloc(sizeof(struct Dcel_Face_T)*nFaces);

	// Check error allocating memory.
	if ((dcel->vertex == NULL) ||
        (dcel->edges == NULL) ||
        (dcel->faces == NULL))
	{
#ifdef LOGGING
		sprintf( log_Text, "Error allocating memory in initialize_DCEL");
		write_Log( log_Text);
#endif
	    ret = FAILURE;
	    printf("Error allocating memory in initialize_DCEL");
	}

	return(ret);
}


int 	read_DCEL( struct DCEL_T *dcel, char *fileName)
{
    int     ret=SUCCESS;    // Return value.
	int		i=0;			// Loop counter.
	FILE 	*fd=NULL;		// File descriptor.
	int		origin=0;
	int		twin=0;
	int		prev=0;
	int		next=0;
	int		edgeFace=0;
	int		nVertex=0, nEdges=0, nFaces=0;
	struct Dcel_Vertex_T	vertex;
	struct Dcel_Face_T	face;

	// Open input file.
	if ((fd = fopen( fileName, "r")) == NULL)
	{
#ifdef LOGGING
		sprintf( log_Text, "Error opening DCEL file: %s\n", fileName);
		write_Log( log_Text);
#endif
		printf("Error opening DCEL file: %s\n", fileName);
		ret = FAILURE;
	}
	else
	{
		// Read number of vertex.
		ret = fscanf( fd, "%d", &nVertex);

		// Read vertex.
		for (i=0; i<nVertex; i++)
		{
			// Read x and y coordinates and edge id.
			ret = fscanf( fd, "%lf", &vertex.vertex.x);
			ret = fscanf( fd, "%lf", &vertex.vertex.y);
			ret = fscanf( fd, "%d", &vertex.origin_Edge);
		}

		// Read number of edges.
		ret = fscanf( fd, "%d", &nEdges);

		// Read edges.
		for (i=0; i<nEdges; i++)
		{
			// Read x and y coordinates and edge id.
			ret = fscanf( fd, "%d", &origin);
			ret = fscanf( fd, "%d", &twin);
			ret = fscanf( fd, "%d", &prev);
			ret = fscanf( fd, "%d", &next);
			ret = fscanf( fd, "%d", &edgeFace);
		}

		// Read number of faces.
		ret = fscanf( fd, "%d", &nFaces);

		// Move to first position of the file.
		fseek( fd, 0, SEEK_SET);

		initialize_DCEL( dcel, nVertex, nEdges, nFaces);

		// Read number of vertex.
		ret = fscanf( fd, "%d", &nVertex);

		// Read vertex.
		for (i=0; i<nVertex; i++)
		{
			// Read x and y coordinates and edge id.
			ret = fscanf( fd, "%lf", &vertex.vertex.x);
			ret = fscanf( fd, "%lf", &vertex.vertex.y);
			ret = fscanf( fd, "%d", &vertex.origin_Edge);

			// Insert vertex.
			insertVertex( dcel, vertex);
		}

		// Read number of edges.
		ret = fscanf( fd, "%d", &nEdges);

		// Read edges.
		for (i=0; i<nEdges; i++)
		{
			// Read x and y coordinates and edge id.
			ret = fscanf( fd, "%d", &origin);
			ret = fscanf( fd, "%d", &twin);
			ret = fscanf( fd, "%d", &prev);
			ret = fscanf( fd, "%d", &next);
			ret = fscanf( fd, "%d", &edgeFace);

			// Insert vertex.
			insertEdge( dcel, origin, twin, prev, next, edgeFace);
		}

		// Read number of faces.
		ret = fscanf( fd, "%d", &nFaces);

		// Read faces.
		for (i=0; i<nFaces; i++)
		{
			// Read x and y coordinates and edge id.
			ret = fscanf( fd, "%d", &face.edge);

			// Insert vertex.
			insertFace( dcel, face.edge);
		}

#ifdef DEBUG
		printf("Read dcel with:\n%d vertex\n%d edges\n%d faces\n",
												dcel->nVertex,
												dcel->nEdges,
												dcel->nFaces);
#endif

		fclose(fd);
	}

	return(ret);
}


int     read_Points_Flat_File( struct DCEL_T *dcel, char *fileName)
{
    int     ret=SUCCESS;        // Return value.
    int     number_Points=0;    // Number of points of set.
	FILE	*fd=NULL;           // File descriptor.

    // Open configuration file.
	if ((fd = fopen( fileName, "r")) == NULL)
	{
	    // Set default number of points.
#ifdef LOGGING
		sprintf( log_Text, "Error %d opening input file: %s\n", errno, fileName);
		write_Log( log_Text);
#endif
	    printf("Error %d opening input file: %s\n", errno, fileName);
		ret = FAILURE;
	}
	else
	{
	    // Read number of points.
	    if (fscanf( fd, "%d", &number_Points) != 1)
	    {
#ifdef LOGGING
			sprintf( log_Text, "Error reading number of points from file: ../data/points_input.txt.\n");
			write_Log( log_Text);
#endif
            printf("Error reading number of points from file: ../data/points_input.txt.");
			ret = FAILURE;
	    }
	    else
	    {
			// Allocate DCEL structure.
			if (initialize_DCEL( dcel, number_Points, number_Points*4*2, (number_Points+2)*2) == FAILURE)
			{
#ifdef LOGGING
				sprintf( log_Text, "Error allocating memory when calling initialize_DCEL\n");
				write_Log( log_Text);
#endif
				printf("Error allocating memory when calling initialize_DCEL");
				ret = FAILURE;
			}
			else
			{
				// Read initial set of points and close input file.
                if (read_Points_DCEL( fd, number_Points, dcel) == FAILURE)
                {
					ret = FAILURE;
                }
			}
	    }

        // Close input file.
        fclose(fd);
	}

	return(ret);
}


void	write_DCEL( struct DCEL_T *dcel, int	type, char *in_FileName)
{
	int		i=0;			// Loop counter.
	int     opened=FALSE;   // Loop control flag.
	int     file_ID=1;     // File id.
	FILE 	*fd=NULL;		// File descriptor.
	char    fileName[1024]; // File name.

	// Initialize loop.
	opened = FALSE;
	while (!opened)
	{
	    // Build file name.
        sprintf( fileName, "%s_%d.txt", in_FileName, file_ID);

        // Check if file exists.
	    if ((fd = fopen( fileName, "r")) == NULL)
	    {
	        opened = TRUE;

	        // Open in writting mode.
	        fd = fopen( fileName, "w");
	    }
	    else
	    {
            // File exists -> close file.
            fclose(fd);
            fd = NULL;

            // Try next file identifier.
            file_ID++;
	    }
	}

	if (type == DCEL_TYPE)
	{
		// Write number of vertex.
		fprintf( fd, "%d\n", dcel->nVertex);

		// Write vertex.
		for (i=0; i<dcel->nVertex; i++)
		{
			// Read x and y coordinates and edge id.
			fprintf( fd, "%f ", dcel->vertex[i].vertex.x);
			fprintf( fd, "%f ", dcel->vertex[i].vertex.y);
			fprintf( fd, "%d\n", dcel->vertex[i].origin_Edge);
		}

		// Write number of edges.
        fprintf( fd, "%d\n", dcel->nEdges);

		// Write edges.
		for (i=0; i<dcel->nEdges; i++)
		{
			// Read x and y coordinates and edge id.
			fprintf( fd, "%d ", dcel->edges[i].origin_Vertex);
			fprintf( fd, "%d ", dcel->edges[i].twin_Edge);
			fprintf( fd, "%d ", dcel->edges[i].previous_Edge);
			fprintf( fd, "%d ", dcel->edges[i].next_Edge);
			fprintf( fd, "%d\n", dcel->edges[i].face);
		}

		// Write number of faces.
		fprintf( fd, "%d\n", dcel->nFaces);

		// Read faces.
		for (i=0; i<dcel->nFaces; i++)
		{
			// Read x and y coordinates and edge id.
			fprintf( fd, "%d\n", dcel->faces[i].edge);
		}
    }
	else
	{
		// Write number of vertex.
		fprintf( fd, "%d ", dcel->nVertex);

		// Write vertex.
		for (i=0; i<dcel->nVertex; i++)
		{
			// Read x and y coordinates and edge id.
			fprintf( fd, "%f ", dcel->vertex[i].vertex.x);
			fprintf( fd, "%f ", dcel->vertex[i].vertex.y);
		}
	}

	// Close file.
    fclose(fd);
}

void	print_DCEL( struct DCEL_T *dcel)
{
    int     i=0;            // Loop counter.

	// Print number of vertex.
    printf("# vertices: %d\n", dcel->nVertex);

	// Write vertex.
	for (i=0; i<dcel->nVertex; i++)
	{
        // Read x and y coordinates and edge id.
		printf("%d\t%f ", i+1, dcel->vertex[i].vertex.x);
		printf("%f ", dcel->vertex[i].vertex.y);
		printf("%d\n", dcel->vertex[i].origin_Edge);
    }

	// Write number of edges.
	printf("# edges: %d\n", dcel->nEdges);

	// Write edges.
	for (i=0; i<dcel->nEdges; i++)
	{
		// Read x and y coordinates and edge id.
		printf("%d\t%d ", i+1, dcel->edges[i].origin_Vertex);
		printf("%d ", dcel->edges[i].twin_Edge);
		printf("%d ", dcel->edges[i].previous_Edge);
		printf("%d ", dcel->edges[i].next_Edge);
		printf("%d\n", dcel->edges[i].face);
    }

	// Write number of faces.
	printf("# faces: %d\n", dcel->nFaces);

	// Read faces.
	for (i=0; i<dcel->nFaces; i++)
	{
        // Read x and y coordinates and edge id.
		printf("%d\t%d\n", i, dcel->faces[i].edge);
    }
}

void	reset_DCEL(struct DCEL_T *dcel)
{
	// Reset vertex.
	dcel->nVertex = 0;
	memset( dcel->vertex, 0, sizeof(struct Dcel_Vertex_T)*dcel->sizeVertex);

	// Reset edges.
	dcel->nEdges = 0;
	memset( dcel->edgeChecked, 0, dcel->sizeEdges*sizeof(int));
	memset( dcel->edges, 0, dcel->sizeEdges*sizeof(struct Dcel_Edge_T));

	// Reset faces.
	dcel->nFaces = 0;
	memset( dcel->faces, 0, sizeof(struct Dcel_Face_T)*dcel->sizeFaces);
}

void	copy_Dcel(struct DCEL_T *in_Dcel, struct DCEL_T *out_dcel)
{
	// Copy vertex.
	out_dcel->nVertex = in_Dcel->nVertex;
	out_dcel->sizeVertex = in_Dcel->sizeVertex;
	memcpy( &out_dcel->vertex, in_Dcel->vertex, sizeof(struct Dcel_Vertex_T)*in_Dcel->nVertex);

	// Copy edges.
	out_dcel->nEdges = in_Dcel->nEdges;
	out_dcel->sizeEdges = in_Dcel->sizeEdges;
	memcpy( &out_dcel->edgeChecked, in_Dcel->edgeChecked, sizeof(int)*in_Dcel->nEdges);
	memcpy( &out_dcel->edges, in_Dcel->edges, sizeof(struct Dcel_Edge_T)*in_Dcel->nEdges);

	// Copy faces.
	out_dcel->nFaces = in_Dcel->nFaces;
	out_dcel->sizeFaces = in_Dcel->sizeFaces;
	memcpy( &out_dcel->faces, in_Dcel->faces, sizeof(struct Dcel_Face_T)*in_Dcel->nFaces);
}

void    check_DCEL_Data(struct DCEL_T *dcel)
{
    struct DCEL_T	check_Dcel;

    // Initialize copy of input DCEL data.
    initialize_DCEL( &check_Dcel, dcel->nVertex, dcel->nEdges, dcel->nFaces);

    // Copy Dcel.
    copy_Dcel( dcel, &check_Dcel);

    // Order in angular order from bottom most point.
	sort( dcel);

    // Deallocate memory.
    finalize_DCEL( &check_Dcel);
}

void	finalize_DCEL(struct DCEL_T *dcel)
{
	// Deallocate vertex memory.
	if (dcel->vertex != NULL)
	{
		free(dcel->vertex);
	}
	dcel->vertex = NULL;

	// Deallocate edges memory.
	if (dcel->edges != NULL)
	{
		free(dcel->edges);
	}
	dcel->edges = NULL;

	// Deallocate edgeChecked memory.
	if (dcel->edgeChecked != NULL)
	{
		free(dcel->edgeChecked);
	}
	dcel->edgeChecked = NULL;

	// Deallocate faces memory.
	if (dcel->faces != NULL)
	{
		free(dcel->faces);
	}
	dcel->faces = NULL;
}






int 	get_Number_Vertex( struct DCEL_T *dcel)
{
    // Return number of vertex in DCEL.
    return(dcel->nVertex);
}



void 	insertVertex( struct DCEL_T *dcel, struct Dcel_Vertex_T vertex)
{
	// Update next vertex.
	dcel->vertex[dcel->nVertex].vertex.x = vertex.vertex.x;
	dcel->vertex[dcel->nVertex].vertex.y = vertex.vertex.y;
	dcel->vertex[dcel->nVertex].origin_Edge = vertex.origin_Edge;

	// Update number of vertex.
	dcel->nVertex++;
}


void 	insert_Vertex_At( struct DCEL_T *dcel, struct Dcel_Vertex_T vertex, int index)
{
	// Update next vertex.
	dcel->vertex[index].vertex.x = vertex.vertex.x;
	dcel->vertex[index].vertex.y = vertex.vertex.y;
	dcel->vertex[index].origin_Edge = vertex.origin_Edge;

	// Update number of vertex.
	dcel->nVertex++;
}

void 	update_Vertex_Edge_At( struct DCEL_T *dcel, int edge_ID, int index)
{
	// Update vertex edge.
	dcel->vertex[index].origin_Edge = edge_ID;
}

void	swap_Vertex( struct DCEL_T *dcel, int index1, int index2)
{
	struct Dcel_Vertex_T vertex;

	vertex = dcel->vertex[index1];
	dcel->vertex[index1] = dcel->vertex[index2];
	dcel->vertex[index2] = vertex;

}

struct Dcel_Vertex_T *get_Vertex( struct DCEL_T *dcel, int index)
{
	// Return vertex.
	return(&dcel->vertex[index]);
}

struct Point_T       *get_Vertex_Point( struct DCEL_T *dcel, int index)
{
	// Return vertex point.
	return(&dcel->vertex[index].vertex);
}

int 	get_Number_Edges( struct DCEL_T *dcel)
{
    // Return number of edges.
    return(dcel->nEdges);
}


void 	insertEdge( struct DCEL_T *dcel, int origin, int twin, int prev, int next, int face)
{
	// Update next edge.
	dcel->edges[dcel->nEdges].origin_Vertex = origin;
	dcel->edges[dcel->nEdges].twin_Edge = twin;
	dcel->edges[dcel->nEdges].previous_Edge = prev;
	dcel->edges[dcel->nEdges].next_Edge = next;
	dcel->edges[dcel->nEdges].face = face;

	// Update number of edges.
	dcel->nEdges++;
}

void    update_Edge( struct DCEL_T *dcel, 	int origin,	int twin, int prev, int next, int face, int index)
{
    // Check if origin point field must be updated.
    if (origin != NO_UPDATE)
    {
        dcel->edges[index].origin_Vertex = origin;
    }

    // Check if twin edge field must be updated.
    if (twin != NO_UPDATE)
    {
        dcel->edges[index].twin_Edge = twin;
    }

    // Check if prev edge field must be updated.
    if (prev != NO_UPDATE)
    {
        dcel->edges[index].previous_Edge = prev;
    }

    // Check if next edge field must be updated.
    if (next != NO_UPDATE)
    {
        dcel->edges[index].next_Edge = next;
    }

    // Check if face field must be updated.
    if (face != NO_UPDATE)
    {
        dcel->edges[index].face = face ;
    }

}


int		are_Twins( struct DCEL_T *dcel, struct Dcel_Edge_T *edge1, struct Dcel_Edge_T *edge2)
{
	int		twins=0;						// Return value.
	struct Dcel_Vertex_T *origin1=NULL;	// Origin of edge 1.
	struct Dcel_Vertex_T *dest1=NULL;	// Destination of edge 1.
	struct Dcel_Vertex_T *origin2=NULL;	// Origin of edge 2.
	struct Dcel_Vertex_T *dest2=NULL;	// Destination of edge 2.

	// Get origin and destination points of edge 1.
	origin1 = get_Vertex( dcel, edge1->origin_Vertex-1);
	dest1 = get_Vertex( dcel, dcel->edges[edge1->next_Edge-1].origin_Vertex-1);

	// Get origin and destination points of edge 2.
	origin2 = get_Vertex( dcel, edge2->origin_Vertex-1);
	dest2 = get_Vertex( dcel, dcel->edges[edge2->next_Edge-1].origin_Vertex-1);

	// Check if both edges share origin and destination.
	twins = 0;
	if (equal_Point( &origin1->vertex, &dest2->vertex) &&
		equal_Point( &origin2->vertex, &dest1->vertex))
	{
		twins = 1;
	}

	return(twins);
}

int		is_External_Edge( struct DCEL_T *dcel, int index)
{
	int is_External=FALSE;			// Return value.

	// Check if this edge or its twin belong face #0.
	if ((dcel->edges[index].face == EXTERNAL_FACE) ||
		(dcel->edges[dcel->edges[index].twin_Edge-1].face == EXTERNAL_FACE))
	{
		// Edge is external.
		is_External = TRUE;
	}

	return(is_External);
}

int     get_Edge_Origin_Vertex( struct DCEL_T *dcel, int edge_Index)
{
    // Return origin of edge at position edge_Index.
    return(dcel->edges[edge_Index].origin_Vertex);
}

struct Dcel_Edge_T *get_Edge( struct DCEL_T *dcel, int index)
{
	// Return edge.
	return(&dcel->edges[index]);
}

void    copy_Edge( struct DCEL_T *dcel, int index, struct Dcel_Edge_T *edge)
{
	// Copy edge information.
	edge->face			= dcel->edges[index].face;
	edge->next_Edge		= dcel->edges[index].next_Edge;
	edge->origin_Vertex = dcel->edges[index].origin_Vertex;
	edge->previous_Edge = dcel->edges[index].previous_Edge;
	edge->twin_Edge		= dcel->edges[index].twin_Edge;
}


int     get_Number_Faces(struct DCEL_T *dcel)
{
    // Return the number of faces.
    return(dcel->nFaces);
}


/***************************************************************************
* Name: insertFace
* IN:	edge_ID			new edge assigned to face
* OUT:		N/A
* IN/OUT:	dcel		dcel data
* RETURN:	N/A
* Description: 	Inserts a new face in the dcel
***************************************************************************/
void 	insertFace( struct DCEL_T *dcel, int edge_ID)
{
	// Update face edge.
	dcel->faces[dcel->nFaces].edge = edge_ID;

	// Update number of faces.
	dcel->nFaces++;
}


/***************************************************************************
* Name: update_Face
* IN:	dcel			dcel data
* 		edge_ID			new edge assigned to face
* 		index			face index
* OUT:		N/A
* RETURN:	N/A
* Description: 	Updates the edge assigned to "index" face.
***************************************************************************/
void	update_Face(struct DCEL_T *dcel, int edge_ID, int index)
{
	// Update edge id of face at index-th position.
	if (edge_ID == INVALID)
	{
		dcel->faces[index].imaginaryFace = INVALID;
	}
	else
	{
		dcel->faces[index].edge = edge_ID;
	}
}


/***************************************************************************
* Name: get_Face
* IN:	dcel			dcel data
* 		index			face index
* OUT:	N/A
* RETURN:	returns the face stored at index position
* Description: 	returns the face stored at index position
***************************************************************************/
struct Dcel_Face_T *get_Face( struct DCEL_T *dcel, int index)
{
	// Return face.
	return(&dcel->faces[index]);
}

//#define DEBUG_GET_FACE_VERTEX
void   get_Face_Vertex( struct DCEL_T *dcel, int face_ID, struct Dcel_Vertex_T *v1, struct Dcel_Vertex_T *v2, struct Dcel_Vertex_T *v3)
{
#ifdef DEBUG_GET_FACE_VERTEX
	printf("Three vertices points are %d %d %d\n", dcel->edges[dcel->faces[face_ID].edge-1].origin_Vertex,
			dcel->edges[dcel->edges[dcel->faces[face_ID].edge-1].next_Edge-1].origin_Vertex,
			dcel->edges[dcel->edges[dcel->faces[face_ID].edge-1].previous_Edge-1].origin_Vertex);
#endif
    // Get first point.
    v1->vertex.x = dcel->vertex[dcel->edges[dcel->faces[face_ID].edge-1].origin_Vertex-1].vertex.x;
    v1->vertex.y = dcel->vertex[dcel->edges[dcel->faces[face_ID].edge-1].origin_Vertex-1].vertex.y;

    // Get second point.
    v2->vertex.x = dcel->vertex[dcel->edges[dcel->edges[dcel->faces[face_ID].edge-1].next_Edge-1].origin_Vertex-1].vertex.x;
    v2->vertex.y = dcel->vertex[dcel->edges[dcel->edges[dcel->faces[face_ID].edge-1].next_Edge-1].origin_Vertex-1].vertex.y;

    // Get third point.
    v3->vertex.x = dcel->vertex[dcel->edges[dcel->edges[dcel->faces[face_ID].edge-1].previous_Edge-1].origin_Vertex-1].vertex.x;
    v3->vertex.y = dcel->vertex[dcel->edges[dcel->edges[dcel->faces[face_ID].edge-1].previous_Edge-1].origin_Vertex-1].vertex.y;
}


/***************************************************************************
* Name: get_Number_Real_Faces
* IN:	dcel			dcel data
* OUT:	pointsSet		array of three Point_T that are the vertex of the triangle
* RETURN:				# real faces in triangulation
* Description: 	returns the number of real faces in a Delaunay triangulation.
* 				Real faces are those triangles that belong to a Delaunay
* 				triangulation that was built incrementally and so it can
* 				contain imaginary faces due to imaginary points taken into
* 				account while building the triangulation.
***************************************************************************/
int     get_Number_Real_Faces(struct DCEL_T *dcel)
{
	int		i=0;					// Loop counter.
	int		index=0;				// Edge index.
	struct Dcel_Face_T  *face=NULL;	// Current face.
	int		nFaces=0;				// Return value.

	// Loop all faces.
	for (i=1; i<dcel->nFaces ;i++)
	{
		// Get i-face.
		face = get_Face( dcel, i);
		index = face->edge-1;

		// Check if it is a real face.
		if (index >= 0)
		{
			nFaces++;
		}
	}

	return(nFaces);
}






int     is_Interior_To_Face( struct DCEL_T *dcel, struct Point_T *p, int face_ID)
{
    int     edge_Index=0;
    int     is_Interior=FALSE;              // Return value.
	struct Dcel_Vertex_T  p1, p2, p3;

    // Get edge index of edge departing from face_ID.
    edge_Index = dcel->faces[face_ID].edge-1;

    // Get points of current face.
    p1 = dcel->vertex[dcel->edges[dcel->edges[edge_Index].previous_Edge-1].origin_Vertex-1];
    p2 = dcel->vertex[dcel->edges[edge_Index].origin_Vertex-1];
    p3 = dcel->vertex[dcel->edges[dcel->edges[edge_Index].next_Edge-1].origin_Vertex-1];

    // Compute if point is interior.
    is_Interior = interior_Triangle( &p1.vertex, &p2.vertex, &p3.vertex, p);

    return(is_Interior);
}

//#define DEBUG_IS_NEGATIVE_ANY_VERTEX

int     is_Negative_Any_Vertex( struct DCEL_T *dcel, int edge_Index)
{
    int     is_Negative=FALSE;      // Return value.

#ifdef DEBUG_IS_NEGATIVE_ANY_VERTEX
    printf("Origin %d Destination %d\n", dcel->edges[edge_Index-1].origin_Vertex,
    					dcel->edges[dcel->edges[edge_Index-1].next_Edge-1].origin_Vertex);
#endif
    // Check if any of the vertex of the triangle is negative.
    if ((dcel->edges[edge_Index-1].origin_Vertex < 0) ||
        (dcel->edges[dcel->edges[edge_Index-1].next_Edge-1].origin_Vertex < 0))
    {
        is_Negative = TRUE;
    }

    return(is_Negative);
}

void    shake_Dcel( struct DCEL_T *dcel)
{
    int     point_Index=0;          // Loop counter.

    // Loop all points.
    for (point_Index=0; point_Index<dcel->nVertex ;point_Index++)
    {
        // Add a delta value to all points in dcel.
        dcel->vertex[point_Index].vertex.x = dcel->vertex[point_Index].vertex.x + 0.1;
        dcel->vertex[point_Index].vertex.y = dcel->vertex[point_Index].vertex.y + 0.1;
    }
}

void  get_Extreme_Point( struct DCEL_T *dcel, int (*f)(struct Point_T *, struct Point_T *), struct Point_T *p)
{
	int	edge_Index=0;			// Edge index.
	int	point_Index=0, first=0;	// Edge index.
	int finished=FALSE;		    // Loop control flag.

	// Initialize loop.
	edge_Index = dcel->faces[EXTERNAL_FACE].edge-1;
	point_Index = dcel->edges[edge_Index].origin_Vertex-1;
	first = point_Index;

	// Get first point.
	p->x = dcel->vertex[point_Index].vertex.x;
	p->y = dcel->vertex[point_Index].vertex.y;

	// Search in convex hull.
	finished = FALSE;
	while (!finished)
	{
		// Get next edge and next origin point.
		edge_Index = dcel->edges[edge_Index].next_Edge-1;
		point_Index = dcel->edges[edge_Index].origin_Vertex-1;

		// Check if all points in convex hull have been checked.
		if (point_Index == first)
		{
			finished = TRUE;
		}
		else
		{
			// Check if new point is better than current.
			if (f( &dcel->vertex[point_Index].vertex, p))
			{
				// Update most point.
				p->x = dcel->vertex[point_Index].vertex.x;
				p->y = dcel->vertex[point_Index].vertex.y;
			}
		}
	}
}

int 	read_Points_DCEL( FILE *fd, int nPoints, struct DCEL_T *dcel)
{
    int     ret=SUCCESS;                // Return value.
	int		i=0;						// Loop counter.
	struct	Dcel_Vertex_T	point;		// Temp point.
	struct	Dcel_Vertex_T	bottom_Most;// Bottom most point.

	// Initialize bottom most.
	bottom_Most.vertex.x = MAX_X_COORD;
	bottom_Most.vertex.y = MAX_Y_COORD;

	// Add elements.
	i=0;
    while ((i<nPoints) && (ret==SUCCESS))
	{
		// Read point.
		if (fscanf( fd, "%lf", &point.vertex.x) != 1)
		{
            ret = FAILURE;
#ifdef LOGGING
			sprintf( log_Text, "Fail reading X coordinate of %dth point", i+1);
			write_Log( log_Text);
#endif
            printf("Fail reading X coordinate of %dth point", i+1);
		}
		else if (fscanf( fd, "%lf", &point.vertex.y) != 1)
		{
            ret = FAILURE;
#ifdef LOGGING
			sprintf( log_Text, "Fail reading Y coordinate of %dth point", i+1);
			write_Log( log_Text);
#endif
            printf("Fail reading Y coordinate of %dth point", i+1);
		}
		else
		{
            // Insert new point.
            point.origin_Edge = -1;
            insertVertex( dcel, point);

            // Update bottom most point.
            if (point.vertex.y < bottom_Most.vertex.y)
            {
                // Save bottom most.
                bottom_Most = point;

                // Swap it with previous vertex that was bottom most.
                swap_Vertex( dcel, 0, dcel->nVertex-1);
            }
            else if (point.vertex.y == bottom_Most.vertex.y)
            {
                if (point.vertex.x < bottom_Most.vertex.x)
                {
                    // Save bottom most.
                    bottom_Most = point;

                    // Swap it with previous vertex that was bottom most.
                    swap_Vertex( dcel, 0, dcel->nVertex-1);
                }
            }

            // Next point.
            i++;
		}
	}

	return(ret);
}

void	generate_Points_DCEL( int nPoints, struct DCEL_T *dcel)
{
	int		i=0;						// Loop counter.
	struct	Dcel_Vertex_T	point;		// Temporary point.
	struct	Dcel_Vertex_T	bottom_Most;// Bottom most point.

	// Initialize bottom most.
	bottom_Most.vertex.x = MAX_X_COORD;
	bottom_Most.vertex.y = MAX_Y_COORD;

	// Generate random set of points.
	// Create seed.
	srand48((int) time(NULL));

	// Add elements.
	for (i=0; i<nPoints ;i++)
	{
        // Generate new point.
		point.vertex.x = (drand48() * (double) MAX_X_COORD) / 10.0;
		point.vertex.y = (drand48() * (double) MAX_Y_COORD) / 10.0;

        // Insert new point.
        point.origin_Edge = -1;
        insertVertex( dcel, point);

        // Update bottom most point.
        if (point.vertex.y < bottom_Most.vertex.y)
        {
            // Save bottom most.
            bottom_Most = point;

            // Swap it with previous vertex that was bottom most.
            swap_Vertex( dcel, 0, dcel->nVertex-1);
        }
        else if (point.vertex.y == bottom_Most.vertex.y)
        {
            if (point.vertex.x < bottom_Most.vertex.x)
            {
                // Save bottom most.
                bottom_Most = point;

                // Swap it with previous vertex that was bottom most.
                swap_Vertex( dcel, 0, dcel->nVertex-1);
            }
        }
    }
}

void	shake_Points_DCEL( struct DCEL_T *dcel)
{
	int		i=0;						// Loop counter.
	int		position=0;					// Random value.

	// Create seed.
	srand((int) time(NULL));

	// Add elements.
	for (i=0; i<dcel->nVertex ;i++)
	{
		// Generate new point.
		position = rand();
		position = position % dcel->nVertex;

		// Swap elements.
		swap_Vertex( dcel, i, position);
	}
}



//*****************************************************************************
//                      PRIVATE FUNCTION BODIES
//*****************************************************************************
void	print_Vertex(struct Dcel_Vertex_T *vertex)
{
	printf("Vertex info:\n");
	print_Point( &vertex->vertex);
	printf("Edge: %d\n\n", vertex->origin_Edge);
}


void	print_Edge(struct Dcel_Edge_T *edge)
{
	printf("Edge info:\n\tOrigin: %d\n\tTwin %d\n\tPrevious %d\n\tNext %d\n\tFace %d\n\n",
												edge->origin_Vertex,
												edge->twin_Edge,
												edge->previous_Edge,
												edge->next_Edge,
												edge->face);
}


void	print_Face(struct Dcel_Face_T *face)
{
	printf("Face info:\n\tEdge: %d\n\n", face->edge);
}


int		get_Vertex_Index( struct DCEL_T *dcel, struct Point_T *point)
{
	int		i=0;			// Loop counter.
	int		found=0;		// Loop control flag.

	// Search point in list.
	i=0;
	found = 0;
	while ((i<dcel->nVertex) && !found)
	{
		// Compare input point and i-point.
		if (equal_Point( point, &dcel->vertex[i].vertex))
		{
			found = 1;
		}
		// Next point.
		else
		{
			i++;
		}
	}

	return(i);
}

int		in_Convex_Hull( struct DCEL_T *dcel, struct Dcel_Edge_T *edge)
{
	int		index=0;
	int		in_Convex=0;						// Return value.
	struct Dcel_Edge_T *prev_Edge=NULL;		// Previous edge.
	struct Dcel_Vertex_T	*v1=NULL;
	struct Dcel_Vertex_T	*v2=NULL;

	// Get three points before.
	index = edge->twin_Edge-1;
	prev_Edge = get_Edge( dcel, index);
	index = prev_Edge->previous_Edge-1;
	prev_Edge = get_Edge( dcel, index);
	index = prev_Edge->previous_Edge-1;
	prev_Edge = get_Edge( dcel, index);

	// Get origin points.
	v1 = get_Vertex( dcel, edge->origin_Vertex-1);
	v2 = get_Vertex( dcel, prev_Edge->origin_Vertex-1);

	// Check if it is the same point -> not in convex hull.
	if (equal_Point( &v1->vertex, &v2->vertex))
	{
		in_Convex = 0;
	}
	else
	{
		in_Convex = 1;
	}

	return(in_Convex);
}


void	set_Edge_Not_Checked( struct DCEL_T *dcel, int index, int	*n)
{
	// If twin not checked set as checked.
	if (dcel->edges[index].face != 0)
	{
		if (dcel->edgeChecked[index])
		{
			(*n)++;
			dcel->edgeChecked[index] = 0;
		}
	}
}


