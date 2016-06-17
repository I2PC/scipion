#include "dcel.h"
#include "defines.h"
#include <errno.h>
#include <float.h>
#ifdef LOGGING
#include "log.h"
#endif
#include <math.h>
#include "point.h"
#include "sorting.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define	DEFAULT_CONVEX_HULL_LEN			20

/**************************************************************************
* Private functions declaration
**************************************************************************/
void	print_Vertex(struct Dcel_Vertex_T *vertex);
void	print_Edge(struct Dcel_Edge_T *edge);
void	print_Face(struct Dcel_Face_T *face);
int		get_Vertex_Index( struct DCEL_T *dcel, struct Point_T *point);
int		in_Convex_Hull( struct DCEL_T *dcel, struct Dcel_Edge_T *edge);
void	copy_Dcel(struct DCEL_T *in_Dcel, struct DCEL_T *out_Dcel);


//#define DEBUG_INITIALIZE_DCEL
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
	dcel->vertex = (struct Dcel_Vertex_T *) calloc(nPoints, sizeof(struct Dcel_Vertex_T));

	// Initialize edges attributes.
	dcel->nEdges = 0;
	if (nEdges != INVALID)
	{
		dcel->sizeEdges = nEdges;
	}
	// Max number of edges is = 3xnPoints - 6.
	// NOTE: It is 6xnPoints because there are two directions per edge.
	else
	{
		dcel->sizeEdges = 6*(nPoints+2) - 6;
	}

	// Allocate array of edges.
	dcel->edgeChecked = (int*) calloc( dcel->sizeEdges, sizeof(int));
	dcel->edges = (struct Dcel_Edge_T *) malloc(sizeof(struct Dcel_Edge_T)*dcel->sizeEdges);

	// Initialize faces attributes. If # vertices and # edges are known then # faces = 2 - #v + #e
	dcel->nFaces = 0;
	if (nFaces == INVALID)
	{
		dcel->sizeFaces = 2 - dcel->sizeVertex + dcel->sizeEdges;
	}
	else
	{
		dcel->sizeFaces = nFaces;
	}

	// Allocate array of vertex.
	dcel->faces = (struct Dcel_Face_T *) calloc( dcel->sizeFaces, sizeof(struct Dcel_Face_T));

#ifdef DEBUG_INITIALIZE_DCEL
	printf("Allocating %d points, %d edges and %d faces\n", dcel->sizeVertex, dcel->sizeEdges, dcel->sizeFaces);
#endif

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
	    printf("Function initialize_DCEL.\tError allocating memory in initialize_DCEL\n");
	}

	return(ret);
}

//#define DEBUG_READ_DCEL
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

#ifdef DEBUG_READ_DCEL
		printf("# vertex %d\n", nVertex);
#endif
		// Read vertex.
		for (i=0; i<nVertex; i++)
		{
			// Read x and y coordinates and edge id.
#ifdef FLOAT_TYPE
			ret = fscanf( fd, "%f", &vertex.vertex.x);
			ret = fscanf( fd, "%f", &vertex.vertex.y);
#else
			ret = fscanf( fd, "%lf", &vertex.vertex.x);
			ret = fscanf( fd, "%lf", &vertex.vertex.y);
#endif
			ret = fscanf( fd, "%d", &vertex.origin_Edge);
		}

		// Read number of edges.
		ret = fscanf( fd, "%d", &nEdges);

#ifdef DEBUG_READ_DCEL
		printf("# edges %d\n", nEdges);
#endif
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

#ifdef DEBUG_READ_DCEL
		printf("# faces %d\n", nFaces);
#endif

		// Move to first position of the file.
		fclose( fd);
		fd = fopen( fileName, "r");

		initialize_DCEL( dcel, nVertex, nEdges, nFaces);

		// Read number of vertex.
		ret = fscanf( fd, "%d", &nVertex);

		// Read vertex.
		for (i=0; i<nVertex; i++)
		{
			// Read x and y coordinates and edge id.
#ifdef FLOAT_TYPE
			ret = fscanf( fd, "%f", &vertex.vertex.x);
			ret = fscanf( fd, "%f", &vertex.vertex.y);
#else
			ret = fscanf( fd, "%lf", &vertex.vertex.x);
			ret = fscanf( fd, "%lf", &vertex.vertex.y);
#endif
			ret = fscanf( fd, "%d", &vertex.origin_Edge);

#ifdef DEBUG_READ_DCEL
			printf("Vertex read (%f,%f) origin edge %d\n", vertex.vertex.x, vertex.vertex.y, vertex.origin_Edge);
#endif

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

#ifdef DEBUG_READ_DCEL
		printf("Read dcel with:\n%d vertex\n%d edges\n%d faces\n",
												dcel->nVertex,
												dcel->nEdges,
												dcel->nFaces);
#endif

		fclose(fd);
	}

	return(ret);
}

//#define DEBUG_READ_FLAT_FILE
int     read_Points_Flat_File( struct DCEL_T *dcel, const char *fileName)
{
    int     ret=SUCCESS;        // Return value.
    int     number_Points=0;    // Number of points of set.
	FILE	*fd=NULL;           // File descriptor.

    // Open file.
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
			sprintf( log_Text, "Error %d reading number of points from file: %s.", errno, fileName);
			write_Log( log_Text);
#endif
            printf("Error %d reading number of points from file: %s.", errno, fileName);
			ret = FAILURE;
	    }
	    else
	    {
#ifdef DEBUG_READ_FLAT_FILE
	    	printf("# points %d\n", number_Points);
#endif
			// Allocate DCEL structure.
			if (initialize_DCEL( dcel, number_Points, INVALID, INVALID) == FAILURE)
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

//#define DEBUG_WRITE_DCEL
int		write_DCEL( struct DCEL_T *dcel, int type, const char *fileName)
{
	int		i=0;			// Loop counter.
	int		ret=SUCCESS;	// Return value.
	FILE 	*fd=NULL;		// File descriptor.

	// Open file.
	if ((fd = fopen( fileName, "w")) == NULL)
	{
	    // Set default number of points.
#ifdef LOGGING
		sprintf( log_Text, "Error %d opening input file: %s\n", errno, fileName);
		write_Log( log_Text);
#endif
	    printf("Error %d opening output file: %s\n", errno, fileName);
		ret = FAILURE;
	}
	else
	{
		if (type == DCEL_TYPE)
		{
#ifdef DEBUG_WRITE_DCEL
			printf("Writing DCEL data to %s file\n", fileName);
#endif

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
#ifdef DEBUG_WRITE_DCEL
			printf("Writing only points to %s file\n", fileName);
#endif
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

	return(ret);
}

void	print_DCEL( struct DCEL_T *dcel)
{
    int     i=0;            // Loop counter.

	// Print number of vertex.
    printf("# vertices: %d\n", dcel->nVertex);

	// Print vertex.
	for (i=0; i<dcel->nVertex; i++)
	{
        // Print x and y coordinates and edge id.
		printf("%d\t%f ", i+1, dcel->vertex[i].vertex.x);
		printf("%f ", dcel->vertex[i].vertex.y);
		printf("%d\n", dcel->vertex[i].origin_Edge);
    }

	// Print number of edges.
	printf("# edges: %d\n", dcel->nEdges);

	// Print edges.
	for (i=0; i<dcel->nEdges; i++)
	{
		// Print x and y coordinates and edge id.
		printf("%d\t%d ", i+1, dcel->edges[i].origin_Vertex);
		printf("%d ", dcel->edges[i].twin_Edge);
		printf("%d ", dcel->edges[i].previous_Edge);
		printf("%d ", dcel->edges[i].next_Edge);
		printf("%d\n", dcel->edges[i].face);
    }

	// Print number of faces.
	printf("# faces: %d\n", dcel->nFaces);

	// Print faces.
	for (i=0; i<dcel->nFaces; i++)
	{
        // Print x and y coordinates and edge id.
		printf("%d\t%d\t%d\n", i, dcel->faces[i].edge,  dcel->faces[i].imaginaryFace);
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


/***************************************************************************
* Name: resize_DCEL
* IN:		resize_Type		type of resize (vertex, edges or faces)
* OUT:		N/A
* IN/OUT:	dcel			dcel data
* RETURN:	SUCCESS if resize successfully. FAILURE i.o.c.
* Description: 	Resize the data defined by "resize_Type" to double its size.
***************************************************************************/
int 	resize_DCEL( struct DCEL_T *dcel, enum Resize_DCEL_T resize_Type)
{
    int             ret=SUCCESS;        		// Return value.

    // Temporary pointers.
	struct	Dcel_Vertex_T	*refVertex;
	int						*refEdgeChecked;
	struct	Dcel_Edge_T	    *refEdges;
	struct	Dcel_Face_T	    *refFaces;

	// Check resize type.
	switch(resize_Type)
	{
		// Resize vertices.
		case RESIZE_POINTS:
		{
			// Get reference to current data.
			refVertex = dcel->vertex;

			// Allocate new array.
			dcel->sizeVertex = dcel->sizeVertex*2;
			dcel->vertex = (struct Dcel_Vertex_T *) calloc(dcel->sizeVertex, sizeof(struct Dcel_Vertex_T));
			if (dcel->vertex != NULL)
			{
				// Copy vertices array.
				memcpy( dcel->vertex, refVertex, sizeof(struct Dcel_Vertex_T)*dcel->nVertex);

				// Free previous array.
				free(refVertex);
			}
			else
			{
#ifdef LOGGING
				sprintf( log_Text, "Error allocating vertex in resize_DCEL.\n");
				write_Log( log_Text);
#endif
				printf("Error allocating vertex in resize_DCEL.\n");
				ret = FAILURE;
			}
			break;
		}
		// Resize edges.
		case RESIZE_EDGES:
		{
			// Get reference to current data.
			refEdgeChecked = dcel->edgeChecked;
			refEdges = dcel->edges;

			// Allocate new arrays.
			dcel->sizeEdges = dcel->sizeEdges*2;
			dcel->edgeChecked = (int *) calloc(dcel->sizeVertex, sizeof(int));
			dcel->edges = (struct Dcel_Edge_T *) calloc(dcel->sizeVertex, sizeof(struct Dcel_Edge_T));
			if ((dcel->edgeChecked != NULL) &&
				(dcel->edges != NULL))
			{
				// Copy edges arrays.
				memcpy( dcel->edgeChecked, refEdgeChecked, sizeof(int)*dcel->nEdges);
				memcpy( dcel->edges, refEdges, sizeof(struct Dcel_Edge_T)*dcel->nEdges);

				// Free previous arrays.
				free(refEdgeChecked);
				free(refEdges);
			}
			else
			{
#ifdef LOGGING
				sprintf( log_Text, "Error allocating edges in resize_DCEL.\n");
				write_Log( log_Text);
#endif
				printf("Error allocating edges in resize_DCEL.\n");
				ret = FAILURE;
			}
			break;
		}
		// Resize faces.
		default:
		{
			// Get reference to current data.
			refFaces = dcel->faces;

			// Allocate new array.
			dcel->sizeFaces = dcel->sizeFaces*2;
			dcel->faces = (struct Dcel_Face_T *) calloc(dcel->sizeFaces, sizeof(struct Dcel_Face_T));
			if (dcel->faces != NULL)
			{
				// Copy faces array.
				memcpy( dcel->faces, refFaces, sizeof(struct Dcel_Face_T)*dcel->nFaces);

				// Free previous array.
				free(refFaces);
			}
			else
			{
#ifdef LOGGING
				sprintf( log_Text, "Error allocating faces in resize_DCEL.\n");
				write_Log( log_Text);
#endif
				printf("Error allocating faces in resize_DCEL.\n");
				ret = FAILURE;
			}
			break;
		}
	}

	return(ret);
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
	if (dcel != NULL)
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
}


/***************************************************************************
* Name: get_Number_Vertex
* IN:		dcel			DCEL data
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	# vertex in DCEL
* Description: 	returns # vertex in DCEL
***************************************************************************/
int 	get_Number_Vertex( struct DCEL_T *dcel)
{
    // Return number of vertex in DCEL.
    return(dcel->nVertex);
}


//#define DEBUG_INSERTPOINT
/***************************************************************************
* Name: insertPoint
* IN:		point			input point
* OUT:		N/A
* IN/OUT:	dcel			DCEL data
* GLOBAL:	N/A
* RETURN:	SUCCESS if point inserted. FAILURE i.o.c.
* Description: Inserts a new point in the DCEL without setting the edge
* 				that departs from it.
***************************************************************************/
int 	insertPoint( struct DCEL_T *dcel, struct Point_T *point)
{
	int		ret=SUCCESS;		// Return value.

	// Check if DCEl vertex is full.
	if (dcel->nVertex == dcel->sizeVertex)
	{
		printf("DCEL vertex is full. Size %d and # elements is %d.\n",
														dcel->sizeVertex,
														dcel->nVertex);
		// Resize vertex array.
		if (resize_DCEL( dcel, RESIZE_POINTS) == FAILURE)
		{
#ifdef LOGGING
			sprintf("Error resizing DCEL in insertPoint.\n");
			write_Log( log_Text);
#endif
			printf("Error resizing DCEL in insertPoint.\n");
			ret = FAILURE;
		}
	}

	if (ret == SUCCESS)
	{
#ifdef DEBUG_INSERTPOINT
		printf("Insert point (%lf,%lf) at position %d. Size %d\n", point->x, point->y,
													dcel->nVertex, dcel->sizeVertex);
#endif

		// Update next vertex.
		dcel->vertex[dcel->nVertex].vertex.x = point->x;
		dcel->vertex[dcel->nVertex].vertex.y = point->y;
		dcel->vertex[dcel->nVertex].origin_Edge = INVALID;

		// Update number of vertex.
		dcel->nVertex++;
	}

	return(ret);
}

//#define DEBUG_INSERTVERTEX
/***************************************************************************
* Name: insertVertex
* IN:		vertex			input vertex
* OUT:		N/A
* IN/OUT:	dcel			DCEL data
* GLOBAL:	N/A
* RETURN:	SUCCESS if point inserted. FAILURE i.oc.
* Description: Inserts a new vertex in the DCEL
***************************************************************************/
int 	insertVertex( struct DCEL_T *dcel, struct Dcel_Vertex_T vertex)
{
	int		ret=SUCCESS;		// Return value.

	// Check if DCEl vertex is full.
	if (dcel->nVertex == dcel->sizeVertex)
	{
		printf("DCEL vertex is full. Size %d and # elements is %d.\n", dcel->sizeVertex,
																dcel->nVertex);
		// Resize vertex array.
		if (resize_DCEL( dcel, RESIZE_POINTS) == FAILURE)
		{
#ifdef LOGGING
			sprintf("Error resizing DCEL in insertPoint.\n");
			write_Log( log_Text);
#endif
			printf("Error resizing DCEL in insertPoint.\n");
			ret = FAILURE;
		}
	}

	if (ret == SUCCESS)
	{
#ifdef DEBUG_INSERTVERTEX
		printf("Insert point (%lf,%lf) at position %d. Size %d\n", vertex.vertex.x,
									vertex.vertex.y, dcel->nVertex, dcel->sizeVertex);
#endif

		// Update next vertex.
		dcel->vertex[dcel->nVertex].vertex.x = vertex.vertex.x;
		dcel->vertex[dcel->nVertex].vertex.y = vertex.vertex.y;
		dcel->vertex[dcel->nVertex].origin_Edge = vertex.origin_Edge;

		// Update number of vertex.
		dcel->nVertex++;
	}

	return(ret);
}


/***************************************************************************
* Name: insert_Vertex_At
* IN:		vertex		vertex data to insert
* 			index		position to insert vertex
* OUT:		N/A
* IN/OUT:	dcel		DCEL data
* RETURN:	SUCCESS if vertex inserted. FAILURE i.o.c.
* Description: 	Inserts "vertex" at "index" position.
***************************************************************************/
int insert_Vertex_At( struct DCEL_T *dcel, struct Dcel_Vertex_T vertex, int index)
{
	int		ret=SUCCESS;		// Return value.

	// Check if index is out of bounds.
	if (index >= dcel->nVertex)
	{
		printf("Error inserting vertex in DCEL. Index %d out of bounds (0,%d).\n",
															index,
															dcel->nVertex-1);
#ifdef LOGGING
		sprintf("Error inserting vertex in DCEL. Index %d out of bounds (0,%d).\n",
				index,
				dcel->nVertex-1);
		write_Log( log_Text);
#endif
		ret = FAILURE;
	}
	else
	{
		// Update next vertex.
		dcel->vertex[index].vertex.x = vertex.vertex.x;
		dcel->vertex[index].vertex.y = vertex.vertex.y;
		dcel->vertex[index].origin_Edge = vertex.origin_Edge;

		// Update number of vertex.
		dcel->nVertex++;
	}

	return(ret);
}

//#define DEBUG_UPDATE_VERTEX_EDGE_AT
/***************************************************************************
* Name: update_Vertex_Edge_At
* IN:		vertex		vertex data to insert
* 			index		position to insert vertex
* OUT:		N/A
* IN/OUT:	dcel		DCEL data
* RETURN:	SUCCESS if vertex inserted. FAILURE i.o.c.
* Description: 	Inserts "vertex" at "index" position.
***************************************************************************/
int update_Vertex_Edge_At( struct DCEL_T *dcel, int edge_ID, int index)
{
	int		ret=SUCCESS;		// Return value.

	// Check if index is out of bounds.
	if (index >= dcel->nVertex)
	{
		printf("Error updating vertex in DCEL. Index %d out of bounds (0,%d).\n",
															index,
															dcel->nVertex-1);
#ifdef LOGGING
		sprintf("Error updating vertex in DCEL. Index %d out of bounds (0,%d).\n",
				index,
				dcel->nVertex-1);
		write_Log( log_Text);
#endif
		ret = FAILURE;
	}
	else
	{
#ifdef DEBUG_UPDATE_VERTEX_EDGE_AT
		printf("Updating vertex %d with edge %d\n", index, edge_ID);
#endif
		// Update vertex edge.
		dcel->vertex[index].origin_Edge = edge_ID;
	}

	return(ret);
}

/***************************************************************************
* Name: swap_Vertex
* IN:		index1		first vertex position
* 			index2		second vertex position
* OUT:		N/A
* IN/OUT:	dcel		DCEL data
* RETURN:	SUCCESS if vertex inserted. FAILURE i.o.c.
* Description: 	Swap vertices at "index1"  and "index2" positions.
***************************************************************************/
int		swap_Vertex( struct DCEL_T *dcel, int index1, int index2)
{
	int		ret=SUCCESS;			// Return value.
	struct Dcel_Vertex_T vertex;	// Temporary variable.

	// Check if index is out of bounds.
	if ((index1 >= dcel->nVertex) ||
		(index2 >= dcel->nVertex))
	{
		printf("Error swapping vertices in DCEL. Index1 %d or Index2 %d out of bounds (0,%d).\n",
															index1,
															index2,
															dcel->nVertex-1);
#ifdef LOGGING
		sprintf("Error swapping vertices in DCEL. Index1 %d or Index2 out of bounds (0,%d).\n",
				index1,
				index2,
				dcel->nVertex-1);
		write_Log( log_Text);
#endif
		ret = FAILURE;
	}
	else
	{
		// Swap vertices.
		vertex = dcel->vertex[index1];
		dcel->vertex[index1] = dcel->vertex[index2];
		dcel->vertex[index2] = vertex;
	}

	return(ret);
}


/***************************************************************************
* Name: get_Vertex
* IN:		dcel		DCEL data
* 			index		vertex position
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	Vertex at index position if "index" into array bounds.
* 			Vertex at 0 position if "index" out of bounds.
* Description: Returns the vertex at "index" position.
***************************************************************************/
struct Dcel_Vertex_T *get_Vertex( struct DCEL_T *dcel, int index)
{
	// Check if index is out of bounds.
	if (index >= dcel->nVertex)
	{
		printf("Error in get_Vertex DCEL. Index %d out of bounds (0,%d).\n",
															index,
															dcel->nVertex-1);
#ifdef LOGGING
		sprintf("Error in get_Vertex DCEL. Index %d out of bounds (0,%d).\n",
				index,
				dcel->nVertex-1);
		write_Log( log_Text);
#endif
		// Return vertex at 0 position.
		return(&dcel->vertex[0]);
	}
	else
	{
		// Return vertex.
		return(&dcel->vertex[index]);
	}
}


/***************************************************************************
* Name: get_Vertex_Point
* IN:		dcel		DCEL data
* 			index		vertex position
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	point at index position if "index" into array bounds.
* 			point at 0 position if "index" out of bounds.
* Description: Returns the point at "index" position.
***************************************************************************/
struct Point_T  *get_Vertex_Point( struct DCEL_T *dcel, int index)
{
	// Check if index is out of bounds.
	if (index >= dcel->nVertex)
	{
		printf("Error in get_Vertex_Point DCEL. Index %d out of bounds (0,%d).\n",
															index,
															dcel->nVertex-1);
#ifdef LOGGING
		sprintf("Error in get_Vertex_Point DCEL. Index %d out of bounds (0,%d).\n",
				index,
				dcel->nVertex-1);
		write_Log( log_Text);
#endif
		// Return point at 0 position.
		return(&dcel->vertex[0].vertex);
	}
	else
	{
		// Return point.
		return(&dcel->vertex[index].vertex);
	}
}


/***************************************************************************
* Name: get_Number_Edges
* IN:		dcel		DCEL data
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	# edges in DCEL
* Description: Returns # edges in DCEL
***************************************************************************/
int 	get_Number_Edges( struct DCEL_T *dcel)
{
    // Return number of edges.
    return(dcel->nEdges);
}

/***************************************************************************
* Name: get_Number_Real_Edges
* IN:		dcel		DCEL data
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	# readl edges in DCEL
* Description: Returns # real edges in DCEL. Real edges has no imaginary
* 				vertex as origin nor destination.
***************************************************************************/
int 	get_Number_Real_Edges( struct DCEL_T *dcel)
{
	int		i=0;				// Loop counter.
	int		nEdges=0;			// Return value.

	// Check all edges.
	for (i=0; i<get_Number_Edges(dcel) ;i++)
	{
		// Check origin and destination has no imaginary vertex.
		if ((dcel->edges[i].origin_Vertex >= 0) &&
			(dcel->edges[dcel->edges[i].twin_Edge-1].origin_Vertex >= 0))
		{
			nEdges++;
		}
	}

	return(nEdges / 2);
}


//#define DEBUG_INSERTEDGE
/***************************************************************************
* Name: insertEdge
* IN:		origin		origin point of the edge
* 			twin		twin edge
* 			prev		previous edge
* 			next		next edge
* 			face		edge face
* OUT:		N/A
* IN/OUT:	dcel		dcel DATA
* RETURN:	SUCCESS if edge inserted. FAILURE i.o.c.
* Description: Inserts an edge in the next position of the edge array.
***************************************************************************/
int insertEdge( struct DCEL_T *dcel, int origin, int twin, int prev, int next, int face)
{
	int		ret=SUCCESS;		// Return value.

	// Check if DCEl edges is full.
	if (dcel->nEdges == dcel->sizeEdges)
	{
		printf("DCEL edges is full. Size %d and # elements is %d.\n",
														dcel->sizeEdges,
														dcel->nEdges);
		// Resize edges array.
		if (resize_DCEL( dcel, RESIZE_EDGES) == FAILURE)
		{
#ifdef LOGGING
			sprintf("Error resizing DCEL in insertEdge.\n");
			write_Log( log_Text);
#endif
			printf("Error resizing DCEL in insertEdge.\n");
			ret = FAILURE;
		}
	}

	if (ret == SUCCESS)
	{
#ifdef DEBUG_INSERTEDGE
		printf("Added edge %d. Origin %d Twin %d Prev %d Next %d Face %d.\n",
															dcel->nEdges + 1,
															origin,
															twin,
															prev,
															next,
															face);
#endif
		// Update next edge.
		dcel->edges[dcel->nEdges].origin_Vertex = origin;
		dcel->edges[dcel->nEdges].twin_Edge = twin;
		dcel->edges[dcel->nEdges].previous_Edge = prev;
		dcel->edges[dcel->nEdges].next_Edge = next;
		dcel->edges[dcel->nEdges].face = face;

		// Update number of edges.
		dcel->nEdges++;
	}

	return(ret);
}

//#define DEBUG_UPDATE_EDGE
/***************************************************************************
* Name: update_Edge
* IN:		origin		origin point of the edge
* 			twin		twin edge
* 			next		next edge
* 			prev		previous edge
* 			face		edge face
* 			index		edge position
* OUT:		N/A
* IN/OUT:	dcel		dcel DATA
* RETURN:	SUCCESS if edge updated. FAILURE i.o.c.
* Description: Updates the edge at "index" position. Fields with INVALID
* 				value are not updated.
***************************************************************************/
int update_Edge( struct DCEL_T *dcel, int origin, int twin, int prev, int next,
																	int face,
																	int index)
{
	int	ret=SUCCESS;			// Return value.

	// Check index is not out of bounds.
	if (index >= dcel->nEdges)
	{
		printf("Trying to update index edge %d out of bounds %d\n",
																index,
																dcel->nEdges);
#ifdef LOGGING
		sprintf("Trying to update index edge %d out of bounds %d\n",
																index,
																dcel->nEdges);
		write_Log( log_Text);
#endif
		ret = FAILURE;
	}
	else
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

#ifdef DEBUG_UPDATE_EDGE
		printf("Updated edge %d. Origin %d Twin %d Prev %d Next %d Face %d.\n",
															index + 1,
															origin,
															twin,
															prev,
															next,
															face);
#endif
	}

	return(ret);
}


/***************************************************************************
* Name: are_Twins
* IN:		dcel		dcel DATA
* 			edge1		first edge
* 			edge2		secomd edge
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	TRUE if edge are twins. FALSE i.o.c.
* Description: Checks if "edge1" and "edge2" are twin in "dcel".
***************************************************************************/
int	are_Twins( struct DCEL_T *dcel, struct Dcel_Edge_T *edge1,
									struct Dcel_Edge_T *edge2)
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


/***************************************************************************
* Name: is_External_Edge
* IN:		dcel		dcel DATA
* 			index		edge index.
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	TRUE if edge at "index" position is in external face. FALSE i.o.c.
* Description: Checks if edge at "index" position is in external face.
***************************************************************************/
int		is_External_Edge( struct DCEL_T *dcel, int index)
{
	int is_External=FALSE;			// Return value.

	// Check index is not out of bounds.
	if (index >= dcel->nEdges)
	{
		printf("Trying to check an edge %d out of bounds %d in is_External_Edge\n",
														index,
														dcel->nEdges);
#ifdef LOGGING
		sprintf("Trying to check an edge %d out of bounds %d in is_External_Edge\n",
				index,
				dcel->nEdges);
		write_Log( log_Text);
#endif
		is_External = FALSE;
	}
	else
	{
		// Check if this edge or its twin belong face #0.
		if ((dcel->edges[index].face == EXTERNAL_FACE) ||
			(dcel->edges[dcel->edges[index].twin_Edge-1].face == EXTERNAL_FACE))
		{
			// Edge is external.
			is_External = TRUE;
		}
	}

	return(is_External);
}


/***************************************************************************
* Name: get_Edge_Origin_Vertex
* IN:		dcel			dcel DATA
* 			edge_Index		edge index.
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	origin vertex of edge at "edge_Index" position in the DCEL
* Description: 	Returns the origin vertex of edge at "edge_Index" position
* 				in the DCEL
***************************************************************************/
int     get_Edge_Origin_Vertex( struct DCEL_T *dcel, int edge_Index)
{
	// Check index is not out of bounds.
	if (edge_Index >= dcel->nEdges)
	{
		printf("Trying to get edge %d out of bounds %d in get_Edge_Origin_Vertex\n",
														edge_Index,
														dcel->nEdges);
#ifdef LOGGING
		sprintf("Trying to check an edge %d out of bounds %d in is_External_Edge\n",
				index,
				dcel->nEdges);
		write_Log( log_Text);
#endif
		return(dcel->edges[0].origin_Vertex);
	}
	else
	{
		// Return origin of edge at position edge_Index.
		return(dcel->edges[edge_Index].origin_Vertex);
	}
}

//#define DEBUG_GET_EDGE_IN_CONVEX_HULL
int 	get_Edge_In_Convex_Hull( struct DCEL_T *dcel, int firstEdge_Index)
{
	int	found=FALSE;		// Loop control flag.
	int	edge_Index=0;		// Edge index.

	// Get first edge index.
	edge_Index = firstEdge_Index;

	// Check all edges in current face.
	while (!found)
	{
#ifdef DEBUG_GET_EDGE_IN_CONVEX_HULL
		printf("Checking convex hull edge index %d. Points are %d and %d\n", edge_Index,
								dcel->edges[edge_Index].origin_Vertex,
								dcel->edges[dcel->edges[edge_Index].twin_Edge-1].origin_Vertex);
#endif
		// Check origin and destination are positive.
		if ((dcel->edges[edge_Index].origin_Vertex >= 0) &&
			(dcel->edges[dcel->edges[edge_Index].twin_Edge-1].origin_Vertex >= 0))
		{
			found = TRUE;
		}
		// Next edge.
		else
		{
			edge_Index = dcel->edges[edge_Index].next_Edge - 1;
		}
	}

#ifdef DEBUG_GET_EDGE_IN_CONVEX_HULL
	printf("Found edge %d.\n", edge_Index+1);
#endif

	return(edge_Index+1);
}

/***************************************************************************
* Name: get_Edge
* IN:		dcel			dcel DATA
* 			index			edge index.
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	edge at "index" position in the DCEL
* Description: 	Returns the edge at "index" position in the DCEL
***************************************************************************/
struct Dcel_Edge_T *get_Edge( struct DCEL_T *dcel, int index)
{
	// Check index is not out of bounds.
	if (index >= dcel->nEdges)
	{
		printf("Trying to get edge %d out of bounds %d in get_Edge\n",
														index,
														dcel->nEdges);
#ifdef LOGGING
		sprintf("Trying to check an edge %d out of bounds %d in get_Edge\n",
														index,
														dcel->nEdges);
		write_Log( log_Text);
#endif
		return(&dcel->edges[0]);
	}
	else
	{
		// Return edge.
		return(&dcel->edges[index]);
	}
}


/***************************************************************************
* Name: copy_Edge
* IN:		dcel			dcel DATA
* 			index			edge index.
* OUT:		edge			copied edge
* IN/OUT:	N/A
* RETURN:	SUCCESS if edge copied. FAILURE i.o.c.
* Description: 	copies edge at "index" position in "edge"
***************************************************************************/
int copy_Edge( struct DCEL_T *dcel, int index, struct Dcel_Edge_T *edge)
{
	int		ret=SUCCESS;			// Return value.

	// Check index is not out of bounds.
	if (index >= dcel->nEdges)
	{
		printf("Trying to get edge %d out of bounds %d in copy_Edge\n",
														index,
														dcel->nEdges);
#ifdef LOGGING
		sprintf("Trying to check an edge %d out of bounds %d in copy_Edge\n",
														index,
														dcel->nEdges);
		write_Log( log_Text);
#endif
		ret = FAILURE;
	}
	else
	{
		// Copy edge information.
		edge->face			= dcel->edges[index].face;
		edge->next_Edge		= dcel->edges[index].next_Edge;
		edge->origin_Vertex = dcel->edges[index].origin_Vertex;
		edge->previous_Edge = dcel->edges[index].previous_Edge;
		edge->twin_Edge		= dcel->edges[index].twin_Edge;
	}

	return(ret);
}


/***************************************************************************
* Name: get_Number_Faces
* IN:		dcel			dcel DATA
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	# faces in the DCEL
* Description: 	returns # faces in the DCEL
***************************************************************************/
int     get_Number_Faces(struct DCEL_T *dcel)
{
    // Return the number of faces.
    return(dcel->nFaces);
}


//#define DEBUG_INSERTFACE
/***************************************************************************
* Name: insertFace
* IN:	edge_ID			new edge assigned to face
* OUT:		N/A
* IN/OUT:	dcel		dcel data
* RETURN:	N/A
* Description: 	Inserts a new face in the dcel
***************************************************************************/
int insertFace( struct DCEL_T *dcel, int edge_ID)
{
	int		ret=SUCCESS;		// Return value.

	// Check if DCEL faces is full.
	if (dcel->nFaces == dcel->sizeFaces)
	{
		printf("DCEL faces array is full. Size %d # elements %d.\n",
														dcel->sizeFaces,
														dcel->nFaces);

		// Resize faces array.
		if (resize_DCEL( dcel, RESIZE_FACES) == FAILURE)
		{
#ifdef LOGGING
			sprintf("Error resizing DCEL in insertFace.\n");
			write_Log( log_Text);
#endif
			printf("Error resizing DCEL in insertFace.\n");
			ret = FAILURE;
		}
	}

	if (ret == SUCCESS)
	{
#ifdef DEBUG_INSERTFACE
		printf("Insert face %d with edge %d\n", dcel->nFaces, edge_ID);
#endif
		// Update face edge.
		dcel->faces[dcel->nFaces].edge = edge_ID;

		// Set face as real.
		dcel->faces[dcel->nFaces].imaginaryFace = VALID;

		// Update number of faces.
		dcel->nFaces++;
	}

	return(ret);
}


//#define DEBUG_UPDATE_FACE
/***************************************************************************
* Name: update_Face
* IN:	dcel			dcel data
* 		edge_ID			new edge assigned to face
* 		index			face index
* OUT:		N/A
* RETURN:	N/A
* Description: 	Updates the edge assigned to "index" face.
***************************************************************************/
int update_Face(struct DCEL_T *dcel, int edge_ID, int index)
{
	int		ret=SUCCESS;		// Return value.

	// Check if DCEL faces is full.
	if (dcel->nFaces == dcel->sizeFaces)
	{
		printf("DCEL faces array is full. Size %d # elements %d.\n",
														dcel->sizeFaces,
														dcel->nFaces);

		// Resize faces array.
		if (resize_DCEL( dcel, RESIZE_FACES) == FAILURE)
		{
#ifdef LOGGING
			sprintf("Error resizing DCEL in update_Face.\n");
			write_Log( log_Text);
#endif
			printf("Error resizing DCEL in update_Face.\n");
			ret = FAILURE;
		}
	}

	if (ret == SUCCESS)
	{
#ifdef DEBUG_UPDATE_FACE
		printf("Updating face %d with edge %d\n", index, edge_ID);
#endif
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

	return(ret);
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
	// Check if DCEL faces is full.
	if (index >= dcel->nFaces)
	{
		printf("Trying to get face %d out of bounds %d in get_Face\n",
														index,
														dcel->nFaces);
#ifdef LOGGING
		sprintf("Trying to get an face %d out of bounds %d in get_Face\n",
														index,
														dcel->nFaces);
		write_Log( log_Text);
#endif
		return(&dcel->faces[0]);
	}
	else
	{
		// Return face.
		return(&dcel->faces[index]);
	}
}

//#define DEBUG_GET_FACE_VERTEX
/***************************************************************************
* Name: get_Face_Vertex
* IN:	dcel			DCEL data
* 		face_ID			face index
* OUT:	v1				first vertex
* 		v2				third vertex
* 		v3				second vertex
* RETURN: SUCCESS if face into face array bounds. FAILURE i.o.c.
* Description: 	returns the three vertices of the face at "face_ID" position.
***************************************************************************/
int   get_Face_Vertex( struct DCEL_T *dcel, int face_ID, struct Dcel_Vertex_T *v1,
														struct Dcel_Vertex_T *v2,
														struct Dcel_Vertex_T *v3)
{
	int		ret=SUCCESS;		// Return value.

	// Check if DCEL faces is full.
	if (face_ID >= dcel->nFaces)
	{
		printf("Trying to get face %d out of bounds %d in get_Face\n",
														face_ID,
														dcel->nFaces);
#ifdef LOGGING
		sprintf("Trying to get an face %d out of bounds %d in get_Face\n",
														face_ID,
														dcel->nFaces);
		write_Log( log_Text);
#endif
		ret = FAILURE;
	}
	else
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

	return(ret);
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
	struct Dcel_Face_T  *face=NULL;	// Current face.
	int		nFaces=0;				// Return value.

	// Loop all faces.
	for (i=1; i<dcel->nFaces ;i++)
	{
		// Get i-face.
		face = get_Face( dcel, i);

		// Check if it is a real face.
		if (face->imaginaryFace == VALID)
		{
			nFaces++;
		}
	}

	return(nFaces);
}



/***************************************************************************
* Name: get_Vertex_Of_Face
* IN:	dcel			DCEL data
* 		face			face id
* OUT:		index1		first vertex id
* 			index2		second vertex id
* 			index3		third vertex id
* IN/OUT:	N/A
* RETURN:	N/A
* Description: 	returns the three vertex id of the "face" face.
***************************************************************************/
void get_Vertex_Of_Face( struct DCEL_T *dcel, int face, int *index1,
														int *index2,
														int *index3)
{
	int edgeIndex=0;			// Edge index.

	// Get index in face.
	edgeIndex = dcel->faces[face].edge-1;

	// Get vertices from face.
	(*index1) = dcel->edges[edgeIndex].origin_Vertex;
	(*index2) = dcel->edges[dcel->edges[edgeIndex].twin_Edge-1].origin_Vertex;
	(*index3) = dcel->edges[dcel->edges[edgeIndex].previous_Edge-1].origin_Vertex;
}


/***************************************************************************
* Name: is_Interior_To_Face
* IN:	dcel			triangulation DCEL
* 		p				input point
* 		face			face id to check
* OUT:	N/A
* RETURN:	True		TRUE if point is interior to face
* Description: check if input point p is interior to "face" triangle.
***************************************************************************/
bool     is_Interior_To_Face( struct DCEL_T *dcel, struct Point_T *p, int face_ID)
{
    int     edge_Index=0;					// Edge index.
    bool    is_Interior=false;              // Return value.
	struct Dcel_Vertex_T  p1, p2, p3;		// Temporary points.

	// Check face is not out of bounds.
	if (face_ID >= dcel->nFaces)
	{
		printf("Checking face %d out of bounds %d in is_Interior_To_Face\n",
														face_ID,
														dcel->nFaces);
#ifdef LOGGING
		sprintf("Checking face %d out of bounds %d in is_Interior_To_Face\n",
														face_ID,
														dcel->nFaces);
		write_Log( log_Text);
#endif
		is_Interior = false;
	}
	else
	{
		// Get edge index of edge departing from face_ID.
		edge_Index = dcel->faces[face_ID].edge-1;

		// Get points of current face.
		p1 = dcel->vertex[dcel->edges[dcel->edges[edge_Index].previous_Edge-1].origin_Vertex-1];
		p2 = dcel->vertex[dcel->edges[edge_Index].origin_Vertex-1];
		p3 = dcel->vertex[dcel->edges[dcel->edges[edge_Index].next_Edge-1].origin_Vertex-1];

		// Compute if point is interior.
		is_Interior = interior_Triangle( &p1.vertex, &p2.vertex, &p3.vertex, p);
	}

    return(is_Interior);
}



//#define DEBUG_IS_NEGATIVE_ANY_VERTEX
/***************************************************************************
* Name: is_Negative_Any_Vertex
* IN:	dcel			DCEL data
* 		edgeID			edge identifier
* OUT:	N/A
* RETURN: TRUE if any vertex in edge "edgeID" is negative. FALSE i.o.c.
* Description: 	checks if any vertex in edge "edgeID" is negative
***************************************************************************/
int     is_Negative_Any_Vertex( struct DCEL_T *dcel, int edgeID)
{
    int     is_Negative=FALSE;      // Return value.

	// Check face is not out of bounds.
	if (edgeID > dcel->nEdges)
	{
		printf("Checking edge %d out of bounds %d in is_Negative_Any_Vertex\n",
														edgeID,
														dcel->nEdges);
#ifdef LOGGING
		sprintf("Checking edge %d out of bounds %d in is_Negative_Any_Vertex\n",
														edgeID,
														dcel->nEdges);
		write_Log( log_Text);
#endif
		is_Negative = FALSE;
	}
	else
	{
#ifdef DEBUG_IS_NEGATIVE_ANY_VERTEX
		printf("Origin %d Destination %d\n", dcel->edges[edgeID-1].origin_Vertex,
											dcel->edges[dcel->edges[edgeID-1].next_Edge-1].origin_Vertex);
#endif
		// Check if any of the vertex of the triangle is negative.
		if ((dcel->edges[edgeID-1].origin_Vertex < 0) ||
			(dcel->edges[dcel->edges[edgeID-1].next_Edge-1].origin_Vertex < 0))
		{
			is_Negative = TRUE;
		}
	}

    return(is_Negative);
}

/***************************************************************************
* Name: printFace
* IN:	dcel			DCEL data
* 		faceID			face identifier
* OUT:	N/A
* RETURN: N/A
* Description: 	prints face information.
***************************************************************************/
void	printFace( struct DCEL_T *dcel, int faceID)
{
	// Print face id.
	printf("Face %d.\n", faceID);

	// Print face edges.
	printf("\tEdges:\t%d\t\t%d\t\t%d.\n", dcel->faces[faceID].edge,
						dcel->edges[dcel->faces[faceID].edge-1].next_Edge,
						dcel->edges[dcel->faces[faceID].edge-1].previous_Edge);

	// Print face vertices.
	printf("\tVertices:\t%d\t\t%d\t\t%d.\n", dcel->edges[dcel->faces[faceID].edge-1].origin_Vertex,
						dcel->edges[dcel->edges[dcel->faces[faceID].edge-1].next_Edge-1].origin_Vertex,
						dcel->edges[dcel->edges[dcel->faces[faceID].edge-1].previous_Edge-1].origin_Vertex);
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
                turn = RIGHT_TURN;
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
				turn = RIGHT_TURN;
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
				turn = RIGHT_TURN;
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
					turn = RIGHT_TURN;
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

void  get_Extreme_Point( struct DCEL_T *dcel, int (*f)(struct Point_T *, struct Point_T *),
																		struct Point_T *p)
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

//#define DEBUG_GET_CONVEX_HULL
/***************************************************************************
* Name: get_Convex_Hull
* IN:		dcel			dcel data
* OUT:		points			vector where convex hull points are returned
* IN/OUT:	length			contains vector length as input and contains
* 							convex hull length as output.
* RETURN:	true if convex hull built. false i.o.c.
* Description: returns in "points" the set of points of the convex hull of
* 				the "dcel" and returns in "length" the number of points.
* 				Returns true if the convex hull was stored in the "points"
* 				vector and true if the vector had not enough space.
***************************************************************************/
bool get_Convex_Hull( struct DCEL_T *dcel, int *length, int *points)
{
	int		nPoints=0;			// # points in convex hull.
	bool	finished=false;		// Loop control flag.
	bool	built=false;		// Return value.
	int		firstIndex=0;		// First index in convex hull.
	int		edgeIndex=0;		// Current index.

	// Check if vector has been allocated.
	if (((*length) > 0) && (dcel->nFaces > 0))
	{
#ifdef DEBUG_GET_CONVEX_HULL
		printf("Initial point %d\n", dcel->edges[0].origin_Vertex);
#endif
		// Insert initial point.
		points[nPoints] = dcel->edges[0].origin_Vertex;
		nPoints++;

		// Get edge departing from 0 point.
		edgeIndex = dcel->vertex[0].origin_Edge - 1;

		// Get edge departing from 0 to MINUS_2.
		finished = false;
		while (!finished)
		{
			if ((dcel->edges[edgeIndex].origin_Vertex == 1) &&
				(dcel->edges[dcel->edges[edgeIndex].twin_Edge-1].origin_Vertex == P_MINUS_2))
			{
				finished = true;
#ifdef DEBUG_GET_CONVEX_HULL
				printf("Found first edge %d. O %d D %d\n", edgeIndex,
						dcel->edges[edgeIndex].origin_Vertex,
						dcel->edges[dcel->edges[edgeIndex].twin_Edge-1].origin_Vertex);
#endif
			}
			else
			{
				// Get next edge departing from 0 point.
				edgeIndex = dcel->edges[dcel->edges[edgeIndex].previous_Edge-1].twin_Edge - 1;
			}
		}

		// Get next edge as it is the first edge in the convex hull.
		edgeIndex = dcel->edges[edgeIndex].previous_Edge-1;
		firstIndex = edgeIndex;

		// Get all convex hull points.
		finished = false;
		while ((!built) && (!finished))
		{
			// Insert next point.
			if (nPoints < (*length))
			{
				points[nPoints] = dcel->edges[edgeIndex].origin_Vertex;
				nPoints++;
			}
#ifdef DEBUG_GET_CONVEX_HULL
			printf("Next point %d\n", dcel->edges[edgeIndex].origin_Vertex);
#endif

			// Get next edge.
			edgeIndex = dcel->edges[edgeIndex].previous_Edge-1;
			edgeIndex = dcel->edges[edgeIndex].twin_Edge-1;
			edgeIndex = dcel->edges[edgeIndex].previous_Edge-1;

			// If point is imaginary then skip edge.
			if (dcel->edges[edgeIndex].origin_Vertex < 0)
			{
				edgeIndex = dcel->edges[edgeIndex].twin_Edge-1;
				edgeIndex = dcel->edges[edgeIndex].previous_Edge-1;
			}

			// Check if first edge reached.
			if (edgeIndex == firstIndex)
			{
				// End loop and update vector length.
				if (nPoints < (*length))
				{
					built = true;
				}
				finished = true;
				(*length) = nPoints - 1;
#ifdef DEBUG_GET_CONVEX_HULL
				printf("Finished\n");
#endif
			}
		}
	}
	else
	{
		built = false;
	}

	return(built);
}

//#define DEBUG_INTERIOR_TO_CONVEX_HULL
/***************************************************************************
* Name: is_Interior_To_Convex_Hull
* IN:		dcel			dcel data
* 			p				point to be checked
* OUT:		error			error while executing function.
* IN/OUT:	N/A
* RETURN:	true if point is interior to convex hull. false i.o.c.
* Description: checks if input point "p" is interior to triangulation in
* 				"dcel" DCEL. To be interior the turn between all edges
* 				and the point must not be a RIGHT turn.
***************************************************************************/
bool is_Interior_To_Convex_Hull( struct DCEL_T *dcel, struct Point_T *p, bool *error)
{
	bool	isInterior=false;					// Return value.
	int		*convex=NULL;						// Set of points in convex hull.
	int		length=DEFAULT_CONVEX_HULL_LEN;		// Convex hull length.
	bool	built=false;						// Convex hull built flag.
	int		index=0;							// Vector index.
	bool	finished=false;						// Loop control flag.

	// Initialize loop.
	built = false;
	(*error) = false;

	// Loop to get convex hull.
	while ((!built) && (!(*error)))
	{
		// Allocate vector.
		if ((convex = (int *)malloc(sizeof(int)*length)) == NULL)
		{
			(*error) = true;
#ifdef LOGGING
			sprintf( log_Text, "Error allocating memory for convex hull vector");
			write_Log( log_Text, "interior_To_Convex_Hull", true);
#else
            printf("Function interior_To_Convex_Hull:\n");
            printf("Error allocating memory for convex hull vector\n");
#endif
		}
		// Check error getting convex hull.
		else if (!get_Convex_Hull( dcel, &length, convex))
		{
			// Deallocate vector.
			free(convex);
		}
		else
		{
			// Convex hull successfully built.
			built = true;
		}
	}

	// Check error while allocating vector.
	if (!(*error))
	{
		// Initialize loop variables.
		index = 0;
		isInterior = true;
		finished = false;

		// Loop to check all edges in convex hull.
		while ((isInterior) && (!finished))
		{
			// Check if all points checked.
			if (index == (length-1))
			{
				finished = true;
			}
			// Next edge.
			else
			{
				//  Check if it is a right turn.
				if (return_Turn( dcel, p, convex[index], convex[index+1]) == RIGHT_TURN)
				{
					// It is not interior.
					isInterior = false;
				}
				else
				{
					// Next point in convex hull.
					index++;
				}
			}
		}

		// Deallocate vector.
		free(convex);
	}

	return(isInterior);
}

//#define DEBUG_READ_POINTS_DCEL
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
#ifdef DEBUG_READ_POINTS_DCEL
    	printf("Reading point %d\n", i+1);
#endif
		// Read point.
		if (fscanf( fd, "%f", &point.vertex.x) != 1)
		{
            ret = FAILURE;
#ifdef LOGGING
			sprintf( log_Text, "Fail reading X coordinate of %dth point", i+1);
			write_Log( log_Text);
#endif
            printf("Fail reading X coordinate of %dth point", i+1);
		}
		else if (fscanf( fd, "%f", &point.vertex.y) != 1)
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
#ifdef DEBUG_READ_POINTS_DCEL
			printf("Point read (%lf,%lf)\n", point.vertex.x, point.vertex.y);
#endif
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


//#define DEBUG_GENERATE_RANDOM
/***************************************************************************
* Name: 	generate_Random_Points_DCEL
* IN:		nPoints					# points to generate.
* 			maxX					max X coordinate
* 			maxY					max Y coordinate
* OUT:		N/A
* IN/OUT:	dcel					structure to store points set.
* RETURN:	N/A
* Description: generate a set of "nPoints" points and store it in "dcel".
***************************************************************************/
void	generate_Random_Points_DCEL( int nPoints, struct DCEL_T *dcel, TYPE maxX, TYPE maxY)
{
	int		i=0;						// Loop counter.
	struct	Dcel_Vertex_T	point;		// Temporary point.
	struct	Dcel_Vertex_T	bottom_Most;// Bottom most point.

	// Initialize bottom most.
	bottom_Most.vertex.x = maxX;
	bottom_Most.vertex.y = maxY;

	// Generate random set of points.
	// Create seed.
	srand48((int) time(NULL));

	// Add elements.
	for (i=0; i<nPoints ;i++)
	{
        // Generate new point.
		point.vertex.x = (drand48() * maxX);
		point.vertex.y = (drand48() * maxY);

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


//#define DEBUG_GENERATE_CLUSTER
/***************************************************************************
* Name: 	generate_Cluster_Points_DCEL
* IN:		nPoints					# points to generate
* 			nClusters				# clusters to generate
* 			radius					cluster radius.
* 			maxX					max X coordinate
* 			maxY					max Y coordinate
* OUT:		N/A
* IN/OUT:	dcel					structure to store points set.
* RETURN:	N/A
* Description: 	generate a set of "nPoints" points grouped in "nClusters"
* 				clusters.
***************************************************************************/
void	generate_Cluster_Points_DCEL( int nPoints, struct DCEL_T *dcel,
							int nClusters, int radius, TYPE maxX, TYPE maxY)
{
	int		i=0, j=0;				// Loop counters.
	int		nElementsCluster=0;		// # points per cluster.
	struct	Dcel_Vertex_T	seed;	// Cluster seed point.
	struct	Dcel_Vertex_T	point;	// Temporary point.

#ifdef DEBUG_GENERATE_CLUSTER
	printf("Generating %d points in %d clusters with radius %d\n", nPoints, nClusters, radius);
#endif

	// Create seed.
	srand48((int) time(NULL));

	// Get # elements per cluster.
	nElementsCluster = nPoints / nClusters;

	// Generate clusters centers.
	for (i=0; i<nClusters ;i++)
	{
		// Generate new point.
		seed.vertex.x = (drand48() * maxX);
		seed.vertex.y = (drand48() * maxY);

        // Insert new point.
		seed.origin_Edge = -1;
        insertVertex( dcel, seed);

#ifdef DEBUG_GENERATE_CLUSTER
        printf("Generating %d cluster with seed (%f, %f)\n", i+1, seed.vertex.x, seed.vertex.y);
#endif

        // Add points center in current seed.
        for (j=0; j<nElementsCluster-1 ;j++)
        {
    		// Generate new point.
        	if (drand48() > 0.5)
        	{
        		point.vertex.x = seed.vertex.x + (drand48() * radius);
        	}
        	else
        	{
        		point.vertex.x = seed.vertex.x - (drand48() * radius);
        	}
        	if (drand48() > 0.5)
        	{
        		point.vertex.y = seed.vertex.y + (drand48() * radius);
        	}
        	else
        	{
        		point.vertex.y = seed.vertex.y - (drand48() * radius);
        	}
    		while (distance( &point.vertex, &seed.vertex) > (float) radius)
    		{
        		// Generate new point.
            	if (drand48() > 0.5)
            	{
            		point.vertex.x = seed.vertex.x + (drand48() * radius);
            	}
            	else
            	{
            		point.vertex.x = seed.vertex.x - (drand48() * radius);
            	}
            	if (drand48() > 0.5)
            	{
            		point.vertex.y = seed.vertex.y + (drand48() * radius);
            	}
            	else
            	{
            		point.vertex.y = seed.vertex.y - (drand48() * radius);
            	}
    		}

    		// Correct points out of bounds.
    		if (point.vertex.x < 0.0)
    		{
    			point.vertex.x = -point.vertex.x;
    		}
    		else if (point.vertex.x > maxX)
    		{
    			point.vertex.x = point.vertex.x - radius;
    		}
    		if (point.vertex.y < 0.0)
    		{
    			point.vertex.y = -point.vertex.y;
    		}
    		else if (point.vertex.y > maxY)
    		{
    			point.vertex.y = point.vertex.y - radius;
    		}

            // Insert new point.
            point.origin_Edge = -1;
            insertVertex( dcel, point);

#ifdef DEBUG_GENERATE_CLUSTER
            printf("Generated %d point (%f, %f)\n", j+1, point.vertex.x, point.vertex.y);
#endif
        }
	}

	// Add elements randomly until nPoints reached.
	for (i=dcel->nVertex; i<nPoints ;i++)
	{
		// Generate new point.
		point.vertex.x = (drand48() * maxX);
		point.vertex.y = (drand48() * maxY);

        // Insert new point.
        point.origin_Edge = -1;
        insertVertex( dcel, point);
	}

	// Clutter set of points.
	clutter( dcel);
}

//#define DEBUG_PRINT_DCEL_STATISTICS
/***************************************************************************
* Name: 	print_Dcel_Statistics
* IN:		fileName				output file.
* 			dcel					DCEL structure data.
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	N/A
* Description: writes DCEL statistics data to an output file.
***************************************************************************/
void 	print_Dcel_Statistics( char *fileName, struct DCEL_T *dcel)
{
	FILE 	*fd=NULL;				// File descriptor.
	int		length=0;				// Convex hull length.
	int		*convex_Hull_Set=NULL;	// Set of points in the convex hull.
	bool	error=false;
	long long int		memory=0;	// Amount of memory use.

	if (dcel->nFaces > 0)
	{
#ifdef DEBUG_PRINT_DCEL_STATISTICS
		printf("DCEL has %d faces\n", dcel->nFaces);
#endif
		// Open input file.
		if ((fd = fopen( fileName, "w")) == NULL)
		{
			printf("Error opening output file: %s\n", fileName);
		}
		else
		{
			// Print # of vertices, edges and faces (also allocated #).
			fprintf( fd, "# vertex %d. Allocated %d\n", dcel->nVertex, dcel->sizeVertex);
			fprintf( fd, "# Real edges %d. Total edges %d. Allocated %d\n", get_Number_Real_Edges( dcel), dcel->nEdges, dcel->sizeEdges);
			fprintf( fd, "# Real faces %d. Total faces %d. Allocated %d\n", get_Number_Real_Faces( dcel), dcel->nFaces, dcel->sizeFaces);

			// Check if Delaunay was built using incremental algorithm.
			if (dcel->incremental)
			{
				// Allocate vector.
				length = dcel->nEdges / CONVEX_HULL_LEN_FACTOR;
				if ((convex_Hull_Set = (int *) calloc( length, sizeof(int))) != NULL)
				{
					// Get convex hull.
					error = false;
					while ((!get_Convex_Hull( dcel, &length, convex_Hull_Set)) && (!error))
					{
						// If error then resize vector.
						free(convex_Hull_Set);
						length = length*2;
						if ((convex_Hull_Set = (int *) calloc( length, sizeof(int))) == NULL)
						{
							error = true;
						}
					}
				}

				fprintf( fd, "Convex length %d\n", length);
				free(convex_Hull_Set);
			}

			// Print memory use.
			memory = sizeof(int) + sizeof(int) + sizeof(struct Dcel_Vertex_T)*dcel->sizeVertex;
			memory += sizeof(int) + sizeof(int) + sizeof(struct Dcel_Edge_T)*dcel->sizeEdges + sizeof(int)*dcel->sizeEdges;
			memory += sizeof(int) + sizeof(int) + sizeof(struct Dcel_Face_T)*dcel->sizeFaces;
			fprintf( fd,"Memory required:\n\t\t%f KB\n\t\t%f MB\n", memory / SIZE_OF_KB, memory / (SIZE_OF_KB*SIZE_OF_KB));

			fclose(fd);
		}
	}
	else
	{
		printf("DCEL has no triangles\n");
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


int	set_Edge_Not_Checked( struct DCEL_T *dcel, int index, int *n)
{
	int		ret=SUCCESS;		// Return value.

	// Check if index out of bounds.
	if (index >= dcel->nEdges)
	{
#ifdef LOGGING
		sprintf("Trying to access index %d out of bounds (0,%d)\n.", index,
				dcel->nEdges-1);
		write_Log( log_Text);
#endif
		printf("Trying to access index %d out of bounds (0,%d)\n.", index,
															dcel->nEdges-1);
		ret = FAILURE;
	}
	else
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

	return(ret);
}


