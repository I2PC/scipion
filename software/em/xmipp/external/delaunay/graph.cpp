#include "defines.h"
#include "graph.h"
#ifdef LOGGING
#include "log.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*****************************************************************************
* Public function bodies
*****************************************************************************/
/*****************************************************************************
* Name: initialize_Graph
* IN:		size			# nodes to allocate
* OUT:		graph			graph data
* IN/OUT:	N/A
* RETURN:	SUCCESS if graph allocated. FAILURE i.o.c
* Description: 	allocates structures to store "size" nodes.
*****************************************************************************/
int  initialize_Graph(struct Graph_T *graph, int  size)
{
    int     ret=SUCCESS;        // Return value.

    // Initialize variables.
    graph->size = size;
    graph->nElements = 0;

    // Allocate memory initializing data.
    graph->nodes = (struct Node_T *) calloc( graph->size, sizeof(struct Node_T));
	graph->face_Node_Assign = (int *) calloc( graph->size, sizeof(int));

	// Check error allocating memory.
	if ((graph->nodes == NULL) ||
		(graph->face_Node_Assign == NULL))
	{
		printf( "Error allocating memory in initialize_Graph\n");
	    ret = FAILURE;
#ifdef LOGGING
		sprintf( log_Text, "Error allocating memory in initialize_Graph\n");
		write_Log( log_Text);
#endif
	}

	return(ret);
}


/***************************************************************************
* Name: get_Graph_Length
* IN:	graph			graph data
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	returns # nodes in the graph
* Description: 	returns # nodes in the graph
***************************************************************************/
int get_Graph_Length(struct Graph_T *graph)
{
    // Return number of elements in graph.
    return(graph->nElements);
}


/***************************************************************************
* Name: get_Graph_Size
* IN:	graph			graph data
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	returns size of the graph
* Description: 	returns size of the graph
***************************************************************************/
int get_Graph_Size(struct Graph_T *graph)
{
    // Return number of elements in graph.
    return(graph->size);
}


/***************************************************************************
* Name: get_nChildren_Node
* IN:	graph			graph data
* 		index			node to check
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	returns # of children in "index" node.
* Description: 	returns # of children in "index" node.
***************************************************************************/
int get_nChildren_Node(struct Graph_T *graph, int index)
{
    // Return number of children of node_at position index.
    return( graph->nodes[index].nChildren);
}


/***************************************************************************
* Name: get_iChildren_Node
* IN:	graph			graph data
* 		node_Index		node to check
* 		child_Index		i child to return
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	returns the "child_Index" child node of the "node_Index" node.
* Description: 	returns the "child_Index" child node of the "node_Index" node.
***************************************************************************/
int get_iChildren_Node(struct Graph_T *graph, int node_Index, int child_Index)
{
    // Return children node identifier at index-th position.
    return( graph->nodes[node_Index].children_Index[child_Index]);
}


/***************************************************************************
* Name: get_Face_Of_Node
* IN:	graph			graph data
* 		index			node index
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	returns the face assigned to node at "index" position
* Description: 	returns the face assigned to node at "index" position
***************************************************************************/
int	get_Face_Of_Node(struct Graph_T *graph, int index)
{
	// Return the face of node at position index.
	return(graph->nodes[index].face_ID);
}


/***************************************************************************
* Name: get_Vertex_Of_Node
* IN:	graph			graph data
* 		node_Index		node index
* OUT:		index1		first vertex id
* 			index2		second vertex id
* 			index3		third vertex id
* IN/OUT:	N/A
* RETURN:	N/A
* Description: 	returns the three vertex id of the "node_Index" node.
***************************************************************************/
void get_Vertex_Of_Node( struct Graph_T *graph, int node_Index, int *index1,
																int *index2,
																int *index3)
{
	// Get point vertexes of node at node_Index position.
	(*index1) = graph->nodes[node_Index].points_Index[0];
	(*index2) = graph->nodes[node_Index].points_Index[1];
	(*index3) = graph->nodes[node_Index].points_Index[2];
}



/***************************************************************************
* Name: is_Graph_Full
* IN:	graph			graph data
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	# nodes in graph
* Description: 	returns the # nodes in graph
***************************************************************************/
int is_Graph_Full(struct Graph_T *graph)
{
    int     is_Full=FALSE;         // Return value.

    // Check if array is full.
    if (graph->nElements == graph->size)
    {
        is_Full = TRUE;
    }

    return(is_Full);
}


//#define DEBUG_INSERT_NODE
/***************************************************************************
* Name: insert_Node
* IN:	node		new node data
* OUT:		N/A
* IN/OUT:	graph			graph data
* RETURN:	N/A
* Description: 	inserts a new node using input "node" data
***************************************************************************/
int insert_Node(struct Graph_T *graph, struct Node_T *node)
{
    int     ret=SUCCESS;

    // Check if array is full.
    if (is_Graph_Full(graph))
    {
        if (resize_Graph( graph) == FAILURE)
        {
#ifdef LOGGING
        	sprintf( log_Text, "Graph is full and error resizing graph.\n");
        	write_Log( log_Text);
#endif
        	printf("Graph is full and error resizing graph.\n");
            ret = FAILURE;
        }
    }

    // Check previous error.
    if (ret == SUCCESS)
    {
        // Set node as a leaf.
        graph->nodes[graph->nElements].nChildren = 0;

        // Update children triangles.
        graph->nodes[graph->nElements].children_Index[0] = INVALID;
        graph->nodes[graph->nElements].children_Index[1] = INVALID;
        graph->nodes[graph->nElements].children_Index[2] = INVALID;

		// Update point indexes fields.
		graph->nodes[graph->nElements].points_Index[0] = node->points_Index[0];
		graph->nodes[graph->nElements].points_Index[1] = node->points_Index[1];
		graph->nodes[graph->nElements].points_Index[2] = node->points_Index[2];

        // Copy face ID.
        graph->nodes[graph->nElements].face_ID = node->face_ID;

		// Set current node as assigned to face_ID.
		graph->face_Node_Assign[node->face_ID] = graph->nElements;

#ifdef DEBUG_INSERT_NODE
        printf("New node %d. Node points are %d,%d,%d. Face %d\n", graph->nElements,
														node->points_Index[0],
														node->points_Index[1],
														node->points_Index[2],
														graph->nodes[graph->nElements].face_ID);
#endif
        // Increase number of elements.
        graph->nElements++;
    }

    return(ret);
}


/***************************************************************************
* Name: get_Node
* IN:	graph			graph data
* 		index			index of the node to update
* OUT:		N/A
* IN/OUT:	N/A
* RETURN:	node at "index" position
* Description: 	returns the node at "index" position
***************************************************************************/
struct Node_T *get_Node(struct Graph_T *graph, int index)
{
    // Return node at position index.
    return(&graph->nodes[index]);
}

//#define DEBUG_UPDATE_NODE
/***************************************************************************
* Name: update_Node
* IN:	index			index of the node to update
* 		nChildren		3 children in new node
* 		node			new node data
* OUT:		N/A
* IN/OUT:	graph			graph data
* RETURN:	N/A
* Description: 	updates the node at "index" position using "nChildren" and
* 				"node" data
***************************************************************************/
void update_Node(struct Graph_T *graph, int index, int  nChildren, struct Node_T *node)
{
	// Check if index is out of bounds
	if ((index > graph->nElements) || (index < 0))
	{
		printf("Index %d out of bounds (0,%d) when updating graph.",
															index,
															graph->nElements-1);
#ifdef LOGGING
		sprintf( log_Text, "Index %d out of bounds (0,%d) when updating graph.",
															index,
															graph->nElements-1);
		write_Log( log_Text);
#endif
	}
	else
	{
		// Update values of node at position index.
		graph->nodes[index].nChildren = nChildren;

		// Update children fields.
		graph->nodes[index].children_Index[0] = node->children_Index[0];
		graph->nodes[index].children_Index[1] = node->children_Index[1];
		graph->nodes[index].children_Index[2] = node->children_Index[2];

		// Update point indexes fields.
		graph->nodes[index].points_Index[0] = node->points_Index[0];
		graph->nodes[index].points_Index[1] = node->points_Index[1];
		graph->nodes[index].points_Index[2] = node->points_Index[2];

		// Update face ID.
		graph->nodes[index].face_ID = node->face_ID;

#ifdef DEBUG_UPDATE_NODE
		printf("Updated node %d. Points are %d, %d and %d. Face %d. # children %d\n", index,
													node->points_Index[0],
													node->points_Index[1],
													node->points_Index[2],
													graph->nodes[index].face_ID,
													nChildren);
#endif
	}
}


/***************************************************************************
* Name: is_Leaf_Node
* IN:	graph			graph data
* 		index			index of the node to check
* OUT:		N/A
* RETURN:	True if any the node is leaf. False i.o.c.
* Description: Checks if any the node is a leaf (has no children).
***************************************************************************/
int is_Leaf_Node(struct Graph_T *graph, int index)
{
    // If no children then the node is a leaf.
    return(graph->nodes[index].nChildren == 0);
}


/***************************************************************************
* Name: get_Node_Assigned
* IN:	graph			graph data
* 		face_ID			face whose node is returned
* OUT:		N/A
* RETURN:	node that stores the "face_ID" face.
* Description: returns the node that stores the "face_ID" face.
***************************************************************************/
int	get_Node_Assigned( struct Graph_T *graph, int face_ID)
{
	// Return the node assigned to face_ID.
	return(graph->face_Node_Assign[face_ID]);
}


/***************************************************************************
* Name: number_Imaginary_Points
* IN:	graph			graph data
* 		nodeIndex		index of the node to check
* OUT:		N/A
* RETURN:	# imaginary points in "nodeIndex" node.
* Description: returns the # imaginary points in "nodeIndex" node.
* 				Imaginary points are used during incremental Delaunay building.
***************************************************************************/
int	number_Imaginary_Points( struct Graph_T *graph, int nodeIndex)
{
	int	i=0;					// Loop counter.
	int	n_Imaginary=0;			// Return value.

	for (i=0; i<N_POINTS ;i++)
	{
		// Check if i point is imaginary.
		if (graph->nodes[nodeIndex].points_Index[i] < 0)
		{
			n_Imaginary++;
		}
	}
	return(n_Imaginary);
}


/***************************************************************************
* Name: is_Any_Point_Imaginary
* IN:	graph			graph data
* 		nodeIndex		index of the node to check
* OUT:		N/A
* RETURN:	True if any of the vertex of the node is imaginary
* Description: Checks if any of the three points is imaginary (negative).
* 				Imaginary points are used during incremental Delaunay building.
***************************************************************************/
int	is_Any_Point_Imaginary( struct Graph_T *graph, int nodeIndex)
{
	int	is_Imaginary=FALSE;			// Return value.

	// Check if any of the three points is imaginary (negative).
	if ((graph->nodes[nodeIndex].points_Index[0] < 0) ||
		(graph->nodes[nodeIndex].points_Index[1] < 0) ||
		(graph->nodes[nodeIndex].points_Index[2] < 0))
	{
		is_Imaginary = TRUE;
	}

	return(is_Imaginary);
}


/***************************************************************************
* Name: finalize_Graph
* IN:		N/A
* OUT:		N/A
* IN/OUT	graph			graph data
* RETURN:	N/A
* Description: deallocates graph data.
***************************************************************************/
void finalize_Graph(struct Graph_T *graph)
{
    // Set invalid values in variables.
    graph->size = -1;
    graph->nElements = 0;

    // Free allocated memory.
    if (graph->nodes != NULL)
    {
        free(graph->nodes);
        graph->nodes = NULL;
    }

	// Free face nodes assignement.
	if (graph->face_Node_Assign != NULL)
	{
		free(graph->face_Node_Assign);
        graph->face_Node_Assign = NULL;
	}
}


/***************************************************************************
* Name: reset_Graph
* IN:		N/A
* OUT:		N/A
* IN/OUT	graph			graph data
* RETURN:	N/A
* Description: sets the graph data to zero.
***************************************************************************/
void reset_Graph(struct Graph_T *graph)
{
    // Reset number of elements.
    graph->nElements = 0;

    // Reset array.
    memset( graph->nodes, 0, sizeof(struct Node_T)*graph->size);
	memset( graph->face_Node_Assign, 0, sizeof(int)*graph->size);
}


/***************************************************************************
* Name: copy_Node
* IN:		in
* OUT:		out
* IN/OUT	N/A
* RETURN:	N/A
* Description: copies "in" node into "out" node.
***************************************************************************/
void copy_Node(struct Node_T *in, struct Node_T *out)
{
    // Copy number of children.
    out->nChildren = in->nChildren;

    // Copy children.
    out->children_Index[0] = in->children_Index[0];
    out->children_Index[1] = in->children_Index[1];
    out->children_Index[2] = in->children_Index[2];

	// Copy face id.
	out->face_ID = in->face_ID;
}


/***************************************************************************
* Name: resize_Graph
* IN:		graph			graph data
* OUT:		N/A
* IN/OUT	N/A
* RETURN:	N/A
* Description: copies "in" node into "out" node.
***************************************************************************/
int resize_Graph(struct Graph_T *graph)
{
    int             ret=SUCCESS;        		// Return value.
	struct Node_T   *nodes_Ref=NULL;        	// Array of elements.
	int				*face_Node_Assign_Ref=NULL; // Faces and nodes assignment reference.

	// Save reference to current graph buffers.
	nodes_Ref = graph->nodes;
	face_Node_Assign_Ref = graph->face_Node_Assign;

    // Allocate new graph.
    graph->nodes = (struct Node_T *) calloc( graph->size*2, sizeof(struct Node_T));
	graph->face_Node_Assign = (int *) calloc( graph->size*2, sizeof(int));

	// Check if memory correctly allocated.
	if ((graph->nodes == NULL) ||
		(graph->face_Node_Assign == NULL))
	{
		// Reset pointers.
		graph->nodes = nodes_Ref;
		graph->face_Node_Assign = face_Node_Assign_Ref;
#ifdef LOGGING
		sprintf( log_Text, "Error allocating memory in resize_Graph.");
		write_Log( log_Text);
#endif
		printf("Error allocating memory in resize_Graph.");
		ret = FAILURE;
	}
	else
	{
		// Copy graph content.
		memcpy( graph->nodes, nodes_Ref, sizeof(struct Node_T)*graph->size);
		memcpy( graph->face_Node_Assign, face_Node_Assign_Ref, sizeof(int)*graph->size);

		// Set graph size as double.
		graph->size = graph->size*2;

		// Deallocate memory.
		free(nodes_Ref);
		free(face_Node_Assign_Ref);
	}

	return(ret);
}


/***************************************************************************
* Name: print_Node
* IN:		node			node data
* OUT:		N/A
* IN/OUT	N/A
* RETURN:	N/A
* Description: prints node type (leaf or branch), children and face.
***************************************************************************/
void print_Node(struct Node_T *node)
{
    int     i=0;            // Loop counter.

    // Check if node is a leaf or an interior node.
    if (node->nChildren == 0)
    {
        printf("Node is a leaf.");
    }
    else
    {
        printf("Node is interior.\n");

        // Print children nodes.
        for (i=0; i<node->nChildren ;i++)
        {
            printf("\tChild %d id is: %d", i+1, node->children_Index[i]);
        }
    }
    printf("\nNode face id: %d\n\n", node->face_ID);
}


/***************************************************************************
* Name: print_Graph
* IN:		graph			graph data
* OUT:		N/A
* IN/OUT	N/A
* RETURN:	N/A
* Description: prints all nodes in graph.
***************************************************************************/
void print_Graph(struct Graph_T *graph)
{
    int     i=0;            // Loop counter.
    int     len=0;          // Graph length.

    // Get graph length.
    len = get_Graph_Length( graph);
    printf("Graph length is: %d.\n", len);

    // Print graph nodes.
    for (i=0; i<len ;i++)
    {
        printf("Node %d.\n", i);

        // Print node at i-th position.
        print_Node( get_Node( graph, i));
    }
}

