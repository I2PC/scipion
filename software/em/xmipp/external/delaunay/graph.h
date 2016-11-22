#ifndef GRAPH_H
#define GRAPH_H
#include <stdio.h>

/**************************************************************************
* Defines
**************************************************************************/
#define     N_POINTS       3


/**************************************************************************
* Data types definition
**************************************************************************/
struct Node_T
{
    int nChildren;                  // Number of triangles inside this triangle.
    int children_Index[3];          // Position of children triangles in array.
	int	points_Index[3];			// Vertexes of node.
    int face_ID;                    // Face id.
};

struct Graph_T
{
    int             size;				// Max number of elements.
    int             nElements;			// Current number of elements.
    struct Node_T   *nodes;				// Array of elements.
	int				*face_Node_Assign;	// Relation between a face and its node.
};



/**************************************************************************
* Public function headers.
**************************************************************************/
int             initialize_Graph(struct Graph_T *graph, int size);
void            finalize_Graph(struct Graph_T *graph);
void            reset_Graph(struct Graph_T *graph);
void            copy_Node(struct Node_T *in, struct Node_T *out);
int             resize_Graph(struct Graph_T *graph);
int             get_Graph_Length(struct Graph_T *graph);
int             get_Graph_Size(struct Graph_T *graph);
int             get_nChildren_Node(struct Graph_T *graph, int index);
int				get_iChildren_Node(struct Graph_T *graph, int node_Index, int child_Index);
int				get_Face_Of_Node(struct Graph_T *graph, int index);
int             is_Graph_Full(struct Graph_T *graph);

// Graph modification interface.
int            	insert_Node(struct Graph_T *graph, struct Node_T *node);
struct Node_T   *get_Node(struct Graph_T *graph, int index);
void            update_Node(struct Graph_T *graph, int index, int  nChildren, struct Node_T *node);
int             is_Leaf_Node(struct Graph_T *graph, int index);
int				get_Node_Assigned( struct Graph_T *graph, int face_ID);
int				is_Any_Point_Imaginary( struct Graph_T *graph, int nodeIndex);
int				number_Imaginary_Points( struct Graph_T *graph, int nodeIndex);

// I/O interface.
void print_Node(struct Node_T *node);
void print_Graph(struct Graph_T *graph);


#endif
