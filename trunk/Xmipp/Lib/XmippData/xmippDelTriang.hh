#ifndef DELTRIANG_H
#define DELTRIANG_H

#include <vector>
using std::vector;
///////////////////////////////Data Types

typedef struct {
   int p1,p2,p3;
} ITRIANGLE;
typedef struct {
   int p1,p2;
} IEDGE;
typedef struct {
   double x,y,z;
} XYZ;




////////////////////////////////  FUNCTIONS

int Triangulate(int ,XYZ *,ITRIANGLE *,int *);

int CircumCircle(double ,double ,
   double ,double ,double ,double ,double ,double ,
   double *,double *,double *);
   
// Returns index of triangle in vector LatTri
//int FindNearest(XYZ &, vector <XYZ> &,vector <ITRIANGLE> & );

#endif

