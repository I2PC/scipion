#ifndef DELTRIANG_H
#define DELTRIANG_H

#include <vector>

typedef struct
{
    int p1, p2, p3;
}
ITRIANGLE;

typedef struct
{
    int p1, p2;
}
IEDGE;

typedef struct
{
    double x, y, z;
}
XYZ;

int Triangulate(int, XYZ*, ITRIANGLE*, int*);

bool CircumCircle(double , double ,
                 double , double , double , double , double , double ,
                 double *, double *, double *);

// Returns index of triangle in vector LatTri
//int FindNearest(XYZ &, std::vector <XYZ> &,std::vector <ITRIANGLE> & );

#endif

