/*

CONDOR 1.06 - COnstrained, Non-linear, Direct, parallel Optimization 
              using trust Region method for high-computing load, 
              noisy functions
Copyright (C) 2004 Frank Vanden Berghen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation version 2
of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

If you want to include this tools in any commercial product, 
you can contact the author at fvandenb@iridia.ulb.ac.be

*/

#ifndef _MATRIXTRIANGLE_H
#define _MATRIXTRIANGLE_H

#include "Vector.h"

class Matrix;

class MatrixTriangle // lower triangular
{
  friend class Matrix;
  protected:
    void destroyCurrentBuffer();
    typedef struct MatrixTriangleDataTag
    {
       int n;
       int ext;
       int ref_count;
       double **p;
    } MatrixTriangleData;
    MatrixTriangleData *d;

  public:

// creation & management of Matrix:
    MatrixTriangle(int _n=0);
    void setSize(int _n);

// allow shallow copy:
    ~MatrixTriangle();
    MatrixTriangle(const MatrixTriangle &A);
    MatrixTriangle& operator=( const MatrixTriangle& A );
    MatrixTriangle clone();
    void copyFrom(MatrixTriangle r);

// accessor method
    inline bool operator==( const MatrixTriangle& A ) { return (A.d==d);}
    inline int nLine() { return d->n; };
    inline double *operator [](int i) { return d->p[i]; };
    inline operator double**() const { return d->p; };
    
// simple math tools:
    void solveInPlace(Vector b);
    void solveTransposInPlace(Vector y);
    //void invert();
    void LINPACK(Vector &u);
    
// default return matrix in case of problem in a function
    static MatrixTriangle emptyMatrixTriangle;
};

#endif
