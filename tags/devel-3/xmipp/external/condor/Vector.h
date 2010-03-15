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

#ifndef _INCLUDE_VECTOR_H
#define _INCLUDE_VECTOR_H

//#include <stdlib.h> // for the declaration of NULL
#include <stdio.h> // for the declaration of FILE*
#include "VectorInt.h"

class Matrix;

class Vector
{
  public: 
    // only use the following method at your own risks!
    void prepareExtend(int new_extention);
    void alloc(int n, int ext);
    typedef struct VectorDataTag
    {
       int n,extention;
       int ref_count;
       double *p;
       char externalData;
    } VectorData;
    VectorData *d;

// creation & management of Vector:
    Vector(int _n=0);
    Vector(int _n, int _ext);
    Vector(int _n, double *dd, char externalData=0);
    Vector(char *filename);
    Vector(char *line, int guess_on_size);
    Vector(Vector a,Vector b);
    
    char *getFromLine(char *line);
    void extend();
    void setSize(int _n);
    void exactshape();
    void print();
    void save(char *filename, char ascii);
    void save(FILE *f, char ascii);
    void appendToMatrixFile(char *saveFileName, char **names=NULL);
    void setExternalData(int _n, double *dd);

// allow shallow copy:
    Vector clone();
    void copyFrom(Vector r, int _n=0);
    Vector( const Vector& P );
	Vector& operator=( const Vector& P );
    void destroyCurrentBuffer();
    ~Vector();

// accessor method
    inline unsigned sz() {return d->n;};
//  inline double &operator [](int i) { return d->p[i]; };
    inline int operator==( const Vector Q) { return d==Q.d; };
    int equals( const Vector Q );
    operator double*() const { if (d) return d->p; else return NULL; };
    //double &operator[]( unsigned i) {return p[i];};
	void setPart(int i, Vector v, int n=0, int ii=0);

// simple math tools:
    double euclidianNorm();
    double L1Norm();
    double LnftyNorm();
    double euclidianDistance(Vector v);
    double L1Distance(Vector v);
    double LnftyDistance(Vector v);

    double square();
    void multiply(double a);
    void multiply(Vector R, double a);
    void zero(int _i=0, int _n=0);
    void set(double dd);
    void shift(int s);
    double scalarProduct(Vector v);
    void oneByOneMutiply(Vector r);
    void oneByOneInvert();
    double mmin();
    double mmax();
    bool isNull();
    Vector operator-( Vector v);
    Vector operator+( Vector v);
    Vector operator-=( Vector v);
    Vector operator+=( Vector v);
    void addInPlace(double a, Vector v); // this+=a*v
    void addInPlace(double a, int i, Matrix m); // this+=a * M(i,:)
    void transposeAndMultiply(Vector vR, Matrix M); // result in vR
    void diagonalizeAndMultiply(Matrix M); 
    void permutIn(Vector vR, VectorInt viP);
    void permutOut(Vector vR, VectorInt viP);

// default return Vector in case of problem in a function
    static Vector emptyVector;
    
};

#endif

