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
#ifndef _INCLUDE_VectorIntINT_H
#define _INCLUDE_VectorIntINT_H

#include <stdlib.h>

class VectorInt
{
  public: 
    typedef struct VectorIntDataTag
    {
       int n,extention;
       int ref_count;
       int *p;
    } VectorIntData;
    VectorIntData *d;

    VectorInt(int _n=0);
    VectorInt(int _n, int _ext);
    VectorInt(int _n, int *d);

// allow shallow copy:
    VectorInt clone();
    void copyFrom(VectorInt r);
    VectorInt( const VectorInt& P );
	VectorInt& operator=( const VectorInt& P );
    void destroyCurrentBuffer();
    ~VectorInt();

    void extend();
    void setSize(int _n);
    void exactshape();
    void print();
    
    // only use the following method at your own risks!
    void prepareExtend(int new_extention);
    
//    int &operator [](int i) { return p[i]; };
    inline int sz() {return d->n;};
    int operator==( const VectorInt& P) {return (d==P.d);}
    operator int*() const { return d->p; };
//    operator unsigned*() const { return (unsigned*)p; };
//    int &operator[]( unsigned i) {return p[i];};
    int equals( const VectorInt& Q );
    
    static VectorInt emptyVectorInt;

private:
    void alloc(int, int);

};

#endif

