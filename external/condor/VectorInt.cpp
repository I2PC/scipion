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

#include <stdio.h>
#include <memory.h>
#include "VectorInt.h"

VectorInt VectorInt::emptyVectorInt;

void VectorInt::alloc(int _n, int _ext)
{
    d=(VectorIntData*)malloc(sizeof(VectorIntData));
    d->n=_n;
    d->extention=_ext;
    d->ref_count=1;

    if (_ext==0) { d->p=NULL; return; };

    d->p=(int*)malloc(_ext*sizeof(int)); 
    if ((d->p)==NULL) { printf("memory allocation error\n"); exit(253); }
}

VectorInt::VectorInt(int _n)
{
    alloc(_n, _n);
    memset(d->p,0,d->extention*sizeof(int));
}

VectorInt::VectorInt(int _n, int _ext)
{
    alloc(_n, _ext);
    memset(d->p,0,d->extention*sizeof(int));
}

VectorInt::VectorInt(int _n, int *_d)
{
    alloc(_n,_n);
    if (_d) memcpy(d->p,_d,_n*sizeof(int));
    else memset(d->p,0,d->extention*sizeof(int));
}

void VectorInt::prepareExtend(int new_extention)
{
	if (d->extention<new_extention)
	{
		d->p=(int*)realloc(d->p,new_extention*sizeof(int));
        memset(d->p+d->extention,0,(new_extention-d->extention)*sizeof(int));
    	d->extention=new_extention;
	};
}

VectorInt::~VectorInt()
{
    destroyCurrentBuffer();
}

void VectorInt::destroyCurrentBuffer()
{
    if (!d) return;
    (d->ref_count) --;
	if (d->ref_count==0)
    {
        if (d->p) free(d->p);
        free(d);
    };
}

VectorInt::VectorInt(const VectorInt &A)
{
    // shallow copy
    d=A.d;
	(d->ref_count)++ ;
}

VectorInt& VectorInt::operator=( const VectorInt& A )
{
    // shallow copy
    if (this != &A)
	{
        destroyCurrentBuffer();
        d=A.d;
		(d->ref_count) ++ ;
	}
	return *this;
}

VectorInt VectorInt::clone()
{
    // a deep copy
    VectorInt r(sz());
    r.copyFrom(*this);
    return r;
}

void VectorInt::copyFrom(VectorInt r)
{
    unsigned n=r.sz();
    if (n==0) return;
    setSize(n);
    memcpy(d->p,r.d->p,n*sizeof(double));
}

void VectorInt::setSize(int _n)
{
    d->n=_n;
    if (_n==0) { free(d->p); d->p=NULL; d->extention=0; return; }
    prepareExtend(_n);
}

void VectorInt::extend()
{
    d->n++;
    if (d->n>d->extention) prepareExtend(d->extention+100);
}

void VectorInt::exactshape()
{
	if (d->extention!=d->n)
	{
		d->p=(int*)realloc(d->p,d->n*sizeof(int));
		d->extention=d->n;
	};
}

int VectorInt::equals( const VectorInt& Q )
{
  if (d->n != Q.d->n) return 0;

  int *cP = d->p, *cQ = Q.d->p;
  int i = d->n;

  while( i-- )
  {
    if (*cP!=*cQ) return 0;
    cP++; cQ++;
  }

  return 1;
}

//ostream& VectorInt::PrintToStream( ostream& out ) const
void VectorInt::print()
{
	printf("[");
	if (!d->n || !d->p) { printf("]\n"); return; }

    int N=d->n,*up=d->p;
	while (--N) printf("%i,",*(up++));
	printf("%i]\n",*up);
}
