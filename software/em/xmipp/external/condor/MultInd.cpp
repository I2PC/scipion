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

//
//	Multiindex
//
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

#define __INSIDE_MULTIND_CPP__
#include "MultInd.h"
#undef __INSIDE_MULTIND_CPP__

#include "tools.h"

unsigned MultInd::maxDim;
unsigned *MultInd::buffer;
MultIndCache cacheMultInd;

MultInd& MultInd::operator=( const MultInd &p )
{
    dim=p.dim; deg=p.deg; next=NULL; 
    lastChangesV=p.lastChangesV; indexesOfCoefInLexOrderV=p.indexesOfCoefInLexOrderV;
    indexV=p.indexV;
    standardInit();
    if (deg==0) memcpy(coeffDeg,p.coeffDeg,dim*sizeof(unsigned));
    return *this;
}

void MultInd::standardInit()
{
    if (deg==0)
    {
        coeffDeg=(unsigned*)malloc(dim*sizeof(unsigned));
        coeffLex=NULL;
    } else
    {
        coeffDeg=buffer;
        coeffLex=buffer+dim;
    };
}

MultInd::MultInd( unsigned _dim, unsigned _deg): 
    dim(_dim), deg(_deg), next(NULL)
{
    standardInit();
    fullInit(); 
    resetCounter();
}

MultInd::MultInd(unsigned d):  dim(d), deg(0), next(NULL) 
{
    standardInit();
    resetCounter();
}

MultInd::~MultInd()
{
    if (deg==0) free(coeffDeg);
}

void MultInd::print()
{
	printf("[");
	if (!dim) { printf("]"); return; }

    unsigned N=dim,*up=coeffDeg;
	while (--N) printf("%i,",*(up++));
	printf("%i]",*up);
}

unsigned MultInd::len()
{
    unsigned l=0, *ccDeg=coeffDeg, j=dim;
    while (j--) l+=*(ccDeg++);
    return l;
}

bool MultInd::operator==( const MultInd& m )
{
    unsigned *p1=(*this), *p2=m, n=dim;
    while (n--)
        if (*(p1++)!=*(p2++)) return false;
    return true;
}

void MultInd::fullInit()
{
    unsigned *ccLex, *ccDeg, degree=deg, n=choose(dim+deg,dim),i,k,sum, d=dim-1;
    int j;

    lastChangesV.setSize(n-1);
    indexesOfCoefInLexOrderV.setSize(n);

    memset(coeffLex+1,0,d*sizeof(int));
	*coeffLex=deg;

    for (i=0; i<n; i++)
    {
        sum=0; ccLex=coeffLex; j=dim;
        while (j--) sum+=*(ccLex++);
        if (sum) k=choose( sum+d, dim ); else k=0;

        resetCounter();
        *coeffDeg=sum;

        while(1)
        {
            ccLex=coeffLex; ccDeg=coeffDeg;
    	    for ( j=d; j>0 ; j--, ccLex++, ccDeg++ ) if (*ccLex != *ccDeg) break;
            if (*ccLex >= *ccDeg) break;
            ++(*this); k++;
        }

        indexesOfCoefInLexOrderV[i]=k;

        if (i==n-1) break;

        // lexical order ++ :
        if (coeffLex[d])
        {
           lastChangesV[i]=d;
           coeffLex[d]--;
        } else
        for (j=d-1; j>=0; j--)
        {
	        if (coeffLex[j])
	        {
	          lastChangesV[i]=j;
	          sum=--coeffLex[j];
	          for (k=0; k<(unsigned)j; k++) sum+=coeffLex[k];
	          coeffLex[++j]=degree-sum;
	          for (k=j+1; k<=d; k++) coeffLex[k]=0;
	          break;
	        }
        }
    }
}

void MultInd::resetCounter()
{
    indexV=0;
    memset(coeffDeg,0,dim*sizeof(unsigned));
}
	
MultInd& MultInd::operator++()
{
	unsigned *cc = coeffDeg;
    int n=dim, pos, i;

	if (!n || !cc) return *this;

	for (pos = n-2; pos >= 0; pos--)
	{
		if (cc[pos])	// Gotcha
		{
		  cc[pos]--;
		  cc[++pos]++;
		  for (i = pos+1; i < n;i++)
		  {
		    cc[pos] += cc[i];
		    cc[i] = 0;
		  }
          indexV++;
		  return *this;
		}
	}

	(*cc)++;
	for ( i = 1; i < n; i++)
	{
		*cc += cc[i];
		cc[i] = 0;
	}

    indexV++;
    return *this;
}

unsigned *MultInd::lastChanges()
{
    if (deg==0) 
    {
        printf("use MultIndCache to instanciate MultInd");
        getchar(); exit(252);
    }
    return (unsigned*)lastChangesV.d->p;
}

unsigned *MultInd::indexesOfCoefInLexOrder()
{
    if (deg==0) 
    {
        printf("use MultIndCache to instanciate MultInd");
        getchar(); exit(252);
    }
    return (unsigned*)indexesOfCoefInLexOrderV.d->p;
}

MultIndCache::MultIndCache(): head(NULL) 
{
    MultInd::maxDim=100;
    MultInd::buffer=(unsigned*)malloc(MultInd::maxDim*2*sizeof(unsigned));
}

MultIndCache::~MultIndCache()
{
    MultInd *d=head, *d1;
    while (d)
    {
        d1=d->next;
        delete d;
        d=d1;
    }
    free(MultInd::buffer);
}

MultInd *MultIndCache::get(unsigned _dim, unsigned _deg )
{
    if (_deg==0)
    {
        printf("use normal constructor of MultiInd");
        getchar(); exit(252);
    }
    if (_dim>MultInd::maxDim)
    {
        free(MultInd::buffer);
        MultInd::maxDim=_dim;
        MultInd::buffer=(unsigned*)malloc(_dim*2*sizeof(unsigned));
    }
    MultInd *d=head;
    while (d)
    {
        if ((_dim==d->dim)&&(_deg==d->deg)) return d;
        d=d->next;
    }

    d=new MultInd(_dim,_deg);
    d->next=head;
    head=d;
    return d;
}
