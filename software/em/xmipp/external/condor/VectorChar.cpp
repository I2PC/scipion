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
#include "VectorChar.h"

#define CHECK(p) if ((p)==NULL) { printf("memory allocation error\n"); exit(253); }

void VectorChar::alloc()
{
    if (extention==0) { p=NULL; return; };
    CHECK(p=(char*)malloc(extention*sizeof(char))); 
}

VectorChar::VectorChar(int _n): np(_n), extention(_n), n(np)
{
    alloc();
    memset(p,0,extention*sizeof(char));
}

VectorChar::VectorChar(int _n, int _ext): np(_n), extention(_ext), n(np)
{
    alloc();
    memset(p,0,extention*sizeof(char));
}

VectorChar::VectorChar(VectorChar *v): np(v->n), extention(v->n), n(np)
{
    alloc();
    memcpy(p,v->p,n*sizeof(char));
}

VectorChar::VectorChar(int _n, char *d): np(_n), extention(_n), n(np) 
{
    alloc();
    if (d) memcpy(p,d,_n*sizeof(char));
    else memset(p,0,extention*sizeof(char));
}

void VectorChar::prepareExtend(int new_extention)
{
	if (extention<new_extention)
	{
		CHECK(p=(char*)realloc(p,new_extention*sizeof(char)));
        memset(p+extention,0,(new_extention-extention)*sizeof(char));
    	extention=new_extention;
	};
}

void VectorChar::setSize(int _n)
{
    np=_n;
    if (_n==0) { free(p); p=NULL; return; }
    prepareExtend(_n);
}

void VectorChar::extend()
{
    np++;
    if (np>extention) prepareExtend(extention+100);
}

void VectorChar::exactshape()
{
	if (extention!=0)
	{
		CHECK(p=(char*)realloc(p,n*sizeof(char)));
		extention=np;
	};
}

VectorChar::~VectorChar()
{
    if (p) free(p);
}

char VectorChar::operator==( const VectorChar& Q )
{
  if (n != Q.n) return 0;

  char *cP = p, *cQ = Q.p;
  int i = n;

  while( i-- )
  {
    if (*cP!=*cQ) return 0;
    cP++; cQ++;
  }

  return 1;
}

VectorChar& VectorChar::operator=( const VectorChar& P )
{
    if (extention<P.n)
    {
        extention=P.n;
        free(p); alloc();
    }
    np=P.n;
    memcpy(p,P.p,n*sizeof(char));
    return *this;
}

VectorChar::VectorChar( const VectorChar& P ): np(P.n), extention(P.n), n(np)
{
    alloc();
    memcpy(p,P.p,n*sizeof(char));
}

//ostream& VectorChar::PrintToStream( ostream& out ) const
void VectorChar::print()
{
	printf("[");
	if (!n || !p) { printf("]"); return; }

    int N=n;
    char *up=p;
	while (--N) printf("%i,",*(up++));
	printf("%i]",*up);
}

void VectorChar::set(char c)
{
    memset(p,c,n);
}
