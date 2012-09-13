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
//	Various tools and auxilary functions
//
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "tools.h"
#ifndef WIN32
#include <unistd.h>
#endif

unsigned long choose( unsigned n, unsigned k )
{
	const unsigned long uupSize = 100;
	static unsigned long uup[uupSize];
	register unsigned long *up;
	static unsigned long Nold = 0;
	static unsigned long Kold = 0;
	register unsigned long l,m;
	register unsigned i,j;

	if ( (n < k) || !n ) return 0;

	if ( (n == k) || !k ) 	// includes n == 1
		return 1;

	if ( k > (n >> 1) )	// Only lower half
		k = n-k;

	if ( (Nold == n) && (k < Kold) )	// We did it last time ...
		return *(uup + k - 1);

	if ( k > uupSize )
	{
		printf( "choose( unsigned, unsigned) :  overflow\n");
		getchar(); exit(-1);
	}

	Nold=n; Kold=k;

	*(up=uup)=2;
	for (i=2; i<n; i++)		// Pascal's triangle
	{
        // todo: remove next line:
		*(up+1)=1;
		l=1;
        m=*(up=uup);
		for (j=0; j<mmin(i,k); j++)
		{
			*up=m+l;
			l=m;
			m=*(++up);
		}
        // todo: remove next line:
		*up=1;
	}

	return *(uup + k - 1);
}
	
unsigned long mysrand; 	
double rand1()
{
    mysrand=1664525*mysrand+1013904223L;
    double r=((double)mysrand)/4294967297.0;
    //if (r>.52)
    //{
    //    printf("whoups\n");
    //}
    return r;

}


void initRandom(int i)
{
    if (i) { mysrand=i; return; }
    mysrand=(unsigned long) (clock());
}

void error(char *s)
{
   printf("Error due to %s.", s);
   getchar(); exit(255);
};

double euclidianNorm(int i, double *xp)
{
// no tested
// same code for the Vector eucilidian norm and for the Matrix Froebenis norm
/*
    double sum=0;
    while (i--) sum+=sqr(*(xp++));
    return sqrt(sum);
*/
    const double SMALL=5.422e-20, BIG=1.304e19/((double)i);
    double s1=0,s2=0,s3=0, x1max=0, x3max=0, xabs;

    while (i--)
    {
        xabs=abs(*(xp++));
        
        if (xabs>BIG)
        {
            if (xabs>x1max)
            {
                s1=1.0+s1*sqr(x1max/xabs);
                x1max=xabs;
                continue;
            }
            s1+=sqr(xabs/x1max);
            continue;
        }
        if (xabs<SMALL)
        {
            if (xabs>x3max)
            {
               s3=1.0+s3*sqr(x3max/xabs);
               x3max=xabs;
               continue;
            }
            if (xabs!=0) s3+=sqr(xabs/x3max);
            continue;
        }
        s2+=sqr(xabs);
    };
    if (s1!=0) return x1max*sqrt(s1+(s2/x1max)/x1max);
    if (s2!=0)
    {
        if (s2>=x3max) return sqrt(s2*(1.0+(x3max/s2)*(x3max*s3)));
        return sqrt(x3max*((s2/x3max)+(x3max*s3)));
    }
    return x3max*sqrt(s3);
}

#define EOL1 13
#define EOL2 10

char isEmpty(const char *line)
{
    line=skipSpaces(line);
    return (*line==0);
}

const char *skipSpaces(const char *t)
{
    while ((*t==' ')||(*t=='\t')) t++;
    return t;
}

char isEmptyline(const char *line)
{
	if (*line==';') return 1;
	line=skipSpaces(line);
	if ((*line==EOL1)||(*line==EOL2)||(*line=='\0')) return 1;
	return 0;
}

char *GetRidOfTheEOL(char *tline)
{
   char *t;
   t=tline=(char*)skipSpaces(tline);
   while ((*tline!=EOL1)&&(*tline!=EOL2)&&(*tline)) tline++;
   *tline='\0';
   return t;
}

char *removeQuotes(char *t)
{
    if ((*t=='\'')||(*t=='"'))
    {
        t[strlen(t)-1]=0;
        return t+1;
    }
    return t;
}

char **getNameTable(const char *line, int *ncolumn)
{
    char *n=(char*)line;
    int j=0,nc;
    if (ncolumn) nc=*ncolumn;
    if (nc==0)
    {
        while (*n) 
        {
            while ((*n)&&(*n!='\t')&&(*n!=' ')&&(*n!=13)&&(*n!=10)) n++;
            nc++;
            while ((*n)&&((*n=='\t')||(*n==' ')||(*n==13)||(*n==10))) n++;
        }
    }
    if (ncolumn) *ncolumn=nc;
    char **names=(char**)malloc(nc*sizeof(char**));
    n=(char*)malloc(strlen(line)+1);
    strcpy(n,line);
    
    for (j=0; j<nc-1; j++)
    {
        names[j]=n;
        while ((*n)&&(*n!='\t')&&(*n!=' ')&&(*n!=13)&&(*n!=10)) n++;
        if (*n) { *n=0; n++; }
        while ((*n)&&((*n=='\t')||(*n==' ')||(*n==13)||(*n==10))) n++;
    }
    names[j]=n;
    return names;
}

char *stringDuplicate(const char *line)
{
    int l=(int)strlen(line);
    // remove quotes:
    if ((*line=='\'')||(*line=='"')) { l-=2; line++; }
    char *t=(char*)malloc(l+1);
	memcpy(t,line,l);
    t[l]=0;
    GetRidOfTheEOL(t);
    return t;
}

char *removeAllEOL(char *t)
{
    char *t2=t;
    while (*t)
    {
        if ((*t=='\r')||(*t=='\n')) *t=' '; 
        t++;
    }
    return t2;
}

const char *skipLine(const char *t)
{
    while ((*t!='\r')&&(*t!='\n')&&(*t)) t++; 
    if (*t=='\n') { if (*(t+1)=='\r') t++; }
    else if (*t=='\r') { if (*(t+1)=='\n') t++; }
    return t+1;
}

char *loadFile(FILE *f)
{
//    FILE *f=fopen(filename,"rb");
    fseek(f,0,SEEK_END);
    int l=ftell(f);
    fseek(f,0,SEEK_SET);
    char *buf=(char*)malloc(l+1);
    fread(buf,l,1,f);
    fclose(f);
    buf[l]=0;
    return buf;
}

void deleteFile(const char *line)
{
#ifdef WIN32
    remove(line);
#else
    unlink(line);
#endif
}

