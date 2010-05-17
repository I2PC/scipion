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
//	Header file for tools
//


#ifndef _MPI_TOOLS_H_
#define _MPI_TOOLS_H_

#include <math.h>

//#define SWAP(a,b) {tempr=(a); (a)=(b); (b)=tempr;}

#define maxDWORD 4294967295 //2^32-1;
#define INF 1.7E+308 
#define EOL 10
#ifndef PI
    #define PI 3.1415926535897932384626433832795
#endif
#define ROOT2 1.41421356
 
inline double condorAbs( const double t1 )
{
	return t1 > 0.0 ? t1 : -t1;
}

inline double sign( const double a )
 // equivalent to sign(1,a)
{
    return a<0?-1:1;
}

inline double sign( const double t1, const double t2 )
{
    if(t2>=0) return condorAbs(t1);
    return -condorAbs(t1);
}

inline double isInt( const double a)
{
    return condorAbs(a-floor( a + 0.5 ))<1e-4;
}

inline double mmin( const double t1, const double t2 )
{
	return t1 < t2 ? t1 : t2;
}

inline unsigned mmin( const unsigned t1, const unsigned t2 )
{
	return t1 < t2 ? t1 : t2;
}

inline int mmin( const int t1, const int t2 )
{
	return t1 < t2 ? t1 : t2;
}

inline double mmax( const double t1, const double t2 )
{
	return t1 > t2 ? t1 : t2;
}

inline unsigned mmax( const unsigned t1, const unsigned t2 )
{
	return t1 > t2 ? t1 : t2;
}

inline int mmax( int t1, int t2 )
{
	return t1 > t2 ? t1 : t2;
}

inline double sqr( const double& t )
{
	return t*t;
}

inline double round (double a)
{
    return (int)(a+.5);
}

unsigned long choose( unsigned n, unsigned k );
double rand1();
void initRandom(int i=0);
double euclidianNorm(int i, double *xp);

#include "Vector.h"
#include "Matrix.h"

//void saveValue(Vector tmp,double valueOF, Matrix data);

char emptyline(const char *line);
char *GetRidOfTheEOL(char *tline);
char *removeQuotes(char *t);
char **getNameTable(const char *line, int *ncolumn);
char *stringDuplicate(const char *l);
char isEmpty(const char *line);
char *removeAllEOL(char *t);
char *loadFile(FILE *f);
const char *skipLine(const char *l);
const char *skipSpaces(const char *l);
void deleteFile(const char *l);

#endif 	/* _MPI_TOOLS_H_ */

