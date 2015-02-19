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
#include "MatrixTriangle.h"

MatrixTriangle::MatrixTriangle(int _n)
{
    d=(MatrixTriangleData*)malloc(sizeof(MatrixTriangleData));
    d->n=_n; d->ext=_n;
    d->ref_count=1;
    if (_n>0)
	{
	    double **t,*t2;
        int i=1;
    	t=d->p=(double**)malloc(_n*sizeof(double*));
    	t2=(double*)malloc((_n+1)*_n/2*sizeof(double));
	    while(_n--)
        {
            *(t++)=t2; t2+=i; i++;
        }
    } else d->p=NULL;
}

void MatrixTriangle::setSize(int _n)
{
    d->n=_n;
    if (_n>d->ext)
    {
        d->ext=_n;
	    double **t,*t2;
        if (!d->p)
        {
    	    t2=(double*)malloc((_n+1)*_n/2*sizeof(double));
            t=d->p=(double**)malloc(_n*sizeof(double));
        } else
        {
    	    t2=(double*)realloc(*d->p,(_n+1)*_n/2*sizeof(double));
            t=d->p=(double**)realloc(d->p,_n*sizeof(double));
        }
        int i=1;
	    while(_n--)
        {
            *(t++)=t2; t2+=i; i++;
        }
    }
}

void MatrixTriangle::solveInPlace(Vector b)
{
    int i,k,n=nLine();
    double **a=(*this), *x=b, sum;

    if ((int)b.sz()!=n)
    {
        printf("error in matrixtriangle solve.\n"); getchar(); exit(254);
    }
    for (i=0; i<n; i++)
    {
        sum=x[i]; k=i;
        while (k--) sum-=a[i][k]*x[k];
        x[i]=sum/a[i][i];
    }
}

void MatrixTriangle::solveTransposInPlace(Vector y)
{
    int n=nLine(),i=n,k;
    double **a=(*this), *x=y, sum;

    while(i--)
    {
        sum=x[i];
        for (k=i+1; k<n; k++) sum-=a[k][i]*x[k];
        x[i]=sum/a[i][i];
    }
}
/*
void MatrixTriangle::invert()
{
    int i,j,k,n=nLine();
    double **a=(*this), sum;
    for (i=0; i<n; i++)
    {
        a[i][i]=1/a[i][i];
        for (j=i+1; j<n; j++) 
        {
            sum=0;
            for (k=i; k<j; k++) sum-=a[j][k]*a[k][i];
            a[j][i]=sum/a[j][j];
        }
    }
}
*/
void MatrixTriangle::LINPACK(Vector &R)
{
    int i,j,n=nLine();
    R.setSize(n);
    double **L=(*this), *w=R, sum;

    for (i=0; i<n; i++)
    {
        if (L[i][i]==0) w[i]=1.0;

        sum=0; j=i-1;
        if (i) while (j--) sum+=L[i][j]*w[j];
        if (((1.0-sum)/L[i][i])>((-1.0-sum)/L[i][i])) w[i]=1.0; else w[i]=-1.0;
    }
    solveTransposInPlace(R);
    R.multiply(1/R.euclidianNorm());
}

MatrixTriangle::~MatrixTriangle()
{
    destroyCurrentBuffer();
}

void MatrixTriangle::destroyCurrentBuffer()
{
    if (!d) return;
    (d->ref_count) --;
	if (d->ref_count==0)
    {
        if (d->p) { free(*d->p); free(d->p); }
        free(d);
    };
}

MatrixTriangle::MatrixTriangle(const MatrixTriangle &A)
{
    // shallow copy
    d=A.d;
	(d->ref_count)++ ;
}

MatrixTriangle& MatrixTriangle::operator=( const MatrixTriangle& A )
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

MatrixTriangle MatrixTriangle::clone()
{
    // a deep copy
    MatrixTriangle r(nLine());
    r.copyFrom(*this);
    return r;
}

void MatrixTriangle::copyFrom(MatrixTriangle r)
{
    int n=r.nLine();
    setSize(n);
    if (n==0) return;
    memcpy(*d->p,*(r.d->p),(n+1)*n/2*sizeof(double));
}

