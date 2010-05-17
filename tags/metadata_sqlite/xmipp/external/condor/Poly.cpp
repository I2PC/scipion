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
//#include <crtdbg.h> 
#include "Vector.h"
#include "MultInd.h"
#include "tools.h"
#include "Poly.h"
#include "IntPoly.h"

const unsigned int Polynomial::NicePrint = 1;
const unsigned int Polynomial::Warning   = 2;
const unsigned int Polynomial::Normalized= 4;  // Use normalized monomials
      unsigned int Polynomial::flags = Polynomial::Warning||Polynomial::NicePrint;

Polynomial Polynomial::emptyPolynomial;

void Polynomial::init(int _dim, int _deg, double *data)
{
    int n;
    d=(PolynomialData*)malloc(sizeof(PolynomialData));
    if (_dim) n=d->n=choose( _dim+_deg, _dim );
    else n=d->n=0;

    d->dim=_dim;
    d->deg=_deg;
    d->ref_count=1;

    if (n==0) { d->coeff=NULL; return; };

    d->coeff=(double*)malloc(n*sizeof(double)); 
    if (d->coeff==NULL) { printf("memory allocation error\n"); getchar(); exit(253); }

    if (data) memcpy(d->coeff, data, d->n*sizeof(double));
    else memset(d->coeff, 0, d->n*sizeof(double));
}
/*
Polynomial::PolyInit( const Polynomial& p )
{ dim = p.dim; deg = p.deg; coeff=((Vector)p.coeff).clone(); }
*/
Polynomial::Polynomial( unsigned Dim, unsigned Deg, double *data ) 
{
    init(Dim,Deg,data);
}

Polynomial::Polynomial( unsigned Dim, double val ) // Constant polynomial
{ 
    init(Dim,0,&val); 
}

Polynomial::Polynomial( MultInd& I ) 
{
    init(I.dim,I.len());
	d->coeff[I.index()] = 1;
}

Polynomial::~Polynomial()
{
    destroyCurrentBuffer();
};

void Polynomial::destroyCurrentBuffer()
{
    if (!d) return;
    (d->ref_count) --;
	if (d->ref_count==0)
    {
        if (d->coeff) free(d->coeff);
        free(d);
    }
}

Polynomial& Polynomial::operator=( const Polynomial& A )
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

Polynomial::Polynomial(const Polynomial &A)
{
    // shallow copy
    d=A.d;
	(d->ref_count)++ ;
}

Polynomial Polynomial::clone()
{
    // a deep copy
    Polynomial m(d->dim,d->deg);
    m.copyFrom(*this);
    return m;
}

void Polynomial::copyFrom(Polynomial m)
{
    if (m.d->dim!=d->dim)
    {
        printf("poly: copyFrom: dim do not agree");
        getchar(); exit(254);
    }

    d->deg = mmax(d->deg,m.d->deg);	// New degree
    unsigned N1=sz(), N2=m.sz();
    if (N1!=N2) 
    {
        d->coeff=(double*)realloc(d->coeff,N2*sizeof(double));
        d->n=m.d->n;
    }
    memcpy((*this),m,N2*sizeof(double));
}

Polynomial Polynomial::operator*( const double t )
{
    int i=sz();
	Polynomial q( d->dim, d->deg );
	double *tq = q.d->coeff, *tp = d->coeff;

	while (i--) *(tq++) = *(tp++) * t;
	return q;
}


Polynomial Polynomial::operator/( const double t )
{
	if (t == 0)
	{
		printf( "op/(Poly,double): Division by zero\n");
		getchar(); exit(-1);
	}
    int i=sz();
	Polynomial q( d->dim, d->deg );
	double *tq = q.d->coeff, *tp = d->coeff;

	while (i--) *(tq++) = *(tp++) / t;
	return q;
}

Polynomial Polynomial::operator+( Polynomial q )
{
	if (d->dim != q.d->dim)
	{
		printf( "Poly::op+ : Different dimension\n");
		getchar(); exit(-1);
	}

	Polynomial r(d->dim,mmax(d->deg,q.d->deg));
	unsigned N1=sz(), N2=q.sz(), Ni=mmin(N1,N2);
	double	*tr = r, *tp = (*this), *tq = q;
    while (Ni--) *(tr++) = *(tp++) + *(tq++);
    if (N1<N2)
    {
        memcpy(tr,tq,(N2-N1)*sizeof(double));
//        N2-=N1; while (N2--) *(tr++)=*(tq++);
    }
    return r;
}
	

Polynomial Polynomial::operator-( Polynomial q )
{
	if (d->dim != q.d->dim)
	{
		printf("Poly::op- : Different dimension\n");
		getchar(); exit(-1);
	}

	Polynomial r(d->dim,mmax(d->deg,q.d->deg));
	unsigned N1=sz(), N2=q.sz(), Ni=mmin(N1,N2);
	double	*tr = r, *tp = (*this), *tq = q;
    while (Ni--) *(tr++) = *(tp++) - *(tq++);
    if (N1<N2)
    {
        N2-=N1; while (N2--) *(tr++)=-(*(tq++));
    }
	return r;
}

Polynomial Polynomial::operator-( void )
{
	unsigned Ni = sz();
    double *tp = (*this);

	if (!Ni || !tp) return *this;	// Take it like it is ...
	
	Polynomial r(d->dim,d->deg);
	double *tq = (r);
	while( Ni-- ) *(tq++) = -(*(tp++));
	return r;
}


Polynomial Polynomial::operator+=( Polynomial p )
{
	if (d->dim != p.d->dim)
	{
		printf("Poly::op+= : Different dimension\n");
		getchar(); exit(-1);
	}

    d->deg = mmax(d->deg,p.d->deg);	// New degree
    unsigned N1=sz(), N2=p.sz(), Ni=mmin(N1,N2);
    if (N1<N2) 
    {
        d->coeff=(double*)realloc(d->coeff,N2*sizeof(double));
        d->n=p.d->n;
    }
	double *tt = (*this),*tp = p;

    while (Ni--) *(tt++) += *(tp++);
    
    if (N1<N2) 
    {
        memcpy(tt,tp,(N2-N1)*sizeof(double));
//        N2-=N1; while (N2--) *(tt++)=*(tp++);
    }

	return *this;
}


Polynomial Polynomial::operator-=( Polynomial p )
{
	if (d->dim != p.d->dim)
	{
		printf( "Poly::op-= : Different dimension\n");
		getchar(); exit(-1);
	}

    d->deg = mmax(d->deg,p.d->deg);	// New degree
    unsigned N1=sz(), N2=p.sz(), Ni=mmin(N1,N2);
    if (N1<N2) 
    {
        d->coeff=(double*)realloc(d->coeff,N2*sizeof(double));
        d->n=p.d->n;
    }
	double *tt = (*this),*tp = p;

    while (Ni--) *(tt++) -= *(tp++);
    
    if (N1<N2) 
    {
        N2-=N1; while (N2--) *(tt++)=-(*(tp++));
    }

	return *this;
}


Polynomial Polynomial::operator*=( const double t )
{
    int i=sz();
	double *tp = (*this);

	while (i--) *(tp++) *=t;
	return *this;
}


Polynomial Polynomial::operator/=( const double t )
{
	if (t == 0)
	{
		printf( "Poly::op/= : Division by zero\n");
		getchar(); exit(-1);
	}

    int i=sz();
	double *tp = (*this);

	while (i--) *(tp++) /=t;
	return *this;
}

int Polynomial::equals( Polynomial q )
{
    if (d==q.d) return 1;
	if ( (d->deg != q.d->deg) || (d->dim != q.d->dim) ) return 0;

	unsigned N = sz();
	double	*tp = (*this),*tq = q;

	while (N--)
		if ( *(tp++) != *(tq++) ) return 0;

	return 1;
}

//ostream& Polynomial::PrintToStream( ostream& out ) const
void Polynomial::print()
{
	MultInd I( d->dim );
    double *tt = (*this);
	unsigned N = sz();
	bool IsFirst=true;
	
	if ( !N || !tt ) { printf("[Void polynomial]\n"); return; }

    if (*tt) { IsFirst=false; printf("%f", *tt); }
    tt++; ++I; 

	for (unsigned i = 1; i < N; i++,tt++,++I)
    {
		if (*tt != 0)
		{
			if (IsFirst)
			{
				if (queryFlag( NicePrint ))
				{
                    if (*tt<0) printf("-");
                    printf("%f x^",condorAbs(*tt)); I.print();
				}
				else
				{
                    printf("+%f x^",*tt); I.print();
				}
				IsFirst = false;
				continue;
			}
			if (queryFlag( NicePrint ))
			{
                if (*tt<0) printf("-"); else printf("+");
                printf("%f x^",condorAbs(*tt)); I.print();
			}
			else
			{
                printf("+%f x^",*tt); I.print();
			}

		}
	}
}

/*
double Polynomial::simpleEval(Vector P)
{
    unsigned i=coeff.sz(),j;
    double *cc=coeff,r=0, r0, *p=(double*)P;
    MultInd I(dim);
    while (i--)
    {
        r0=*(cc++); j=dim;
        while (j--) r0*=pow(p[j],I[j]);
        r+=r0; I++;
    }
    return r;
}
*/

double Polynomial::shiftedEval( Vector Point, double minusVal)
{
    double tmp1=d->coeff[0], tmp2;
    d->coeff[0]-=minusVal;
    tmp2=(*this)(Point);
    d->coeff[0]=tmp1;
    return tmp2;
}

// Evaluation oprator
// According to Pena, Sauer, "On the multivariate Horner scheme",
//   SIAM J. Numer. Anal., to appear

double Polynomial::operator()( Vector Point ) 
{
// I didn't notice any difference in precision:
//  return simpleEval(P);

  unsigned dim=d->dim, deg=d->deg;
  double r,r0;                             // no static here because of the 2 threads !
  double rbuf[100];	// That should suffice // no static here because of the 2 threads !
  double *rbufp = rbuf;
  unsigned lsize = 100;
  double *rptr;
  int i,j;

  if (Point==Vector::emptyVector) return *d->coeff;

  if ( dim != (unsigned)Point.sz() )
  {
    printf( "Polynomial::operator()( Vector& ) : Improper size\n");
    getchar(); exit(-1);
  }

  if ( !sz() )
  {
    if ( queryFlag( Warning ) )
    {
      printf( "Polynomial::operator()( Vector& ) : evaluating void polynomial\n");
    }
    return 0;
  }

  if ( dim > lsize )	// Someone must be crazy !!!
  {
    if ( queryFlag( Warning ) )
    {
	printf( "Polynomial::operator()( Vector& ) : Warning -> 100 variables\n");
    }
    if ((rbufp != rbuf) && rbufp) delete rbufp;

    lsize=dim;
    rbufp = (double*)malloc(lsize*sizeof(double));	// So be it ...

    if ( !rbufp )
    {
      printf( "Polynomial::operator()( Vector& ) : Cannot allocate <rbufp>\n");
      getchar(); exit( -1 );
    }
  }

  if (deg==0) return *d->coeff;

  // Initialize
  MultInd *mic=cacheMultInd.get( dim, deg );
  unsigned *nextI=mic->indexesOfCoefInLexOrder(),
           *lcI=mic->lastChanges();
  double *cc = (*this), *P=Point;
  unsigned nxt, lc;

  // Empty buffer (all registers = 0)
  memset(rbufp,0,dim*sizeof(double));

  r0=cc[*(nextI++)];
  i=sz()-1;
  while (i--)
  {
    nxt= *(nextI++);
    lc = *(lcI++);

    r=r0; rptr=rbufp+lc; j=dim-lc;
    while (j--) { r+=*rptr; *(rptr++)=0; }
    rbufp[lc]=P[lc]*r; 
    r0=cc[nxt];
  }
  r=r0; rptr=rbufp; i=(int)dim;
  while (i--) r+=*(rptr++);

  return r;
}

Polynomial Polynomial::derivate(int i)
{
    unsigned dim=d->dim, deg=d->deg;
    if (deg<1) return Polynomial(dim,0.0);

    Polynomial r(dim, deg-1);
	MultInd I( dim );
	MultInd J( dim );
    double *tS=(*this), *tD=r;
	unsigned j=sz(), k, *cc, sum,
             *allExpo=(unsigned*)I, *expo=allExpo+i, *firstOfJ=(unsigned*)J;

    while (j--) 
    { 
        if (*expo)
        {
            (*expo)--;

            sum=0; cc=allExpo; k=dim;
            while (k--) sum+=*(cc++);
            if (sum) k=choose( sum-1+dim, dim ); else k=0;
            J.resetCounter(); *firstOfJ=sum;
            while (!(J==I)) { k++; J++; }

            (*expo)++;
            tD[k]=(*tS) * (double)*expo;
        }
        tS++;
        I++;
    }
    return r;
}

void Polynomial::gradient(Vector P, Vector G)
{
    unsigned i=d->dim;
    G.setSize(i);
    double *r=G;
    if (P.equals(Vector::emptyVector))
    {
        memcpy(r,d->coeff+1,i*sizeof(double));
        return;
    }
    while (i--) r[i]=(derivate(i))(P);
}

void Polynomial::gradientHessian(Vector P, Vector G, Matrix H)
{
    unsigned dim=d->dim;
    G.setSize(dim);
    H.setSize(dim,dim);
    double *r=G, **h=H;
    unsigned i,j;

    if (d->deg==2)
    {        
        double *c=d->coeff+1;
        memcpy(r,c,dim*sizeof(double));
        c+=dim;
        for (i=0; i<dim; i++)
        {
            h[i][i]=2* *(c++);
            for (j=i+1; j<dim; j++)
                h[i][j]=h[j][i]=*(c++);
        }
        if (P.equals(Vector::emptyVector)) return;
        G+=H.multiply(P);
        return;
    }

    Polynomial *tmp=new Polynomial[dim], a;
    i=dim;
    while (i--)
    {
        tmp[i]=derivate(i);
        r[i]=(tmp[i])(P);
    }

    i=dim;
    while (i--)
    {
        j=i+1;
        while (j--)
        {
            a=tmp[i].derivate(j);
            h[i][j]=h[j][i]=a(P);
        }
    }

//    _CrtCheckMemory(); 

    delete []tmp;
}

void Polynomial::translate(Vector translation)
{
    if (d->deg>2)
    {
        printf("Translation only for polynomial of degree lower than 3.\n");
        getchar(); exit(255);
    }
    d->coeff[0]=(*this)(translation);
    if (d->deg==1) return;
    int dim=d->dim;
    Vector G(dim);
    Matrix H(dim,dim);
    gradientHessian(translation, G, H);
    memcpy(((double*)d->coeff)+1, (double*)G, dim*sizeof(double));
}

void Polynomial::save(char *name)
{
    FILE *f=fopen(name,"wb");
    fwrite(&d->dim, sizeof(int),1, f);
    fwrite(&d->deg, sizeof(int),1, f);
    fwrite(d->coeff, d->n*sizeof(double),1, f);
    fclose(f);
}

Polynomial::Polynomial(char *name)
{
    unsigned _dim,_deg;
    FILE *f=fopen(name,"rb");
    fread(&_dim, sizeof(int),1, f);
    fread(&_deg, sizeof(int),1, f);
    init(_dim,_deg);
    fread(d->coeff, d->n*sizeof(double),1, f);
    fclose(f);
}
