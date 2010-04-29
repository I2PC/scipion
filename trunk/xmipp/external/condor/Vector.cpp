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
#include <string.h> // for memmove: microsoft bug
#include "Vector.h"
#include "Matrix.h"
#include "tools.h"

Vector Vector::emptyVector;

void Vector::alloc(int _n, int _ext)
{
    d=(VectorData*)malloc(sizeof(VectorData));
    d->n=_n;
    d->extention=_ext;
    d->ref_count=1;

    if (_ext==0) { d->p=NULL; return; };

    d->p=(double*)malloc(_ext*sizeof(double)); 
    if (d->p==NULL) { printf("memory allocation error\n"); getchar(); exit(253); }
}

Vector::Vector(int n)
{
    alloc(n,n);
    zero();
};

Vector::Vector(int n, int ext)
{
    alloc(n,ext); 
    zero();
};

Vector::Vector(int n, double *dd, char _exte)
{
    alloc(n,n);
    if (dd) 
    {
        if (_exte) { d->externalData=1; d->p=dd; }
        else memcpy(d->p,dd,n*sizeof(double));
    }
    else zero();
}

Vector::Vector(Vector a,Vector b)
{
    int n=a.sz()+b.sz();
    alloc(n,n);
    memcpy(d->p       ,a,a.sz()*sizeof(double));
    memcpy(d->p+a.sz(),b,b.sz()*sizeof(double));
}

void Vector::setExternalData(int _n, double *dd)
{
    if ((d->extention==_n)||(!d->extention))
    {
        d->n=_n; d->extention=_n; d->externalData=1; d->p=dd; 
    } else
    {
        printf("do not use this function ('setExternalData'): it's too dangerous.\n");
        getchar(); exit(255);
    }
}


void Vector::zero(int i, int _n)
{
    if (_n==0) _n=d->n-i;
    if (d->p) memset(d->p+i,0,_n*sizeof(double));
}

void Vector::prepareExtend(int new_extention)
{
	if (d->extention<new_extention)
	{
		d->p=(double*)realloc(d->p,new_extention*sizeof(double));
        if (d->p==NULL) { printf("memory allocation error\n"); getchar(); exit(253); }

        // not really necessary (fill with zero's):
        memset(d->p+d->extention,0,(new_extention-d->extention)*sizeof(double));
    	d->extention=new_extention;
	};
};

void Vector::setSize(int _n)
{
    d->n=_n;
    if (_n==0) { if (d->p) free(d->p); d->p=NULL; d->extention=0; return; }
    prepareExtend(_n);
}

void Vector::extend()
{
    d->n++;
    if (d->n>d->extention) prepareExtend(d->extention+100);
}

void Vector::exactshape()
{
	if (d->extention!=d->n)
	{
		d->p=(double*)realloc(d->p,d->n*sizeof(double));
        if (d->p==NULL) { printf("memory allocation error\n"); getchar(); exit(253); }
		d->extention=d->n;
	};
};

int Vector::equals( Vector Q )
{
  if (Q.d==d) return 1;
  if (Q.d==emptyVector.d)
  {
      double *cP=d->p;
      int i=sz();
      while (i--) if (*(cP++)) return 0;
      return 1;
  }

  if (sz() != Q.sz()) return 0;

  double *cP = d->p, *cQ = Q.d->p;
  int i = sz();

  while( i-- )
  {
    if (*cP!=*cQ) return 0;
    cP++; cQ++;
  }

  return 1;
}

//ostream& Vector::PrintToStream( ostream& out ) const
void Vector::print()
{
    int N=sz();
	printf("[");
	if (!N || !d->p) { printf("]\n"); return; }

    double *up=d->p;
	while (--N) printf("%f,",*(up++));
	printf("%f]\n",*up);
}

Vector::~Vector()
{
    destroyCurrentBuffer();
};

void Vector::destroyCurrentBuffer()
{
    if (!d) return;
    (d->ref_count) --;
	if (d->ref_count==0)
    {
        if ((d->p)&&(!d->externalData)) free(d->p);
        free(d);
    };
}

Vector::Vector(const Vector &A)
{
    // shallow copy
    d=A.d;
	(d->ref_count)++ ;
}

Vector& Vector::operator=( const Vector& A )
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

Vector Vector::clone()
{
    // a deep copy
    Vector r(sz());
    r.copyFrom(*this);
    return r;
}

void Vector::copyFrom(Vector r, int _n)
{
    if (_n==0) _n=r.sz();
    setSize(_n);
    if (_n) memcpy(d->p,r.d->p,_n*sizeof(double));
}

double Vector::euclidianNorm()
{
    return ::euclidianNorm(sz(), d->p);
}

double Vector::L1Norm()
{
    if (sz()==0) return 0;
    double *x=d->p, sum=0;
    int ni=sz();
    while (ni--) sum+=condorAbs(*(x++));
    return sum;
}

double Vector::square()
{
    double *xp=d->p, sum=0;
    int ni=sz();
    while (ni--) sum+=sqr(*(xp++));
    return sum;
}

double Vector::euclidianDistance(Vector v)
{
    Vector t=(*this)-v;
    return ::euclidianNorm(sz(), t.d->p);
/*
    double *xp1=d->p, *xp2=v.d->p, sum=0;
    int ni=sz();
    while (ni--) sum+=sqr(*(xp1++)-*(xp2++));
    return sqrt(sum);
*/
}

double Vector::LnftyDistance(Vector v)
{
    double *xp1=d->p, *xp2=v.d->p, sum=-1.0;
    int ni=sz();
    while (ni--) sum=::mmax(sum, condorAbs(*(xp1++)-*(xp2++)));
    return sum;
}

double Vector::LnftyNorm()
{
    double *xp1=d->p, sum=-1.0;
    int ni=sz();
    while (ni--) sum=::mmax(sum, condorAbs(*(xp1++)));
    return sum;
}

double Vector::L1Distance(Vector v)
{
    if (sz()==0) return 0;
    double *xp1=d->p, *xp2=v.d->p, sum=0;
    int ni=sz();
    while (ni--) sum+=condorAbs(*(xp1++)-*(xp2++));
    return sum;
}

void Vector::multiply(double a)
{
    double *xp=d->p;
    int ni=sz();
    while (ni--) *(xp++)*=a;
}

void Vector::multiply(Vector R, double a)
{
    int ni=sz();
    R.setSize(ni);
    double *xs=d->p, *xd=R;
    while (ni--) *(xd++)=a * (*(xs++));
}


void Vector::transposeAndMultiply(Vector vR, Matrix M)
{
    if ((int)sz()!=M.nLine())
    {
        printf("error in V^t * M.\n"); getchar(); exit(254);
    }
    int n=sz(), szr=M.nColumn(), i;
    vR.setSize(szr);
    double sum, *dv=(*this), **dm=M, *dd=vR;

    while (szr--)
    {
        sum=0.0;
        i=n;
        while (i--) sum+=dv[i]*dm[i][szr];
        dd[szr]=sum;
    }
}

void Vector::diagonalizeAndMultiply(Matrix M)
{
    int i,j,nl=M.nLine(),nc=M.nColumn();
    if ((int)sz()!=nl)
    {
        printf("(matrix_diagonal * matrix) error");
        getchar(); exit(249);
    }
    double **p1=M,*p2=(*this);

    for (i=0; i<nl; i++) 
        for (j=0; j<nc; j++)
            p1[i][j]*=p2[i];

}

double Vector::scalarProduct(Vector v)
{
    double *xp1=d->p, *xp2=v.d->p, sum=0;
    int ni=sz();
    while (ni--) { sum+=*(xp1++) * *(xp2++); };
    return sum;
}

double Vector::mmin()
{
    if (sz()==0) return 0;
    double *xp=d->p, m=INF;
    int ni=sz();
    while (ni--) m=::mmin(m,*(xp++));
    return m;
}

double Vector::mmax()
{
    if (sz()==0) return 0;
    double *xp=d->p, m=-INF;
    int ni=sz();
    while (ni--) m=::mmax(m,*(xp++));
    return m;
}

bool Vector::isNull()
{
    double *xp=d->p;
    int ni=sz();
    while (ni--) if (*(xp++)!=0) return false;
    return true;
}

Vector Vector::operator-( Vector v)
{
    int ni=sz();
    Vector r(sz());
    double *xp1=r.d->p, *xp2=v.d->p, *xp3=d->p;
    while (ni--) 
        *(xp1++)+=*(xp3++)-*(xp2++);
    return r;
}

void Vector::oneByOneMutiply(Vector rescaling)
{
    int i=::mmin(sz(),rescaling.sz());
    double *xb=(*this), *r=rescaling;
    while (i--) xb[i]*=r[i];
}

void Vector::oneByOneInvert()
{
    int i=sz();
    double *xb=(*this);
    while (i--) xb[i]=1/xb[i];
}

Vector Vector::operator+( Vector v)
{
    int ni=sz();
    Vector r(sz());
    double *xp1=r.d->p, *xp2=v.d->p, *xp3=d->p;
    while (ni--) 
        *(xp1++)+=*(xp3++)+*(xp2++);
    return r;
}

Vector Vector::operator-=( Vector v)
{
    int ni=sz();
    double *xp1=d->p, *xp2=v.d->p;
    while (ni--) *(xp1++)-=*(xp2++);
    return *this;
}

Vector Vector::operator+=( Vector v)
{
    int ni=sz();
    double *xp1=d->p, *xp2=v.d->p;
    while (ni--) *(xp1++)+=*(xp2++);
    return *this;
}

void Vector ::addInPlace(double a, Vector v)
{
    int ni=sz();
    double *xp1=d->p, *xp2=v.d->p;
    if (a==1.0)  while (ni--) *(xp1++)+=     *(xp2++);
    else         while (ni--) *(xp1++)+=a * (*(xp2++));
}

void Vector::addInPlace(double a, int i, Matrix m)
{
    int ni=sz();
    double *xp1=d->p, *xp2=((double**)m)[i];
    while (ni--) *(xp1++)+=a * (*(xp2++));
}

Vector::Vector(char *filename)
{
    unsigned _n;
    FILE *f=fopen(filename,"rb");
    fread(&_n, sizeof(int),1, f);
    alloc(_n,_n);
    fread(d->p, d->n*sizeof(double),1, f);
    fclose(f);
}

void Vector::save(char *filename, char ascii)
{
    FILE *f;
    if (ascii) f=fopen(filename,"w"); else f=fopen(filename,"wb");
    if (f==NULL)
    {
		printf("Cannot write to '%s'\n",filename);
		exit(255);
    }
    save(f,ascii);
    fclose(f);
}

void Vector::save(FILE *f, char ascii)
{
    char *header="CONDORVBv1.0";
    if (ascii)
    {
        unsigned i;
        double *pp=d->p;
        if (ascii!=2) fprintf(f,"CONDORVAv1.0\n");
        if (sz())
        {
            for (i=0; i<sz()-1; i++) fprintf(f,"%.16e\t",pp[i]);
            fprintf(f,"%.16e\n",pp[i]);
        }
        return;
    }
    fwrite(header,sizeof(char),13,f);
    fwrite(&d->n, sizeof(int),1, f);
    fwrite(d->p, d->n*sizeof(double),1, f);
}


void Vector::setPart(int i, Vector v, int n, int ii)
{
	if (n==0) n=v.sz()-ii;
    n=::mmin((int)n,(int)sz()-i);
    memcpy(d->p+i, ((double*)v)+ii, n*sizeof(double));
}

void Vector::set(double dd)
{
    double *p=(*this);
    if (!p) return;
    int n=sz();
    while (n--) *(p++)=dd;
}

void Vector::shift(int s)
{
    int n=sz();
    if (!n) return;
    double *ps=(*this), *pd=ps; // pointer source / destination
    if (s==0) return;
    if (s>0) { n-=s; pd+=s; } 
    else { n+=s; ps+=s; }
    memmove(pd,ps,n*sizeof(double));
}

void Vector::permutIn(Vector vR, VectorInt viP)
{
    int i,n=sz(), *ii=viP;
    if (!n) return;
    if (n!=viP.sz())
    {
        printf("error in permutation IN: sizes don't agree.\n"); getchar(); exit(255);
    }
    vR.setSize(n);
    double *ps=(*this), *pd=vR; // pointer source / destination
    for (i=0; i<n; i++)
//        *(pd++)=ps[ii[i]];
        pd[ii[i]]=*(ps++);
}

void Vector::permutOut(Vector vR, VectorInt viP)
{
    int i,n=sz(), *ii=viP;
    if (!n) return;
    if (n!=viP.sz())
    {
        printf("error in permutation IN: sizes don't agree.\n"); getchar(); exit(255);
    }
    vR.setSize(n);
    double *ps=(*this), *pd=vR; // pointer source / destination
    for (i=0; i<n; i++)
//        pd[ii[i]]=*(ps++);
        *(pd++)=ps[ii[i]];
}

#define EOL1 13
#define EOL2 10

char isANumber(char c)
{
 return (((c>='0')&&(c<='9'))||
         (c=='.')||
         (c=='e')||
         (c=='E')||
         (c=='+')||
         (c=='-'));
}
Vector::Vector(char *line, int gn)
{
	char *tline=line,*oldtline=NULL;

    if (gn==0) 
    {
	    while ((*tline!=EOL1)&&(*tline!=EOL2))
	    {
		    while ((*tline==' ')|| 
			    (*tline=='\t'))tline++;
		    if ((*tline==EOL1)||(*tline==EOL2)||(*tline==0)) break;
		    while (isANumber(*tline)) tline++;
		    gn++;
            if (oldtline==tline)
            {
                alloc(gn-1,gn-1);
                return;
                //tline[10]=0;
                //printf("Error in simulation output file. The full line is:\n"
                //       "%s\n"
                //       "There is an error here:\n"
                //       "%s\n",line,tline);
            }
            oldtline=tline;
	    };
    };

	if (gn==0) { alloc(0,0); return; };
    alloc(gn,gn);
    getFromLine(line);
};

char *Vector::getFromLine(char *line)
{
    double *dp=d->p;
    int n=sz(),k;
	char *tline=line, *oldtline; 
	for (k=0; k<n; k++) 
    {
		while ((*tline==' ')|| 
			    (*tline=='\t'))tline++;
		if ((*tline==EOL1)||(*tline==EOL2)) 
        {
            setSize(k);
	        return tline;
        }
        oldtline=tline;
        while(isANumber(*tline)) tline++;
        if (!isANumber(*(tline-1)))
        {
            setSize(k);
            return tline;
            //tline[10]=0;
            //printf( "Error in simulation output file. The full line is:\n"
            //        "%s\n"
            //        "There is an error here:\n"
            //        "%s\n",line,tline);
        }
        if (oldtline==tline)
        {
            setSize(k);
	        return tline;
        };
        if (*tline) { *tline='\0'; tline++; }
        dp[k]=atof(oldtline);
    }
    return tline;
}

void Vector::appendToMatrixFile(char *saveFileName, char **names)
{
    FILE *f=fopen(saveFileName,"rb+");
    if (f==NULL)
    {
        Matrix t(*this);
        if (names) t.setColNames(names);
        t.save(saveFileName,0);
        return;
    }
    int nc=sz();
    fseek(f,0,SEEK_END);
    unsigned l=ftell(f);
    if (l==0)
    {
        Matrix t(*this);
        if (names) t.setColNames(names);
        t.save(saveFileName,0);
        return;
    }
    fwrite(d->p,sizeof(double)*nc,1,f);
    unsigned nlfile;
    fseek(f,13,SEEK_SET);
    fread(&nlfile,sizeof(int),1,f); 
    nlfile++;
    fseek(f,13,SEEK_SET);
    fwrite(&nlfile, sizeof(unsigned), 1, f);
    fclose(f);
}
