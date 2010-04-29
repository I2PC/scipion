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

#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include "tools.h"

Matrix Matrix::emptyMatrix;

void Matrix::init(int _nLine, int _nColumn, int _extLine, int _extColumn,MatrixData* md)
{
    if (md==NULL)
    {
        d=(MatrixData*)malloc(sizeof(MatrixData));
        d->columnNames=NULL;
        d->ref_count=1;
    }
    else d=md;
    d->nLine=_nLine;
    d->nColumn=_nColumn;
    d->extLine=_extLine;
    d->extColumn=_extColumn;

    if ((_extLine>0)&&(_extColumn>0))
    {
        double **t,*t2;
        t=d->p=(double**)malloc(_extLine*sizeof(double*));
        t2=(double*)malloc(_extColumn*_extLine*sizeof(double));
        while(_extLine--)
        {
            *(t++)=t2;
            t2+=_extColumn;
        }
    }
    else d->p=NULL;
}

Matrix::Matrix(int _nLine,int _nColumn)
{
    init(_nLine,_nColumn,_nLine, _nColumn);
};

void Matrix::diagonal(double dd)
{
    zero();
    double *p=*d->p;
    int n=nLine(), i=d->extColumn+1;
    while (n--)
    {
        *p=dd;
        p+=i;
    }

}

Matrix::Matrix(Vector a)
{
    int nl=1, nc=a.sz();
    init(nl,nc,nl,nc);
    memcpy(*d->p,a,a.sz()*sizeof(double));
}

Matrix::Matrix(Vector a, Vector b)  // a * b^T
{
    int nl=a.sz(), nc=b.sz(), i,j;
    double *pa=a, *pb=b;
    init(nl,nc,nl,nc);
    double **pp=d->p;
    for (i=0; i<nl; i++)
        for (j=0; j<nc; j++)
            pp[i][j]=pa[i]*pb[j];
}

Matrix::Matrix(int _nLine,int _nColumn,int _extLine,int _extColumn)
{
    init(_nLine,_nColumn,_extLine,_extColumn);
};

Matrix::Matrix(const char *filename, char ascii)
{
    int _nLine,_nColumn;
    char c[13];
    FILE *f=fopen(filename,"rb");
    if (f==NULL)
    {
        printf("file not found.\n");
        exit(255);
    }
    if (ascii)
    {
        char line[30000];
        char *r=fgets(line,30000,f);
        if (r==NULL)
        {
            init(0,0,0,0);
            return;
        }
        if (line[7]!='A')
        {
            printf("not a ASCII matrix.\n");
            exit(255);
        }
        fgets(line,30000,f);
        _nLine=atol(line);
        fgets(line,30000,f);
        _nColumn=atol(line);
        init(_nLine,_nColumn,_nLine,_nColumn);
        fgets(line,30000,f);
        setColNames(getNameTable(line,&_nColumn));
        Vector tt(_nColumn);
        int i;
        for (i=0; i<_nLine; i++)
        {
            fgets(line,30000,f);
            tt.getFromLine(line);
            setLine(i,tt);
        }
        return;
    }
    if (fread(c, 13, sizeof(char), f)==0)
    {
        init(0,0,0,0);
        return;
    }
    if (c[7]!='B')
    {
        printf("not a binary matrix.\n");
        exit(255);
    }
    fread(&_nLine, sizeof(unsigned), 1, f);
    fread(&_nColumn, sizeof(unsigned), 1, f);
    init(_nLine,_nColumn,_nLine,_nColumn);
    int i,j=0;
    fread(&i, sizeof(int), 1, f);
    if (i)
    {
        char **names=(char**)malloc(_nColumn*sizeof(char**)),
               *n=(char*)malloc(i);
        fread(n,i,1,f);
        for (j=0; j<_nColumn-1; j++)
        {
            names[j]=n;
            while (*(n++));
        }
        names[j]=n;
        setColNames(names);
        free(*names);
        free(names);
    }
    if (d->nColumn*d->nLine)
        fread(*d->p,sizeof(double)*d->nColumn*d->nLine,1,f);
    fclose(f);
}

void Matrix::extendLine()
{
    d->nLine++;
    if (d->nLine>d->extLine) setExtSize(d->nLine+9,d->extColumn);
}

void Matrix::setNLine(int _nLine)
{
    d->nLine=_nLine;
    if (_nLine>d->extLine) setExtSize(_nLine,d->extColumn);
}

void Matrix::extendColumn()
{
    d->nColumn++;
    if (d->nColumn>d->extColumn) setExtSize(d->extLine,d->nColumn+9);
}

void Matrix::setNColumn(int _nColumn)
{
    d->nColumn=_nColumn;
    if (_nColumn>d->extColumn) setExtSize(d->extLine,_nColumn);
}

void Matrix::setSize(int _nLine,int _nColumn)
{
    d->nLine=_nLine;
    d->nColumn=_nColumn;
    if ((_nLine>d->extLine)||(_nColumn>d->extColumn)) setExtSize(_nLine,_nColumn);
}

void Matrix::setExtSize(int _extLine, int _extColumn)
{
    int ec=d->extColumn;
    if ((ec==0)||(d->extLine==0))
    {
        init(d->nLine,d->nColumn,_extLine,_extColumn,d);
        return;
    }
    if (_extColumn>ec)
    {
        int nc=d->nColumn,i;
        double *tmp,*tmp2,**tmp3=d->p,*oldBuffer=*tmp3;

        if (d->extLine<_extLine)
            tmp3=d->p=(double**)realloc(tmp3,_extLine*sizeof(double*));
        else _extLine=d->extLine;

        tmp2=tmp=(double *)malloc(_extLine*_extColumn*sizeof(double));
        if (tmp==NULL)
        {
            printf("memory allocation error");
            getchar();
            exit(255);
        }

        i=_extLine;
        while (i--)
        {
            *(tmp3++)=tmp2;
            tmp2+=_extColumn;
        };

        if ((nc)&&(d->nLine))
        {
            tmp2=oldBuffer;
            i=d->nLine;
            nc*=sizeof(double);
            while(i--)
            {
                memmove(tmp,tmp2,nc);
                tmp+=_extColumn;
                tmp2+=ec;
            };
            free(oldBuffer);
        };
        d->extLine=_extLine;
        d->extColumn=_extColumn;
        return;
    }
    if (_extLine>d->extLine)
    {
        int i;
        double *tmp,**tmp3;
        tmp=(double *)realloc(*d->p,_extLine*ec*sizeof(double));
        if (tmp==NULL)
        {
            printf("memory allocation error");
            getchar();
            exit(255);
        }
        free(d->p);
        tmp3=d->p=(double **)malloc(_extLine*sizeof(double*));
        i=_extLine;
        while (i--)
        {
            *(tmp3++)=tmp;
            tmp+=ec;
        };
        d->extLine=_extLine;
    }
}

void Matrix::save(char *filename,char ascii)
{
    FILE *f;
    if (ascii)
    {
        f=fopen(filename,"w");
        if (f==NULL)
        {
            printf("cannot save ascii Matrix into file '%s'.\n",filename);
            exit(255);
        }
    }
    else
    {
        f=fopen(filename,"wb");
        if (f==NULL)
        {
            printf("cannot save binary Matrix into file '%s'.\n",filename);
            exit(255);
        }
    }
    save(f,ascii);
    fclose(f);
}

void Matrix::save(FILE *f,char ascii)
{
    char *cc="CONDORMBv1.0";
    double **p=(d->p);
    int i,j;
    if (ascii)
    {
        if (ascii<2) fprintf(f,"CONDORMAv1.0\n%i\n%i\n",d->nLine,d->nColumn);
        if (ascii<3)
        {
            if (d->columnNames)
            {
                for (i=0; i<d->nColumn-1; i++)
                    fprintf(f,"%s\t",d->columnNames[i]);
                fprintf(f,"%s\n",d->columnNames[i]);
            }
            else fprintf(f,"null\n");
        }
        for (i=0; i<d->nLine; i++)
        {
            for (j=0; j<d->nColumn-1; j++)
                fprintf(f,"%.16e\t",p[i][j]);
            fprintf(f,"%.16e\n",p[i][j]);
        }
    }
    else
    {
        fwrite(cc, sizeof(char), 13, f);
        fwrite(&d->nLine, sizeof(unsigned), 1, f);
        fwrite(&d->nColumn, sizeof(unsigned), 1, f);
        if (d->columnNames)
        {
            j=0;
            for (i=0; i<d->nColumn; i++) j+=(int)strlen(d->columnNames[i])+1;
            fwrite(&j, sizeof(int), 1, f);
            for (i=0; i<d->nColumn; i++)
                fwrite(d->columnNames[i],strlen(d->columnNames[i])+1,1,f);
        }
        else
        {
            j=0;
            fwrite(&j, sizeof(int), 1, f);
        }
        for (i=0; i<d->nLine; i++)
            fwrite(p[i],sizeof(double)*d->nColumn,1,f);
    };
}

void Matrix::updateSave(char *saveFileName)
{
    FILE *f=fopen(saveFileName,"rb+");
    if (f==NULL)
    {
        save(saveFileName,0);
        return;
    }
    fseek(f,0,SEEK_END);
    long l=ftell(f);
    if (l==0)
    {
        save(saveFileName,0);
        return;
    }
    int nc=d->nColumn, nlfile, nl=d->nLine, i;
    fseek(f,13,SEEK_SET);
    fread(&nlfile,sizeof(int),1,f);
    fseek(f,13,SEEK_SET);
    fwrite(&d->nLine, sizeof(unsigned), 1, f);
    fseek(f,0,SEEK_END);
    double **p=d->p;
    for (i=nlfile; i<nl; i++)
        fwrite(p[i],sizeof(double)*nc,1,f);
    fclose(f);
}

void Matrix::exactshape()
{
    int i, nc=d->nColumn, ec=d->extColumn, nl=d->nLine, el=d->extLine;
    double *tmp,*tmp2,**tmp3;

    if ((nc==ec)&&(nl==el)) return;

    if (nc!=ec)
    {
        i=nl;
        tmp=tmp2=*d->p;
        while(i--)
        {
            memmove(tmp,tmp2,nc*sizeof(double));
            tmp+=nc;
            tmp2+=ec;
        };
    }

    tmp=(double *)realloc(*d->p,nl*nc*sizeof(double));
    if (tmp==NULL)
    {
        printf("memory allocation error");
        getchar();
        exit(255);
    }
    if (tmp!=*d->p)
    {
        tmp3=d->p=(double **)realloc(d->p,nl*sizeof(double*));
        i=nl;
        while (i--)
        {
            *(tmp3++)=tmp;
            tmp+=nc;
        };
    }
    else d->p=(double **)realloc(d->p,nl*sizeof(double*));

    d->extLine=nl;
    d->extColumn=nc;
};


void Matrix::print()
{
    double **p=d->p;
    int i,j;

    printf("[");
    for (i=0; i<d->nLine; i++)
    {
        for (j=0; j<d->nColumn; j++)
            if (p[i][j]>=0.0) printf(" %2.3f ",p[i][j]);
            else printf("%2.3f ",p[i][j]);
        if (i==d->nLine-1) printf("]\n");
        else printf(";\n");
    }
    fflush(0);
}

Matrix::~Matrix()
{
    destroyCurrentBuffer();
};

void Matrix::destroyCurrentBuffer()
{
    if (!d) return;
    (d->ref_count) --;
    if (d->ref_count==0)
    {
        if (d->columnNames)
        {
            free(*d->columnNames);
            free(d->columnNames);
        }
        if (d->p)
        {
            free(*d->p);
            free(d->p);
        }
        free(d);
    }
}

Matrix& Matrix::operator=( const Matrix& A )
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

Matrix::Matrix(const Matrix &A)
{
    // shallow copy
    d=A.d;
    (d->ref_count)++ ;
}

Matrix Matrix::clone()
{
    // a deep copy
    Matrix m(nLine(),nColumn());
    m.copyFrom(*this);
    return m;
}

void Matrix::copyFrom(Matrix m)
{
    int nl=m.nLine(),nc=m.nColumn(), ec=m.d->extColumn;
    if ((nl!=nLine())||(nc!=nColumn()))
    {
        printf("Matrix: copyFrom: size do not agree");
        getchar();
        exit(254);
    }
    if (ec==nc)
    {
        memcpy(*d->p,*m.d->p,nc*nl*sizeof(double));
        return;
    }
    double *pD=*d->p,*pS=*m.d->p;
    while(nl--)
    {
        memcpy(pD,pS,nc*sizeof(double));
        pD+=nc;
        pS+=ec;
    };
}

void Matrix::transposeInPlace()
{
    int nl=nLine(),nc=nColumn(),i,j;
    if (nl==nc)
    {
        double **p=(*this),t;
        for (i=0; i<nl; i++)
            for (j=0; j<i; j++)
            {
                t=p[i][j];
                p[i][j]=p[j][i];
                p[j][i]=t;
            }
        return;
    }
    Matrix temp=clone();
    setSize(nc,nl);
    double **sp=temp, **dp=(*this);
    i=nl;
    while (i--)
    {
        j=nc;
        while (j--) dp[j][i]=sp[i][j];
    }
}

void Matrix::transpose(Matrix temp)
{
    int nl=nLine(),nc=nColumn(),i,j;
    temp.setSize(nc,nl);
    double **sp=temp, **dp=(*this);
    i=nl;
    while (i--)
    {
        j=nc;
        while (j--) sp[j][i]=dp[i][j];
    }
}

Matrix Matrix::transpose()
{
    Matrix temp(nColumn(),nLine());
    transpose(temp);
    return temp;
}

//Matrix Matrix::deepCopy()
//{
//    Matrix cop(this); // contructor of class matrix
//    return cop;    // copy of class Matrix in return Variable
//                   // destruction of instance cop.
//};

void Matrix::zero()
{
    memset(*d->p,0,nLine()*d->extColumn*sizeof(double));
};

void Matrix::multiplyByDiagonalMatrix(Vector v)
{
    int i,j,nc=nColumn(),nl=nLine();
    if ((int)v.sz()!=nc)
    {
        printf("(matrix * matrix_diagonal) error");
        getchar();
        exit(249);
    }
    double **p1=(*this),*p2=v;

    for (i=0; i<nl; i++)
        for (j=0; j<nc; j++)
            p1[i][j]*=p2[j];
}

void Matrix::multiply(Matrix R, Matrix Bplus)
{
    if (Bplus.nLine()!=nColumn())
    {
        printf("(matrix * matrix) error");
        getchar();
        exit(249);
    }
    int i,j,k, nl=nLine(), nc=Bplus.nColumn(), n=nColumn();
    R.setSize(nl,nc);
    double sum,**p1=(*this),**p2=Bplus,**pr=R;

    for (i=0; i<nl; i++)
        for (j=0; j<nc; j++)
        {
            sum=0;
            for (k=0; k<n; k++) sum+=p1[i][k]*p2[k][j];
            pr[i][j]=sum;
        }
}

void Matrix::transposeAndMultiply(Matrix R, Matrix Bplus)
{
    if (Bplus.nLine()!=nLine())
    {
        printf("(matrix^t * matrix) error");
        getchar();
        exit(249);
    }
    int i,j,k, nl=nColumn(), nc=Bplus.nColumn(), n=nLine();
    R.setSize(nl,nc);
    double sum,**p1=(*this),**p2=Bplus,**pr=R;

    for (i=0; i<nl; i++)
        for (j=0; j<nc; j++)
        {
            sum=0;
            for (k=0; k<n; k++) sum+=p1[k][i]*p2[k][j];
            pr[i][j]=sum;
        }
}

void Matrix::multiplyByTranspose(Matrix R, Matrix Bplus)
{
    if (Bplus.nColumn()!=nColumn())
    {
        printf("(matrix * matrix^t) error");
        getchar();
        exit(249);
    }
    int i,j,k, nl=nLine(), nc=Bplus.nLine(), n=nColumn();
    R.setSize(nl,nc);
    double sum,**p1=(*this),**p2=Bplus,**pr=R;

    for (i=0; i<nl; i++)
        for (j=0; j<nc; j++)
        {
            sum=0;
            for (k=0; k<n; k++) sum+=p1[i][k]*p2[j][k];
            pr[i][j]=sum;
        }
}

Matrix Matrix ::multiply(Matrix Bplus)
{
    Matrix R(nLine(),Bplus.nColumn());
    multiply(R,Bplus);
    return R;
}

void Matrix::multiplyInPlace(double dd)
{
    int i,j, nl=nLine(), nc=nColumn();
    double **p1=(*this);

    for (i=0; i<nl; i++)
        for (j=0; j<nc; j++)
            p1[i][j]*=dd;
}

void Matrix::multiply(Vector rv, Vector v)
{
    int i,j, nl=nLine(), nc=nColumn();
    rv.setSize(nl);
    if (nc!=(int)v.sz())
    {
        printf("matrix multiply error");
        getchar();
        exit(250);
    };
    double **p=(*this), *x=v, *r=rv, sum;

    for (i=0; i<nl; i++)
    {
        sum=0;
        j=nc;
        while (j--) sum+=p[i][j]*x[j];
        r[i]=sum;
    }
}

void Matrix::transposeAndMultiply(Vector rv, Vector v)
{
    int i,j, nc=nLine(), nl=nColumn();
    rv.setSize(nl);
    if (nc!=(int)v.sz())
    {
        printf("matrix multiply error");
        getchar();
        exit(250);
    };
    double **p=(*this), *x=v, *r=rv, sum;

    for (i=0; i<nl; i++)
    {
        sum=0;
        j=nc;
        while (j--) sum+=p[j][i]*x[j];
        r[i]=sum;
    }
}

Vector Matrix::multiply(Vector v)
{
    Vector r(nLine());
    multiply(r,v);
    return r;
}

double Matrix::scalarProduct(int nl, Vector v)
{
    double *x1=v, *x2=d->p[nl], sum=0;
    int n=v.sz();
    while (n--)
    {
        sum+=*(x1++) * *(x2++);
    };
    return sum;
}

void Matrix::addInPlace(Matrix B)
{
    if ((B.nLine()!=nLine())||
        (B.nColumn()!=nColumn()))
    {
        printf("matrix addition error");
        getchar();
        exit(250);
    };

    int i,j, nl=nLine(), nc=nColumn();
    double **p1=(*this),**p2=B;

    for (i=0; i<nl; i++)
        for (j=0; j<nc; j++)
            p1[i][j]+=p2[i][j];
}

void Matrix::addMultiplyInPlace(double d, Matrix B)
{
    if ((B.nLine()!=nLine())||
        (B.nColumn()!=nColumn()))
    {
        printf("matrix addition error");
        getchar();
        exit(250);
    };

    int i,j, nl=nLine(), nc=nColumn();
    double **p1=(*this),**p2=B;

    for (i=0; i<nl; i++)
        for (j=0; j<nc; j++)
            p1[i][j]+=d*p2[i][j];
}

//inline double sqr(double a){return a*a;};

#ifndef NOMATRIXTRIANGLE
MatrixTriangle MatrixTriangle::emptyMatrixTriangle(0);

Matrix::Matrix(MatrixTriangle A, char bTranspose)
{
    int n=A.nLine(),i,j;
    init(n,n,n,n);
    double **pD=(*this), **pS=A;

    if (bTranspose)
    {
        for (i=0; i<n; i++)
            for (j=0; j<n; j++)
                if (j>=i) pD[i][j]=pS[j][i];
                else pD[i][j]=0;
    }
    else
    {
        for (i=0; i<n; i++)
            for (j=0; j<n; j++)
                if (j<=i) pD[i][j]=pS[i][j];
                else pD[i][j]=0;
    }
}

bool Matrix::cholesky(MatrixTriangle matL, double lambda, double *lambdaCorrection)
// factorize (*this)+lambda.I into L.L^t
{
    double s,s2;
    int i,j,k,n=nLine();
    matL.setSize(n);

    double **A=(*this), **L_=matL;
    if (lambdaCorrection) *lambdaCorrection=0;

    for (i=0; i<n; i++)
    {
        s2=A[i][i]+lambda;
        k=i;
        while ( k-- ) s2-=sqr(L_[i][k]);
        if (s2<=0)
        {
            if (lambdaCorrection)
            {
                // lambdaCorrection
                n=i+1;
                Vector X(n); // zero everywhere
                double *x=X, sum;
                x[i]=1.0;
                while(i--)
                {
                    sum=0.0;
                    for (k=i+1; k<n; k++) sum-=L_[k][i]*x[k];
                    x[i]=sum/L_[i][i];
                }
                *lambdaCorrection=-s2/X.euclidianNorm();
            }
            return false;
        }
        L_[i][i] = s2 = sqrt(s2);

        for (j=i+1; j<n; j++)
        {
            s=A[i][j];
            k=i;
            while (k--) s-=L_[j][k]*L_[i][k];
            L_[j][i]=s/s2;
        }
    }
    return true;
}

void Matrix::choleskySolveInPlace(Vector b)
{
    MatrixTriangle M(nLine());
    if (!cholesky(M))
    {
        b.setSize(0); // no cholesky decomposition => return emptyVector
        return;
    }
    M.solveInPlace(b);
    M.solveTransposInPlace(b);
}

void Matrix::QR(Matrix Q, MatrixTriangle Rt, VectorInt vPermutation)
{
//      QR factorization of the transpose of (*this)
//      beware!! :
//          1. (*this) is destroyed during the process.
//          2. Rt contains the tranpose of R (get an easy manipulation matrix using:
//                  Matrix R(Rt,1); ).
//          3. use of permutation IS tested
//
//
//      subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
//      integer m,n,lda,lipvt
//      integer ipvt(lipvt)
//      logical pivot
//
//
//      double precision a(lda,n),rdiag(n),acnorm(n),wa(n)
//c     **********
//c
//c     subroutine qrfac
//c
//c     this subroutine uses householder transformations with column
//c     pivoting (optional) to compute a qr factorization of the
//c     m by n matrix a. that is, qrfac determines an orthogonal
//c     matrix q, a permutation matrix p, and an upper trapezoidal
//c     matrix r with diagonal elements of nonincreasing magnitude,
//c     such that a*p = q*r. the householder transformation for
//c     column k, k = 1,2,...,min(m,n), is of the form
//c
//c                           t
//c           i - (1/u(k))*u*u
//c
//c     where u has zeros in the first k-1 positions. the form of
//c     this transformation and the method of pivoting first
//c     appeared in the corresponding linpack subroutine.
//c
//c     the subroutine statement is
//c
//c       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
//c
//c     where
//c
//c       m is a positive integer input variable set to the number
//c         of rows of a.
//c
//c       n is a positive integer input variable set to the number
//c         of columns of a.
//c
//c       a is an m by n array. on input a contains the matrix for

//c         which the qr factorization is to be computed. on output
//c         the strict upper trapezoidal part of a contains the strict
//c         upper trapezoidal part of r, and the lower trapezoidal
//c         part of a contains a factored form of q (the non-trivial
//c         elements of the u vectors described above).
//c
//c       lda is a positive integer input variable not less than m
//c         which specifies the leading dimension of the array a.
//c
//c       pivot is a logical input variable. if pivot is set true,
//c         then column pivoting is enforced. if pivot is set false,
//c         then no column pivoting is done.
//c
//c       ipvt is an integer output array of length lipvt. ipvt
//c         defines the permutation matrix p such that a*p = q*r.
//c         column j of p is column ipvt(j) of the identity matrix.
//c         if pivot is false, ipvt is not referenced.
//c
//c       lipvt is a positive integer input variable. if pivot is false,
//c         then lipvt may be as small as 1. if pivot is true, then
//c         lipvt must be at least n.
//c
//c       rdiag is an output array of length n which contains the
//c         diagonal elements of r.
//c
//c       wa is a work array of length n. if pivot is false, then wa
//c         can coincide with rdiag.
    char pivot=!(vPermutation==VectorInt::emptyVectorInt);
    int i,j,k,kmax,minmn;
    double ajnorm,sum,temp;
    // data one,p05,zero /1.0d0,5.0d-2,0.0d0/

    const double epsmch = 1e-20; // machine precision
    int nc=nColumn(), nl=nLine();

    if (nl>nc)
    {
        printf("QR factorisation of A^t is currently not possible when number of lines is greater than number of columns.\n");
        getchar();
        exit(255);
    }

    Vector vWA(nl), vRDiag;
    int *ipvt;
    double *wa=vWA, *rdiag, **a=*this;

    if (pivot)
    {
        vPermutation.setSize(nl);
        ipvt=vPermutation;
        vRDiag.setSize(nl);
        rdiag=vRDiag;
    }
    else rdiag=wa;

//c
//c     compute the initial line norms and initialize several arrays.
//c
    for (j=0; j<nl; j++)
    {
        rdiag[j]=wa[j]=euclidianNorm(j);
        if (pivot) ipvt[j]=j;
    }
//c
//c     reduce a to r with householder transformations.
//c
    minmn=mmin(nl,nc);
    for (j=0; j<minmn; j++)
    {
        if (pivot)
        {
//c
//c        bring the line of largest norm into the pivot position.
//c
            kmax=j;
            for (k=j+1; k<nl; k++)
                if (rdiag[k]>rdiag[kmax]) kmax=k;

            if (kmax!=j)
            {
                for (i=0; i<nc; i++)
                {
                    temp = a[j][i];
                    a[j][i] = a[kmax][i];
                    a[kmax][i] = temp;
                }
                rdiag[kmax] = rdiag[j];
                wa[kmax] = wa[j];
                k = ipvt[j];
                ipvt[j] = ipvt[kmax];
                ipvt[kmax] = k;
            }
        }
//c
//c        compute the householder transformation to reduce the
//c        j-th line of a to a multiple of the j-th unit vector.
//c

//        ajnorm = enorm(nl-j+1,a(j,j))
        ajnorm=::euclidianNorm(nc-j, &a[j][j]);

        if (ajnorm==0.0)
        {
            rdiag[j]=0.0;
            continue;
        }
        if (a[j][j]<0.0) ajnorm = -ajnorm;
        for (i=j; i<nc; i++) a[j][i]=a[j][i]/ajnorm;
        a[j][j]+=1.0;

//c
//c        apply the transformation to the remaining lines
//c        and update the norms.
//c
        if (j>=nc)
        {
            rdiag[j] = -ajnorm;
            continue;
        }

        for (k = j+1; k<nl; k++)
        {
            sum=0.0;
            for (i=j; i<nc; i++) sum=sum+a[j][i]*a[k][i];

            temp = sum/a[j][j];
            for (i=j; i<nc; i++) a[k][i]=a[k][i]-temp*a[j][i];

            if ((!pivot)||(rdiag[k]==0.0)) continue;

            temp = a[k][j]/rdiag[k];
            rdiag[k] *= sqrt(mmax(0.0,1.0-temp*temp));

            if (0.05*sqr(rdiag[k]/wa[k])> epsmch) continue;

            //rdiag(k) = enorm(nl-j,a(jp1,k))
            rdiag[k]=::euclidianNorm(nc-j, &a[k][j+1]);
            wa[k] = rdiag[k];
        }
        rdiag[j] = -ajnorm;
    }
//c
//c     last card of subroutine qrfac.
//c
    if (!(Rt==MatrixTriangle::emptyMatrixTriangle))
    {
        Rt.setSize(minmn);
        double **r=Rt;
        for (i=0; i<minmn; i++)
        {
            r[i][i]=rdiag[i];
            for (j=i+1; j<minmn; j++)
                r[j][i]=a[j][i];
        }
    }

    if (!(Q==Matrix::emptyMatrix))
    {
        Q.setSize(nc,nc);
        double **q=Q;
        Q.diagonal(1.0);
        for (j=nl-1; j>=0; j--)
        {
            if (a[j][j]==0.0) continue;
            for (k=j; k<nc; k++)
            {
                sum=0.0;
                for (i=j; i<nc; i++) sum=sum+a[j][i]*q[i][k];

                temp = sum/a[j][j];
                for (i=j; i<nc; i++) q[i][k]=q[i][k]-temp*a[j][i];
            }
        }
    }
}

#endif

void Matrix::addUnityInPlace(double dd)
{
    int nn=d->extColumn+1, i=nLine();
    double *a=*d->p;
    while (i--)
    {
        (*a)+=dd;
        a+=nn;
    };
}

double Matrix::frobeniusNorm()
{
// no tested
// same code for the Vector eucilidian norm
    /*
        double sum=0, *a=*p;
        int i=nLine()*nColumn();
        while (i--) sum+=sqr(*(a++));
        return sqrt(sum);
    */
    return ::euclidianNorm(nLine()*nColumn(),*d->p);
}

double Matrix::LnftyNorm()
{
// not tested
    double m=0, sum;
    int j,nl=nLine(), nc=nColumn();
    double **a=(*this), *xp;
    while (nl--)
    {
        sum=0;
        j=nc;
        xp=*(a++);
        while (j--) sum+=condorAbs(*(xp++));
        m=::mmax(m,sum);
    }
    return m;
}

Vector Matrix::getMaxColumn()
{
    double **a=(*this), sum, maxSum=0;
    int i=nColumn(),j,k=0, nl=nLine();
    while (i--)
    {
        sum=0;
        j=nl;
        while(j--) sum+=sqr(a[j][i]);
        if (sum>maxSum)
        {
            maxSum=sum;
            k=i;
        }
    }
    Vector rr(nl);
    double *r=rr;
    j=nl;
//    while (j--) *(r++)=a[j][k];
    while (j--) r[j]=a[j][k];
    return rr;
}

Vector Matrix::getLine(int i, int n, int startc)
{
    if (n==0) n=nColumn()-startc;
    Vector r(n,d->p[i]+startc);
    return r;
}

void Matrix::getLine(int i, Vector r, int n, int startc)
{
    if (n==0) n=nColumn()-startc;
    r.setSize(n);
    memcpy((double*)r, d->p[i]+startc, n*sizeof(double));
}

Vector Matrix::getColumn(int i, int n)
{
    if (n==0) n=nLine();
    Vector r(n);
    double **d=*this, *rr=r;
    while (n--) rr[n]=d[n][i];
    return r;
}

void Matrix::getColumn(int i, Vector r, int n)
{
    if (n==0) n=nLine();
    r.setSize(n);
    double **d=*this, *rr=r;
    while (n--) rr[n]=d[n][i];
}

void Matrix::setLine(int i, Vector v, int n)
{
    if (n==0) n=nColumn();
    memcpy(d->p[i], (double*)v, n*sizeof(double));
}

void Matrix::setLines(int indexDest, Matrix Source, int indexSource, int number)
{
    if (!Source.nLine()) return;
    double **dest=(*this), **sour=Source;
    int snl=d->nColumn*sizeof(double);
    if (number==0) number=Source.nLine()-indexSource;
    while (number--) memcpy(dest[indexDest+number], sour[indexSource+number], snl);
}

double Matrix::euclidianNorm(int i)
{
    return ::euclidianNorm(nColumn(), d->p[i]);
}

void Matrix::getSubMatrix(Matrix R, int startL, int startC, int nl, int nc)
{
    if (nl==0) nl=  nLine()-startL;
    else nl=mmin(nl,  nLine()-startL);
    if (nc==0) nc=nColumn()-startC;
    else nc=mmin(nc,nColumn()-startC);
    R.setSize(nl,nc);
    double **sd=(*this), **dd=R;
    while (nl--)
        memcpy(dd[nl], sd[nl+startL]+startC, nc*sizeof(double));
}

void Matrix::swapLines(int i, int j)
{
    if (i==j) return;
    int n=nColumn();
    double *p1=d->p[i], *p2=d->p[j], t;
    while (n--)
    {
        t=p1[n];
        p1[n]=p2[n];
        p2[n]=t;
    }
}

int Matrix::lineIndex(Vector r, int nn)
{
    if (nn==0) nn=mmin((int)nColumn(),(int)r.sz())*sizeof(double);
    else nn*=sizeof(double);

    int i=nLine();
    double **dp=d->p, *dp2=r;
    while (i--)
        if (memcmp(dp[i],dp2,nn)==0) break;
    return i;
}

void Matrix::setColNames(char **c, int nc)
{
    if (c==NULL)
    {
        if (d->columnNames)
        {
            free(*d->columnNames);
            free(d->columnNames);
        }
        return;
    }
    int l=0,i;
    if (nc==0) nc=d->nColumn;
    d->columnNames=(char**)malloc(nc*sizeof(char*));
    for (i=0; i<nc; i++) l+=(int)strlen(c[i])+1;
    char *t1=(char*)malloc(l);
    for (i=0; i<nc; i++)
    {
        d->columnNames[i]=t1;
        strcpy(t1,c[i]);
        t1+=(int)strlen(t1)+1;
    }
}

void Matrix::merge(Matrix m,int eliminateDoubles)
{
    int nc=nColumn(), nlm=m.nLine();
    if (nlm==0) return;
    if (nc!=m.nColumn())
    {
        printf("Merging: Size do not agree.\n");
        exit(255);
    }
    int nl=nLine(),i,j;
    nc*=sizeof(double);
    double *pdi;
    for (i=0; i<nlm; i++)
    {
        pdi=((double**)m)[i];
        for (j=0; j<nl; j++)
            if (memcmp(d->p[j],pdi,nc)==0) break;
        if (j!=nl) continue;
        extendLine();
        memcpy(d->p[nl],pdi,nc);
        nl++;
    }
}

void Matrix::append(Vector tmp)
{
    int nl=nLine(), nc=nColumn(), mdim=tmp.sz();
    if (nc==0) setSize(1,mdim);
    else extendLine();
    setLine(nl,tmp,mdim);
}

/*
int Matrix::solve(Vector vB)
{
    double t;
    int i, j, k, l, info=0;
    int nl=nLine(), nc=nColumns();

    // gaussian elimination with partial pivoting
    if ( nl>1 )
    {
        for ( k=0; k<nl-1 ; k++ )
        {
            // find l = pivot index
            l=k; maxp=condorAbs(x[k][k]);
            for (i=k+1; i<nl; i++)
                if (condorAbs(x[i][k])>maxp) { maxp=condorAbs(x[i][k]); l=i; }

            jpvt[k] = l;
            // zero pivot implies this column already triangularized
            if ( maxp==0.0 ) info = k;
            else
            {
                // interchange if necessary
                if ( l!=k )
                {
                    for (i=k; i<nc; i++)
                    {
                        t=x[l][i];
                        x[l][i]=x[k][k];
                        x[k][k]=t;
                    }
                    t=b[l]; b[l]=b[k]; b[k]=t;
                }
                // compute multipliers
                maxp=-1.0/maxp;
                for (i=k+1; i<nc; j++ ) x[k][i]*=maxp;

                // row elimination
                for (j=k+1; j<nl; j++ )
                {
                    t=x[k][j];
                    for (i=k+1; i<nc; i++) x[j][i] += t*x[k][i];
                }
            }
        }
    }
    if ( x[nl-1][nl-1]==0.0 ) info=nl-1;
    return;
}
*/
