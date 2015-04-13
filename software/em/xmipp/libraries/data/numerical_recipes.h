/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
/*****************************************************************************/
/* Variable and prototype definitions for the Numerical Core                 */
/*****************************************************************************/
#ifndef _NUMERICAL_HH
#define _NUMERICAL_HH

#include <math.h>
#include "xmipp_memory.h"
#include "xmipp_macros.h"

//@defgroup NumericalRecipes Functions from the Numerical Recipes
//@ingroup DataLibrary
//@{

// Utilities --------------------------------------------------------------
void nrerror(const char error_text[]);

// Random numbers ---------------------------------------------------------
double ran1(int *idum);                                 // Uniform random
double gasdev(int *idum);                               // Gaussian random
double tdev(double nu, int *idum);                      // t-student random

// Cumulative distribution functions and Kolmogorov-Smirnov test
void   ksone(double data[], int n, double(*func)(double), double * d, double * prob);  // Chapter 13.5
double probks(double alam);  // Chapter 13.5

// FFT ---------------------------------------------------------------------
void four1(double data[], int nn, int isign);           // Complex FFT 1D
void realft(double data[], int n, int isign);           // Real FFT 1D
void fourn(double data[], int nn[], int ndim, int isign);  // Complex FFT 2D,3D,...

// Sorting -----------------------------------------------------------------
void indexx(int n, double arrin[], int indx[]);         // Sorting indexes
void qcksrt(int n, double arr[]);                       // Sorting

// Bessel functions --------------------------------------------------------
double bessj0(double x);
double bessj3_5(double x);
double bessj1_5(double x);

double bessi0(double x);
double bessi1(double x);
double bessi0_5(double x);
double bessi1_5(double x);
double bessi2(double x);
double bessi3(double x);
double bessi2_5(double x);
double bessi3_5(double x);
double bessi4(double x);

// Special functions -------------------------------------------------------
double gammln(double xx);
double gammp(double a, double x);
double betacf(double a, double b, double x);
double betai(double a, double b, double x);
inline double sinc(double x)
{
    if (fabs(x)<0.0001)
        return 1;
    else
    {
        double arg=PI*x;
        return sin(arg)/arg;
    }
}
inline size_t fact(int num)
{
    size_t value = 1;
    if (num !=1 && num != 0)
    {
        for (int i = num; i > 0; --i)
            value *= i;
    }
    return value;

}

inline size_t binom(int n, int k)
{
    size_t factor=1;
    for (int i = n; i > (n-k); --i)
        factor *= i;
    return factor/fact(k);
}

// Singular value descomposition of matrix a (numerical recipes, chapter 2-6 for details)
void svdcmp(double *a, int m, int n, double *w, double *v);
void svbksb(double *u, double *w, double *v, int m, int n, double *b, double *x);

// -------------------------------------------------------------------------
void convlv(double *data, int n, double *respns, int m, int isign, double *ans);
void realft(double *data, int n, int isign);
void twofft(double *data1, double *data2, double *fft1, double *fft2, int n);
void savgol(double *c, int np, int nl, int nr, int ld, int m);
void four1(double *data, int nn, int isign);

// Optimization ------------------------------------------------------------
void powell(double *p, double *xi, int n, double ftol, int &iter,
            double &fret, double(*func)(double *, void *), void *prm,
            bool show);

void amebsa(double **p, double y[], int ndim, double pb[], double *yb,
            double ftol, double (*funk)(double []), int *iter, double temptr);

// These two routines have been taken from
// http://users.utu.fi/vesoik/userdocs/programs/libpet
// and they implement an algorithm of Lawson-Hanson of
// nonnegative least squares
int nnls(double *a, int m, int n, double *b, double *x,
         double *rnorm, double *w, double *zz, int *index);
int nnlsWght(int N, int M, double *A, double *b, double *weight);

// CFSQP -------------------------------------------------------------------
// These routines are from
// http://www.aemtechnology.com/aemdesign/downloadfsqp.htm
// They implement the CFSQP algorithm
/* Declare and initialize user-accessible flag indicating    */
/* whether x sent to user functions has been changed within  */
/* CFSQP.            */
// Gradients - Finite Difference
void    grobfd(int, int, double *, double *, void(*)(int, int,
               double *, double *, void *), void *);
void    grcnfd(int, int, double *, double *, void(*)(int, int,
               double *, double *, void *), void *);

// CFSQP
void    cfsqp(int, int, int, int, int, int, int, int, int, int *, int, int,
              int, int *, double, double, double, double, double *,
              double *, double *, double *, double *, double *,
              void(*)(int, int, double *, double *, void *),
              void(*)(int, int, double *, double *, void *),
              void(*)(int, int, double *, double *,
                      void(*)(int, int, double *, double *, void *), void *),
              void(*)(int, int, double *, double *,
                      void(*)(int, int, double *, double *, void *), void *),
              void *);

// Wavelets ----------------------------------------------------------------
void wt1(double a[], unsigned long n, int isign,
         void(*wtstep)(double [], unsigned long, int));
void wtn(double a[], unsigned long nn[], int ndim, int isign,
         void(*wtstep)(double [], unsigned long, int));
void pwtset(int n);
void pwt(double a[], unsigned long n, int isign);

// Working with matrices ---------------------------------------------------
// LU decomposition
#define TINY 1.0e-20;
/* Chapter 2 Section 3: LU DECOMPOSITION */
template <class T>
void ludcmp(T *a, int n, int *indx, T *d)
{
    int i, imax=0, j, k;
    T big, dum, sum, temp;
    T *vv;

    ask_Tvector(vv, 1, n);
    *d = (T)1.0;
    for (i = 1;i <= n;i++)
    {
        big = (T)0.0;
        for (j = 1;j <= n;j++)
            if ((temp = (T)fabs((double)a[i*n+j])) > big)
                big = temp;
        if (big == (T)0.0)
            nrerror("Singular matrix in routine LUDCMP");
        vv[i] = (T)1.0 / big;
    }
    for (j = 1;j <= n;j++)
    {
        for (i = 1;i < j;i++)
        {
            sum = a[i*n+j];
            for (k = 1;k < i;k++)
                sum -= a[i*n+k] * a[k*n+j];
            a[i*n+j] = sum;
        }
        big = (T)0.0;
        for (i = j;i <= n;i++)
        {
            sum = a[i*n+j];
            for (k = 1;k < j;k++)
                sum -= a[i*n+k] * a[k*n+j];
            a[i*n+j] = sum;
            if ((dum = vv[i] * (T)fabs((double)sum)) >= big)
            {
                big = dum;
                imax = i;
            }
        }
        if (j != imax)
        {
            for (k = 1;k <= n;k++)
            {
                dum = a[imax*n+k];
                a[imax*n+k] = a[j*n+k];
                a[j*n+k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j*n+j] == 0.0)
            a[j*n+j] = (T) TINY;
        if (j != n)
        {
            dum = (T)1.0 / (a[j*n+j]);
            for (i = j + 1;i <= n;i++)
                a[i*n+j] *= dum;
        }
    }
    free_Tvector(vv, 1, n);
}
#undef TINY

// Solve Ax=b
/* Chapter 2 Section 3: LU BACKWARD-FORWARD SUBSTITUTION */
template <class T>
void lubksb(T *a, int n, int *indx, T b[])
{
    int i, ii = 0, ip, j;
    T sum;

    for (i = 1;i <= n;i++)
    {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii)
            for (j = ii;j <= i - 1;j++)
                sum -= a[i*n+j] * b[j];
        else if (sum)
            ii = i;
        b[i] = sum;
    }
    for (i = n;i >= 1;i--)
    {
        sum = b[i];
        for (j = i + 1;j <= n;j++)
            sum -= a[i*n+j] * b[j];
        b[i] = sum / a[i*n+i];
    }
}

/* Chapter 2, Section 1. Gauss-Jordan equation system resolution ----------- */
// Solve Ax=b (b=matrix)
template <class T>
void gaussj(T *a, int n, T *b, int m)
{
    T temp;
    int *indxc, *indxr, *ipiv;
    int i, icol=0, irow=0, j, k, l, ll;
    T big, dum;
    double pivinv;

    ask_Tvector(indxc, 1, n);
    ask_Tvector(indxr, 1, n);
    ask_Tvector(ipiv, 1, n);
    for (j = 1;j <= n;j++)
        ipiv[j] = 0;
    for (i = 1;i <= n;i++)
    {
        big = (T)0;
        for (j = 1;j <= n;j++)
            if (ipiv[j] != 1)
                for (k = 1;k <= n;k++)
                {
                    if (ipiv[k] == 0)
                    {
                        if (fabs((double)a[j*n+k]) >= (double) big)
                        {
                            big = ABS(a[j*n+k]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] > 1)
                        nrerror("GAUSSJ: Singular Matrix-1");
                }
        ++(ipiv[icol]);
        if (irow != icol)
        {
            for (l = 1;l <= n;l++)
                SWAP(a[irow*n+l], a[icol*n+l], temp)
                for (l = 1;l <= m;l++)
                    SWAP(b[irow*n+l], b[icol*n+l], temp)
                }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol*n+icol] == 0.0)
            nrerror("GAUSSJ: Singular Matrix-2");
        pivinv = 1.0f / a[icol*n+icol];
        a[icol*n+icol] = (T)1;
        for (l = 1;l <= n;l++)
            a[icol*n+l] = (T)(pivinv * a[icol*n+l]);
        for (l = 1;l <= m;l++)
            b[icol*n+l] = (T)(pivinv * b[icol*n+l]);
        for (ll = 1;ll <= n;ll++)
            if (ll != icol)
            {
                dum = a[ll*n+icol];
                a[ll*n+icol] = (T)0;
                for (l = 1;l <= n;l++)
                    a[ll*n+l] -= a[icol*n+l] * dum;
                for (l = 1;l <= m;l++)
                    b[ll*n+l] -= b[icol*n+l] * dum;
            }
    }
    for (l = n;l >= 1;l--)
    {
        if (indxr[l] != indxc[l])
            for (k = 1;k <= n;k++)
                SWAP(a[k*n+indxr[l]], a[k*n+indxc[l]], temp);
    }
    free_Tvector(ipiv, 1, n);
    free_Tvector(indxr, 1, n);
    free_Tvector(indxc, 1, n);
}

// Cholesky factorization
void choldc(double *a, int n, double *p);
// Cholesky backsubstitution
void cholsl(double *a, int n, double *p, double *b, double *x);
// Polynomial interpolation
void polint(double *xa, double *ya, int n, double x, double &y, double &dy);
//@}

#endif
