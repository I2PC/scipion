/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/
/*****************************************************************************/
/* Variable and prototype definitions for the Numerical Core                 */
/*****************************************************************************/
#ifndef _NUMERICAL_HH
#   define _NUMERICAL_HH

// Random numbers
double ran1(int *idum);                                 // Uniform random
double gasdev(int *idum);                               // Gaussian random

// Working with matrices
template <class T> void ludcmp(T **a, int n, int *indx, T *d);  
                                                       // LU decomposition
template <class T> void lubksb(T **a, int n, int *indx,T b[]);  
                                                       // Solve Ax=b
template <class T> void gaussj(T **a, int n, T **b, int m);
                                                       // Solve Ax=b (b=matrix)

// FFT
void four1(double data[],int nn,int isign);             // Complex FFT 1D
void realft(double data[],int n,int isign);             // Real FFT 1D
void fourn(double data[],int nn[],int ndim,int isign);  // Complex FFT 2D,3D,...

// Sorting
void indexx(int n, double arrin[], int indx[]);         // Sorting indexes
void qcksrt(int n, double arr[]);                       // Sorting

// Bessel functions
double bessj0(double x);
double bessj3_5(double x);

double bessi0(double x);
double bessi1(double x);
double bessi0_5(double x);
double bessi1_5(double x);
double bessi2(double x);
double bessi2_5(double x);
double bessi3_5(double x);

// Special functions
double gammln(double xx);
double betacf(double a, double b, double x);
double betai(double a, double b, double x);

// Eigen decomposition
void jacobi(double **a, int n, double *d, double **v, int *nrot);
// These functions are needed two obtain eigenvalues of real non-symetric
// matrices 
void balanc(double **a, int n);
void elmhes(double **a, int n);
void hqr(double **a, int n,double *wr,double *wi);

// Singular value descomposition of matrix a (numerical recipes, chapter 2-6 for details)
void svdcmp(double **a, int m,int n,double *w, double **v);
void svbksb(double **u,double *w,double **v, int m,int n,double *b,double *x);

// 
void convlv(double *data,int n,double *respns,int m,int isign,double *ans);
void realft(double *data,int n,int isign);
void twofft(double *data1,double *data2,double *fft1,double *fft2,int n);
void savgol(double *c, int np, int nl, int nr, int ld, int m);
void four1(double *data,int nn,int isign);

// Optimization
void powell(double *p, double **xi, int n, double ftol, int &iter,
   double &fret, double (*func)(double *), bool show);

// Wavelets
void wt1(double a[], unsigned long n, int isign,
	void (*wtstep)(double [], unsigned long, int));
void wtn(double a[], unsigned long nn[], int ndim, int isign,
	void (*wtstep)(double [], unsigned long, int));
void pwtset(int n);
void pwt(double a[], unsigned long n, int isign);

#endif
