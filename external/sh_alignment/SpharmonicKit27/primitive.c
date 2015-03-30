/***************************************************************************
  **************************************************************************
  
                Spherical Harmonic Transform Kit 2.7
  
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
  
   Copyright 1997-2003  Sean Moore, Dennis Healy,
                        Dan Rockmore, Peter Kostelec
  
  
   Copyright 2004  Peter Kostelec, Dan Rockmore


     SpharmonicKit is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
  
     SpharmonicKit is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
  
     You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
  
  
   Commercial use is absolutely prohibited.
  
   See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/


/*

  some "primitive" functions that are used in cospmls.c
  and newmathx.c

  */

#include <math.h>
#include <string.h>  /* to declare memcpy */

#ifndef PI
#define PI 3.14159265358979
#endif

/************************************************************************/
/* Recurrence coefficients */
/************************************************************************/
/* Recurrence coefficents for L2-normed associated Legendre
   recurrence.  When using these coeffs, make sure that
   inital Pmm function is also L2-normed */
/* l represents degree, m is the order */

double L2_an(int m,
	     int l)
{
  return (sqrt((((double) (2*l+3))/((double) (2*l+1))) *
	       (((double) (l-m+1))/((double) (l+m+1)))) *
	  (((double) (2*l+1))/((double) (l-m+1))));

}

/* note - if input l is zero, need to return 0 */
double L2_cn(int m,
	     int l) 
{
  if (l != 0) {
    return (-1.0 *
	  sqrt((((double) (2*l+3))/((double) (2*l-1))) *
	       (((double) (l-m+1))/((double) (l+m+1))) *
	       (((double) (l-m))/((double) (l+m)))) *
	  (((double) (l+m))/((double) (l-m+1))));
  }
  else
    return 0.0;

}

/* when using the reverse recurrence, instead of calling
   1/L2_cn_tr(m,l), let me just define the function ...
   it might be more stable */

double L2_cn_inv(int m,
		 int l)
{
  double dl, dm;

  dl = (double) l;
  dm = (double) m;

  return ( -(1.0 + (1. - 2. * dm)/(dm + dl)) *
	   sqrt( ((-1. + 2.*dl)/(3. + 2.*dl)) *
		 ((dl + dl*dl + dm + 2.*dl*dm + dm*dm)/
		  (dl + dl*dl - dm - 2.*dl*dm + dm*dm)) )
	   );

}

/* when using the reverse recurrence, instead of calling
   -L2_an(m,l)/L2_cn(m,l), let me just define the
   function ... it might be more stable */

double L2_ancn(int m,
	       int l)
{
  double dl, dm;

  dl = (double) l;
  dm = (double) m;

  return( sqrt( 4.0 + ( (4.0 * dm * dm - 1.0)/
			(dl * dl - dm * dm) ) ) );
}

/************************************************************************/
/* vector arithmetic operations */
/************************************************************************/
/* does result = data1 + data2 */
/* result and data are vectors of length n */

void vec_add(double *data1,
	     double *data2,
	     double *result,
	     int n)
{
  int k;


  for (k = 0; k < n % 4; ++k)
    result[k] = data1[k] + data2[k];

  for ( ; k < n ; k += 4)
    {
      result[k] = data1[k] + data2[k];
      result[k + 1] = data1[k + 1] + data2[k + 1];
      result[k + 2] = data1[k + 2] + data2[k + 2];
      result[k + 3] = data1[k + 3] + data2[k + 3];
    }
}
/************************************************************************/
/************************************************************************/
/*
   vec_mul(scalar,data1,result,n) multiplies the vector 'data1' by
   'scalar' and returns in result 
*/
void vec_mul(double scalar,
	     double *data1,
	     double *result,
	     int n)
{
   int k;


   for( k = 0; k < n % 4; ++k)
     result[k] = scalar * data1[k];

   for( ; k < n; k +=4)
     {
       result[k] = scalar * data1[k];
       result[k + 1] = scalar * data1[k + 1];
       result[k + 2] = scalar * data1[k + 2];
       result[k + 3] = scalar * data1[k + 3];
     }

}
/************************************************************************/
/* point-by-point multiplication of vectors */

void vec_pt_mul(double *data1,
		double *data2,
		double *result,
		int n)
{
   int k;

  
  for(k = 0; k < n % 4; ++k)
    result[k] = data1[k] * data2[k];
  
  for( ; k < n; k +=4)
    {
      result[k] = data1[k] * data2[k];
      result[k + 1] = data1[k + 1] * data2[k + 1];
      result[k + 2] = data1[k + 2] * data2[k + 2];
      result[k + 3] = data1[k + 3] * data2[k + 3];
    }
 
}


/************************************************************************/
/* returns an array of the angular arguments of n Chebyshev nodes */
/* eval_pts points to a double array of length n */

void ArcCosEvalPts(int n,
		   double *eval_pts)
{
    int i;
    double twoN;

    twoN = (double) (2 * n);

   for (i=0; i<n; i++)
     eval_pts[i] = (( 2.0*((double)i)+1.0 ) * PI) / twoN;

}
/************************************************************************/
/* returns an array of n Chebyshev nodes */

void EvalPts( int n,
	      double *eval_pts)
{
    int i;
    double twoN;

    twoN = (double) (2*n);

   for (i=0; i<n; i++)
     eval_pts[i] = cos((( 2.0*((double)i)+1.0 ) * PI) / twoN);

}

/************************************************************************/
/* L2 normed Pmm.  Expects input to be the order m, an array of
 evaluation points arguments of length n, and a result vector of length n */
/* The norming constant can be found in Sean's PhD thesis */
/* This has been tested and stably computes Pmm functions thru bw=512 */

void Pmm_L2( int m,
	     double *eval_pts,
	     int n,
	     double *result)
{
  int i;
  double md, id, mcons;

  id = (double) 0.0;
  md = (double) m;
  mcons = sqrt(md + 0.5);

  for (i=0; i<m; i++) {
    mcons *= sqrt((md-(id/2.0))/(md-id));
    id += 1.0;
  }
  if (m != 0 )
    mcons *= pow(2.0,-md/2.0);
  if ((m % 2) != 0) mcons *= -1.0;

  for (i=0; i<n; i++) 
    result[i] = mcons * pow(sin(eval_pts[i]),((double) m));

}

/************************************************************************/
/************************************************************************/
/* This piece of code synthesizes a function which is the weighted sum of 
   associated Legendre functions.  The coeffs array should contain
   bw - m coefficients ordered from zeroth degree to bw-1, and eval_pts
   should be an array of the arguments (arccos) of the desired 
   evaluation points of length 2*bw.  Answer placed
   in result (and has length 2*bw).

   workspace needs to be of size 16 * bw

   workspace needs to be of size 14 * bw
   */
/************************************************************************/
void P_eval(int m,
	    double *coeffs,
	    double *eval_args,
	    double *result,
	    double *workspace,
	    int bw)
{
    double *prev, *prevprev, *temp1, *temp2, *temp3, *temp4, *x_i;
    int i, j, n;
    double splat;

    prevprev = workspace;
    prev = prevprev + (2*bw);
    temp1 = prev + (2*bw);
    temp2 = temp1 + (2*bw);
    temp3 = temp2 + (2*bw);
    temp4 = temp3 + (2*bw);
    x_i = temp4 + (2*bw);

    n = 2*bw;

    /* now get the evaluation nodes */
    EvalPts(n,x_i);

    /*   for(i=0;i<n;i++)
      fprintf(stderr,"in P_eval evalpts[%d] = %lf\n", i, x_i[i]);
      */   
    for (i=0; i<n; i++) 
      prevprev[i] = 0.0;

    if (m == 0) {
	for (i=0; i<n; i++) {
	  /* prev[i] = 0.707106781186547; sqrt(1/2) */
	   prev[i] = 1.0;
	    /* now mult by first coeff and add to result */
	    result[i] = coeffs[0] * prev[i];
	}
    }
    else {
	Pmm_L2(m, eval_args, n, prev);
	splat = coeffs[0];
	for (i=0; i<n; i++)
	  result[i] = splat * prev[i];
    }

    for (i=0; i<bw-m-1; i++) {
	vec_mul(L2_cn(m,m+i),prevprev,temp1,n);
	vec_pt_mul(prev, x_i, temp2, n);
	vec_mul(L2_an(m,m+i), temp2, temp3, n);
	vec_add(temp3, temp1, temp4, n); /* temp4 now contains P(m,m+i+1) */
	/* now add weighted P(m,m+i+1) to the result */
	splat = coeffs[i+1];
	for (j=0; j<n; j++)
	  result[j] += splat * temp4[j];
	memcpy(prevprev, prev, sizeof(double) * n);
	memcpy(prev, temp4, sizeof(double) * n);
    }

}

