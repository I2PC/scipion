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


/*************************************************************************/

/* Source code to synthesize functions using a naive method
   based on recurrence.  This is slow but does not require any
   precomputed functions, and is also stable. 
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "csecond.h"
#include "weights.h"

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

double ns_L2_an(int m, int l)
{
  return (sqrt((((double) (2*l+3))/((double) (2*l+1))) *
	       (((double) (l-m+1))/((double) (l+m+1)))) *
	  (((double) (2*l+1))/((double) (l-m+1))));
}

double ns_L2_bn()
{
  return ((double) 0.0);
}

/* note - if input l is zero, need to return 0 */
double ns_L2_cn(int m, int l) 
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
/************************************************************************/
/* vector arithmetic operations */
/************************************************************************/
/* does result = data1 + data2 */
/* result and data are vectors of length n */

void ns_vec_add(double *data1,
		double *data2,
		double *result,
		int n)
{

  int k;

  for(k = 0; k < n % 4; ++k)
    {
      result[k] = data1[k] + data2[k];
    }

  for( ; k < n; k += 4)
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
   vec_mul(scalar,data1,n) multiplies the vector 'data1' by
   'scalar' and returns in result 
*/
void ns_vec_mul(double scalar,
		double *data1,
		double *result,
		int n)
{
   int k;


   for( k = 0; k < n % 4; ++k)
     result[k] = scalar * data1[k];

   for( ; k < n; k += 4)
     {
      result[k] = scalar * data1[k];
      result[k + 1] = scalar * data1[k + 1];
      result[k + 2] = scalar * data1[k + 2];
      result[k + 3] = scalar * data1[k + 3];
    }
}


/************************************************************************/
/* point-by-point multiplication of vectors */

void ns_vec_pt_mul(const double *data1,
		   const double *data2,
		   double *result,
		   int n)
{
   int k;

   for (k = 0; k < n % 4; ++k)
     result[k] = data1[k] * data2[k];

   for ( ; k < n; k += 4)
     {
       result[k] = data1[k] * data2[k];
       result[k + 1] = data1[k + 1] * data2[k + 1];
       result[k + 2] = data1[k + 2] * data2[k + 2];
       result[k + 3] = data1[k + 3] * data2[k + 3];
     }

}


/************************************************************************/
/************************************************************************/
/* returns an array of the angular arguments of n Chebyshev nodes */
/* eval_pts points to a double array of length n */

void ns_ArcCosEvalPts(int n,
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

void ns_EvalPts(int n,
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

void ns_Pmm_L2(int m,
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
/* This is the procedure that synthesizes a function from a list
   of coefficients of a Legendre series.  Function is synthesized
   at the (2*bw) Chebyshev nodes.
   
   bw - bandwidth
   m - order
   coeffs - a pointer to double array of size (bw-m).  First coefficient is
            coefficient for Pmm
   result - a pointer to double array of size (2*bw) and containing the
            synthesized function
   workspace - a pointer to double array of size (32*bw)
   

*/

void Naive_Synthesize(int bw,
		      int m,
		      double *coeffs,
		      double *result,
		      double *workspace)
{
    double *prev, *prevprev, *temp1, *temp2, *temp3, *temp4, *x_i, *eval_args;
    double t1;
    int i, j;

    prevprev = workspace;
    prev = prevprev + (2*bw);
    temp1 = prev + (2*bw);
    temp2 = temp1 + (2*bw);
    temp3 = temp2 + (2*bw);
    temp4 = temp3 + (2*bw);
    x_i = temp4 + (2*bw);
    eval_args = x_i + (2*bw);


    /* get the evaluation nodes */
    ns_EvalPts(2*bw,x_i);
    ns_ArcCosEvalPts(2*bw,eval_args);
    

    /* set initial values of first two Pmls */
    for (i=0; i<2*bw; i++) 
      prevprev[i] = 0.0;
    if (m == 0) {
	for (i=0; i<2*bw; i++) {
	    prev[i] = 1.0;
	}
    }
    else 
      ns_Pmm_L2(m, eval_args, 2*bw, prev);

    /* this may be useful for generating Gmls later but
       removed for now

    if ((m % 2) == 1) { 
	for (i=0; i<2*bw; i++)
	  prev[i] /= sin(eval_args[i]);
    }
    */

    /* make sure result is zeroed out */
    for (i=0; i<(2*bw); i++)
      result[i] = 0.0;

    /* add in Pmm contribution */

    if (coeffs[0] != 0.0)
      {
	t1 = coeffs[0];
	for (j=0; j<(2*bw); j++)
	  result[j] += t1 * prev[j];
      }

    /* now generate remaining pmls while synthesizing function */

    for (i=0; i<bw-m-1; i++) {
	ns_vec_mul(ns_L2_cn(m,m+i),prevprev,temp1,2*bw);
	ns_vec_pt_mul(prev, x_i, temp2, 2*bw);
	ns_vec_mul(ns_L2_an(m,m+i), temp2, temp3, 2*bw);
	ns_vec_add(temp3, temp1, temp4, 2*bw); /* temp4 now contains P(m,m+i+1) */

	/* add in contribution */
	if (coeffs[i+1] != 0.0)
	  {
	    t1 = coeffs[i+1];
	    for (j=0; j<(2*bw); j++)
	      result[j] += (t1 * temp4[j]);
	  }

	/* now update Pi and P(i+1) */
	memcpy(prevprev, prev, (size_t) sizeof(double) * 2 * bw);
	memcpy(prev, temp4, (size_t) sizeof(double) * 2 * bw);

    }
  }

/************************************************************************/
/* Naive Analysis function - used to get a measure of the error inherent
   in the synthesis of functions using three-term recurrence.  Essentially,
   this is a naive transform function

*/
/************************************************************************/
/*
   bw - bandwidth
   m - order
   data - a pointer to double array of size (2*bw) containing a synthesized
          function.  
   result - a pointer to double array of size (bw-m) and containing the
            computed Legendre coefficients, starting with the Pmm
	    coefficient.
   workspace - a pointer to double array of size (32*bw)
*/


void Naive_Analysis(int bw,
		    int m,
		    double *data,
		    double *result,
		    double *workspace)
{
    double *prev, *prevprev, *temp1, *temp2, *temp3, *temp4, *x_i, *eval_args;
    double *wdata;
    const double *weight_vec;
    double t1;
    int i, j;

    prevprev = workspace;
    prev = prevprev + (2*bw);
    temp1 = prev + (2*bw);
    temp2 = temp1 + (2*bw);
    temp3 = temp2 + (2*bw);
    temp4 = temp3 + (2*bw);
    x_i = temp4 + (2*bw);
    eval_args = x_i + (2*bw);
    wdata = eval_args + (2*bw);


    /* get the evaluation nodes */
    ns_EvalPts(2*bw,x_i);
    ns_ArcCosEvalPts(2*bw,eval_args);
    

    /* set initial values of first two Pmls */
    for (i=0; i<2*bw; i++) 
      prevprev[i] = 0.0;
    if (m == 0) {
	for (i=0; i<2*bw; i++) {
	    prev[i] = 0.5;
	}
    }
    else 
      ns_Pmm_L2(m, eval_args, 2*bw, prev);

    /* make sure result is zeroed out */
    for (i=0; i<(bw-m); i++)
      result[i] = 0.0;

    /* apply quadrature weights */
    weight_vec = get_weights(bw);

    for (i=0; i<(2*bw); i++)
      wdata[i] = data[i] * weight_vec[i];


    /* compute Pmm coefficient */
    t1 = 0.0;
    for (j=0; j<(2*bw); j++)
      t1 += wdata[j] * prev[j];
    result[0] = t1;

    /* now generate remaining pmls while computing coefficients */

    for (i=0; i<bw-m-1; i++) {
	ns_vec_mul(ns_L2_cn(m,m+i),prevprev,temp1,2*bw);
	ns_vec_pt_mul(prev, x_i, temp2, 2*bw);
	ns_vec_mul(ns_L2_an(m,m+i), temp2, temp3, 2*bw);
	ns_vec_add(temp3, temp1, temp4, 2*bw); /* temp4 now contains P(m,m+i+1) */
	
	/* compute this coefficient */
	t1 = 0.0;
	for (j=0; j<(2*bw); j++)
	  t1 += wdata[j] * temp4[j];
	result[i+1] = t1;

	/* now update Pi and P(i+1) */
	/***
	for (j=0; j<2*bw; j++) {
	    prevprev[j] = prev[j];
	    prev[j] = temp4[j];
	}
	***/
	memcpy( prevprev, prev, sizeof(double) * 2 * bw );
	memcpy( prev, temp4, sizeof(double) * 2 * bw );
    }
  }


/************************************************************************/
/* Naive Analysis function for timing - used to get a measure of the 
   naive transformation time.

*/
/************************************************************************/
/*
   bw - bandwidth
   m - order
   data - a pointer to double array of size (2*bw) containing a synthesized
          function.  
   result - a pointer to double array of size (bw-m) and containing the
            computed Legendre coefficients, starting with the Pmm
	    coefficient.
   workspace - a pointer to double array of size (32*bw)

   timing - 1 to turn on timing info

   runtime - a double slot for writing time

   loops - number of timing loops

*/


void Naive_Analysis_Timing(double *data,
			   int bw,
			   int m,
			   double *result, 
			   int timing,
			   double *runtime,
			   int loops,
			   double *workspace)
{
    double *prev, *prevprev, *temp1, *temp2, *temp3, *temp4, *x_i, *eval_args;
    double *wdata;
    const double *weight_vec;
    int i, j, l;

    double total_time, tstart, tstop;

    prevprev = workspace;
    prev = prevprev + (2*bw);
    temp1 = prev + (2*bw);
    temp2 = temp1 + (2*bw);
    temp3 = temp2 + (2*bw);
    temp4 = temp3 + (2*bw);
    x_i = temp4 + (2*bw);
    eval_args = x_i + (2*bw);
    wdata = eval_args + (2*bw);


    /* get the evaluation nodes */
    ns_EvalPts(2*bw,x_i);
    ns_ArcCosEvalPts(2*bw,eval_args);
    


  /* start the timer
    if (timing)
      gettimeofday(&ftstart,0); */
    if (timing)
      tstart = csecond();

  /* main timing loop */
  for (l=0; l< loops; l++)
    {
      /* set initial values of first two Pmls */
      for (i=0; i<2*bw; i++) 
	prevprev[i] = 0.0;
      if (m == 0) {
	for (i=0; i<2*bw; i++) {
	  prev[i] = 0.5;
	}
      }
      else 
	ns_Pmm_L2(m, eval_args, 2*bw, prev);

      /* make sure result is zeroed out */
      for (i=0; i<(bw-m); i++)
	result[i] = 0.0;

      /* apply quadrature weights */
      weight_vec = get_weights(bw);

      for (i=0; i<(2*bw); i++)
	wdata[i] = data[i] * weight_vec[i];

      /* compute Pmm coefficient */

      for (j=0; j<(2*bw); j++)
	  result[0] += wdata[j] * prev[j];

      /* now generate remaining pmls while computing coefficients */

      for (i=0; i<bw-m-1; i++) {
	ns_vec_mul(ns_L2_cn(m,m+i),prevprev,temp1,2*bw);
	ns_vec_pt_mul(prev, x_i, temp2, 2*bw);
	ns_vec_mul(ns_L2_an(m,m+i), temp2, temp3, 2*bw);
	ns_vec_add(temp3, temp1, temp4, 2*bw); /* temp4 now contains P(m,m+i+1) */
	
	/* compute this coefficient */
	for (j=0; j<(2*bw); j++)
	  result[i+1] += wdata[j] * temp4[j];

	/* now update Pi and P(i+1) */

	memcpy(prevprev, prev, sizeof(double) * 2 * bw);
	memcpy(prev, temp4, sizeof(double) * 2 * bw);

      }
    } /* closes main timing loop */

  if (timing)
    {

      tstop = csecond();
      total_time = tstop - tstart;
      *runtime = total_time;

      fprintf(stdout,"\n");
      fprintf(stdout,"Program: Naive Legendre Transform \n");
      fprintf(stdout,"m = %d\n", m);
      fprintf(stdout,"Bandwidth = %d\n", bw);
#ifndef WALLCLOCK
      fprintf(stdout,"Total elapsed cpu time: %f seconds.\n\n", total_time); 
#else
      fprintf(stdout,"Total elapsed wall time: %f seconds.\n\n", total_time); 
#endif

    }
  else
    *runtime = 0.0;
}



/************************************************************************/
/*
   bw - bandwidth
   m - order
   data - a pointer to double array of size (2*bw) containing a synthesized
          function.  
   result - a pointer to double array of size (bw-m) and containing the
            computed Legendre coefficients, starting with the Pmm
	    coefficient.
   workspace - a pointer to double array of size (32*bw)

   timing - 1 to turn on timing info

   runtime - a double slot for writing time

   loops - number of timing loops


   Just like the above routine EXCEPT I precompute the
   legendre functions before turning on the stopwatch.

*/


void Naive_Analysis_TimingX(double *data,
			    int bw,
			    int m,
			    double *result, 
			    int timing,
			    double *runtime,
			    int loops,
			    double *workspace)
{
  double *prev, *prevprev, *temp1, *temp2, *temp3, *temp4, *x_i, *eval_args;
  double *wdata;
  const double *weight_vec;
  int i, j, l;
  
  double *storeplm, *storeplm_ptr;
  double result0, result1, result2, result3;

  double total_time, tstart, tstop;
  
  prevprev = workspace;
  prev = prevprev + (2*bw);
  temp1 = prev + (2*bw);
  temp2 = temp1 + (2*bw);
  temp3 = temp2 + (2*bw);
  temp4 = temp3 + (2*bw);
  x_i = temp4 + (2*bw);
  eval_args = x_i + (2*bw);
  wdata = eval_args + (2*bw);
  
  storeplm = (double *) malloc(sizeof(double) * 2 * bw *
			       (bw - m));

  storeplm_ptr = storeplm;

  /* get the evaluation nodes */
  ns_EvalPts(2*bw,x_i);
  ns_ArcCosEvalPts(2*bw,eval_args);
  
  /* set initial values of first two Pmls */
  for (i=0; i<2*bw; i++) 
    prevprev[i] = 0.0;
  if (m == 0)
    for (i=0; i<2*bw; i++)
      prev[i] = 0.5;
  else 
    ns_Pmm_L2(m, eval_args, 2*bw, prev);

  memcpy(storeplm, prev, (size_t) sizeof(double) * 2 * bw);

  for(i = 0; i < bw - m - 1; i++)
    {
      ns_vec_mul(ns_L2_cn(m,m+i),prevprev,temp1,2*bw);
      ns_vec_pt_mul(prev, x_i, temp2, 2*bw);
      ns_vec_mul(ns_L2_an(m,m+i), temp2, temp3, 2*bw);
      ns_vec_add(temp3, temp1, temp4, 2*bw); /* temp4 now contains P(m,m+i+1) */
      
      storeplm += (2 * bw);
      memcpy(storeplm, temp4, (size_t) sizeof(double) * 2 * bw);
      memcpy(prevprev, prev, (size_t) sizeof(double) * 2 * bw);
      memcpy(prev, temp4, (size_t) sizeof(double) * 2 * bw);
    }

  storeplm = storeplm_ptr;

  /* start the timer */
  if (timing)
    tstart = csecond();
  
  /* main timing loop */
  for (l=0; l< loops; l++)
    {
      
      /* make sure result is zeroed out */
      for (i=0; i<(bw-m); i++)
	result[i] = 0.0;
      
      /* apply quadrature weights */
      weight_vec = get_weights(bw);
      
      ns_vec_pt_mul(data, weight_vec, wdata, 2 * bw);
      
      storeplm = storeplm_ptr;

      for (i = 0; i < bw - m; i++)
	{
	  result0 = 0.0; result1 = 0.0;
	  result2 = 0.0; result3 = 0.0;

	  for(j = 0; j < (2 * bw) % 4; ++j)
	    result0 += wdata[j] * storeplm[j];
	  for( ; j < (2 * bw); j += 4)
	    {
	      result0 += wdata[j] * storeplm[j];
	      result1 += wdata[j + 1] * storeplm[j + 1];
	      result2 += wdata[j + 2] * storeplm[j + 2];
	      result3 += wdata[j + 3] * storeplm[j + 3];
	    }
	  result[i] = result0 + result1 + result2 + result3;

	  storeplm += (2 * bw);
	}

    } /* closes main timing loop */
  
  if (timing)
    {
      
      tstop = csecond();
      total_time = tstop - tstart;
      *runtime = total_time;
      
      fprintf(stdout,"\n");
      fprintf(stdout,"Program: Naive Legendre Transform \n");
      fprintf(stdout,"m = %d\n", m);
      fprintf(stdout,"Bandwidth = %d\n", bw);
#ifndef WALLCLOCK
      fprintf(stdout,"Total elapsed cpu time: %f seconds.\n\n", total_time); 
#else
      fprintf(stdout,"Total elapsed wall time: %f seconds.\n\n", total_time); 
#endif

    }
  else
    *runtime = 0.0;

  storeplm = storeplm_ptr;

  free(storeplm);
}





/************************************************************************/
/*
   bw - bandwidth
   m - order
   data - a pointer to double array of size (2*bw) containing a synthesized
          function.
   plmtable - a pointer to a double array of size (2*bw*(bw-m));
	      contains the precomputed plms

   result - a pointer to double array of size (bw-m) and containing the
            computed Legendre coefficients, starting with the Pmm
	    coefficient.

   workspace - array of size 2 * bw;

   A minimal, hacked version of the above routine, except that
   as input ones of the arguments is the precomputed table of
   necessary associated Legendre functions.

*/


void Naive_AnalysisX(double *data,
		     int bw,
		     int m,
		     double *result,
		     double *plmtable,
		     double *workspace)
{
  int i, j;
  const double *weight_vec;
  double result0, result1, result2, result3;
  register double *wdata;

  wdata = workspace;

  /* make sure result is zeroed out */
  for (i=0; i<(bw-m); i++)
    result[i] = 0.0;



  /* apply quadrature weights */
  weight_vec = get_weights(bw);
      
  /*  ns_vec_pt_mul(data, weight_vec, wdata, 2 * bw); */

  for(i = 0; i < 2 * bw; i++)
    wdata[i] = data[i] * weight_vec[i];

  for (i = 0; i < bw - m; i++)
	{
	  result0 = 0.0; result1 = 0.0;
	  result2 = 0.0; result3 = 0.0;
	  
	  for(j = 0; j < (2 * bw) % 4; ++j)
	    result0 += wdata[j] * plmtable[j];
	  for( ; j < (2 * bw); j += 4)
	    {
	      result0 += wdata[j] * plmtable[j];
	      result1 += wdata[j + 1] * plmtable[j + 1];
	      result2 += wdata[j + 2] * plmtable[j + 2];
	      result3 += wdata[j + 3] * plmtable[j + 3];
	    }
	  result[i] = result0 + result1 + result2 + result3;

	  plmtable += (2 * bw);
	}

}


/************************************************************************/
/* This is the procedure that synthesizes a function from a list
   of coefficients of a Legendre series.  Function is synthesized
   at the (2*bw) Chebyshev nodes. Associated Legendre functions are
   assumed to be precomputed.
   
   bw - bandwidth
   m - order
   plmtable - precomputed associated Legendre functions
   coeffs - a pointer to double array of size (bw-m).  First coefficient is
            coefficient for Pmm
   result - a pointer to double array of size (2*bw) and containing the
            synthesized function

*/

void Naive_SynthesizeX(double *coeffs,
		       int bw,
		       int m,
		       double *result,
		       double *plmtable)
{
    int i, j;
    double tmpcoef;

    /* make sure result is zeroed out
    for (i=0; i<(2*bw); i++)
      result[i] = 0.0;
      */


    /* add in Pmm contribution */

    tmpcoef = coeffs[0];
    if (tmpcoef != 0.0)
      {
	if (m == 0)
	  for (j=0; j<(2*bw); j++)
	    result[j] = tmpcoef;
	else
	  for (j=0; j<(2*bw); j++)
	    result[j] = tmpcoef * plmtable[j];
      }

    plmtable += ( 2 * bw );

    /* now generate remaining pmls while synthesizing function */
    if (m == 0)
      {
	for (i=0; i<bw-m-1; i++) {
	  /* add in contribution  */
	  tmpcoef = coeffs[i+1] * 2.0;
	  if (tmpcoef != 0.0)
	    {
	      for (j=0; j<(2*bw); j++)
		result[j] += (tmpcoef * plmtable[j]);
	      plmtable += (2 * bw);
	    }
	}
      }
    else
      {
	for (i=0; i<bw-m-1; i++) {
	  /* add in contribution  */
	  tmpcoef = coeffs[i+1];
	  if (tmpcoef != 0.0)
	    {
	      for (j=0; j<(2*bw); j++)
		result[j] += (tmpcoef * plmtable[j]);
	      plmtable += (2 * bw);
	    }
	}
      }



}
