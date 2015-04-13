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


/*  FFT routines for use in fast convolution */

#include <stdio.h>
#include <string.h>  /* for memcpy */

#include "indextables.h"
#include "permroots.h"


#define compmult(a,b,c,d,e,f) (e) = ((a)*(c))-((b)*(d)); (f) = ((a)*(d))+((b)*(c))
#define compadd(a,b,c,d,e,f) (e) = (a) + (c); (f) = (b) + (d)


/************************************************************************/
/* table lookup to compute logs of powers of 2 */
static int ilog2fast(int value)
{
  switch (value)
    {
    case 1: return 0;
    case 2: return 1;
    case 4: return 2;
    case 8: return 3;
    case 16: return 4;
    case 32: return 5;
    case 64: return 6;
    case 128: return 7;
    case 256: return 8;
    case 512: return 9;
    case 1024: return 10;
    case 2048: return 11;
    case 4096: return 12;
    case 8192: return 13;
    case 16384: return 14;
    case 32768: return 15;
    case 65536: return 16;
    default: return -1;
    }
}
/************************************************************************/

/************************************************************************/
/* This is an FFT procedure which has been modified to do the evaluation
   at K nodes given N coefficients.  This is used to speed up convolution
   by zero-striped sequences.  Expects K to be larger than N */
/* arguments

   reals, imags - input double arrays each of size N
   r_vals, i_vals - output double arrays each of size K
   N is the number of coefficients, K the number of samples 

   Assumes that coefficients are ordered from 0 to N-1.
   NOTA BENE!
   When the evaluation is complete, the data is in bit-reverse
   order.  
   Since the purpose of this is to speed up convolution,
   once the series is evaluated it should be left in bit reversed order.
   In this case, brflag should be set to 0.  If you set brflag to
   1, then the results are bit-reversed 

   This evaluation is modeled on polynomial division tree.

   workspace needs to be double array of size 2 * K

*/

void FFTEval( double *reals, double *imags,
	      double *r_vals,double *i_vals,
	      int N,
	      int K,
	      double *workspace,
	      int brflag)
{
  const double *rootptr;
  double *rcptr, *icptr, *rrptr, *irptr;
  double rtemp, itemp;
  int rem2, j, k, l, toggle;

  double tmproot0, tmproot1;

  memcpy(workspace, reals, sizeof(double) * N);
  memcpy(workspace+K, imags, sizeof(double) * N);

  rootptr = r4096; /* point to roots of unity in BR order */
  rem2 = N/2; /* current size of remainders */
  rcptr = workspace;  /* points to real part of current poly coeff */
  icptr = workspace+K;  /* points to imag part of current poly coeff */
  rrptr = r_vals; /* points to real part of current remainder coeff */
  irptr = i_vals; /* points to imag part of current remainder coeff */
  toggle = 0;

  for (l=0; l < ilog2fast(N); l++)
    {
      for (k=0; k<K; k+=rem2) 
	{
	  tmproot0 = rootptr[0]; tmproot1 = rootptr[1];
	  for (j=0; j<rem2; j++) 
	    {
	      compmult(tmproot0,tmproot1,rcptr[rem2+j],icptr[rem2+j],
		       rtemp,itemp);
	      compadd(rtemp,itemp,rcptr[j],icptr[j],
		      rrptr[j],irptr[j]);
	    }
	  rrptr += rem2;
	  irptr += rem2;
	  rootptr += 2;
	  toggle++;
	  if ((toggle % 2) == 0) 
	    {
	      /* check l==0 here. If l==0, then we need to
		 fake it into thinking that there are multiple
		 copies of the input coeffs */
	      if (l != 0) 
		{ 
		  rcptr += 2*rem2;
		  icptr += 2*rem2;
		}
	      else 
		{
		  rcptr = workspace;
		  icptr = workspace+K;
		}
	    }
	} /* closes k loop - done with one level */
    
      if ((l % 2) == 0) 
	{
	  rcptr = r_vals;
	  icptr = i_vals;
	  rrptr = workspace;
	  irptr = workspace + K;
	}
      else 
	{
	  rcptr = workspace;
	  icptr = workspace + K;
	  rrptr = r_vals;
	  irptr = i_vals;
	}
    
      rem2 /= 2;
      rootptr = r4096;
    }
  if ((l % 2) == 0) 
    {
      memcpy(r_vals, workspace, sizeof(double) * K);
      memcpy(i_vals, workspace+K, sizeof(double) * K);
    }
  
  if (brflag == 1) { /* bitreverse the data */
    bitreverse(r_vals,K, workspace);
    bitreverse(i_vals,K, workspace);
  }
  
}

/************************************************************************/
/************************************************************************/
/* This is an Inverse FFT procedure which has been modified to do the 
   interpolation to K coefficients given N samples of a function.
   Expects K <= N, all numbers are powers of 2 */
/* arguments

   reals, imags - input double arrays each of size N
   r_vals, i_vals - output double arrays each of size N - not K -
                    since this double as workspace.

   N is the number of samples, K the number of coefficents

   ALERT!!  ALERT!!  ALERT!!  ALERT!!  ALERT!!  ALERT!!  ALERT!!  
   If brflag == 0, assumes that input samples are in 
   bit-reverse order !!!!!!!!!! If brflag == 1, then
   bit-reverse is performed.
   NOTA BENE!

   Returns a coefficient list ordered from 0 to K-1

   workspace needs to be a double array of size 4 * N 

*/

void FFTInterp( double *reals, double *imags,
		double *r_vals, double *i_vals,
		int N,
		int K,
		double *workspace,
		int brflag )
{
  const double *rootptr;
  double *rhptr, *ihptr, *rlptr, *ilptr;
  double rtemp0, itemp0, rtemp1, itemp1, dn;
  int degree, j, k, l, toggle, upperlim;
  double rptr0, rptr1, rptr2, rptr3;
  double tmprval, tmpival;

  dn = (double) N;

  /* normalize and permute if necessary */
  for (j=0; j<N; j++)
    {
      workspace[j] = reals[j]/dn;
      workspace[j+N] = imags[j]/dn;
    }
  if (brflag == 1)
    {
      bitreverse(workspace, N, workspace+(2*N));
      bitreverse(workspace+N, N, workspace+(2*N));
    }

  rootptr = r4096; /* point to roots of unity in BR order */
  degree = 1; /* next order of polynomials - stride length */
  rlptr = workspace;  /* points to real part of lower order polys */
  ilptr = workspace+N;  /* points to imag part of lower order polys */
  rhptr = r_vals; /* points to real part of higher order polys */
  ihptr = i_vals; /* points to imag part of higher order polys */
  toggle = 0;

  for (l=1; l <= ilog2fast(K); l++)
    {
      upperlim = N/(2*degree);
      for (k=0; k < upperlim; k++)
	{
	  for (j=0; j < degree; j++)
	    {
	      compadd(rlptr[j],ilptr[j],rlptr[j+degree],ilptr[j+degree],
		      rhptr[j],ihptr[j]);
	    }
	  rhptr += degree;
	  ihptr += degree;

	  rptr0 = rootptr[0]; rptr1 = rootptr[1];
	  rptr2 = rootptr[2]; rptr3 = rootptr[3];

	  for (j=0; j < degree; j++)
	    {

	      compmult(rptr0,-rptr1,rlptr[j],ilptr[j],
		       rtemp0,itemp0);
	      compmult(rptr2,-rptr3,rlptr[j+degree],ilptr[j+degree],
		       rtemp1,itemp1);
	      compadd(rtemp0,itemp0,rtemp1,itemp1,
		      rhptr[j],ihptr[j]);

	    }

	  rhptr += degree;
	  ihptr += degree;
	  rlptr += (2*degree);
	  ilptr += (2*degree);
	  rootptr += 4;
	} /* end of k loop */
    
      if ((l % 2) != 0)
	{
	  rlptr = r_vals;
	  ilptr = i_vals;
	  rhptr = workspace;
	  ihptr = workspace + N;
	}
      else
	{
	  rhptr = r_vals;
	  ihptr = i_vals;
	  rlptr = workspace;
	  ilptr = workspace + N;
	}
      degree *= 2;
      rootptr = r4096;
      toggle++;
    
    } /* end of l loop */
  
  if ((toggle % 2) == 0)
    {
      memcpy(r_vals, workspace, sizeof(double) * N);
      memcpy(i_vals, workspace+N, sizeof(double) * N);
    }
  
  for (j=0; j<K; j++)
    {
      tmprval = r_vals[j]; tmpival = i_vals[j];
      for (k=K; k<N; k+=K)
	{
	  tmprval += r_vals[j+k];
	  tmpival += i_vals[j+k];
	}
      r_vals[j] = tmprval; i_vals[j] = tmpival;
    }
  
}
  


/************************************************************************/
