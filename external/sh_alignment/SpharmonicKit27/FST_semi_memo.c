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


/********************************************************************

  FST_semi_memo.c - routines to perform convolutions on the
  2-sphere using a combination of semi-naive and naive algorithms.

  ASSUMES THAT ALL PRECOMPUTED-DATA IS IN MEMORY, AND NOT TO BE
  READ FROM THE DISK.

  The primary functions in this package are

  1) FST_semi_memo() - computes the spherical harmonic expansion.
  2) InvFST_semi_memo() - computes the inverse spherical harmonic transform.
  3) FZT_semi_memo() - computes the zonal harmonic transform.
  4) TransMult() - Multiplies harmonic coefficients using Driscoll-Healy
                    result.  Dual of convolution in "time" domain.
  5) Conv2Sphere_semi_memo() - Convolves two functins defined on the 2-sphere,
                          using seminaive transform

  For descriptions on calling these functions, see the documentation
  preceding each function.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cospmls.h"
#include "fft_grids.h"
#include "naive_synthesis.h"
#include "primitive.h"
#include "primitive_FST.h"
#include "seminaive.h"

#define compmult(a,b,c,d,e,f) (e) = ((a)*(c))-((b)*(d)); (f) = ((a)*(d))+((b)*(c))

/************************************************************************/

#ifndef FFTPACK

/************************************************************************/
/* performs a spherical harmonic transform using the semi-naive
   and naive algorithms */
/* size is the dimension of the input array (size x size) and it is
   expected that size=2*bw.  The inputs rdata and idata are expected
   to be pointers to size x size arrays.  rcoeffs and icoeffs are expected
   to be pointers to bw x bw arrays, and will contain the harmonic
   coefficients in a "linearized" form.


   spharmonic_pml_table should be a (double **) pointer to
   the result of a call to Spharmonic_Pml_Table.  Because this
   table is re-used in the inverse transform, and because for
   timing purposes the computation of the table is not included,
   it is passed in as an argument.  Also, at some point this
   code may be used as par of a series of convolutions, so
   reducing repetitive computation is prioritized.

   spharmonic_pml_table will be an array of (double *) pointers
   the array being of length TableSize(m,bw)

   workspace needs to be a double pointer to an array of size
   (8 * bw^2) + (29 * bw).

   cutoff -> what order to switch from semi-naive to naive
             algorithm.

*/

/* 
   Output Ordering of coeffs f(m,l) is
   f(0,0) f(0,1) f(0,2) ... f(0,bw-1)
          f(1,1) f(1,2) ... f(1,bw-1)
          etc.
                 f(bw-2,bw-2), f(bw-2,bw-1)
		               f(bw-1,bw-1)
			       f(-(bw-1),bw-1)
		 f(-(bw-2),bw-2) f(-(bw-2),bw-1)
	  etc.
	          f(-2,2) ... f(-2,bw-1)
	  f(-1,1) f(-1,2) ... f(-1,bw-1)
    
   This only requires an array of size (bw*bw).  If zero-padding
   is used to make the indexing nice, then you need a an
   (2bw-1) * bw array - but that is not done here.
   Because of the amount of space necessary for doing
   large transforms, it is important not to use any
   more than necessary.

*/

/*      dataformat =0 -> samples are complex, =1 -> samples real */

void FST_semi_memo(double *rdata, double *idata,
		   double *rcoeffs, double *icoeffs, 
		   int size, double **seminaive_naive_table,
		   double *workspace,
		   int dataformat,
		   int cutoff)
{
  int bw, m, i, j;
  double *rres, *ires;
  double *rdataptr, *idataptr;
  double *fltres, *scratchpad;
  double *eval_pts;
  int tmpindex[2];
  double pow_one;
  double *cos_even;
  int l, dummy ;
  double tmpA, tmpB ;

  bw = size/2;

  /* assign space */

  cos_even = (double *) malloc(sizeof(double) * bw);

  rres = workspace;  /* needs (size * size) = (4 * bw^2) */
  ires = rres + (size * size); /* needs (size * size) = (4 * bw^2) */ 
  fltres = ires + (size * size); /* needs bw  */
  eval_pts = fltres + bw; /* needs (2*bw)  */
  scratchpad = eval_pts + (2*bw); /* needs (24 * bw)  */
 

  /* total workspace is (8 * bw^2) + (29 * bw) */

  /* do the FFTs along phi */
  grid_fourier(rdata, idata, rres, ires, size, scratchpad);
  
  /* point to start of output data buffers */
  rdataptr = rcoeffs;
  idataptr = icoeffs;

  for (m=0; m<bw; m++) {
    /*
      fprintf(stderr,"m = %d\n",m);
      */
    
    /*** test to see if before cutoff or after ***/
    if (m < cutoff){
      
      /* do the real part */
      SemiNaiveReduced(rres+(m*size), 
		       bw, 
		       m, 
		       fltres, 
		       seminaive_naive_table[m],
		       scratchpad,
		       cos_even);
      
      /* now load real part of coefficients into output space */  
      memcpy(rdataptr, fltres, sizeof(double) * (bw - m));
      
      rdataptr += bw-m;
      
      /* do imaginary part */
      SemiNaiveReduced(ires+(m*size),
		       bw, 
		       m, 
		       fltres, 
		       seminaive_naive_table[m],
		       scratchpad,
		       cos_even);

      /* now load imaginary part of coefficients into output space */  
      memcpy(idataptr, fltres, sizeof(double) * (bw - m));
      
      idataptr += bw-m;
      
    }
    else{
      /* do real part */
      
      Naive_AnalysisX(rres+(m*size),
		      bw,
		      m,
		      fltres,
		      seminaive_naive_table[m],
		      workspace);
      memcpy(rdataptr, fltres, sizeof(double) * (bw - m));
      rdataptr += bw-m;
      
      /* do imaginary part */
      Naive_AnalysisX(ires+(m*size),
		      bw,
		      m,
		      fltres,
		      seminaive_naive_table[m],
		      workspace);
      memcpy(idataptr, fltres, sizeof(double) * (bw - m));
      idataptr += bw-m;
    }
    
    
  }

  /*** now do upper coefficients ****/
  
  /* now if the data is real, we don't have to compute the
     coefficients whose order is less than 0, i.e. since
     the data is real, we know that
     f-hat(l,-m) = (-1)^m * conjugate(f-hat(l,m)),
     so use that to get the rest of the coefficients
     
     dataformat =0 -> samples are complex, =1 -> samples real
     
     */
  
  if( dataformat == 0 ){
    
    /* note that m is greater than bw here, but this is for
       purposes of indexing the input data arrays.  
       The "true" value of m as a parameter for Pml is
       size - m  */
    /*   fprintf(stderr,"\n now the higher order terms \n\n"); */
    
    for (m=bw+1; m<size; m++) {
      /*
	fprintf(stderr,"m = %d\n",-(size-m));
	*/
      
      /* do real part */
      SemiNaiveReduced(rres+(m*size), 
		       bw, 
		       size-m, 
		       fltres, 
		       seminaive_naive_table[size-m],
		       scratchpad,
		       cos_even);
      
      /* now load real part of coefficients into output space */  
      if ((m % 2) != 0) {
	for (i=0; i<m-bw; i++)
	  rdataptr[i] = -fltres[i];
      }
      else {
	memcpy(rdataptr, fltres, sizeof(double) * (m - bw));
      }
      rdataptr += m-bw;
      
      /* do imaginary part */
      SemiNaiveReduced(ires+(m*size), 
		       bw, 
		       size-m, 
		       fltres, 
		       seminaive_naive_table[size-m],
		       scratchpad,
		       cos_even);
      
      /* now load real part of coefficients into output space */  
      if ((m % 2) != 0) {
	for (i=0; i<m-bw; i++)
	  idataptr[i] = -fltres[i];
      }
      else {
	memcpy(idataptr, fltres, sizeof(double) * (m - bw));
      }
      idataptr += m-bw;
      
    }
  }
  else        /**** if the data is real ****/
    {            
      pow_one = 1.0;
      for(i = 1; i < bw; i++){
	pow_one *= -1.0;
	for( j = i; j < bw; j++){	
	  seanindex2(i, j, bw, tmpindex);
	  rcoeffs[tmpindex[1]] =
	    pow_one * rcoeffs[tmpindex[0]];
	  icoeffs[tmpindex[1]] =
	    -1.0 * pow_one * icoeffs[tmpindex[0]];
	}
      }
    }
  
  free(cos_even);

  /*****************

  New in 2.6: Need to massage the coefficients one more time,
              to get that right factor in there (how embarrassing)

  ******************/

  tmpA = 2. * sqrt( PI );
  tmpB = sqrt( 2. * PI );

  for(m=0;m<bw;m++)
    for(l=m;l<bw;l++)
      {
	dummy = seanindex(m,l,bw);
	if ( m == 0 )
	  {
	    rcoeffs[dummy] *= tmpA ;
	    icoeffs[dummy] *= tmpA ;
	  }
	else
	  {
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	  }
	/* now for the negative-order coefficients */
	if ( m != 0 )
	  {
	    dummy = seanindex(-m,l,bw);
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	    
	  }
      }

  return;
}

/************************************************************************/
/* Inverse spherical harmonic transform.  Inputs rcoeffs and icoeffs
   are harmonic coefficients stored in (bw * bw) arrays in the order
   spec'ed above.  rdata and idata are (size x size) arrays with
   the transformed result.  size is expected to be 2 * bw.
   transpose_spharmonic_pml_table should be the (double **) result of a call
   to Transpose_Spharmonic_Pml_Table()

   workspace is (8 * bw^2) + (32 * bw)

*/

/*      dataformat =0 -> samples are complex, =1 -> samples real */

void InvFST_semi_memo(double *rcoeffs, double *icoeffs, 
		      double *rdata, double *idata, 
		      int size, 
		      double **transpose_seminaive_naive_table,
		      double *workspace,
		      int dataformat,
		      int cutoff)
{
  int bw, m, i, n;
  double *rdataptr, *idataptr;
  double *rfourdata, *ifourdata;
  double *rinvfltres, *iminvfltres, *scratchpad;
  double *sin_values, *eval_pts;
  int l, dummy ;
  double tmpA, tmpB ;
  FILE *fp;

  bw = size/2;

  /* allocate space */

  rfourdata = workspace;                  /* needs (size * size) */
  ifourdata = rfourdata + (size * size);  /* needs (size * size) */
  rinvfltres = ifourdata + (size * size); /* needs (2 * bw) */
  iminvfltres = rinvfltres + (2 * bw);    /* needs (2 * bw) */
  sin_values = iminvfltres + (2 * bw);    /* needs (2 * bw) */
  eval_pts = sin_values + (2 * bw);       /* needs (2 * bw) */
  scratchpad = eval_pts + (2 * bw);       /* needs (24 * bw) */
  
  /* total workspace = (8 * bw^2) + (32 * bw) */

  /* load up the sin_values array */
  n = 2*bw;

  ArcCosEvalPts(n, eval_pts);
  for (i=0; i<n; i++)
    sin_values[i] = sin(eval_pts[i]);

  /**********************

  New in 2.6: Need to massage the coefficients, to get that
              right factor in there (how embarrassing)

  ***********************/

  tmpA = 1./(2. * sqrt(PI) );
  tmpB = 1./sqrt(2. * PI ) ;

  for(m=0;m<bw;m++)
    for(l=m;l<bw;l++)
      {
	dummy = seanindex(m,l,bw);
	if ( m == 0 )
	  {
	    rcoeffs[dummy] *= tmpA ;
	    icoeffs[dummy] *= tmpA ;
	  }
	else
	  {
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	  }
	/* now for the negative-order coefficients */
	if ( m != 0 )
	  {
	    dummy = seanindex(-m,l,bw);
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	  }
      }

  /* Now do all of the inverse Legendre transforms */
  rdataptr = rcoeffs;
  idataptr = icoeffs;
  for (m=0; m<bw; m++)
    {
      /*
	fprintf(stderr,"m = %d\n",m);
	*/
    
      if(m < cutoff)
	{
      
	  /* do real part first */ 
	  InvSemiNaiveReduced(rdataptr,
			      bw,
			      m,
			      rinvfltres,
			      transpose_seminaive_naive_table[m],
			      sin_values,
			      scratchpad);
      
	  /* now do imaginary part */
      
	  InvSemiNaiveReduced(idataptr,
			      bw,
			      m,
			      iminvfltres,
			      transpose_seminaive_naive_table[m],
			      sin_values,
			      scratchpad);
      
	  /* will store normal, then tranpose before doing inverse fft */
	  memcpy(rfourdata+(m*size), rinvfltres, sizeof(double) * size);
	  memcpy(ifourdata+(m*size), iminvfltres, sizeof(double) * size);
      
	  /* move to next set of coeffs */
      
	  rdataptr += bw-m;
	  idataptr += bw-m;
      
	}
      else
	{
	  /* first do the real part */
	
	  Naive_SynthesizeX(rdataptr,
			    bw,
			    m,
			    rinvfltres,
			    transpose_seminaive_naive_table[m]);

	  /* now do the imaginary */	
	  Naive_SynthesizeX(idataptr,
			    bw,
			    m,
			    iminvfltres,
			    transpose_seminaive_naive_table[m]);

	  /* will store normal, then tranpose before doing inverse fft    */
	  memcpy(rfourdata+(m*size), rinvfltres, sizeof(double) * size);
	  memcpy(ifourdata+(m*size), iminvfltres, sizeof(double) * size);
 	
	  /* move to next set of coeffs */
	
	  rdataptr += bw-m;
	  idataptr += bw-m;
	
	}
    
    } /* closes m loop */
  
  
  /* now fill in zero values where m = bw (from problem definition) */
  memset(rfourdata + (bw * size), 0, sizeof(double) * size);
  memset(ifourdata + (bw * size), 0, sizeof(double) * size);
  
  /* now if the data is real, we don't have to compute the
     coefficients whose order is less than 0, i.e. since
     the data is real, we know that
     invf-hat(l,-m) = conjugate(invf-hat(l,m)),
     so use that to get the rest of the real data
     
     dataformat =0 -> samples are complex, =1 -> samples real
     
     */

  if(dataformat == 0){
    
    /* now do negative m values */
    
    for (m=bw+1; m<size; m++) 
      {
	/*
	  fprintf(stderr,"m = %d\n",-(size-m));
	  */
	
	/* do real part first */
	
	InvSemiNaiveReduced(rdataptr,
			    bw,
			    size - m,
			    rinvfltres,
			    transpose_seminaive_naive_table[size - m],
			    sin_values,
			    scratchpad);
	
	/* now do imaginary part */
	
	InvSemiNaiveReduced(idataptr,
			    bw,
			    size - m,
			    iminvfltres,
			    transpose_seminaive_naive_table[size - m],
			    sin_values,
			    scratchpad);
	

	/* will store normal, then tranpose before doing inverse fft    */
	if ((m % 2) != 0)
	  for(i=0; i< size; i++){
	    rinvfltres[i] = -rinvfltres[i];
	    iminvfltres[i] = -iminvfltres[i];
	  }
	
	memcpy(rfourdata + (m*size), rinvfltres, sizeof(double) * size);
	memcpy(ifourdata + (m*size), iminvfltres, sizeof(double) * size);
	
	
	/* move to next set of coeffs */
	rdataptr += bw-(size-m);
	idataptr += bw-(size-m);
	
      } /* closes m loop */
  }
  else {
    for(m = bw + 1; m < size; m++){
      
      memcpy(rfourdata+(m*size), rfourdata+((size-m)*size),
	     sizeof(double) * size);
      memcpy(ifourdata+(m*size), ifourdata+((size-m)*size),
	     sizeof(double) * size);
      for(i = 0; i < size; i++)
	ifourdata[(m*size)+i] *= -1.0;
    }
    
    
  }
  
  /** now transpose **/
  transpose(rfourdata, size);
  transpose(ifourdata, size);
  
  /* now do inverse fourier grid computation */
  grid_invfourier(rfourdata, ifourdata, rdata, idata, size, scratchpad);
  
  /* amscray */
  
}

/************************************************************************/
/* Zonal Harmonic transform using seminaive algorithm - used in convolutions */
/* rdata and idata should be pointers to size x size arrays.
   rres and ires should be pointers to double arrays of size bw.
   size = 2 * bw
   cos_pml_table contains Legendre coefficients of P(0,l) functions
   and is result of CosPmlTableGen for m = 0;
   FZT_semi only computes spherical harmonics for m=0.

   dataformat =0 -> samples are complex, =1 -> samples real
   
   workspace needed is (12 * bw) 

*/
void FZT_semi_memo(double *rdata, double *idata,
		   double *rres, double *ires,
		   int size,
		   double *cos_pml_table,
		   double *workspace,
		   int dataformat)
{

  int i, j, bw;
  double *r0, *i0, dsize;
  double tmpreal, tmpimag;
  double *scratchpad;
  double *cos_even;
  int l, dummy ;
  double tmpA ;

  bw = size/2;

  /* assign memory */
  r0 = workspace;             /* needs (2 * bw) */
  i0 = r0 + (2 * bw);         /* needs (2 * bw) */
  cos_even = i0 + (2 * bw);   /* needs (1 * bw) */
  scratchpad = cos_even + (1 * bw); /* needs (8 * bw) */
                              /* total workspace = 13*bw */


  dsize = (double) size;
  /* compute the m = 0 components */

  for (i=0; i<size; i++) {
    tmpreal = 0.0;
    tmpimag = 0.0;
    for(j=0; j<size; j++) {
      tmpreal += rdata[(i*size)+j];
      tmpimag += idata[(i*size)+j];
    }
    /* normalize */
    r0[i] = tmpreal/dsize;
    i0[i] = tmpimag/dsize;
  }

  /* do the real part */
  SemiNaiveReduced(r0, 
		   bw, 
		   0,
		   rres, 
		   cos_pml_table,
		   scratchpad,
		   cos_even);
   
  if(dataformat == 0)   /* do imaginary part */
    SemiNaiveReduced(i0,
		     bw, 
		     0, 
		     ires, 
		     cos_pml_table,
		     scratchpad,
		     cos_even);
  else                 /* otherwise set coefficients = 0 */
    memset(ires, 0, sizeof(double) * size);

  /*****************

  New in 2.6: Need to massage the coefficients one more time,
              to get that right factor in there (how embarrassing)

  ******************/

  tmpA = 2. * sqrt( PI );

  for(l=0;l<bw;l++)
    {
      dummy = seanindex(0,l,bw);

      rres[dummy] *= tmpA ;
      ires[dummy] *= tmpA ;
    }

}
/************************************************************************/
/* multiplies harmonic coefficients of a function and a filter */
/* See convolution theorem of Driscoll and Healy */
/* datacoeffs should be output of an FST, filtercoeffs the 
   output of an FZT.  There should be (bw * bw) datacoeffs,
   and bw filtercoeffs.
   rres and ires should point to arrays
   of dimension bw * bw. size parameter is 2*bw  
*/

void TransMult(double *rdatacoeffs, double *idatacoeffs, 
	       double *rfiltercoeffs, double *ifiltercoeffs, 
	       double *rres, double *ires,
	       int size)
{

  int m, l, bw;
  double *rdptr, *idptr, *rrptr, *irptr;

  bw = size/2;

  rdptr = rdatacoeffs;
  idptr = idatacoeffs;
  rrptr = rres;
  irptr = ires;

  for (m=0; m<bw; m++) {
    for (l=m; l<bw; l++) {
      compmult(rfiltercoeffs[l], ifiltercoeffs[l],
	       rdptr[l-m], idptr[l-m],
	       rrptr[l-m], irptr[l-m]);
      rrptr[l-m] *= sqrt(4*M_PI/(2*l+1));
      irptr[l-m] *= sqrt(4*M_PI/(2*l+1));

    }
    rdptr += bw-m; idptr += bw-m;
    rrptr += bw-m; irptr += bw-m;
  }
  for (m=bw+1; m<size; m++) {
    for (l=size-m; l<bw; l++){
      compmult(rfiltercoeffs[l], ifiltercoeffs[l],
	       rdptr[l-size+m], idptr[l-size+m],
	       rrptr[l-size+m], irptr[l-size+m]);

      rrptr[l-size+m] *= sqrt(4*M_PI/(2*l+1));
      irptr[l-size+m] *= sqrt(4*M_PI/(2*l+1));

    }
    rdptr += m-bw; idptr += m-bw;
    rrptr += m-bw; irptr += m-bw;
  }

}
/************************************************************************/
/* Here's the big banana
   Convolves two functions defined on the 2-sphere.
   Uses seminaive algorithms for spherical harmonic transforms
   size = 2*bw
   Inputs:
   rdata, idata - (size * size) arrays containing real and
                  imaginary parts of sampled function.
   rfilter, ifilter - (size * size) arrays containing real and
                      imaginary parts of sampled filter function.
   rres, ires - (size * size) arrays containing real and
                  imaginary parts of result function.
   size - should be 2 * bw

   Suggestion - if you want to do multiple convolutions,
   don't keep allocating and freeing space with every call,
   or keep recomputing the spharmonic_pml tables.
   Allocate workspace once before you call this function, then
   just set up pointers as first step of this procedure rather
   than mallocing.  And do the same with the FST, FZT, and InvFST functions.

   Memory requirements for Conv2Sphere

   Need space for spharmonic tables and local workspace and
   scratchpad space for FST_semi

   2 * (Spharmonic_TableSize(bw))  +
   8 * (bw*bw)  +
   2 * bw +
   8 * (bw*bw) + 32 * bw
   

   ASSUMPTIONS:
   1. data is strictly REAL
   2. will do semi-naive algorithm for ALL orders

*/
void Conv2Sphere_semi_memo(double *rdata, double *idata, 
			   double *rfilter, double *ifilter, 
			   double *rres, double *ires, 
			   int size,
			   double *workspace)
{
  int bw, spharmonic_bound;
  double *frres, *fires, *filtrres, *filtires, *trres, *tires;
  double **spharmonic_pml_table, **transpose_spharmonic_pml_table;
  double *spharmonic_result_space, *transpose_spharmonic_result_space;
  double *scratchpad;

  bw = size/2;

  /* assign space */

  spharmonic_bound = Spharmonic_TableSize(bw);

  spharmonic_result_space = workspace;
  transpose_spharmonic_result_space = spharmonic_result_space + 
                                      spharmonic_bound;
  frres = transpose_spharmonic_result_space +
          spharmonic_bound;   /* needs (bw*bw) */
  fires = frres + (bw*bw);    /* needs (bw*bw) */
  trres = fires + (bw*bw);    /* needs (bw*bw) */
  tires = trres + (bw*bw);    /* needs (bw*bw) */
  filtrres = tires + (bw*bw); /* needs bw */
  filtires = filtrres + bw;   /* needs bw */
  scratchpad = filtires + bw; /* needs (8 * bw^2) + (32 * bw) */

  spharmonic_pml_table = 
    Spharmonic_Pml_Table(bw,
			 spharmonic_result_space,
			 scratchpad);

  transpose_spharmonic_pml_table = 
    Transpose_Spharmonic_Pml_Table(spharmonic_pml_table, 
				   bw,
				   transpose_spharmonic_result_space,
				   scratchpad);
  FST_semi_memo(rdata, idata, 
		frres, fires, 
		size, 
		spharmonic_pml_table,
		scratchpad,
		1,
		bw);
  
  FZT_semi_memo(rfilter, ifilter, 
		filtrres, filtires, 
		size, 
		spharmonic_pml_table[0],
		scratchpad,
		1);

  TransMult(frres, fires, filtrres, filtires, trres, tires, size);

  InvFST_semi_memo(trres, tires, 
		   rres, ires, 
		   size,
		   transpose_spharmonic_pml_table,
		   scratchpad,
		   1,
		   bw);

  /***
    have to free the memory that was allocated in
    Spharmonic_Pml_Table() and
    Transpose_Spharmonic_Pml_Table()
    ***/

  free(spharmonic_pml_table);

  free(transpose_spharmonic_pml_table);
}

#else /* will use FFTPACK */

/*********************************************/
/*********************************************/
/*********************************************/
/*********************************************/
/*********************************************/




/************************************************************************/
/* performs a spherical harmonic transform using the semi-naive
   and naive algorithms */
/* size is the dimension of the input array (size x size) and it is
   expected that size=2*bw.  The inputs rdata and idata are expected
   to be pointers to size x size arrays.  rcoeffs and icoeffs are expected
   to be pointers to bw x bw arrays, and will contain the harmonic
   coefficients in a "linearized" form.


   spharmonic_pml_table should be a (double **) pointer to
   the result of a call to Spharmonic_Pml_Table.  Because this
   table is re-used in the inverse transform, and because for
   timing purposes the computation of the table is not included,
   it is passed in as an argument.  Also, at some point this
   code may be used as par of a series of convolutions, so
   reducing repetitive computation is prioritized.

   spharmonic_pml_table will be an array of (double *) pointers
   the array being of length TableSize(m,bw)

   workspace needs to be a double pointer to an array of size
   (8 * bw^2) + (29 * bw).

   cutoff -> what order to switch from semi-naive to naive
             algorithm.

*/

/* 
   Output Ordering of coeffs f(m,l) is
   f(0,0) f(0,1) f(0,2) ... f(0,bw-1)
          f(1,1) f(1,2) ... f(1,bw-1)
          etc.
                 f(bw-2,bw-2), f(bw-2,bw-1)
		               f(bw-1,bw-1)
			       f(-(bw-1),bw-1)
		 f(-(bw-2),bw-2) f(-(bw-2),bw-1)
	  etc.
	          f(-2,2) ... f(-2,bw-1)
	  f(-1,1) f(-1,2) ... f(-1,bw-1)
    
   This only requires an array of size (bw*bw).  If zero-padding
   is used to make the indexing nice, then you need a an
   (2bw-1) * bw array - but that is not done here.
   Because of the amount of space necessary for doing
   large transforms, it is important not to use any
   more than necessary.

*/

/*      dataformat =0 -> samples are complex, =1 -> samples real */

/* 	  input data is strictly real  */
/****************************************************/

void FST_semi_memo(double *rdata,
		   double *rcoeffs,
		   double *icoeffs, 
		   int size,
		   double **seminaive_naive_table,
		   double *workspace,
		   double *wSave,
		   double *CoswSave,
		   int cutoff)
{
  int bw, m, i, j;
  double *rres, *ires;
  double *rdataptr, *idataptr;
  double *scratchpad;
  int tmpindex[2];
  double pow_one;
  double *cos_even;
  double *tmprdata_ptr;
  int l, dummy ;
  double tmpA, tmpB ;

  bw = size/2;

  /* assign space */
  rres = workspace;  /* needs (size * size) = (4 * bw^2) */
  cos_even = rres + (size * size); /* needs bw  */
  scratchpad = cos_even + (2*bw); /* needs (24 * bw)  ?? */
 
  /* total workspace is (4 * bw^2) + (29 * bw) */

  /* first save input data: don't want to write over it! */
  memcpy(rres, rdata, sizeof(double) * size * size);
  tmprdata_ptr = rres; /* to be used for advancing */
  
  /* do the FFTs along phi */
  grid_fourier_FFTPACK(rres, size, wSave);
  
  /* point to start of output data buffers */
  rdataptr = rcoeffs;
  idataptr = icoeffs;
  
  for (m=0; m<bw; m++) {
    /*
      fprintf(stderr,"m = %d\n",m);
      */
    
    /*** test to see if before cutoff or after ***/
    if (m < cutoff){
      
      /* do the real part */
      SemiNaiveReduced(tmprdata_ptr, 
		       bw, 
		       m, 
		       rdataptr, 
		       seminaive_naive_table[m],
		       scratchpad,
		       cos_even,
		       CoswSave );
      
      rdataptr += bw-m;

      tmprdata_ptr += size;
      
      /* do ONLY IF m != 0; otherwise set imaginary coeffs = 0 */
      if (m != 0)
	{
	  /* do imaginary part */
	  SemiNaiveReduced(tmprdata_ptr, 
			   bw, 
			   m, 
			   idataptr, 
			   seminaive_naive_table[m],
			   scratchpad,
			   cos_even,
			   CoswSave );
      
	  idataptr += bw-m;
	  tmprdata_ptr += size;
	}
      else
	{
	  memset(idataptr, 0, sizeof(double) * bw);
	  idataptr += bw;
	}
      
    }
    else{
      /* do real part */
      
      Naive_AnalysisX(tmprdata_ptr,
		      bw,
		      m,
		      rdataptr,
		      seminaive_naive_table[m],
		      workspace);

      rdataptr += bw-m;
      tmprdata_ptr += size;
  
      /* do ONLY IF m != 0; otherwise set imaginary coeffs = 0 */
      if (m != 0)
	{
	  /* do imaginary part */
	  Naive_AnalysisX(tmprdata_ptr,
			  bw,
			  m,
			  idataptr,
			  seminaive_naive_table[m],
			  workspace);

	  idataptr += bw-m;
	  tmprdata_ptr += size;
	}
      else
	{
	  memset(idataptr, 0, sizeof(double) * bw);
	  idataptr += bw;
	}
	
    }
    
  }
  
  /*** now do upper coefficients ****/
  
  pow_one = 1.0;
  for(i = 1; i < bw; i++)
    {
      pow_one *= -1.0;
      for( j = i; j < bw; j++)
	{	
	  seanindex2(i, j, bw, tmpindex);
	  rcoeffs[tmpindex[1]] = pow_one * rcoeffs[tmpindex[0]];
	  icoeffs[tmpindex[1]] = -1.0 * pow_one * icoeffs[tmpindex[0]];
	}
    }

  /*****************

  New in 2.6: Need to massage the coefficients one more time,
              to get that right factor in there (how embarrassing)

  ******************/

  tmpA = 2. * sqrt( PI );
  tmpB = sqrt( 2. * PI );

  for(m=0;m<bw;m++)
    for(l=m;l<bw;l++)
      {
	dummy = seanindex(m,l,bw);
	if ( m == 0 )
	  {
	    rcoeffs[dummy] *= tmpA ;
	    icoeffs[dummy] *= tmpA ;
	  }
	else
	  {
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	  }
	/* now for the negative-order coefficients */
	if ( m != 0 )
	  {
	    dummy = seanindex(-m,l,bw);
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	    
	  }
      }

  
}

/************************************************************************/
/* Inverse spherical harmonic transform.  Inputs rcoeffs and icoeffs
   are harmonic coefficients stored in (bw * bw) arrays in the order
   spec'ed above.  rdata and idata are (size x size) arrays with
   the transformed result.  size is expected to be 2 * bw.
   transpose_spharmonic_pml_table should be the (double **) result of a call
   to Transpose_Spharmonic_Pml_Table()

   workspace is (20 * bw)

*/


/* input data is strictly real  */

void InvFST_semi_memo(double *rcoeffs,
		      double *icoeffs, 
		      double *rdata,
		      int size, 
		      double **transpose_seminaive_naive_table,
		      double *workspace,
		      double *wSave,
		      int cutoff)
{
  int bw, m, i, n;
  double *rdataptr, *idataptr;
  double *scratchpad;
  double *sin_values, *eval_pts;
  double *tmp_ptr;
  int l, dummy ;
  double tmpA, tmpB ;
  FILE *fp;


  bw = size/2;

  /* allocate space */

  sin_values = workspace;                  /* needs (2 * bw) */
  eval_pts = sin_values + (2 * bw);       /* needs (2 * bw) */
  scratchpad = eval_pts + (2 * bw);       /* needs (16 * bw) */

  /* load up the sin_values array */
  n = 2*bw;
  ArcCosEvalPts(n, eval_pts);
  for (i=0; i<n; i++)
    sin_values[i] = sin(eval_pts[i]);

  /*****************

  New in 2.6: Need to massage the coefficients one more time,
              to get that right factor in there (how embarrassing)

  ******************/

  tmpA = 1./(2. * sqrt(PI) );
  tmpB = 1./sqrt(2. * PI ) ;

  for(m=0;m<bw;m++)
    for(l=m;l<bw;l++)
      {
	dummy = seanindex(m,l,bw);
	if ( m == 0 )
	  {
	    rcoeffs[dummy] *= tmpA ;
	    icoeffs[dummy] *= tmpA ;
	  }
	else
	  {
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	  }
	/* now for the negative-order coefficients */
	if ( m != 0 )
	  {
	    dummy = seanindex(-m,l,bw);
	    rcoeffs[dummy] *= tmpB ;
	    icoeffs[dummy] *= tmpB ;
	  }
      }


  /* Now do all of the inverse Legendre transforms */

  tmp_ptr = rdata;

  rdataptr = rcoeffs;
  idataptr = icoeffs;
  for (m=0; m<bw; m++) {
    /*
      fprintf(stderr,"m = %d\n",m);
      */
    
    if(m < cutoff){
      
      /* do real part first */ 
      
      InvSemiNaiveReduced(rdataptr,
			  bw,
			  m,
			  tmp_ptr,
			  transpose_seminaive_naive_table[m],
			  sin_values,
			  scratchpad);
      tmp_ptr += size;

      /* move to next set of coeffs */
      rdataptr += bw - m;

      if ( m != 0)
	{
	  /* now do imaginary part */
	  InvSemiNaiveReduced(idataptr,
			      bw,
			      m,
			      tmp_ptr,
			      transpose_seminaive_naive_table[m],
			      sin_values,
			      scratchpad);
	  tmp_ptr += size;
	  /* move to next set of coefs */
	  idataptr += bw - m;
	}
      else
	{
	  /* move to next set of coeffs */
	  idataptr += bw - m;
	}
    }
    else
      {
	/* first do the real part */
	
	Naive_SynthesizeX(rdataptr,
			  bw,
			  m,
			  tmp_ptr,
			  transpose_seminaive_naive_table[m]);
	tmp_ptr += size;
	/* move to next set of coeffs */
	rdataptr += bw - m;


	if ( m!= 0 )
	  {
	    
	    /* now do the imaginary */	
	    Naive_SynthesizeX(idataptr,
			      bw,
			      m,
			      tmp_ptr,
			      transpose_seminaive_naive_table[m]);
	    tmp_ptr += size;
	    /* move to next set of coefs */
	    idataptr += bw - m;
	  }
	else
	  {
	    /* move to next set of coeffs */
	    idataptr += bw - m;
	  }
      }
    
  } /* closes m loop */
  
  
  /* now fill in zero values where m = bw (from problem definition) */
  memset(tmp_ptr, 0, sizeof(double) * size);

  /** now transpose **/
  transpose(rdata, size);
  
  /* now do inverse fourier grid computation */
  grid_invfourier_FFTPACK(rdata, size, wSave);
  
  /* amscray */
  
}

/************************************************************************/
/* Zonal Harmonic transform using seminaive algorithm - used in convolutions */
/* rdata and idata should be pointers to size x size arrays.
   rres and ires should be pointers to double arrays of size bw.
   size = 2 * bw
   cos_pml_table contains Legendre coefficients of P(0,l) functions
   and is result of CosPmlTableGen for m = 0;
   FZT_semi only computes spherical harmonics for m=0.

   dataformat =0 -> samples are complex, =1 -> samples real
   
   workspace needed is (12 * bw) 

   input data assumed real

*/

void FZT_semi_memo(double *rdata,
		   double *rres, double *ires,
		   int size,
		   double *cos_pml_table,
		   double *workspace,
		   int dataformat,
		   double *wSave)
{

  int i, j, bw;
  double *r0, *i0, dsize;
  double tmpreal, tmpimag;
  double *scratchpad;
  double *cos_even;
  int l, dummy ;
  double tmpA ;

  bw = size/2;

  /* assign memory */
  r0 = workspace;             /* needs (2 * bw) */
  i0 = r0 + (2 * bw);         /* needs (2 * bw) */
  cos_even = i0 + (2 * bw);   /* needs (1 * bw) */
  scratchpad = cos_even + (1 * bw); /* needs (8 * bw) */
                              /* total workspace = 13*bw */


  dsize = (double) size;
  /* compute the m = 0 components */

  for (i=0; i<size; i++) {
    tmpreal = 0.0;
    for(j=0; j<size; j++) {
      tmpreal += rdata[(i*size)+j];
    }
    /* normalize */
    r0[i] = tmpreal/dsize;
  }

  /* do the real part */
  SemiNaiveReduced(r0, 
		   bw, 
		   0,
		   rres, 
		   cos_pml_table,
		   scratchpad,
		   cos_even,
		   wSave);

  /* set imag coefficients = 0 */
  memset(ires, 0, sizeof(double) * size);
  
  /*****************

  New in 2.6: Need to massage the coefficients one more time,
              to get that right factor in there (how embarrassing)

  ******************/

  tmpA = 2. * sqrt( PI );

  for(l=0;l<bw;l++)
    {
      dummy = seanindex(0,l,bw);

      rres[dummy] *= tmpA ;
    }

}
/************************************************************************/
/* multiplies harmonic coefficients of a function and a filter */
/* See convolution theorem of Driscoll and Healy */
/* datacoeffs should be output of an FST, filtercoeffs the 
   output of an FZT.  There should be (bw * bw) datacoeffs,
   and bw filtercoeffs.
   rres and ires should point to arrays
   of dimension bw * bw. size parameter is 2*bw  
*/

void TransMult(double *rdatacoeffs, double *idatacoeffs, 
	       double *rfiltercoeffs, double *ifiltercoeffs, 
	       double *rres, double *ires,
	       int size)
{

  int m, l, bw;
  double *rdptr, *idptr, *rrptr, *irptr;

  bw = size/2;

  rdptr = rdatacoeffs;
  idptr = idatacoeffs;
  rrptr = rres;
  irptr = ires;

  for (m=0; m<bw; m++) {
    for (l=m; l<bw; l++) {
      compmult(rfiltercoeffs[l], ifiltercoeffs[l],
	       rdptr[l-m], idptr[l-m],
	       rrptr[l-m], irptr[l-m]);

      rrptr[l-m] *= sqrt(4*M_PI/(2*l+1));
      irptr[l-m] *= sqrt(4*M_PI/(2*l+1));

    }
    rdptr += bw-m; idptr += bw-m;
    rrptr += bw-m; irptr += bw-m;
  }
  for (m=bw+1; m<size; m++) {
    for (l=size-m; l<bw; l++){
      compmult(rfiltercoeffs[l], ifiltercoeffs[l],
	       rdptr[l-size+m], idptr[l-size+m],
	       rrptr[l-size+m], irptr[l-size+m]);

      rrptr[l-size+m] *= sqrt(4*M_PI/(2*l+1));
      irptr[l-size+m] *= sqrt(4*M_PI/(2*l+1));

    }
    rdptr += m-bw; idptr += m-bw;
    rrptr += m-bw; irptr += m-bw;
  }

}
/************************************************************************/
/* Here's the big banana
   Convolves two functions defined on the 2-sphere.
   Uses seminaive algorithms for spherical harmonic transforms
   size = 2*bw
   Inputs:
   rdata, idata - (size * size) arrays containing real and
                  imaginary parts of sampled function.
   rfilter, ifilter - (size * size) arrays containing real and
                      imaginary parts of sampled filter function.
   rres, ires - (size * size) arrays containing real and
                  imaginary parts of result function.
   size - should be 2 * bw

   Suggestion - if you want to do multiple convolutions,
   don't keep allocating and freeing space with every call,
   or keep recomputing the spharmonic_pml tables.
   Allocate workspace once before you call this function, then
   just set up pointers as first step of this procedure rather
   than mallocing.  And do the same with the FST, FZT, and InvFST functions.

   Memory requirements for Conv2Sphere

   Need space for spharmonic tables and local workspace and
   scratchpad space for FST_semi

   2 * (Spharmonic_TableSize(bw))  +
   8 * (bw*bw)  +
   2 * bw +
   8 * (bw*bw) + 32 * bw
   

   ASSUMPTIONS:
   1. data is strictly REAL
   2. will do semi-naive algorithm for ALL orders

*/
void Conv2Sphere_semi_memo(double *rdata,
			   double *rfilter,
			   double *rres,
			   int size,
			   double *workspace,
			   double *wSave,
			   double *CoswSave )
{
  int bw, spharmonic_bound;
  double *frres, *fires, *filtrres, *filtires, *trres, *tires;
  double **spharmonic_pml_table, **transpose_spharmonic_pml_table;
  double *spharmonic_result_space, *transpose_spharmonic_result_space;
  double *scratchpad;

  bw = size/2;

  /* assign space */

  spharmonic_bound = Spharmonic_TableSize(bw);

  spharmonic_result_space = workspace;
  transpose_spharmonic_result_space = spharmonic_result_space + 
                                      spharmonic_bound;
  frres = transpose_spharmonic_result_space +
          spharmonic_bound;   /* needs (bw*bw) */
  fires = frres + (bw*bw);    /* needs (bw*bw) */
  trres = fires + (bw*bw);    /* needs (bw*bw) */
  tires = trres + (bw*bw);    /* needs (bw*bw) */
  filtrres = tires + (bw*bw); /* needs bw */
  filtires = filtrres + bw;   /* needs bw */
  scratchpad = filtires + bw; /* needs (8 * bw^2) + (32 * bw) */

  spharmonic_pml_table = 
    Spharmonic_Pml_Table(bw,
			 spharmonic_result_space,
			 scratchpad);

  transpose_spharmonic_pml_table = 
    Transpose_Spharmonic_Pml_Table(spharmonic_pml_table, 
				   bw,
				   transpose_spharmonic_result_space,
				   scratchpad);

  FST_semi_memo(rdata,
		frres, fires, 
		size, 
		spharmonic_pml_table,
		scratchpad,
		wSave,
		CoswSave,
		bw);
  
  FZT_semi_memo(rfilter,
		filtrres, filtires, 
		size, 
		spharmonic_pml_table[0],
		scratchpad,
		1,
		CoswSave );
  TransMult(frres, fires, filtrres, filtires, trres, tires, size);

  InvFST_semi_memo(trres, tires, 
		   rres, 
		   size,
		   transpose_spharmonic_pml_table,
		   scratchpad,
		   wSave,
		   bw);

  /***
    have to free the memory that was allocated in
    Spharmonic_Pml_Table() and
    Transpose_Spharmonic_Pml_Table()
    ***/

  free(spharmonic_pml_table);

  free(transpose_spharmonic_pml_table);
}

#endif
