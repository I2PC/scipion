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

/************************************************************************

  If FFTPACK *is* defined, use the FFTPACK-based dcts. Small
  C-interfaces have been written, to act as front ends to the
  functions which live in the fftpack library that I'll be linking
  to.

  NOTE THAT even if FFTPACK is defined, I still need the inverse
  DCT routine ExpIFCT. That is because, in the inverse discrete
  Legendre transform, I need to take N coefficients and inverse
  DCT them to produce 2*N samples. I don't know how to do that
  with FFTPACK's inverse dct. If I have to do an inverse DCT whose
  input and output are the same lengths, *then* I can use FFTPACK's
  routine.

  The "standard" (i.e. non-FFTPACK) forward and inverse dct
  routines defined in this file are

  forward dct: kFCT, kFCTX
  inverse dct: ExpIFCT

  When using FFTPACK, the routines are

  forward dct: DCTf    ("f" for forward)
  inverse dct: DCTb    ("b" for backward)

  If you have your own optimized dct routines, just modify
  the definitions of DCTf and DCTb accordingly (and still
  compile with the -DFFTPACK flag and change the library
  you're linking to in the Makefile).

  NOTE: DCTf and DCTb use slightly modified versions of
        the fftpack functions cosqb ( in DCTf ) and cosqf (in DCTb ).
	Before compiling the fftpack library, these functions
	were modified in order to scale the coefficients the same way
	as they are in kFCT, kFCTX, and ExpIFCT. Originally,
	the transforms cosqb and cosqf were not orthogonal.
	THEY ARE NOW.

	Exactly how was the scaling modified in the original,
	unmodified fftpack routines cosqb, cosqf ?

	Let X be an input array of length N. Let Y be the output
	of kFCT (or kFCTX - they're identical), Z be the
	output of the original, unmodified fftpack routine cosqb.
	Let SCALE = 1 / (2 * N) .Then (using C and not Fortran indexing)

	Y[0] = SCALE * Z[0] / 2
	Y[i] = SCALE * Z[i]    for i = 1, 2, ..., N-1

	The original fftpack library routine cosqb was further
	modified to return only the first p <= N many
	coefficients. This required modifying the original,
	unmodified fftpack-source file cosqb.f .


	Let X be an input array of length N. Let Y be the output
	of ExpIFCT. Now let me multiply the first coefficient of
	X by 2 ( so X[0] *= 2 ) and plug this modified X into
	the original, unmodified fftpack-routine cosqf.
	Let Z be the resulting output of cosqf. Let SCALE = 1/2.
	Then (using C and not Fortran indexing)

	Y[i] = SCALE * Z[i]    for i = 0, 1, 2, ..., N-1


	All this scaling is now done within the source files
	cosqf1.f and cosqb1.f . It's a little hokey, but it
	works.


  IMPORTANT NOTE: DCTf and DCTb both write the output
                  over the input array!!!

************************************************************************/

/************************************************

  Some general background on the routines kFCT, kFCTX and ExpIFCT ...


  This is new C cource code for doing fast cosine transforms, and
  it contains variations for speeding up FCTs and convolutions
  for certain structures appearing in Legendre transforms.
  Many of these commands have both a real-valued version and
  a complex-valued version to save time when performing 
  spherical harmonic transforms.
  It also is written with memory management in mind.  
  
  This cosine transform is based on the Steidl and Tasche algorithm
  with some modifications.  See Sean's Ph.D. thesis for refs.
  

  NOTE: I will try to explicitly unwrap the ChebyDivRem routine
  as called in ExpIFCT: the while-loop will be for divsize > 4.
  Then I'll explicitly type in for the cases of divsize = 4 and
  divsize = 2. This, I believe, significantly reduces the overhead
  involved in calling ChebyDivRem so often. Thank you Doug and Sumit !


  Note about structures: some functions use the lowhigh structure,
  defined to be
  
  struct lowhigh{ double low; double high; }
  
  The reason for this is to "interweave" elements. That is, if
  I have a vector A of length 2*n, then I'll create an array with
  n-many lowhigh elements which will look like
  
  {{A[0],A[n]},{A[1],A[n+1]},{A[2],A[n+2]},...,{A[n-1],A[2*n]}}
  
  Having the elements arranged this way increases in a small
  (but real) way the efficiency of the forward fast cosine transform.
  So, for example, kFCTX is the "structure" version of kFCT. The
  structure is well suited for the transform, given how this
  particular algorithm is implemented. Basically, I believe,
  using this structure cuts down on cache-thrashing.
  
  General note about flags - throughout this code flags are often
  used as input arguments to tell the routine something about
  the input information.  The convention adopted throughout is
  that when a flag is set to 0, that means that no action needs
  to be performed on the pertinent data, and 1 means that some action
  must be taken.
  

  ************************************************/



#include <math.h>
#include <string.h>

#include "OURperms.h"
#include "OURmods.h"
#ifdef FFTPACK
#include "fftpack.h"
#endif
#include "newFCT.h"


/************************************************************************

  Utility functions for kFCT, kFCTX and ExpIFCT

  *********************************************************************/

/************************************************************************

  reverse the elements of a (double) vector of length n
  i.e. if a = [0 1 2 3] then what's returned is
          b = [3 2 1 0]

	  input -> vector to reverse
	  n     -> length of input (must be even !)
	  output -> where to store result

************************************************************************/

static void reverse_it(double *input,
		       double *output,
		       int n)
{
  int i, len;

  len = n / 2;

  if (len >= 4)
    {
      for(i = 0; i < (len % 4); ++i)
	{
	  output[i] = input[n - 1 - i];
	  output[n - 1 - i] = input[i];
	}
      for( ; i < len ; i += 4)
	{
	  output[i] = input[n - i - 1];
	  output[i + 1] = input[n - i - 2];
	  output[i + 2] = input[n - i - 3];
	  output[i + 3] = input[n - i - 4];

	  output[n - i - 1] = input[i];
	  output[n - i - 2] = input[i + 1];
	  output[n - i - 3] = input[i + 2];
	  output[n - i - 4] = input[i + 3];
	}
    }
  else if (len == 2)
    {
      output[0] = input[n - 1];
      output[1] = input[n - 2];
      output[n - 1] = input[0];
      output[n - 2] = input[1];
    }
  else if (len == 3)
    {
      output[0] = input[n - 1];
      output[1] = input[n - 2];
      output[2] = input[n - 3];
      output[n - 1] = input[0];
      output[n - 2] = input[1];
      output[n - 3] = input[2];
    }
  else
    {
      output[0] = input[n - 1];
      output[n - 1] = input[0];
    }

}

/************************************************************************/
/* permutation management code */
/* OUR permutations - a self-inverting permutation used for
   performing fast cosine transforms.  See Sean's Ph.D. thesis */
/************************************************************************/
/* permutes data and puts the result in result */
/* data and result should be arrays of size n */

static void OURpermute(double *data,
		       double *result,
		       int n)
{
  int i;
  const int *pptr;

  pptr = get_perm(n);
  
  for (i=0; i<n; i++)
    result[i] = data[pptr[i]];
  
}
/************************************************************************/
/* permutes data and puts the result in result */
/* data is size n, result array of structures of size n/2 */
/* just like above except designed to work with structures */

static void OURpermuteX(double *data,
			struct lowhigh *result,
			int n)
{
  register int i, dummy;
  const int *pptr;

  dummy = n / 2;

  pptr = get_perm(n);

  for ( i = 0 ; i < dummy; i ++)
    {
      result->low = data[*pptr];
      result->high = data[*(pptr+dummy)];
      result++; pptr++;
    }
}

/************************************************************************/
/* used by kFCT to add up polynomial coefficients once k level is reached */
/* data is double array of size n, result is double array of size k */

static void PartitionAdd(double *data,
			 double *result,
			 int n,
			 int k)
{
  int i, j;
  double tmp;

  for(i = 0; i < k ; i++)
    {
      tmp = 0.0;
      for(j = i; j < n; j += k)
	tmp += data[j];
      result[i] = tmp;
    }
}

/************************************************************************/
/* used by kFCT to add up polynomial coefficients once k level is reached */
/* data is double array of size n, result is double array of size k */
/* just like above except used for lowhigh structures ...
 ASSUMING k = (n/2) */

static void PartitionAddX(struct lowhigh *data,
			  double *result,
			  int n,
			  int k)
{
  int i, j;
  double tmp0, tmp1;
  int dummy;

  dummy = n / 2;

  for (i = 0 ; i < k ; i++)
    {
      tmp0 = 0.0 ; tmp1 = 0.0;
      for(j = i ; j < dummy ; j += k)
	{
	  tmp0 += data[j].low;
	  tmp1 += data[j].high;
	}
      result[i] = ( tmp0 + tmp1 );
    }
}
/************************************************************************/

/************************************************************************

  Ok, if FFTPACK *is* defined, use DCTf, DCTb. Otherwise use
  kFCT, kFCTX

  *********************************************************************/

#ifdef FFTPACK


/****

  If will use the fftpack routines cosqf and cosqb, first
  have to precompute an array of values that both routines
  use. This interface function will call the fftpack routine
  
       cosqi

  This function will allocate memory to store the precomputed
  values and then will return a double pointer to that array.

  NOTE: The memory allocated in this function is freed in
        the routine that calls precomp_dct!!!

  
  size = length of sequence to be transformed
  
  wSave = the array created by precomp_dct, of length
          3 * size + 15, which contains the precomputed
	  values. precomp_dct will return a double pointer
	  to this array.


  ****/

double *precomp_dct( int size )
{
  double *wSave;

  /* allocate space */
  wSave = (double *) malloc( sizeof(double) * (3 * size + 15) );

  /* now precompute data */
  cosqi_( &size, wSave );

  return wSave;

}


/***********************************

  This is the forward DCT routine DCTf. It uses the fftpack library
  function cosqb (NOT cosqf).

  array = double array of samples, of length n

  p = how many coefficients desired.

  CoswSave = array of precomputed data that cosqb requires (produced
             by a call to the function cosqi ).

  OUTPUT IS WRITTEN OVER INPUT !!!

  *******************************/

void DCTf ( double *array,
	    int n,
	    int p,
	    double *CoswSave )
{
  
  cosqb_( &n, &p, array, CoswSave );

}

/***********************************

  This is the inverse DCT routine DCTb. It uses the fftpack library
  function cosqf (NOT cosqb).

  array = double array of coefficients, of length n;
          n many samples will be produced

  CoswSave = array of precomputed data that cosqf requires (produced
             by a call to the function cosqi ).

  
  OUTPUT IS WRITTEN OVER INPUT !!!

  *******************************/

void DCTb( double *array ,
	   int n,
	   double *CoswSave )
{
  cosqf_( &n , array , CoswSave );

}

/************************************************************************/
/************************************************************************/

#else   /* FFTPACK is *not* defined, use the "local" functions kFCT, kFCTX */

/************************************************************************/
/************************************************************************/

/************************************************************************/
/* OK - here is the forward transform code.  The main modification here
   is that given n samples of data, kFCT will compute only the first
   k <= n coefficients, where n and k are assumed to be powers of 2.
   The basic idea is that the algorithm keeps interpolating a polynomial
   of degree 2d-1 from two polynomials of degree d-1.  Initially, lowpoly
   contains polynomials of degree 0 - the data samples.  From these
   values, polynomials of degree 1 are interpolated and stored in highpoly.
   This transform is (essentially) the transpose of the matrix formulation
   of ExpIFCT.  Series is returned with coefficients ordered from low
   degree to high degree

   data - double array of size n
   result - double array of size k
   workspace - double array of size 2*n
   n - number of samples
   k = number of coefficients desired 
   permflag - 0 if data in OUR permuted order, 1 if data needs to 
              be permuted */

void kFCT(double *data,
	  double *result,
	  double *workspace,
	  int n,
	  int k,
	  int permflag)
{
    double *lowpoly, *highpoly, *lowptr, *highptr;
    const double *modptr;
    double *temp, dn2;
    int p2, p3;
    unsigned int i, j, levelsize, highsize, lowsize, loopcntr;
    double modptr0;
    double e0, e1, f0, f1;

    lowpoly = workspace; 
    lowptr = lowpoly;
    highpoly = workspace + n; 
    highptr = highpoly;
    
    levelsize = n;
    modptr = get_mods(levelsize);

    if (permflag == 1)
      OURpermute(data,lowpoly,n);
    else {
      memcpy(lowpoly, data, sizeof(double) * n);
    }

    /***** in what follows, am using the fact
      that modptr[i+1] = -modptr[i] for i even ****/

    /* do first level of interpolation  */
    for (i=0; i<n; i+=2) {
	*highptr++ = *lowptr + *(lowptr+1);
	*highptr++ = *modptr * (*(lowptr+1) - *lowptr);
	modptr+=2;lowptr+=2;
    }

    lowptr = highpoly;
    highptr = lowpoly;
    loopcntr = 0;

    /* now do higher-order interpolations */
    if (k > 2)
      {
	/* reset pointers to logical polynomial workspace */	
	levelsize /= 2;
	modptr = get_mods(levelsize);
	highsize = 4;
	lowsize = 2;

	/* now set counter pointers to appropriate matrix locations */
	p2 = lowsize - 1;
	p3 = 1;
	
	/* main loop */
	while (highsize <= k)
	  {
	    for (j=0; j<n; j+=highsize)
	      {
		/* do loworder coeffs first */
		for (i=0; i<lowsize; i++)
		  *highptr++ = lowptr[i] + lowptr[i+lowsize];

		/* now do higher order coeffs */
		
		modptr0 = modptr[0];
		*highptr++ = modptr0 * (lowptr[lowsize] - lowptr[0]);
		modptr0 *= 2.0;
		
		e0 = modptr0 * ( lowptr[p3+lowsize] - lowptr[p3] );
		e1 = lowptr[p2+lowsize] + lowptr[p2];
		*highptr++ = e0 - e1;
		p2--;
		p3++;
		
		for ( i = 0 ; i<(lowsize-2); i += 2)
		  {
		    e0 = -lowptr[p3] + lowptr[p3+lowsize];
		    f0 = -lowptr[p3+1] + lowptr[p3+lowsize+1];
		    f1 = lowptr[p2-1] + lowptr[p2+lowsize-1];
		    e1 = lowptr[p2] + lowptr[p2+lowsize];		    
		    *highptr++ = modptr0 * e0 - e1;
		    *highptr++ = modptr0 * f0 - f1;
		    p2 -= 2;
		    p3 += 2;
		  }

		lowptr += highsize;
		modptr += 2;
		p2 = lowsize - 1;
		p3 = 1;
	      }

	    levelsize /= 2;
	    modptr = get_mods(levelsize);
	    highsize *= 2;
	    lowsize *= 2;
	    p2 = lowsize - 1;
	    p3 = 1;
	    
	    if ((loopcntr % 2) == 0)
	      {
		lowptr = lowpoly;
		highptr = highpoly;
	      }
	    else
	      {
		lowptr = highpoly;
		highptr = lowpoly;
	      }
	    loopcntr++;
	  }
      }
    
    /* find where the top level polynomials are */
    if (!(loopcntr % 2))
      temp = highpoly;
    else
      temp = lowpoly;
    
    if (k < n) 
      PartitionAdd(temp,result,n,k);
    else
      memcpy(result, temp, sizeof(double) * n);
    
    /* normalize */
    dn2 = 2.0/((double) n);
    result[0] /= ((double) n);
    for(i=1; i<k; i++)
      result[i] *= dn2;
        
    /* later */
}

/************************************************************************/
/* OK - here is the forward transform code.  The main modification here
   is that given n samples of data, kFCT will compute only the first
   k <= n coefficients, where n and k are assumed to be powers of 2.
   The basic idea is that the algorithm keeps interpolating a polynomial
   of degree 2d-1 from two polynomials of degree d-1.  Initially, lowpoly
   contains polynomials of degree 0 - the data samples.  From these
   values, polynomials of degree 1 are interpoltaed and stored in highpoly.
   This transform is (essentially) the transpose of the matrix formulation
   of ExpIFCT.  Series is returned with coefficients ordered from low
   degree to high degree

   data - double array of size n
   result - double array of size k
   workspace - double array of size 2*n
   n - number of samples
   k = number of coefficients desired 
   permflag - 0 if data in OUR permuted order, 1 if data needs to 
              be permuted
   lowpoly, highpoly - lowhigh arrays of length n / 2

 */

/*** just like kFCT EXCEPT will use structures ***/

void kFCTX( double *data,
	    double *result,
	    double *workspace,
	    int n,
	    int k,
	    int permflag,
	    struct lowhigh *lowpoly,
	    struct lowhigh *highpoly)
{
  double dn2;
  int p2, p3;
  register int i, j;
  int levelsize, highsize, lowsize, loopcntr;

  int dummy_int, dummy_int2;

  double e0, e1, f0, f1;
  double e0x, e1x, f0x, f1x;

  struct lowhigh *lowptr, *highptr, *temp;
  const double *modptrlow, *modptrhigh;
  register double modptrlow0, modptrhigh0;
  struct lowhigh *tmphigh;
  struct lowhigh *tmplowa, *tmplowb;

  lowptr = lowpoly;
  highptr = highpoly;
  
  levelsize = n;

  modptrlow = get_mods(levelsize);
  modptrhigh = modptrlow + (levelsize/2);
  
  if (permflag == 1)
    {
      OURpermuteX(data,lowpoly,n);
    }
  else
    {
      dummy_int = n / 2;
      for(i = 0; i < dummy_int ; i++)
	{
	  lowpoly[i].low = data[i];
	  lowpoly[i].high = data[i + dummy_int];
	}
    }
  
  /* do first level of interpolation  */
  dummy_int = n / 2;

  for ( i = 0 ; i < dummy_int ; i += 2)
    {
      (*highptr).low = (*lowptr).low + (*(lowptr+1)).low;
      (*highptr).high = (*lowptr).high + (*(lowptr+1)).high;
      highptr++;
    
      (*highptr).low = -(*modptrlow) *
	( (*lowptr).low - (*(lowptr+1)).low) ;
      (*highptr).high = -(*modptrhigh) *
	( (*lowptr).high - (*(lowptr+1)).high) ;
      highptr++;
      modptrlow += 2; modptrhigh += 2;
      lowptr += 2;

    }

  /* lowptr = highpoly; */
  lowptr = highpoly;
  highptr = lowpoly;
  loopcntr = 0;
  
  /* now do higher-order interpolations */
  if (k > 2) {
    
    /* reset pointers to logical polynomial workspace */
    levelsize /= 2;
    modptrlow = get_mods(levelsize);
    modptrhigh = modptrlow + (levelsize/2);
    highsize = 4;
    lowsize = 2;
    
    /* now set counter pointers to appropriate matrix locations */
    p2 = lowsize - 1;
    p3 = 1;
    
    /* main loop */
    while (highsize <= k) {
      dummy_int = n / 2;      

      for (j=0; j< dummy_int; j+=highsize) {
	/* do loworder coeffs first */

	tmplowa = lowptr;
	tmphigh = highptr;
	i = 0;
	do
	  {

	    tmphigh->low =
	      tmplowa->low + (tmplowa+lowsize)->low;
	    (tmphigh) ->high =
	      (tmplowa) ->high + (tmplowa+lowsize)->high;

	    tmplowa ++;
	    tmphigh ++;

	    tmphigh->low =
	      tmplowa->low + (tmplowa+lowsize)->low;
	    (tmphigh) ->high =
	      (tmplowa) ->high + (tmplowa+lowsize)->high;

	    tmplowa ++; tmphigh ++;
	    i += 2;

	  } while (i < lowsize );

	highptr += lowsize;
	
	/* now do higher order coeffs */
	
	modptrlow0 = *modptrlow;
	highptr->low = modptrlow0 * ((lowptr+lowsize)->low) -
	  modptrlow0 * (lowptr->low) ;

	modptrhigh0 = *modptrhigh;
	highptr->high = modptrhigh0 * ((lowptr+lowsize)->high) -
	  modptrhigh0 * (lowptr->high);

	modptrlow0 *= 2.0;
	modptrhigh0 *= 2.0;

	tmplowa = lowptr+p3;
	tmplowb = lowptr+p2;
	tmphigh = highptr+1;

	e0 = (tmplowa + lowsize)->low - tmplowa->low;
	e1 = (tmplowb + lowsize)->low + tmplowb->low;
	tmphigh->low = modptrlow0 * e0 - e1;
	
	f0 = (tmplowa + lowsize)->high - tmplowa->high;
	f1 = (tmplowb + lowsize)->high + tmplowb->high;
	tmphigh->high = modptrhigh0 * f0 - f1;
	
	tmplowa++; tmplowb--; tmphigh++;

	dummy_int2 = lowsize - 1;

	for( i = 1 ; i < dummy_int2 ; i += 2 )
	  {
	    e0 = (tmplowa + lowsize)->low - tmplowa->low;
	    e0x = (tmplowa + lowsize)->high - (tmplowa)->high;
	    e1 = (tmplowb + lowsize)->low + tmplowb->low;
	    e1x = (tmplowb + lowsize)->high + (tmplowb)->high;

	    tmphigh->low = modptrlow0 * e0 - e1;
	    (tmphigh) ->high = modptrhigh0 * e0x - e1x;

	    tmphigh++;
	    tmplowa++; tmplowb--;

	    f0 = (tmplowa + lowsize)->low - tmplowa->low;
	    f0x = (tmplowa + lowsize)->high - (tmplowa)->high;

	    f1 = (tmplowb + lowsize)->low + tmplowb->low;
	    f1x = (tmplowb + lowsize)->high + (tmplowb)->high;

	    tmphigh->low = modptrlow0 * f0 - f1;
	    (tmphigh)->high = modptrhigh0 * f0x - f1x;

	    tmplowa ++ ; tmplowb --; tmphigh ++;
	  }

	highptr += lowsize;
	lowptr += highsize;
	modptrlow += 2;
	modptrhigh += 2;
	p2 = lowsize - 1;
	p3 = 1;

      }
      
      levelsize /= 2;
      modptrlow = get_mods(levelsize);
      modptrhigh = modptrlow + (levelsize/2);
      highsize *= 2;
      lowsize *= 2;
      p2 = lowsize - 1;
      p3 = 1;
      
      if ((loopcntr % 2) == 0) {
	lowptr = lowpoly;
	highptr = highpoly;
      }
      else {
	lowptr = highpoly;
	highptr = lowpoly;
      }
      loopcntr++;
    }
  }
  
  /* find where the top level polynomials are */
  if (!(loopcntr % 2))
    temp = highpoly;
  else
    temp = lowpoly;
  
  if (k < n) 
    PartitionAddX(temp,result,n,k);
  else
    for (i=0; i<n; i++){
      result[i] = temp[i].low;
      result[i+(n/2)] = temp[i].high;
    }

  /* normalize  */
  dn2 = 2.0/((double) n);
  result[0] /= ((double) n);
  for(i=1; i<k; i++)
    result[i] *= dn2;

  /* later */
}


#endif   /* FFTPACK */


/***********************************************************************

  Now for the inverse DCT routine, ExpIFCT, which is used regardless
  of FFTPACK being defined or not.

  **********************************************************************/

/***********************************************************************
  Chebyshev division program.  This is specifically coded for dividing 
  a chebyshev series of degree (n-1) by a modulus of the form
  T_n/2 + m(0)T_0, which returns a remainder sequence of degree
  (n/2) - 1.  Expects dividend coefficient sequence to be ordered from high
  degree to low degree.  Modulus should be a double value, m(0).  size is
  just n, the degree + 1 of the dividend.  The remainder is returned 
  with coefficients ordered from high degree to low degree.
  
  This function is used in the inverse DCT routine ExpIFCT.
  

  dividend - double array of size n
  remres - double array of size n/2
  m0 - supermodulus zero-degree coefficient value 
  
  *******************************************************************/

/****
  THE ORIGINAL ... easiest to understand how the algorithm actually
  works ... the code uses more efficient implementations of this
  original version.

  ****/

static void ChebyDivRemOrig(double *dividend,
			    double m0,
			    double *remres,
			    int n)
{
  int rlim, dd;

  rlim = n/2;

  memcpy(remres, dividend + rlim, sizeof(double) * rlim );
  for (dd=0; dd<rlim-1; dd++) {
    remres[rlim-1-dd-1] -= dividend[dd];
    remres[dd] -= (2.0 * m0 * dividend[dd]);
  }
  remres[rlim-1] -= m0 * dividend[rlim-1];
}

/***

  The version that is actually used in the code.

  ***/

/****************************************************/

static void ChebyDivRem(double *dividend,
			double m0,
			double *remres,
			int n)
{
  int rlim, dd;
  int tmpfloor;
  double tm0;
  int halfrlim;

  rlim = n/2;
  tm0 = 2.0 * m0;
  tmpfloor = (rlim - 1)/2;
  halfrlim = rlim/2;

  for(dd = 0; dd < tmpfloor; dd++)
    {
      remres[dd] =
	-tm0 * dividend[dd] - dividend[rlim-dd-2] + dividend[rlim+dd];
      remres[dd+halfrlim] =
	-dividend[halfrlim-2-dd] - tm0 * dividend[dd+halfrlim]
	+ dividend[rlim+dd+halfrlim];
    }
  remres[tmpfloor] = dividend[tmpfloor+rlim] -
    dividend[tmpfloor] - (tm0 * dividend[tmpfloor]);

  remres[rlim-1] = dividend[n-1] - (m0 * dividend[rlim-1]);

}

/****************************************************/

/*****************************************************************
  Fast Expanded Inverse Chebyshev transform

  Expects input data to be a Chebyshev series, with low-order 
  coefficients occuring first.
  Evaluates a Chebyshev series of length k at n data points,
  where n >= k.  Assumes power of 2

  data - double array of length k
  result - double array of length n
  workspace - double array of length (2*n)
  n - number of samples (evaluation points) desired
  k - number of coefficients 
  permflag = if 0, then don't permute the result.  If 1,
             permute data using OUR permutation

  *********************************************************************/

void ExpIFCT( double *data,
	      double *result,
	      double *workspace,
	      int n,
	      int k,
	      int permflag)
{
  int j, loopcntr;
  int levelsize, divsize;
  double *dividend, *remres;
  const double *modptr;
  double *divptr, *remptr;
  double divptr0, divptr1, divptr2, divptr3;
  double tmpmod, tmpmod2;

  /* assign workspace locations */
  remres = workspace;
  dividend = workspace + n;
  
  /* first reverse the Chebyshev series to order from high to low */
  reverse_it(data, dividend, k);
  
  /* compute initial values of variables */
  /* levelsize is the number of supermoduli at first level of division */
  /* divsize is the size of each divisor */
  /* modptr points to list of supermods at this level */
  
  levelsize = 2 * (n/k);
  divsize = k;
  modptr = get_mods(levelsize);
  divptr = dividend;
  remptr = remres;
  loopcntr = 0;
  
  /* now do divisions */
  for(j = 0; j < levelsize; j += 2)
    {
      ChebyDivRem(divptr,modptr[j],remptr,divsize);
      remptr += (divsize>>1);
      ChebyDivRem(divptr,modptr[j+1],remptr,divsize);
      remptr += (divsize>>1);
    }
  
  /* now reset everything */
  divptr = remres;
  remptr = dividend;
 
  loopcntr++;
  divsize /= 2;
  levelsize *= 2;
  modptr = get_mods(levelsize);

  /*** now loopcntr = 1 ***/
  while (divsize > 4)
    {

      for(j = 0 ; j < levelsize; j += 2)
	{
	  ChebyDivRem(divptr,modptr[j],remptr,divsize);
	  remptr += (divsize>>1);
	  ChebyDivRem(divptr,modptr[j+1],remptr,divsize);
	  remptr += (divsize>>1);
	  divptr += divsize;
	}
      
      /* now reset everything */
      if ( (loopcntr % 2) == 0 )
	{
	  divptr = remres;
	  remptr = dividend;
	}
      else 
	{
	  divptr = dividend;
	  remptr = remres;
	}
      
      loopcntr++;
      divsize /= 2;
      levelsize *= 2;
      modptr = get_mods(levelsize);
    }


  /***** in what follows, am using the fact
    that modptr[i+1] = -modptr[i] for i even ****/
  
  if(divsize == 4)
    {
      for (j=0; j<levelsize; j+=2)
	{

	  divptr0 = divptr[0]; divptr1 = divptr[1];
	  divptr2 = divptr[2]; divptr3 = divptr[3];
	  tmpmod = modptr[j];
	  tmpmod2 = 2.0 * modptr[j];
	  remptr[0] = divptr2 - divptr0 - (tmpmod2 * divptr0);
	  remptr[1] = (divptr3 - (tmpmod * divptr1));
	  remptr[2] = divptr2 - divptr0 + (tmpmod2 * divptr0);
	  remptr[3] = (divptr3 + (tmpmod * divptr1));
	  
	  remptr += 4;
	  divptr += 4; 
	}
      
      /* now reset everything */
      if ( (loopcntr % 2) == 0)
	{
	  divptr = remres;
	  remptr = dividend;
	}
      else
	{
	  divptr = dividend;
	  remptr = remres;
	}      
      loopcntr++;
      divsize /= 2;
      levelsize *= 2;
      modptr = get_mods(levelsize);
    }
  
  if(divsize == 2)
    {
      for( j = 0; j < levelsize; j += 2)
	{
	  divptr0 = divptr[0]; divptr1 = divptr[1];
	  
	  remptr[0] = divptr1 - (modptr[j] * divptr0);
	  remptr[1] = divptr1 + (modptr[j] * divptr0);
	  
	  remptr += 2;
	  divptr += divsize;
	}
      
      /* now reset everything */
      if ( (loopcntr % 2) == 0 )
	{
	  divptr = remres;
	  remptr = dividend;
	}
      else
	{
	  divptr = dividend;
	  remptr = remres;
	}
    }
  
  /******** *********/
  /* divptr now points to result - copy into result space 
     after checking permute flag */
  if ( permflag == 0 ) {
    memcpy(result, divptr, sizeof(double) * n);
  }
  else 
    OURpermute(divptr,result,n);
  
  /* later */
  
}

/************************************************************************/

