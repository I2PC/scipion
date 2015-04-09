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

  fft_grids.c - routines to perform 1-d ffts on grids, expected
                to be used in spherical transforms!!!


  Compilation flag information: 

  If compiled with -DFFTPACK (i.e. defined in the Makefile), then
  will use the fftpack-based routines. Otherwise, will use the
  (slower) ffts defined in FFTcode.c.


  The FFTcode-based routines in this file are:

  grid_fourier
  grid_invfourier


  The fftpack-based routines in this file are:

  grid_fourier_FFTPACK     ( calls the fftpack routine rfftf )
  grid_invfourier_FFTPACK  ( calls the fftpack routine rfftb )

  These _FFTPACK function definitions can be easily modified
  to allow use of other optimized fft routines, if so desired.
  Just be careful that scaling of coefficients and samples
  is correct.
*/

#include <math.h>
#include <string.h>

#include "primitive_FST.h"

#ifndef FFTPACK
#include "FFTcode.h"
#else
#include "fftpack.h"
#endif


/************************************************************************/
/************************************************************************/

/* compile the appropriate fft routine, depending on whether
   or not FFTPACK is being used */


#ifndef FFTPACK        /* use the routines in FFTcode.c */

/************************************************************************
  Computes the fourier transform of each row of the grid.  This
  is NOT the same as a 2-D Fourier transform.

  Used by FST_semi procedure

  Since this will be input to an associated legendre transform,
  the lines of longitude, or zones, are loaded into the rows
  in a transposed fashion for easy access by the Legendre
  transform rpocedure.  The grid is expected to
  be size * size, which is probably (2*bw) * (2*bw).
  
  realgrid, imaggrid - (size x size) arrays of real and imag
                       input
  rmatrix, imatrix - (size x size) arrays of real and imag
                     output
  size = 2 * bw
  
  workspace - double pointer to array of (6 * size) = (24 * bw)

  *********************************************************/


void grid_fourier(double *realgrid,
		  double *imaggrid,
		  double *rmatrix,
		  double *imatrix,
		  int size,
		  double *workspace)
{

  double *rout, *iout, *scratchpad;
  int i;

  /* need to assign workspace - need 6*size total */

  rout = workspace; /* needs size space */
  iout = rout + size; /* needs size space */
  scratchpad = iout + size; /* needs 4 * size space */

  if(1)
    {
      for (i=0; i<size; i++) 
	{
	  FFTInterp(realgrid+(i*size), imaggrid+(i*size),
		    rmatrix+(i*size),
		    imatrix+(i*size),
		    size, size, scratchpad, 1);
	}
    }
  else
    {
      for (i=0; i<size; i++) 
	{
	  FFTInterp(realgrid+(i*size), imaggrid+(i*size),
		    rout,
		    iout,
		    size, size, scratchpad, 1);

	  memcpy(rmatrix+(i*size),rout, sizeof(double) * size);
	  memcpy(imatrix+(i*size),iout, sizeof(double) * size);

	}
    }

  
  /* now transpose the results */
  transpose(rmatrix,size);
  transpose(imatrix,size);

}

/***********************************************************************

  Same as above except for inverse Fourier transform is used
  used by InvFST_semi procedure 

  workspace = (24 * bw)

  **********************************************************************/

void grid_invfourier(double *realgrid, double *imaggrid,
		     double *rmatrix, double *imatrix,
		     int size, double *workspace)
{

  double *rout, *iout, *scratchpad;
  int i;

  /* need to assign workspace - need 6*size total */
  
  rout = workspace; /* needs size space */
  iout = rout + size; /* needs size space */
  scratchpad = iout + size; /* needs 4 * size space */
  
  for (i=0; i<size; i++) {
    FFTEval(realgrid+(i*size), imaggrid+(i*size),
	    rmatrix+(i*size),
	    imatrix+(i*size), size, size, scratchpad, 1);
  }
  
}



/******************************************/

#else    /* if FFTPACK *is* defined */

/******************************************/

/****

  If will use the fftpack routines rfftf and rfftb, first
  have to precompute an array of values that both routines
  use. This interface function will call the fftpack routine
  
       rffti

  This function will allocate memory to store the precomputed
  values and then will return a double pointer to that array.

  NOTE: The memory allocated in this function is freed in
        the routine that calls precomp_fft!!!

  
  size = length of sequence to be transformed
  
  wSave = the array created by precomp_fft, of length
          2 * size + 15, which contains the precomputed
	  values. precomp_fft will return a double pointer
	  to this array.


  ****/

double *precomp_fft( int size )
{
  double *wSave;

  /* allocate space */
  wSave = (double *) malloc( sizeof(double) * (2 * size + 15) );

  /* now precompute data */
  rffti_( &size, wSave );

  return wSave;

}


/****

  Computes the fourier transform of each row of the grid.  This
  is NOT the same as a 2-D Fourier transform.

  Basically, this function does the same as the above forward
  routine EXCEPT that it uses FFTPACK

  NOTE: Writes over the input data!!!

  size = 2 * bw

  wSave = double array, holds the precomputed data,
          has dimensions 2 * size + 15
  
  ***/

void grid_fourier_FFTPACK(double *grid,
			  int size,
			  double *wSave)
{
  int i;
  double tmp;

  /* normalizing factor */
  tmp = 1.0 / ((double) size);

  for (i = 0 ; i < size ; i ++)
    {
      rfftf_(&size, grid+(i*size), wSave);
    }

  /* normalize result */
  for ( i = 0 ; i < size * size ; i ++ )
    grid[i] *= tmp;

  /* now transpose the result */
  transpose(grid, size);

}


/****

  Same as above except for inverse Fourier transform is used
  used by InvFST_semi procedure 

  Basically, this function does the same as the above inverse
  routine EXCEPT that it uses FFTPACK

  NOTE: Writes over the input data!!!

  wSave = double array, holds the precomputed data,
          has dimensions 2 * size + 15
  
  ***/

void grid_invfourier_FFTPACK( double *grid,
			      int size,
			      double *wSave )
{

  int i;

  for (i=0; i<size; i++)
    {
      rfftb_( &size, grid+(i*size), wSave);
    }

}
#endif
