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


/**********************************************************

  Just some primitive (i.e. utility) functions that the
  various flavours (hybrid, seminaive) of spherical transforms
  share ...


  ********************************************************/

/*****************************************************************

  Given bandwidth bw, seanindex(m,l,bw) will give the position of the
  coefficient f-hat(m,l) in the one-row array that Sean stores the spherical
  coefficients. This is needed to help preserve the symmetry that the
  coefficients have: (l = degree, m = order, and abs(m) <= l)

  f-hat(l,-m) = (-1)^m * conjugate( f-hat(l,m) )

  Thanks for your help Mark!

  ******************************************************************/

int seanindex(int m,
	      int l,
	      int bw)
{     
  int bigL;

  bigL = bw - 1;

  if( m >= 0 )
    return( m * ( bigL + 1 ) - ( ( m * (m - 1) ) /2 ) + ( l - m ) );
  else
    return( ( ( bigL * ( bigL + 3 ) ) /2 ) + 1 +
	    ( ( bigL + m ) * ( bigL + m + 1 ) / 2 ) + ( l - abs( m ) ) );
}


/*****************************************************************
  just like seanindex(m,l,bw) but returns the array coordinates
  for (l,m) AND (l,-m)

  ASSUMING THE M IS GREATER THAN 0 !!!

  this is used in the FST_semi routine

  loc is a 2 element integer array

  ******************************************************************/

void seanindex2(int m,
		int l,
		int bw,
		int *loc)
{     
  int bigL;
  
  bigL = bw - 1;
  
  /* first index for (l,m) */
  loc[0] = m * ( bigL + 1 ) - ( ( m * (m - 1) ) /2 ) + ( l - m );
  
  /* second index for (l,-m) */
  loc[1] = ( ( bigL * ( bigL + 3 ) ) /2 ) + 1 +
    ( ( bigL - m ) * ( bigL - m + 1 ) / 2 ) + ( l -  m ) ;

}


/****************************************************

  just a function to transpose a square array IN PLACE !!!

  array = array to transpose
  size = dimension of array (assuming the array is square, size * size)

  **************************************************/

void transpose(double *array,
	       int size)
{
  register int i, j;
  double t1, t2, t3, t4;

  for(i = 0; i < size; i += 2)
    {
      t1 = array[(i * size) + i + 1];
      array[(i * size) + i + 1] = array[((i + 1) * size) + i];
      array[((i + 1) * size) + i] = t1;
      for(j = (i + 2); j < size; j += 2)
	{
	  t1 = array[(i*size)+j]; t2 = array[(i*size)+j+1];
	  t3 = array[((i+1)*size)+j]; t4 = array[((i+1)*size)+j+1];
	  array[(i*size)+j] = array[(j*size)+i];
	  array[(i*size)+j+1] = array[((j+1)*size)+i];
	  array[((i+1)*size)+j] = array[(j*size)+i+1];
	  array[((i+1)*size)+j+1] = array[((j+1)*size)+i+1];
	  array[(j*size)+i] = t1;
	  array[((j+1)*size)+i] = t2;
	  array[(j*size)+i+1] = t3;
	  array[((j+1)*size)+i+1] = t4;
	}
    }
}

