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


/* externally available commands from newFCT.c */

#ifndef _NEWFCT_H
#define _NEWFCT_H

/* structure for making things parallel */
#ifndef _STRUCT_LOWHIGH
#define _STRUCT_LOWHIGH

struct lowhigh{
  double low;
  double high;
} ;

#endif /* _STRUCT_LOWHIGH */

extern void ExpIFCT( double * ,
		     double * ,
		     double * ,
		     int ,
		     int ,
		     int );

#ifndef FFTPACK   /* use local dcts */

extern void kFCT(double *,
		 double *,
		 double *,
		 int ,
		 int ,
		 int );

extern void kFCTX( double *,
		   double *,
		   double *,
		   int ,
		   int ,
		   int ,
		   struct lowhigh *,
		   struct lowhigh *);

#else   /* use FFTPACK-based routines */

extern double *precomp_dct( int );

extern void DCTf( double *,
		  int ,
		  int ,
		  double * );

extern void DCTb( double *,
		  int ,
		  double * );

#endif /* FFTPACK */

#endif /* _NEWFCT_H */
