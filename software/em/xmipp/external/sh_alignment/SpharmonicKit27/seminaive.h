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


#ifndef _SEMINAIVE_H
#define _SEMINAIVE_H

extern void SemiNaive( double *,
		       int , 
		       int , 
		       int , 
		       double *,
		       double *,
		       int ,
		       double *,
		       int , 
		       double * ) ;

extern void SemiNaiveX( double *,
			int , 
			int , 
			int , 
			double *,
			double * ,
			double * ,
			int ,
			double * ,
			int ) ;

#ifndef FFTPACK
extern void SemiNaiveReduced( double * ,
			      int , 
			      int , 
			      double * , 
			      double * ,
			      double * ,
			      double * ) ;
#else
extern void SemiNaiveReduced( double * ,
			      int , 
			      int , 
			      double * , 
			      double * ,
			      double * ,
			      double * ,
			      double * ) ;
#endif /* FFTPACK */


extern void SemiNaiveReducedX( double * ,
			       int , 
			       int ,
			       int ,
			       double * ,
			       double * ,
			       double * ) ;

extern void InvSemiNaiveReduced( double * ,
				 int , 
				 int , 
				 double * , 
				 double * ,
				 double * ,
				 double * ) ;

#endif /* _SEMINAIVE_H */
