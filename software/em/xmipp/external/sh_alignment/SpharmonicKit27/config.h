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


/*********************************************************************

  Define PRECOMP_DIR to be the directory where the precomputed data lives.

  DON'T FORGET that last backslash !!! (I.e. "/usr/local/" instead of
  "/usr/local" )

  For bandwidth bw = X , the seminaive algorithm's precomputed data
  files will be called

  Semi_bwX.dat     (* for the forward spherical transform *)
  InvSemi_bwX.dat  (* for the inverse spherical transform *)

  The hybrid algorithm's precomputed data will be called

  Hybrid_bwX.dat   (* for the forward spherical transform *)


  Example: if bw = 256 and PRECOMP_DIR = "/usr/local/data/" then
  it is assumed that the files
  
  Semi_bw256.dat
  InvSemi_bw256.dat
  Hybrid_bw256.dat

  will be written to and read from /usr/local/data/

  The default value of PRECOMP_DIR is the directory where the
  executables live.

  ********************************************************************/

#ifndef _CONFIG_H
#define _CONFIG_H

#define PRECOMP_DIR "./" 

#endif /* _CONFIG_H */
