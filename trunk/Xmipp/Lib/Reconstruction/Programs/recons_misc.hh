/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or   
 * (at your option) any later version.                                 
 *                                                                     
 * This program is distributed in the hope that it will be useful,     
 * but WITHOUT ANY WARRANTY; without even the implied warranty of      
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
 * GNU General Public License for more details.                        
 *                                                                     
 * You should have received a copy of the GNU General Public License   
 * along with this program; if not, write to the Free Software         
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA            
 * 02111-1307  USA                                                     
 *                                                                     
 *  All comments concerning this program package may be sent to the    
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/
#ifndef _RECONS_MISC_HH
#  define _RECONS_MISC_HH

#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippSelFiles.hh>
#include "../symmetries.hh"
#include "../projection.hh"

/**@name Reconstruction Miscellanea */
//@{

/** Reconstruction information.
   This structure contains information for all projections which are
   going to participate in the reconstruction.
   This structure has also got
   information for the symmetry implementation. If there is any symmetry
   then an entry in this table is created using the same projection name
   but different symmetry matrices (only the matrix index is annotated
   in this structure). The Euler angles stored for the symmetrized image
   are the final ones, ie, the original Euler angles symmetrized according
   to the symmetry matrix. Then the symmetry identificator kept in this
   structure is only used to keep some track of what matrix was used to
   symmetrize.
*/
struct Recons_info {
   /// Projection filename
   FileName fn_proj;
   /// CTF filename
   FileName fn_ctf;
   /// Rotational angle
   float  rot;
   /// Tilting angle
   float  tilt;
   /// Psi angle
   float  psi;
   /** Symmetry number.
       This number express to which symmetry matrix this projection
       is related to (-1: without symmetry, 0: using symmetry matrix 0,
       1: using symmetry matrix 1 ...) */
   int    sym;
}; 

/** Build from a Selection File and a Symmetry List. 
    The result is stored in the Recons_info array which should point
    to NULL when it is not initialized. */
void build_recons_info(SelFile &selfile, SelFile &selctf, const FileName &fn_ctf,
   const SymList &SL, Recons_info * &IMG_Inf);
//@}

#endif
