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
#ifndef _PROG_IDR_ART_HH
#  define _PROG_IDR_ART_HH

#include "Basic_art.hh"
#include "Prog_art.hh"
#include "Prog_FourierFilter.hh"

/**@name IDR */
//@{
/* IDR Parameters ---------------------------------------------------------- */
/** IDR Parameters. */
class Prog_IDR_ART_Parameters {
public:
   /// Selfile with experimental images
   FileName fn_exp;
   /// Reference volume
   FileName fn_vol;
   /// Selfile with CTFs
   FileName fn_ctf;
   /// File output rootname
   FileName fn_root;
   /// Relaxation factor
   double mu;
   /// Adjust gray levels
   bool adjust_gray_levels;

   /// Side Info: CTF
   FourierMask          ctf;
   /// Side info: Original set of images
   SelFile              SF_exp;
   /// Side info: New set of images
   SelFile              SF_new;
   /// Side info: set of CTFs (maybe it is not necessary)
   SelFile              SF_ctf;
   /// Side Info: volume
   VolumeXmipp          V;
public:
   /// Read parameters from the command line
   void read(int argc, char **argv);
   
   /** Produce Side Information
       Exceptions are thrown if the number of images in the different
       selfiles do not match */
   void produce_side_info();

   /// Show parameters
   void show();
   
   /// Usage
   void Usage();
   
   /** IDR correction.
       This function process all images in the selfile.
   */
   void IDR_correction();
};
//@}
#endif
