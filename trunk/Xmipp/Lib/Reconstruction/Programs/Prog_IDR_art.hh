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

#ifdef _HAVE_VTK

#include "Basic_art.hh"
#include "Prog_art.hh"
#include "Prog_FourierFilter.hh"

/**@name IDR with ART */
//@{
/* IDR Parameters ---------------------------------------------------------- */
/** IDR Parameters. */
class Prog_IDR_ART_Parameters {
public:
   /// ART Basic parameters
   Basic_ART_Parameters *art_prm;
   /// Number of IDR iterations
   int                  idr_iterations;
   /// Dont rewrite projections
   bool                 dont_rewrite;
   /// Only reproject
   FileName             fn_blob_volume;
   /// List of relaxation parameters: left border
   matrix1D<double>     mu0_list;
   /// List of relaxation parameters: right border
   matrix1D<double>     muF_list;
   /// List of relaxation parameters
   matrix1D<double>     mu_list;
   /// CTF filename
   FileName             fn_ctf;
   /// Iteration number
   int                  it;

   /// Side Info: CTF
   FourierMask          ctf;
   /// Side info: Original set of images
   SelFile              SF_original;
   /// Side info: Current set of images
   SelFile              SF_current;
   /// Side info: set of CTFs (maybe it is not necessary)
   SelFile              SF_ctf;
   /// Side info: multiple CTF mode
   bool                 multiple_CTFs;
   /// Side info: Root name of the original set of images
   FileName             fn_root;
   /// Side Info: Blob volume
   GridVolume           vol_blobs;
public:
   /// Empty creator
   Prog_IDR_ART_Parameters(): art_prm(NULL) {};

   /// Read parameters from file
   void read(const FileName &fn);
   
   /** Produce Side Information
       Exceptions are thrown if the number of images in the different
       selfiles do not match */
   void produce_side_info() _THROW;

   /// Show parameters
   void show();
   
   /// Usage
   void Usage();
   
   /** IDR relaxation parameter for iteration n.
       An exception is thrown if no mus are provided. */
   double mu(int n) _THROW {
      int imax=XSIZE(mu_list);
      if (imax==0)
         REPORT_ERROR(1,"IDR_ART: There are no mus\n");
      if (n>=imax) return mu_list(imax-1);
      else         return mu_list(n);
   }

   /** IDR correction.
       Given the reconstructed blob volume and the iteration step
       this routine change the current set of images to a new set
       rewriting or not the current set of images.
   */
   void IDR_correction(GridVolume &vol_blobs, int it);
};

/** Core of the IDR-ART routine.
    This is the routine which does everything except the TRUE IDR correction. */
void Basic_ROUT_IDR_Art(Prog_IDR_ART_Parameters &prm, VolumeXmipp &vol_recons);
//@}
#endif

#endif
