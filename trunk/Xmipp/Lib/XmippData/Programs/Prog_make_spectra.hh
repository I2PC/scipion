/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2002)
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


#ifndef _PROG_MAKE_SPECTRA
   #define _PROG_MAKE_SPECTRA

#include "../xmippRotationalSpectrum.hh"
#include "../xmippProgs.hh"
#include "../xmippImages.hh"
#include <vector>

/**@name Make spectra */
//@{
class Prog_make_spectra_prm: public Prog_parameters {
public:
   /// Output filename
   FileName fn_out;

   /// Rotational spectrum
   Rotational_Spectrum rot_spt;
public:
   // Set of harmonics
   vector< matrix1D<double> > Harmonics;
   // Set of images
   vector<FileName> Img_name;
public:
   /// Empty constructor
   Prog_make_spectra_prm();

   /// Read parameters from command line
   void read(int argc, char **argv) _THROW;
   
   /** Produce side info. */
   void produce_side_info() _THROW;

   /// Show parameters. This function calls show_specific
   void show();

   /// Show specific
   void show_specific();

   /// Usage. This function calls usage_specific
   void usage();

   /// Show specific parameters
   void usage_specific();

   /// Process image
   void process_img(ImageXmipp &img);

   /// Finish processing
   void finish_processing();
};
//@}
#endif
