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


#ifndef _PROG_DENOISING
   #define _PROG_DENOISING

#include "../xmippProgs.hh"

/**@name Denoising */
//@{
/** Denoising parameters */
class Denoising_parameters: public Prog_parameters {
public:
   typedef enum {REMOVE_SCALE,
         SOFT_THRESHOLDING,
	 BAYESIAN,
	 ADAPTIVE_SOFT,
	 CENTRAL,
	 SHAH} Denoising_type;

   /** Wavelet type.
       Valid types DAUB4, DAUB12, DAUB20 */
   string DWT_type;
   
   /** Denoising type.
       Valid types are REMOVE_SCALE, SOFT_THRESHOLDING,
	  BAYESIAN, ADAPTIVE_SOFT, CENTRAL, SHAH. */
   Denoising_type denoising_type;

   /** Scale to which the denoising is applied.
       It is used by remove scale, adaptive soft */
   int    scale;

   /** Threshold for soft thresholding*/
   double threshold;

   /** Radius for central */
   int    R;

   /** Shah number of outer iterations */
   int    Shah_outer;

   /** Shah number of inner iterations */
   int    Shah_inner;

   /** Shah number of refinement iterations */
   int    Shah_refinement;

   /** Shah weight.
       w0=data matching (=0) \\
       w1=1st derivative smooth (=50)\\
       w2=edge strength (=50)\\
       w3=edge smoothness (=0.02)*/
   matrix1D<double> Shah_weight;

   /** Produce Shah edge instead of Shah smooth. */
   bool   Shah_edge;
   
   /** Adjust range in Shah */
   bool   Shah_adjust_range;
public:
   /// Empty constructor
   Denoising_parameters();

   /// Read parameters from command line
   void read(int argc, char **argv) _THROW;
   
   /** Produce side info.
       The DWT type is translated and set */
   void produce_side_info() _THROW;

   /// Show parameters. This function calls show_specific
   void show();

   /// Show specific
   void show_specific();

   /// Usage. This function calls usage_specific
   void usage();

   /// Show specific parameters
   void usage_specific();

   /// Denoise an image
   void denoise(matrix2D<double> &img) _THROW;
   
   /// Denoise a volume
   void denoise(matrix3D<double> &vol) _THROW;
};
//@}
#endif
