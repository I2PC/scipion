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

#ifndef _PROG_SSNR
   #define _PROG_SSNR

#include <iostream>
#include <XmippData/xmippVolumes.hh>
#include <XmippData/xmippSelFiles.hh>

/**@name Spectral Signal to Noise Ratio */
//@{
/** SSNR parameters. */
class Prog_SSNR_prm {
public:
   /// Signal reconstructed volume
   FileName fn_S;
   
   /// Noise reconstructed volume
   FileName fn_N;
   
   /// Selfile with all the experimental images
   FileName fn_Ssel;
   
   /// Selfile with all the noise images
   FileName fn_Nsel;

   /// Filename of the Volumetric SSNR, used only for radial averaging
   FileName fn_VSSNR;

   /// Ringwidth
   double ring_width;   

   /// Sampling rate
   double Tm;

   /** Output filename.
       If empty, SSNR is inserted before the extension in fn_S */
   FileName fn_out;

   /** Generate SSNR images.*/
   bool generate_images;
   
   /** Min_power: Threshold for not dividing */
   double min_power;
public:
   /* Side info -------------------------------------------------------- */
   // Signal volume
   VolumeXmipp S;
   
   // Noise volume
   VolumeXmipp N;
   
   // Selfile with all experimental images
   SelFile SF_S;
   
   // Selfile with all noisy images
   SelFile SF_N;

   // SSNR3D for the radial_avg
   VolumeXmipp VSSNR;

public:
   /// Read parameters from command line
   void read(int argc, char **argv);
   
   /// Show parameters
   friend ostream & operator << (ostream &out, const Prog_SSNR_prm &prm);
   
   /// Usage
   void usage() const;
   
   /// Produce side Info
   void produce_side_info() _THROW;

   /** Estimate SSNR 1D.
    The columns of output are the following:
    Column 0: sample number in Fourier Space,
    Column 1: corresponding frequency in continuous freq (1/A),
    Column 2: corrected SSNR1D
    Column 3: Uncorrected Signal SSNR1D=S21D/N21D
    Column 4: S21D,
    Column 5: N21D
    Column 6: Uncorrected Noise SSNR1D=S21D/N21D
    Column 7: S21D,
    Column 8: N21D
    */
   void Estimate_SSNR_1D(matrix2D<double> &output);

   /** Estimate SSNR 2D.
       Generate images with the particular SSNR. The output filename
       is used as a rootname */
   void Estimate_SSNR_2D();

   /** Radial average of a Volumetric SSNR.
       The Volumetric SSNR is stored as 10*log10(VSSNR+1). To perform
       a radial average that is consistent with the one produced
       by the 1D estimation the +1 must be properly eliminated.

       The columns of output are the following:
       Column 0: sample number in Fourier Space,
       Column 1: corresponding frequency in continuous freq (1/A),
       Column 2: corrected radial_avg
   */
   void Radial_average(matrix2D<double> &output);
};

/** Compute the Uncorrected SSNR 1D.
    The input volume is projected in the same directions
    as the Selfile. Then for each image the SSNR1D is computed
    and finally all SSNR1D are averaged.
    
    The treat_as_noise field is useful for the computation of alpha
    
    The columns of output are the following:
    Column 0: sample number in Fourier Space,
    Column 1: corresponding frequency in continuous freq (1/A),
    Column 2: Uncorrected SSNR1D=S21D/N21D,
    Column 3: S21D,
    Column 4: N21D
    */
    void Compute_Uncorrected_SSNR_1D(
       matrix3D<double> &V, SelFile &SF,
       double ring_width, double Tm,
       matrix2D<double> &output, bool treat_as_noise);

/** Perform all the work.
    For the meaning of the output matrix look at the documentation
    of the function Estimate_SSNR_1D of the class Prog_SSNR_prm. */
    void ROUT_SSNR(Prog_SSNR_prm &prm, matrix2D<double> &output);
//@}
#endif
