/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.uam.es (2002)
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
#ifndef _PROG_ALIGN2D
   #define _PROG_ALIGN2D

#include <XmippData/xmippFFT.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippDocFiles.hh>
#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippImages.hh>
#include <vector>

/**@name Align2D */
//@{
/** Align2D parameters. */
class Prog_align2d_prm {
public:
  /** Filename selection file containing the images */
  FileName fn_sel;
  /** Filename reference image */
  FileName fn_ref;
  /** Filename output reference image */
  FileName fn_ave;
  /**  Filename output document file */
  FileName fn_doc;
  /**  Filename output extension */
  FileName oext;

  /** Integer inner radius for rotational alignment */
  int Ri;
  /** Integer outer radius for rotational alignment */
  int Ro;
  /** Integer number of iterations to perform */
  int Niter;
  /** Integer number of images to be aligned */
  int n_images;

  /** Float maximum allowed shift (discard images that shift more in last iteration)*/
  float max_shift;
  /** Float maximum allowed rotational change (discard images that rotate more in last iteration)*/
  float max_rot;
  /** Float resolution limit for low-pass filter [Angstrom]*/
  float resol;
  /** Float sampling rate (i.e. pixel size) [Angstrom] */
  float sam;

  /** Boolean to apply low-pass filter to images */
  bool do_filter;
  /** Boolean to perform rotational alignment */
  bool do_rot;
  /** Boolean to perform translational alignment */
  bool do_trans;

  /** Boolean to perform complete-search (psi and translations) alignment */
  bool do_complete;
  /** Float psi-interval in complete search */
  float psi_interval;

public:
  // SelFile images
  SelFile SF;
  // Stack of input images
  vector<ImageXmipp>  images;  
  // Stack of optimal correlations for all images
  vector<double>  corr;  
  // Boolean for successful alignment of image
  vector<bool>  success;
  // Image holding current reference
  ImageXmipp Iref;

public:
  /// Read argument from command line
  void read(int argc, char **argv) _THROW;
  
  /// Show
  void show();

  /// Usage
  void usage();
  
  /// Rotational alignment of an image
  bool align_rot(ImageXmipp &img, const matrix2D<double> &Mref, 
		 const float &max_rot, const float &Rin, const float &Rout, const double &outside=0.) _THROW;

   /// Translational alignment of an image
  bool align_trans(ImageXmipp &img, const matrix2D<double> &Mref, const float &max_shift, const double &outside=0.) _THROW;

  /// Alignment by complete search of rotations and translations
  bool align_complete_search(ImageXmipp &img, const matrix2D<double> &Mref, 
                       const float &max_shift, const float &max_rot, const float &psi_interval, 
		       const float &Rin, const float &Rout, const double &outside=0.) _THROW;

  /// Piramidal combination of images to construct a reference
  void do_pspc() _THROW;

  /// Alignment of all images by iterative refinement 
  void refinement() _THROW;

  /// Calculate optimal correlation for in document file
  void calc_correlation(const matrix2D<double> &Mref, const float &Rin, const float &Rout) _THROW;

  /// Main routine
  void align2d() _THROW;

};				    
//@}
#endif
