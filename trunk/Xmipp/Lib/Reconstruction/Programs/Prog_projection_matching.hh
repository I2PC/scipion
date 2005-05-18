/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.uam.es (2004)
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
#include <XmippData/xmippFFT.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippDocFiles.hh>
#include <XmippData/xmippImages.hh>
#include <XmippData/xmippFilters.hh>
#include <XmippData/xmippMasks.hh>
#include <Reconstruction/projection.hh>
#include <Reconstruction/directions.hh>
#include <Reconstruction/symmetries.hh>
#include <malloc.h>

#define FOR_ALL_DIRECTIONS() for (int dirno=0;dirno<nr_dir; dirno++)
#define FOR_ALL_ROTATIONS() for (int ipsi=0; ipsi<(nr_psi); ipsi++ ) 

/**@name projection_matching */
//@{
/** projection_matching parameters. */
class Prog_projection_matching_prm {
public:

  /** Filenames reference selfile/image, fraction docfile & output rootname */
  FileName fn_vol,fn_root,fn_sym,fn_ref;
  /** Selfile with experimental images */
  SelFile SF;
  /** Vector with reference library projections */
  vector<matrix2D<double> > ref_img;
  vector<matrix2D<double> >::iterator idirno;
  /** Vectors for standard deviation and mean of reference library projections */
  double * ref_stddev, * ref_mean;
  /** Vector with reference library angles */
  double * ref_rot, * ref_tilt;
  /** Number of projection directions */
  int nr_dir;
  /** Number of steps to sample in-plane rotation in 90 degrees */
  int nr_psi;
  /** dimension of the images */
  int dim;
  /** Verbose level:
      1: gives progress bar (=default)
      0: gives no output to screen at all */
  int verb;
  /** Flag whether to store optimal transformations in the image headers */
  bool modify_header;
  /** Flag whether to output reference projections, selfile and docfile */
  bool output_refs;
  /** Angular sampling rate */
  double sampling;
  /** Maximum allowed shift */
  double max_shift;
  /** Maximum allowed angular search ranges */
  double rot_range, tilt_range, psi_range;
  /** Mask for shifts */
  matrix2D<int> shiftmask;

public:
  /// Read arguments from command line
  void read(int argc, char **argv) _THROW;
  
  /// Show
  void show();

  /// Usage
  void usage();

  /// Extended Usage
  void extended_usage();

  /** Make shiftmask and calculate nr_psi */
  void produce_Side_info() _THROW;

  /** Actual projection matching for one image */
  void PM_process_one_image(matrix2D<double> &Mexp,
			    float &img_rot, float &img_tilt, float &img_psi, 
			    int &opt_dirno, double &opt_psi,
			    double &opt_xoff, double &opt_yoff, 
			    double &maxCC) _THROW;

  /** Loop over all images */
  void PM_loop_over_all_images(SelFile &SF, DocFile &DFo, double &sumCC) _THROW;


};				    
//@}
