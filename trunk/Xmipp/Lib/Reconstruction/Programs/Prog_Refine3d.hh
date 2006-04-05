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
#include <XmippData/xmippVolumes.hh>
#include <XmippData/xmippFilters.hh>
#include <XmippData/xmippMasks.hh>
#include "Reconstruction/grids.hh"
#include "Reconstruction/symmetries.hh"
#include "Reconstruction/blobs.hh"
#include "Reconstruction/projection.hh"
#include "Reconstruction/directions.hh"
#include <Reconstruction/Programs/Prog_art.hh>
#include <Reconstruction/Programs/Prog_WBP.hh>
#include <Reconstruction/Programs/Prog_MLalign2D.hh> 
#include <Reconstruction/Programs/Prog_symmetrize.hh>
#include <vector>


/**@name Refine3d */
//@{
/** Refine3d parameters. */
class Prog_Refine3d_prm {

public:
  // Filename for reference volume, symmetry file and output rootname
  FileName fn_sel, fn_vol, fn_sym, fn_root, fn_solv, fn_iter;
  // Selfile with experimental images and reference volumes
  SelFile SF, SFvol;
  // Number of volumes to refine
  int Nvols;
  // vector with integers which projections are valid for which volume
  vector<int> eachvol_start, eachvol_end;
  // Iteration numbers
  int istart, Niter;
  // Verbosity flag
  int verb;
  // Convergence check
  double eps;
  // Angular sampling interval (degree)
  double angular;
  // Low-pass filter digital frequency
  double lowpass;
  // Radius for masking of volume 
  double mask_radius;
  // Initial estimate for maximum resolution
  double ini_maxres;
   /// File handler for the history file
  ofstream fh_hist;
  // Use WBP instead of WLS-ART for reconstruction in ML
  bool do_wbp;
  // For user-provided tilt range
  double tilt_range0, tilt_rangeF;

public:

  /// Read additional arguments for 3D-process from command line
  void read(int argc, char **argv) ;

  /// Usage
  void usage();

  /// Show
  void show();

  /// Project the reference volume in evenly sampled directions
  void project_reference_volume(SelFile &SFlib, int rank=0) ;

  /// reconstruction by (weighted ART) or WBP
  void reconstruction(int argc, char **argv, 
		      int iter, int volno);

  /// After reconstruction update reference volume selfile
  void remake_SFvol(int iter, bool rewrite=false) ;

  /// Maksing, filtering etc. of the volume
  void post_process_volumes(int argc, char **argv) ;

  /// Convergency check
  bool check_convergence(int iter) ;

};
//@}
