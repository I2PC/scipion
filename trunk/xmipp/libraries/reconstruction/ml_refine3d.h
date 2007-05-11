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
#include <data/fft.h>
#include <data/args.h>
#include <data/funcs.h>
#include <data/selfile.h>
#include <data/docfile.h>
#include <data/image.h>
#include <data/volume.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/morphology.h>
#include <data/volume_segment.h>

#include "grids.h"
#include "symmetries.h"
#include "blobs.h"
#include "projection.h"
#include "directions.h"
#include "reconstruct_art.h"
#include "reconstruct_wbp.h"
#include "ml_align2d.h"
#include "symmetrize.h"

#include <vector>


/**@name Refine3d */
//@{
/** Refine3d parameters. */
class Prog_Refine3d_prm {

public:
  // Filename for reference volume, symmetry file and output rootname
  FileName fn_sel, fn_vol, fn_sym, fn_root, fn_solv, fn_iter, fn_symmask;
  // Selfile with reference volumes
  SelFile SFvol;
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
   /// File handler for the history file
  ofstream fh_hist;
  // Use WBP instead of WLS-ART for reconstruction in ML
  bool do_wbp;
  // Low-pass filter digital frequency
  double lowpass;
  // For user-provided tilt range
  double tilt_range0, tilt_rangeF;
  // Do not use a starting volume in wlsART reconstruction
  bool wlsart_no_start;
  // Threshold for flooding-like solvent mask
  double threshold_solvent;
  // Do Wang/Terwilliger-like probabilistic solvent flattening
  bool do_prob_solvent;
  // Do separate_object-like selection of largest connected volume in
  // binarized solvent mask
  bool do_deblob_solvent;
  // Dilate (binarized) solvent mask?
  int dilate_solvent;
  // Use MLF mode
  bool fourier_mode;

public:

  /// Read additional arguments for 3D-process from command line
  void read(int &argc, char ** &argv) ;

  /// ML Usage
  void usage();

  /// MLF Usage
  void MLF_usage();

  /// Additional options
  void extended_usage();

  /// Show
  void show();

  /// Project the reference volume in evenly sampled directions
  void project_reference_volume(SelFile &SFlib, int rank=0) ;

  /// (For mpi-version only:) calculate noise averages and write to disc
  void make_noise_images(vector<ImageXmipp> &Iref) ;

  /// reconstruction by (weighted ART) or WBP
  void reconstruction(int argc, char **argv,
		      int iter, int volno, int noise=0);

 /// Calculate 3D SSNR according to Unser ea. (2005)
  void calculate_3DSSNR(matrix1D<double> &spectral_signal, int iter);

  /// After reconstruction update reference volume selfile
  void remake_SFvol(int iter, bool rewrite=false, bool include_noise=false) ;

  /// Merge MLalign2D classification selfiles into volume classes
  void concatenate_selfiles(int iter);

  /// Maksing, filtering etc. of the volume
  void post_process_volumes(int argc, char **argv) ;

  /// Convergency check
  bool check_convergence(int iter) ;

};
//@}
