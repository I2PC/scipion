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
#include <XmippData/xmippProjection.hh>
#include <Reconstruction/projection.hh>
#include <Reconstruction/directions.hh>
#include <Reconstruction/Programs/Prog_symmetrize.hh>

typedef struct Column {double zero; double one; double two; double count;} column;

/**@name WBP */
//@{
/** WBP parameters. */
class Prog_WBP_prm {
public:
  /** Filenames */
  FileName fn_out, fn_sym, fn_sel;
  /** SelFile containing all projections */
  SelFile SF;
  /** If true: apply shifts upon reading the images (default) */
  bool apply_shifts;
  /** Lower threshold for the filter */
  double threshold;
  /** Diameter for reconstruction */
  int diameter;
  /** verbosity flag */
  int verb;
  /** dimensions of the images */
  int dim;
  /** Number of elements in matrix array */
  int no_mats;
  /** columns of matrices*/
  column * mat_g, * mat_f;
  /** Angular sampling for projection directions of arbitrary geometry filter */
  double sampling;
  /** Flag whether to use all experimental projection directions instead of
      sampled projection directions for arbitrary geometry filter */
  bool do_all_matrices;

public:
  /// Read arguments from command line
  void read(int argc, char **argv) _THROW;

  /// Show
  void show();

  /// Usage
  void usage();

  /// Produce side info: fill arrays with relevant transformation matrices
  void produce_Side_info() ;

  /// Fill array with transformation matrices needed for arbitrary geometry filter
  void get_all_matrices(SelFile &SF) ;

  /// Fill array with transformation matrices for representative 
  /// evenly sampled projection directions
  void get_sampled_matrices(SelFile &SF) ;

  // Simple (i.e. unfiltered) backprojection of a single image
  void simple_backprojection(Projection &img, VolumeXmipp &vol, 
			     int diameter) ;

  // Calculate the filter for arbitrary tilt geometry in 2D and apply 
  void apply_2Dfilter_arbitrary_geometry(SelFile &SF, VolumeXmipp &vol) ;

};
//@}
