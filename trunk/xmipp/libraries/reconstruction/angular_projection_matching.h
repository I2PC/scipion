/***************************************************************************
 *
 * Authors:    Sjors Scheres            scheres@cnb.uam.es (2004)
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
#include <data/filters.h>
#include <data/mask.h>

#include "projection.h"
#include "directions.h"
#include "symmetries.h"

#define FOR_ALL_DIRECTIONS() for (int dirno=0;dirno<nr_dir; dirno++)
#define FOR_ALL_ROTATIONS() for (int ipsi=0; ipsi<nr_psi; ipsi++ )

/**@defgroup ProjectionMatching angular_projection_matching (Discrete angular assignment by projection matching)
   @ingroup ReconsLibraryPrograms */
//@{
/** projection_matching parameters. */
class Prog_projection_matching_prm {
public:

  /** Filenames reference selfile/image, fraction docfile & output rootname */
  FileName fn_vol,fn_root,fn_sym,fn_ref,fn_ang,fn_control;
  /** Selfile with experimental images */
  SelFile SF;
  /** Vector with reference library projections */
  std::vector<Matrix2D<double> > ref_img;
  std::vector<Matrix2D<double> >::iterator idirno;
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
  /** Create Projections
      1 yes
      0 no
   **/
  int create_proyections;
   
  /** Flag whether to store optimal transformations in the image headers */
  bool modify_header;
  /** Flag whether to output reference projections, selfile and docfile */
  bool output_refs;
  /** Angular sampling rate */
  double sampling;
  /** Maximum allowed shift */
  double max_shift;
    /** For user-provided tilt range */
  double tilt_range0, tilt_rangeF;
  /** Maximum allowed angular search ranges for rot and tilt */
  double ang_search;
  /** Mask for shifts */
  Matrix2D<int> rotmask;
  /** Number of white pixels in rotmask */
  int nr_pixels_rotmask;
  /** Inner and outer radii to limit the rotational search */
  double Ri, Ro;
    /** Flag to write class averages and class selfiles (one class for
     * each projection direction) */
    bool output_classes;
    /** Vector with all running class averages */
    std::vector<ImageXmipp> class_avgs;
    std::vector<SelFile> class_selfiles;

public:
  /// Read arguments from command line
  virtual void read(int argc, char **argv);

  /// Show
  virtual void show();

  /// Usage
  virtual void usage();

  /// Extended Usage
  void extended_usage();

  /** Make shiftmask and calculate nr_psi */
  void produce_Side_info();

  /** Actual projection matching for one image */
  void PM_process_one_image(Matrix2D<double> &Mexp,
			    float &img_rot, float &img_tilt, float &img_psi,
			    int &opt_dirno, double &opt_psi,
			    double &opt_xoff, double &opt_yoff,
			    double &maxCC, double &Zscore);

  /** Loop over all images */
  void PM_loop_over_all_images(SelFile &SF, DocFile &DFo, double &sumCC);


    /** Write to disc all class averages and selfiles */
    void write_classes();
};				
//@}
