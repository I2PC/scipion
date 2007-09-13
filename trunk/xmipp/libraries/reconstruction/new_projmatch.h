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
#include <data/gridding.h>
#include <data/polar.h>

#include "projection.h"
#include "symmetries.h"
#include "sampling.h"

/**@name projection_matching */
//@{
/** projection_matching parameters. */
class Prog_new_projection_matching_prm {
public:

  /** Filenames reference selfile/image, fraction docfile & output rootname */
  FileName fn_img,fn_vol,fn_root;
  /** Selfile with experimental images */
  SelFile SF;
  /** The reference volume */
  VolumeXmipp vol;
  /** Vector with reference FTs of polar rings */
  vector<Polar<complex<double> > > fP_ref;
  /** Vector with reference images */
  vector<Matrix2D<double> > proj_ref;
  /** vector with stddevs for all reference projections */
  vector<double> stddev_ref;
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
    /** For user-provided tilt range */
  double tilt_range0, tilt_rangeF;
  /** Maximum allowed angular search ranges for rot and tilt */
  double ang_search;
  /** Inner and outer radii to limit the rotational search */
  int Ri, Ro;
  /** Flag to write class averages and class selfiles (one class for
   * each projection direction) */
  bool output_classes;
  /** Vector with all running class averages */
  vector<ImageXmipp> class_avgs;
  vector<SelFile> class_selfiles;

   /** Symmetry. One of the 17 possible symmetries in
      single particle electron microscopy.
      See details at url
      Possible values are: c1, ci, cs, cn, cnv, cnh, sn,
      dn, dnv, dnh, t, td, th, o, oh, i, ih */
   string symmetry;
   /** For infinite groups symmetry order*/
   int sym_order;
   /** sampling object */
   XmippSampling mysampling;
   /** One common kb object for all images! */
   KaiserBessel kb;
   /** Reproject reference projections for translational search
    * instead of storing them all in memory */
   bool do_reproject;


public:
  /// Read arguments from command line
  void read(int argc, char **argv);

  /// Show
  void show();

  /// Usage
  void usage();

  /// Extended Usage
  void extended_usage();

  /** Make shiftmask and calculate nr_psi */
  void produce_Side_info();

  /** Rotational alignment using polar coordinates 
   *  The input image is assumed to be in FTs of polar rings 
   */
  void rotationally_align_one_image(const Matrix2D<double> &img,
				    const int &samplenr, int &opt_samplenr,  
				    double &opt_psi, double &opt_flip, double &maxcorr);

  /** Translational alignment using cartesian coordinates 
   *  The optimal direction is re-projected from the volume
   */
  void translationally_align_one_image(const Matrix2D<double> &img,
				       const int &samplenr, const double &psi, const double &opt_flip,
				       double &opt_xoff, double &opt_yoff, double &maxcorr);

  /** Loop over all images */
  void PM_loop_over_all_images(SelFile &SF, DocFile &DFo, double &sumCC);


    /** Write to disc all class averages and selfiles */
    void write_classes();
};				
//@}
