/***************************************************************************
 *
 * Authors:     Roberto Marabini roberto@cnb.uam.es
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
#include <data/projection.h>
#include <data/filters.h>

#include <reconstruction/directions.h>
#include <reconstruction/symmetrize.h>

/**@defgroup WBP reconstruct_wbp (Weighted Back Projection)
   @ingroup ReconsLibraryPrograms */
//@{
typedef struct Column
{
    double zero;
    double one;
    double two;
    double count;
}
column;

/** WBP parameters. */
class Prog_WBP_prm
{
public:
    /** Filenames */
    FileName fn_out, fn_sym, fn_sel, fn_doc, fn_control;
    /** SelFile containing all projections */
    SelFile SF;
    /** DocFile containing all angles */
    DocFile DF;
    /** Column numbers in the docfile */
    int col_rot, col_tilt, col_psi, col_xoff, col_yoff, col_flip, col_weight;
    /** Lower threshold for the filter */
    double threshold;
    /** Counter for how many times the threshold was not reached */
    int count_thr;
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
    /** Flag whether to use the weights in the image headers */
    bool do_weights;
    /** Symmetry list for symmetric volumes */
    SymList SL;

public:
    /// Read arguments from command line
    void read(int argc, char **argv);

    /// Show
    void show();

    /// Usage
    void usage();

    /// Produce side info: fill arrays with relevant transformation matrices
    void produce_Side_info() ;

    /// Get angles (either from reading the header or from a docfile)
    void get_angles_for_image(FileName fn, double &rot, double &tilt, double &psi,
                              double &xoff, double &yoff, double &flip, double &weight);

    /// Fill array with transformation matrices needed for arbitrary geometry filter
    void get_all_matrices(SelFile &SF) ;

    /// Fill array with transformation matrices for representative
    /// evenly sampled projection directions
    void get_sampled_matrices(SelFile &SF) ;

    // Simple (i.e. unfiltered) backprojection of a single image
    void simple_backprojection(Projection &img, VolumeXmipp &vol,
                               int diameter) ;

    // Calculate the filter and apply it to a projection
    void filter_one_image(Projection &proj);

    // Calculate the filter for arbitrary tilt geometry in 2D and apply
    void apply_2Dfilter_arbitrary_geometry(SelFile &SF, VolumeXmipp &vol) ;

};
//@}
