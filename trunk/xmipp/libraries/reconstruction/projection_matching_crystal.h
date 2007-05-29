/***************************************************************************
 *
 * Authors:    Roberto Marabini
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

#include <reconstruction/projection.h>
#include <reconstruction/directions.h>
#include <reconstruction/symmetries.h>

#include <malloc.h>

#define FOR_ALL_DIRECTIONS() for (int dirno=0;dirno<nr_dir; dirno++)
#define FOR_ALL_ROTATIONS() for (int ipsi=0; ipsi<nr_psi; ipsi++ )

/**@name projection_matching_crystal */
//@{
/** projection_matching_crystal parameters. */
class Prog_projection_matching_crystal_prm
{
public:

    /** Selfile with reference images */
    SelFile SFref;
    /** Selfile with experimental images */
    SelFile SFexp;
    /** dimension of the images */
    int dim;
    /** Filename output rootname */
    FileName fn_root;
    /** Maximum distance to exp angles and shifts plus number
       of samples in psi and shift*/
    double psi_sampling, psi_distance;
    double rot_distance;
    double tilt_distance;
    double shift_sampling, shift_distance ;
    double scale_sampling, scale_distance ;
    /** Number of projection directions */
    int nr_dir;
    /** Vector with reference library projections */
    vector<Matrix2D<double> > ref_img;
    /** Vectors for standard deviation and mean of reference library projections */
    double * ref_stddev, * ref_mean;
    /** Vector with reference library angles */
    double * ref_rot, * ref_tilt;
    /** Vectors to stored valid shifts */
    vector <Matrix1D<double> > shift_vector;
    /** Maximum allowed shift */
    double max_shift;
    /** Flag whether to store optimal transformations in the image headers */
    bool modify_header;
#ifdef NEVERDEFINE
    vector<Matrix2D<double> >::iterator idirno;
    /** Number of steps to sample in-plane rotation in 90 degrees */
    int nr_psi;
    /** Verbose level:
        1: gives progress bar (=default)
        0: gives no output to screen at all */
    int verb;
    /** Flag whether to output reference projections, selfile and docfile */
    bool output_refs;
    /** Flag whether to input files are in Fourier, so far it is
    implemented only with the -ref flag */
    bool fourier_input;
    /** Angular sampling rate */
    double sampling;
    /** Mask for shifts */
    Matrix2D<int> rotmask;
    /** Number of white pixels in rotmask */
    int nr_pixels_rotmask;
    /** Inner and outer radii to limit the rotational search */
    double Ri, Ro;




    /// Extended Usage
    void extended_usage();



#endif

public:
    /** Actual projection matching for one image */
    void PM_process_one_image(Matrix2D<double> &Mexp,
                              float img_rot, float img_tilt, float img_psi,
                              float scale,
                              int &opt_dirno, double &opt_psi, double &opt_scale,
                              double &opt_xoff, double &opt_yoff,
                              double &maxCC);
    /// Read arguments from command line
    void read(int argc, char **argv);
    /// Show
    void show();
    /// Usage
    void usage();
    /** Make shiftmask and calculate nr_psi */
    void produce_Side_info();
    /** Loop over all images */
    void PM_loop_over_all_images(DocFile &DFo, double &sumCC);

};
//@}
