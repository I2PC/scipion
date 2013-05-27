/***************************************************************************
 *
 * Authors:     Roberto Marabini roberto@cnb.csic.es
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/xmipp_fft.h>
#include <data/args.h>
#include <data/xmipp_funcs.h>
#include <data/metadata_extension.h>
#include <data/xmipp_image.h>
#include <data/projection.h>
#include <data/filters.h>

#include <reconstruction/recons.h>

#include <reconstruction/directions.h>
#include <reconstruction/symmetrize.h>

/**@defgroup WBP reconstruct_wbp (Weighted Back Projection)
   @ingroup ReconsLibrary */
//@{
typedef struct
{
    double x; // Projection direction (x,y,z)
    double y;
    double z;
    double count;
}
WBPInfo;

/** WBP parameters. */
class ProgRecWbp: public ProgReconsBase
{
public:
    /** Filenames */
    FileName fn_out, fn_sym, fn_sel;
    /** SelFile containing all projections and angles */
    MetaData SF;
    /** Lower threshold for the filter */
    double threshold;
    /** Counter for how many times the threshold was not reached */
    int count_thr;
    /** Diameter for reconstruction */
    int diameter;
    /** verbosity flag */
    //int verb;
    /** dimensions of the images */
    size_t dim;
    /** Number of elements in matrix array */
    int no_mats;
    /** columns of matrices*/
    WBPInfo * mat_g, * mat_f;
    /** Angular sampling for projection directions of arbitrary geometry filter */
    double sampling;
    /** Flag whether to use all experimental projection directions instead of
        sampled projection directions for arbitrary geometry filter */
    bool do_all_matrices;
    /** Flag whether to use the weights in the image headers */
    bool do_weights;
    /** Symmetry list for symmetric volumes */
    SymList SL;
    /// Time bar variables
    size_t time_bar_step, time_bar_size, time_bar_done;
    /// Iterator over input metadata
    MDIterator * iter;
    /// Reconstructed volume
    Image<double> reconstructedVolume;
public:

    ProgRecWbp();
    ~ProgRecWbp();

    /// Read arguments from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /** Do the job */
    void run();

    /** Show progress */
    virtual void showProgress();

    /** Finish processing */
    virtual void finishProcessing();

    /** Set IO for a new reconstruction*/
    void setIO(const FileName &fn_in, const FileName &fn_out);

    /// Fill arrays with relevant transformation matrices
    virtual void produceSideInfo() ;

    /// Get 1 image to process
    virtual bool getImageToProcess(size_t &objId, size_t &objIndex);

    /// Get angles (either from reading the header or from a docfile)
    void getAnglesForImage(size_t id, double &rot, double &tilt, double &psi,
                              double &xoff, double &yoff, double &flip, double &weight);

    /// Fill array with transformation matrices needed for arbitrary geometry filter
    void getAllMatrices(MetaData &SF) ;

    /// Fill array with transformation matrices for representative
    /// evenly sampled projection directions
    void getSampledMatrices(MetaData &SF) ;

    // Simple (i.e. unfiltered) backprojection of a single image
    void simpleBackprojection(Projection &img, MultidimArray<double> &vol,
                               int diameter) ;

    // Calculate the filter and apply it to a projection
    void filterOneImage(Projection &proj, Tabsinc &TSINC);

    // Calculate the filter for arbitrary tilt geometry in 2D and apply
    void apply2DFilterArbitraryGeometry() ;
};
//@}
