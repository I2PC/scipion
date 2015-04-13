/***************************************************************************
 *
 * Authors:    Sjors Scheres            scheres@cnb.csic.es (2004)
 * Authors:    Roberto Marabini            roberto@cnb.csic.es (2012)
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

#ifndef _angular_projection_matching_H
#define _angular_projection_matching_H

#include <data/args.h>
#include <data/xmipp_funcs.h>
#include <data/metadata.h>
#include <data/xmipp_image.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/polar.h>
#include <data/xmipp_fftw.h>
#include <data/xmipp_threads.h>
#include <pthread.h>

#include <data/projection.h>
#include <data/symmetries.h>
#include <data/sampling.h>
#include <data/ctf.h>

#include <data/xmipp_program.h>

#define MY_OUPUT_SIZE 10

class ProgAngularProjectionMatching;

// Thread declaration
void * threadRotationallyAlignOneImage( void * data );

// This structure is needed to pass parameters to threadRotationallyAlignOneImage
typedef struct{
    size_t thread_id;
    ProgAngularProjectionMatching *prm;
    MultidimArray<double> *img;
    size_t this_image;
    int opt_refno;
    double opt_psi;
    bool opt_flip;
    double maxcorr;
} structThreadRotationallyAlignOneImage ;

/**@defgroup angular_projection_matching new_projmatch (Discrete angular assignment using a new projection matching)
   @ingroup ReconsLibrary */
//@{
/** projection_matching parameters. */
class ProgAngularProjectionMatching : public XmippProgram
{

public:

    /** Filenames */
    FileName fn_exp, fn_ref, fn_out, fn_ctf;
    /** Docfile with experimental images */
    MetaData DFexp;
    /** Docfile with results */
    MetaData DFo;
    /** dimension of the images and padded images */
    size_t dim, paddim;
    /** Padding factor (only for applying CTF to references) */
    double pad;
    /** Maximum allowed shift */
    double max_shift;
    /** Inner and outer radii to limit the rotational search */
    int Ri, Ro;
    /** Available memory for storage of all references (in Gb) */
    double avail_memory;
    /** Maximum number of references to store in memory */
    int max_nr_refs_in_memory;
    /** Maximum number of references that can be stored in memory 
    The difference between this and the previous varible is that
    the previous one is the MIN(max_nr_refs_in_memory, number of references)
    */
    int max_nr_imgs_in_memory;
    /** Total number of references */
    int total_nr_refs;
    /** Counter for current filling of memory with references */
    int counter_refs_in_memory;
    /** Pointers for reference retrieval */
    std::vector<int> pointer_allrefs2refsinmem;   //order after removing redundant
    std::vector<int> pointer_refsinmem2allrefs;   //order in memory
    /** Vector to assign reference number to stack positions*/
    std::vector <size_t> convert_refno_to_stack_position;
    /** Array containing the images ids in metadata */
    std::vector<size_t> ids;
    /** Array with Polars of references and of translated images and their mirrors */
    Polar<std::complex<double> >   *fP_ref, *fP_img, *fPm_img;
    /** Array with reference images */
    MultidimArray<double> *proj_ref;
    /** Global plans for fftw transformers of all polar rings */
    Polar_fftw_plans global_plans;
    /** vector with stddevs for all reference projections */
    double *stddev_ref, *stddev_img;
    /** sampling object */
    Sampling mysampling;
    /** Flag whether to loop from low to high or from high to low
     * through the references */
    bool loop_forward_refs;
    /** 5D-search: maximum offsets (+/- pixels) */
    int search5d_shift;
    /** 5D-search: offset step (pixels) */
    int search5d_step;
    /** 5D-search: actual displacement vectors */
    std::vector<int> search5d_xoff, search5d_yoff;
    /** CTF image */
    MultidimArray<double> Mctf;
    /** Are the experimental images phase flipped? */
    bool phase_flipped;
    /** Threads */
    int threads;
    /** Number of translations in 5D search */
    size_t nr_trans;
    /** Thread barrier */
    barrier_t thread_barrier;

    /** scale params */
    bool do_scale;
    WriteModeMetaData do_overwrite;
    double scale_step;
    double scale_nsteps;

    /** Progress bar */
    int progress_bar_step;
public:

    /// Read arguments from command line
    virtual void readParams();

    /// Define arguments accepted
    virtual void defineParams();

    /** Show. */
    void show();

    /** Run. */
    void run();

    /** Make shiftmask and calculate nr_psi */
    virtual void produceSideInfo();

    /** Rotational alignment using polar coordinates
     *  The input image is assumed to be in FTs of polar rings
     */
    void rotationallyAlignOneImage(Matrix2D<double> &img, int imgno, int &opt_samplenr,
    		double &opt_psi, bool &opt_flip, double &maxcorr);

    /** Translational alignment using cartesian coordinates
     *  The optimal direction is re-projected from the volume
     */
    void translationallyAlignOneImage(MultidimArray<double> &img,
    		const int &samplenr, const double &psi, const bool &opt_flip,
    		double &opt_xoff, double &opt_yoff, double &maxcorr);

    /** Translational alignment using cartesian coordinates
         *  The optimal direction is re-projected from the volume
         */
    void scaleAlignOneImage(MultidimArray<double> &img,
        		const int &samplenr, const double &psi, const bool &opt_flip,
        		const double &opt_xoff, const double &opt_yoff, const double &old_scale, double &opt_scale, double &maxcorr);

    /** Get pointer to the current reference image
      If this image wasn't stored in memory yet, read it from disc and
      store FT of the polar transform as well as the original image */
    void getCurrentReference(int refno, Polar_fftw_plans &local_plans);

    /** Get images to process.
     * This function will return the id's of images to process.
     * It will be specially useful for MPI case when images will be distributed
     * by the master node.
     * Return false if there aren't more images to process
     */
    virtual void processAllImages();

    /** Loop over all images */
    void processSomeImages(const std::vector<size_t> &imagesToProcess);

    /** Read current image into memory and translate accoring to
      previous optimal Xoff and Yoff */
    void getCurrentImage(size_t imgid, Image<double> &img);

    /** Write out results to disk
     * This function should be override in MPI class, only master should write.
     */
    virtual void writeOutputFiles();

    /** destroy and clean **/
    void destroyAndClean();
};				
//@}
#endif
