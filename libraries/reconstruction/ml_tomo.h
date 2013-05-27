/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2004)
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

#ifndef _MLTOMO_H
#define _MLTOMO_H

#include <data/xmipp_fftw.h>
#include <data/xmipp_fft.h>
#include <data/args.h>
#include <data/xmipp_funcs.h>
#include <data/metadata_extension.h>
#include <data/xmipp_image.h>
#include <data/geometry.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/ctf.h>
#include <data/sampling.h>
#include <data/symmetries.h>
#include "symmetrize.h"
#include <data/xmipp_threads.h>
#include <vector>
#include <data/xmipp_program.h>

#define SIGNIFICANT_WEIGHT_LOW 1e-8
#undef SMALLANGLE
#define SMALLANGLE 2.75
#define MLTOMO_DATALINELENGTH 10
#define MLTOMO_BLOCKSIZE 10

#define FN_ITER_VOL(iter, base, refno) formatString("%s_it%d_%s%d.vol", fn_root.c_str(), iter, base, refno+1)
#define FN_ITER_MD(iter) formatString("%s_it%d.sel", fn_root.c_str(), iter)

class ProgMLTomo;

// Thread declaration
void * threadMLTomoExpectationSingleImage( void * data );

// This structure is needed to pass parameters to threadMLTomoExpectationSingleImage
typedef struct
{
    int thread_id;
    int thread_num;
    ProgMLTomo *prm;
    MetaData *MDimg;
    int *iter;
    double *wsum_sigma_noise;
    double *wsum_sigma_offset;
    double *sumfracweight;
    double *LL;
    std::vector<MultidimArray<double> > *wsumimgs;
    std::vector<MultidimArray<double> > *wsumweds;
    std::vector<Image<double> > *Iref;
    std::vector<MultidimArray<double> > *docfiledata;
    MultidimArray<double> * sumw;
    std::vector<size_t> * imgs_id;
    ThreadTaskDistributor * distributor;
}
structThreadExpectationSingleImage ;

/**@defgroup ml_tomo Maximum likelihood for tomograms
   @ingroup ReconsLibrary */
//@{
/** ml_tomo parameters. */
class ProgMLTomo: public XmippProgram
{
public:
    /** Filenames reference selfile/image, fraction MetaData & output rootname */
    FileName fn_sel, fn_ref, fn_root, fn_frac, fn_sym, fn_missing, fn_doc, fn_mask;
    /** Command line */
    std::string cline;
    /** Sigma value for expected pixel noise */
    double sigma_noise;
    /** sigma-value for origin offsets */
    double sigma_offset;
    /** Vector containing estimated fraction for each model */
    MultidimArray<double> alpha_k;
    /** Flag whether to fix estimates for model fractions */
    bool fix_fractions;
    /** Flag whether to fix estimate for sigma of origin offset */
    bool fix_sigma_offset;
    /** Flag whether to fix estimate for sigma of noise */
    bool fix_sigma_noise;
    /** Starting iteration */
    int istart;
    /** Number of iterations to be performed */
    int Niter, Niter2;
    /** dimension of the images */
    int oridim, dim, dim3, hdim;
    double ddim3;
    /** Number of reference images */
    int nr_ref;
    /** Keep angles from MetaData in generation of random subset averages */
    bool do_keep_angles;
    /** Total number of experimental images */
    size_t nr_images_global;
    /** Number of experimental images assigned to an MPI node */
    size_t nr_images_local;
    /** Index of the images to work an MPI node */
    size_t myFirstImg, myLastImg;
    /** Store data before write to md */
    MultidimArray<double > docfiledata;
    /** Sum of squared amplitudes of the references */
    std::vector<double> A2, corrA2;
    /** Stopping criterium */
    double eps;
    /** SelFile images, references and missingregions */
    MetaData MDimg, MDref, MDmissing;
    /** Flag whether md_contains_angles */
    bool mdimg_contains_angles;
    /** Vector for images to hold references (new & old) */
    std::vector < Image<double> > Iref, Iold, Iwed;
    /** Matrices for calculating PDF of (in-plane) translations */
    MultidimArray<double> P_phi, Mr2;
    /** Flag for generation of initial models from random subsets */
    bool do_generate_refs;
    /** Vector to assign reference number to stack positions*/
    std::vector <size_t> convert_refno_to_stack_position;

    /** Optimal refno and angno from previous iteration */
    std::vector<int> imgs_optrefno, imgs_optangno;
    std::vector<double> imgs_trymindiff, imgs_optpsi;
    std::vector<Matrix1D<double> > imgs_optoffsets;
    /** Number for missing data structure group */
    std::vector<int> imgs_missno;
    /** For initial guess of mindiff */
    double trymindiff_factor;
    /** Local search angular distance */
    double ang_search;
    /** Also limit psi search range */
    bool do_limit_psirange;
    /** Prohibit translations larger than this value */
    double limit_trans;
    /** Perturb angular sampling */
    bool do_perturb;
    /** Low-pass filter at FSC=0.5 resolution in each step */
    bool do_filter, do_ini_filter;
    /** Flag is true if symmetry > c1 */
    bool do_sym;

    /// Internally store all missing wedges or re-compute on the fly?

    /** Use missing wedges, cones or pyramids */
    bool do_missing;//, do_wedge, do_pyramid, do_cone;
    /** Do imputaion-like algorithm? */
    bool do_impute;
    /** Do maximum-likelihood algorithm? */
    bool do_ml;
    /** Threshold for no-imputation algorithm */
    double noimp_threshold;
    /** Maximum resolution (dig.freq.) */
    double maxres, scale_factor;
    MultidimArray<double> fourier_mask, real_mask, real_omask;
    MultidimArray<unsigned char> fourier_imask;

    /** Do not search any rotation/translation, only classify */
    bool dont_align;
    /** Do not search any rotation, only search translations and classify */
    bool dont_rotate;

    /** Classify only in defined mask */
    bool do_mask;
    /** Volume with the mask */
    Image<double> Imask;
    /** Number of non-zero voxels in the mask */
    double nr_mask_pixels;
    /** Just apply optimal rotations and use weighted-averaging to get class averages */
    bool do_only_average;

    // Missing data information
    enum MissingType { MISSING_WEDGE_Y, MISSING_WEDGE_X, MISSING_PYRAMID, MISSING_CONE };
    struct MissingInfo
    {
        MissingType type;
        double thy0, thyF, thx0, thxF;
        double tg0_y, tgF_y, tg0_x, tgF_x;
        bool do_limit_x, do_limit_y, do_cone;
        double nr_pixels;
    };
    /** Number of missing data structures */
    int nr_miss;
    /** vector to store all missing info */
    std::vector<MissingInfo> all_missing_info;

    // Angular sampling information
    struct AnglesInfo
    {
        size_t direction;
        double rot, rot_ori;
        double tilt, tilt_ori;
        double psi, psi_ori;
        Matrix2D<double> A, A_ori;
    };
    typedef std::vector<AnglesInfo> AnglesInfoVector;

    /// Missing data information

    /** Angular sampling  */
    double angular_sampling, psi_sampling;
    /** Vector with all angle combinations */
    AnglesInfoVector all_angle_info;
    /** Number of angle combinations */
    int nr_ang;
    /** Pixel size */
    double pixel_size;

    /** Regularization parameters */
    double reg0, regF, reg_current;
    int reg_steps;

    /** Switch off SMALL_ANGLE addition (for phantoms) */
    bool no_SMALLANGLE;

    // sampling object
    Sampling mysampling;
    // For user-provided tilt range
    double tilt_range0, tilt_rangeF;
    // Symmetry setup
    int symmetry, sym_order;

    /** Threads */
    int threads;

    /** FFTW objects */
    FourierTransformer transformer;

    /** debug flag */
    int debug;

    /** the vector with the images id in MetaData */
    std::vector<size_t> imgs_id;

public:
    /// Define the arguments accepted
    void defineParams();

    /// Read arguments from command line
    void readParams();

    /// Main body of the program
    void run();

    /// Show
    void show();

    /// Setup lots of stuff
    virtual void produceSideInfo();

    /// Generate initial references from random subset averages
    virtual void generateInitialReferences();

    /// Set the number of images, this function is useful only for MPI
    virtual void setNumberOfLocalImages();

    /** Read reference images in memory & set offset vectors
        (This produceSideInfo is Selfile-dependent!) */
    virtual void produceSideInfo2(int nr_vols = 1);

    /// Calculate Angular sampling
    void perturbAngularSampling();

    /// Calculate probability density distribution for in-plane transformations
    void calculatePdfTranslations();

    /** Read missing info from metadata
     * this function will read the metadata with the missing info and set the
     * number of missing regions
     */
    void readMissingInfo();

    /// Get binary missing wedge (or pyramid)
    void getMissingRegion(MultidimArray<unsigned char> &Mmeasured,
                          const Matrix2D<double> &A,
                          const int missno);

    void maskSphericalAverageOutside(MultidimArray<double> &Min);

    // Resize a volume, based on the max_resol
    // if down_scale=true: go from oridim to dim
    // if down_scale=false: go from dim to oridim
    void reScaleVolume(MultidimArray<double> &Min, bool down_scale=true);

    void postProcessVolume(Image<double> &Vin, double resolution = -1.);

    /// Fill vector of matrices with all rotations of reference
    void precalculateA2(std::vector< Image<double> > &Iref);

    /// ML-integration over all hidden parameters
    void expectationSingleImage(MultidimArray<double> &Mimg, int imgno, const int missno, double old_rot,
                                std::vector<Image<double> > &Iref,
                                std::vector<MultidimArray<double> > &wsumimgs,
                                std::vector<MultidimArray<double> > &wsumweds,
                                double &wsum_sigma_noise, double &wsum_sigma_offset,
                                MultidimArray<double> &sumw, double &LL, double &dLL,
                                double &fracweight, double &sumfracweight, double &trymindiff,
                                int &opt_refno, int &opt_angno, Matrix1D<double> &opt_offsets);

    /// Maximum constrained correlation search over all hidden parameters
    void maxConstrainedCorrSingleImage(MultidimArray<double> &Mimg, int imgno, int missno, double old_rot,
                                       std::vector<Image<double> > &Iref,
                                       std::vector<MultidimArray<double> > &wsumimgs,
                                       std::vector<MultidimArray<double> > &wsumweds,
                                       MultidimArray<double> &sumw, double &maxCC, double &sumCC,
                                       int &opt_refno, int &opt_angno, Matrix1D<double> &opt_offsets);

    /// Integrate over all experimental images
    virtual void expectation(MetaData &MDimg, std::vector< Image<double> > &Iref, int iter,
                     double &LL, double &sumfracweight,
                     std::vector<MultidimArray<double> > &wsumimgs,
                     std::vector<MultidimArray<double> > &wsumweds,
                     double &wsum_sigma_noise, double &wsum_sigma_offset,
                     MultidimArray<double> &sumw);

    /// Update all model parameters
    void maximization(std::vector<MultidimArray<double> > &wsumimgs,
                      std::vector<MultidimArray<double> > &wsumweds,
                      double &wsum_sigma_noise, double &wsum_sigma_offset,
                      MultidimArray<double> &sumw, double &sumfracweight,
                      double &sumw_allrefs, std::vector<MultidimArray<double> > &fsc, int iter);

    /// Calculate resolution by FSC
    void calculateFsc(MultidimArray<double> &M1, MultidimArray<double> &M2,
                      MultidimArray<double> &W1, MultidimArray<double> &W2,
                      MultidimArray< double >& freq, MultidimArray< double >& fsc,
                      double &resolution);

    /// Apply regularization
    void regularize(int iter);

    /// check convergence
    bool checkConvergence(std::vector<double> &conv);

    ///Add info of some processed images to later write to files
    virtual void addPartialDocfileData(const MultidimArray<double> &data, size_t first, size_t last);

    /// Write out reference images, selfile and logfile
    virtual void writeOutputFiles(const int iter,
                          std::vector<MultidimArray<double> > &wsumweds,
                          double &sumw_allrefs, double &LL, double &avefracweight,
                          std::vector<double> &conv, std::vector<MultidimArray<double> > &fsc);

};
//@}
#endif
