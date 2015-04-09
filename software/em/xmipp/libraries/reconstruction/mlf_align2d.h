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

#ifndef MLFALIGN2D_H
#define MLFALIGN2D_H

#include "ml2d.h"
#include <numeric>

/**@defgroup MLFalign2D mlf_align2d (Maximum likelihood in 2D in Fourier space)
   @ingroup ReconsLibrary */
//@{

#define FOR_ALL_MODELS() for (int refno=0;refno<model.n_ref; refno++)
#define FOR_ALL_ROTATIONS() for (size_t ipsi=0; ipsi<nr_psi; ipsi++ )
#define FOR_ALL_FLIPS() for (size_t iflip=0; iflip<nr_flip; iflip++)
#define FOR_ALL_LIMITED_TRANSLATIONS() for (size_t itrans=0; itrans<nr_trans; itrans++)
#define FOR_ALL_DEFOCUS_GROUPS() for (size_t ifocus=0; ifocus<nr_focus; ifocus++)
#define FOR_ALL_DIGITAL_FREQS() for (size_t irr = 0; irr < hdim; irr++)
#define FOR_ALL_POINTS() for (size_t ipoint = 0; ipoint < nr_points_2d; ++ipoint)

//Helper macros to define the elements of vectors:
//Vsig, Vctf and Vsnr indexed by ifocus and irr
#define VSNR_ITEM dAi(Vsnr[ifocus], irr)
#define VCTF_ITEM dAi(Vctf[ifocus], irr)
#define VDEC_ITEM dAi(Vdec[ifocus], irr)
#define VSIG_ITEM dAi(Vsig[ifocus], irr)

#define SIGNIFICANT_WEIGHT_LOW 1e-8
#define SMALLVALUE 1e-4
#define HISTMIN -6.
#define HISTMAX 6.
#define HISTSTEPS 120

/** Some filename convetions for output files */
#define FN_EXTRA(file) formatString("%sextra/%s", fn_root.c_str(), file)
#define FN_NOISE_IMG_MD FN_EXTRA("noise_images.xmd")
#define FN_NOISE_IMG    FN_EXTRA("noise_images.stk")
#define FN_CREF_IMG     FN_EXTRA("cref_classes.stk")
#define FN_CREF_IMG_MD  FN_EXTRA("cref_classes.xmd")

#define FN_ITER_BASE(iter) getIterExtraPath(fn_root, iter)
#define FN_REF(base, refno) formatString("%06d@%sclasses.stk", (refno), (base).c_str())
#define FN_VSIG(base, ifocus, ext) ((nr_focus > 1) ? formatString("ctf%06d@%s%s", ((ifocus) + 1), (base).c_str(), (ext)) : ((base) + "_ctf" + (ext)))


/** MLFalign2D parameters. */
class ProgMLF2D: public ML2DBaseProgram
{
public:

    /** sigma-value for origin offsets */
    double sigma_offset;
    /** Vector containing estimated fraction for each model */
    std::vector<double> alpha_k;
    /** Vector containing estimated fraction for mirror of each model */
    std::vector<double> mirror_fraction;
    /** Number of steps to sample in-plane rotation in 90 degrees */
    size_t max_nr_psi;
    /** Vary psi and translational sampling with resolution */
    bool do_variable_psi, do_variable_trans;
    /** Vector for images to hold references (new & old) */
    std::vector< Image<double> > Ictf;
    /** Limit translational searches */
    bool limit_trans;
    /** Number of limited translations */
    size_t nr_trans;
    /** Number for which limited translation is zero */
    size_t zero_trans;
    /** Limited search range for origin offsets */
    int search_shift;
    /** Number of subdirectories to keep for unique offsets filenames */
    int offsets_keepdir;

    /** CTFDat file for all images */
    //CTFDat ctfdat;
    /** Flag whether to include CTFs in the image formation model */
    bool do_ctf_correction;
    /** Pixel size in Angstroms */
    double sampling;
    /** Vector with number of images per defocuss group */
    std::vector<int> count_defocus;
    /** Flag whether the phases of the experimental images are flipped already */
    bool phase_flipped;
    /** Matrix with resolution shell at each Fourier pixel */
    MultidimArray<size_t> Mresol_int;
    /** Vectors with sigma2 (for each defocus group) */
    std::vector<MultidimArray<double> > Vsig, Vctf, Vdec;
    /** Multiplicative factor for SSNR */
    double reduce_snr;
    /** number of defocus groups */
    size_t nr_focus;
    /** Overall low and high resolution cutoffs for fourier-mode (in Fourier pixels) */
    double lowres_limit, highres_limit, ini_highres_limit;
    /** Do not multiply signal with CTF in the first iteration */
    bool first_iter_noctf;
    /** Divide by CTF (until first zero) instead of wiener filter */
    bool do_divide_ctf;
    /** Include all frequencies in the refinement */
    bool do_include_allfreqs;
    /** Fix high-resolution limit */
    double fix_high;
    /** Pointers to the 2D matrices (in FourierTransformHalf format) */
    std::vector<int> pointer_2d, pointer_i, pointer_j;
    size_t nr_points_prob, nr_points_2d, dnr_points_2d;
    /** Current highest resolution shell */
    size_t current_highres_limit;

    /// IN DEVELOPMENT

    /// USe t-distribution instead of normal one
    /** Use t-student distribution instead of normal one */
    bool do_student;
    /** Degrees of freedom for the t-student distribution */
    double df, df2;
    /** Perform sigma-trick for faster convergence (Mclachlan&Peel, p. 228)*/
    bool do_student_sigma_trick;

    /// Statistical analysis of the noise distributions
    /** Perform Kolmogorov-Smirnov test on noise distribution */
    bool do_kstest;
    /** Iteration at which to write out histograms */
    int iter_write_histograms;
    /** Average histogram */
    Histogram1D sumhist;
    std::vector<Histogram1D > resolhist;

    std::vector<double> refs_avgscale;

    int rank, size, nr_vols; //mpi rank, size and 3d number of vols

    /** Taken from expectation and/or maximization */
    double LL, sumcorr, wsum_sigma_offset, sumw_allrefs;
    std::vector<MultidimArray<double> > wsum_Mref, wsum_ctfMref;
    std::vector< std::vector<double> > Mwsum_sigma2;
    std::vector<double> sumw, sumw2, sumwsc, sumwsc2, sumw_mirror, sumw_defocus;
    std::vector<double> Fref, Fwsum_imgs, Fwsum_ctfimgs;


    /** Constructor
     * nr_vols: Used for 3D case, if 0, then doML3D = false
     * rank and size are for MPI use
     */
    ProgMLF2D(int nr_vols = 0, int rank = 0, int size = 1);

    /// Params definition
    void defineParams();

    /// Read arguments from command line
    void readParams();


public:
    /// Show
    void show(bool ML3D = false);

    /// Setup lots of stuff
    virtual void produceSideInfo();

    /// Read reference images in memory & set offset vectors
    /// (This produceSideInfo is Selfile-dependent!)
    virtual void produceSideInfo2();

    /// Some redefinitions of the basic and additional params
    virtual void defineBasicParams(XmippProgram * prog);
    virtual void defineAdditionalParams(XmippProgram * prog, const char * sectionLine);

    /// Calculate initial sigma2 from average power spectrum of the
    /// experimental images
    void estimateInitialNoiseSpectra();

    /// Calculate Wiener filter for defocus series as defined by Frank
    /// (2nd ed. formula 2.32b on p.60)
    void updateWienerFilters(const MultidimArray<double> &spectral_signal,
                             std::vector<double> &sumw_defocus, int iter);

    /// Vary in-plane and translational sampling rates with resolution
    void setCurrentSamplingRates(double current_probres_limit);

    /// Generate initial references from random subset averages
    void generateInitialReferences();


    /// Calculate probability density distribution for in-plane transformations
    void calculatePdfInplane();

    // Append a FourierTransform (in half format!) to a vector
    void appendFTtoVector(const MultidimArray<std::complex<double> > &Fin,
                          std::vector<double> &out);

    // get a FourierTransform (in half format!) from a vector
    void getFTfromVector(const std::vector<double> &in,
                         const int start_point,
                         MultidimArray<std::complex<double> > &out,
                         bool only_real = false);

    /// Fill vector of matrices with all rotations of reference
    void rotateReference(std::vector<double> &out);

    /// Apply reverse rotations to all matrices in vector and fill new matrix with their sum
    void reverseRotateReference(const std::vector<double> &in,
                                std::vector<MultidimArray<double > > &out);

    /// Calculate which references have projection directions close to
    /// phi and theta
    void preselectDirections(float &phi, float &theta,
                             std::vector<double> &pdf_directions);

    /// Pre-calculate which model and phi have significant probabilities
    /// without taking translations into account!
    void preselect_significant_model_phi(MultidimArray<double> &Mimg, std::vector<double> &offsets,
                                         std::vector <std::vector< MultidimArray<double > > > &Mref,
                                         MultidimArray<int> &Msignificant,
                                         std::vector<double > &pdf_directions);

    // Calculate the FT of a translated matrix using a phase shift in
    // Fourier space
    void fourierTranslate2D(const std::vector<double> &in,
                            MultidimArray<double> &trans,
                            std::vector<double> &out,
                            int point_start = 0);

    // If not determined yet: search optimal offsets using maxCC
    // Then for all optimal translations, calculate all translated FTs
    // for each of the flipped variants
    void calculateFourierOffsets(const MultidimArray<double> &Mimg,
                                 const std::vector<double > &offsets,
                                 std::vector<double>  &out,
                                 MultidimArray<int> &Moffsets,
                                 MultidimArray<int> &Moffsets_mirror);

    /// Perform expectation step for a single image
    void processOneImage(const MultidimArray<double> &Mimg,
                         size_t focus, bool apply_ctf,
                         double &fracweight,  double &maxweight2, double &sum_refw2,
                         double &opt_scale, size_t &opt_refno, double &opt_psi,
                         size_t &opt_ipsi, size_t &opt_iflip,
                         MultidimArray<double> &opt_offsets,
                         std::vector<double> &opt_offsets_ref,
                         std::vector<double > &pdf_directions,
                         bool do_kstest, bool write_histograms,
                         FileName fn_img, double &KSprob);

    /// Perform Kolmogorov-Smirnov test
    double performKSTest(MultidimArray<double> &Mimg,  const int focus, bool apply_ctf,
                         FileName &fn_img, bool write_histogram,
                         std::vector <std::vector< MultidimArray<std::complex<double> > > > &Fref,
                         double &opt_scale, int &opt_refno, int &opt_ipsi, int &opt_iflip,
                         MultidimArray<double> &opt_offsets);

    /// One iteration of the process
    void iteration();
    /// Integrate over all experimental images
    void expectation();
    /// Update all model parameters (maximization step)
    void maximization();
    /// Redefine endIteration to update Wiener filters
    void endIteration();

    /// Write out reference images, selfile and logfile
    virtual void writeOutputFiles(const ModelML2D &model, OutputType outputType);

    /// Write noise estimation per resolution
    void writeNoiseFile(const FileName &fn_base, int ifocus);

    virtual void addPartialDocfileData(const MultidimArray<double> &data, size_t first, size_t last);

};
//@}

#endif
