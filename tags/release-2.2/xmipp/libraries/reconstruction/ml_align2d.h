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

#ifndef _MLALIGN3D_H
#define _MLALIGN3D_H

// This program will run 3x faster with FFTW routines than with Xmipp's own FFT routines
// However, in some cases the FFTW may not be available (although it comes with Xmipp)
// In that case, comment out the next line to use Xmipp's routines
#define HAVE_FFTW
#ifdef HAVE_FFTW
#include <data/fftw.h>
#endif
#include <data/fft.h>
#include <data/args.h>
#include <data/funcs.h>
#include <data/selfile.h>
#include <data/docfile.h>
#include <data/image.h>
#include <data/geometry.h>
#include <data/filters.h>
#include <data/mask.h>

#include "ctf.h"

#include <vector>

#define FOR_ALL_MODELS() for (int refno=0;refno<n_ref; refno++)
#define FOR_ALL_ROTATIONS() for (int ipsi=0; ipsi<nr_psi; ipsi++ )
#define FOR_ALL_FLIPS() for (int iflip=0; iflip<nr_flip; iflip++)
#define FOR_ALL_LIMITED_TRANSLATIONS() for (int itrans=0; itrans<nr_trans; itrans++)
#define SIGNIFICANT_WEIGHT_LOW 1e-8
#define SMALLVALUE 1e-4
#define SMALLANGLE 1.75

/**@defgroup MLalign2D ml_align2d (Maximum likelihood in 2D)
   @ingroup ReconsLibraryPrograms */
//@{
/** MLalign2D parameters. */
class Prog_MLalign2D_prm
{
public:
    /** Filenames reference selfile/image, fraction docfile & output rootname */
    FileName fn_sel, fn_ref, fn_root, fn_frac, fn_sig, fn_doc, fn_oext, fn_scratch, fn_control;
    /** If true: use maximum cross-correlation instead of maximum likelihood target */
    bool maxCC_rather_than_ML;
    /** Command line */
    std::string cline;
    /** Sigma value for expected pixel noise */
    double sigma_noise;
    /** sigma-value for origin offsets */
    double sigma_offset;
    /** Vector containing estimated fraction for each model */
    std::vector<double> alpha_k;
    /** Vector containing estimated fraction for mirror of each model */
    std::vector<double> mirror_fraction;
    /** Flag for checking mirror images of all references */
    bool do_mirror;
    /** Flag whether to fix estimates for model fractions */
    bool fix_fractions;
    /** Flag whether to fix estimate for sigma of origin offset */
    bool fix_sigma_offset;
    /** Flag whether to fix estimate for sigma of noise */
    bool fix_sigma_noise;
    /** Starting iteration */
    int istart;
    /** Number of iterations to be performed */
    int Niter;
    /** dimension of the images */
    int dim, dim2, hdim;
    /** Number of steps to sample in-plane rotation in 90 degrees */
    int nr_psi;
    /** Number of operations in "flip-array" (depending on do_mirror) */
    int nr_flip;
    /** Sampling rate for in-plane rotation */
    float psi_step;
    /** Total degrees in FOR_ALL_ROTATIONS */
    double psi_max;
    /** Total number of no-mirror rotations in FOR_ALL_FLIPS */
    int nr_nomirror_flips;
    /** Number of reference images */
    int n_ref;
    /** Total number of experimental images */
    int nr_exp_images;
    /** Sum of squared amplitudes of the references */
    std::vector<double> A2;
    /** Verbose level:
        1: gives progress bar (=default)
        0: gives no output to screen at all */
    int verb;
    /** Stopping criterium */
    double eps;
    /** Write out document file with orientations & models that give
        max. probability for each image */
    bool write_docfile;
    /** Write out selection files according to model assignments */
    bool write_selfiles;
    /** Write out images and selfile after each iteration */
    bool write_intermediate;
    /** SelFile images (working, test and reference set) */
    SelFile SF, SFr;
    /** vector for flipping (i.e. 90/180-degree rotations) matrices */
    std::vector<Matrix2D<double> > F;
    /** Vector for images to hold references (new & old) */
    std::vector < ImageXmippT<double> > Iref, Iold;
    /** Matrices for calculating PDF of (in-plane) translations */
    Matrix2D<double> P_phi, Mr2;
    /** Fast mode */
    bool fast_mode;
    /** Fast mode */
    double C_fast;
    /** Very crude origin pixel artifact correction */
    bool do_esthetics;
    /** Maximum shift to be trusted */
    double max_shift;
    /** Limit translational searches */
    bool limit_trans;
    /** Number of limited translations */
    int nr_trans;
    /** Number for which limited translation is zero */
    int zero_trans;
    /** Offsets for limited translations */
    std::vector<Matrix1D<double> > Vtrans;
    /** Start all optimal offsets from zero values */
    bool zero_offsets;
    /** Limited search range for origin offsets */
    double search_shift;
    /** Limit orientational searches */
    bool limit_rot;
    /** Limited search range for projection directions */
    double search_rot;
    /** Save memory options */
    bool save_mem1, save_mem2, save_mem3;
    /** Vectors to store old phi and theta for all images */
    std::vector<float> imgs_oldphi, imgs_oldtheta;
    /** Number of subdirectories to keep for unique offsets filenames */
    int offsets_keepdir;
    /** Flag for using ML3D */
    bool do_ML3D;
    /** Flag for generation of initial models from random subsets */
    bool do_generate_refs;
    /** Flag whether to write offsets to disc */
    bool do_write_offsets;
    /** Vector to store optimal origin offsets (if not written to disc) */
    std::vector<std::vector<double> > imgs_offsets;

    /// Students t-distribution
    /** Use t-student distribution instead of normal one */
    bool do_student;
    /** Degrees of freedom for the t-student distribution */
    double df, df2;
    /** Do sigma-division trick in student-t*/
    bool do_student_sigma_trick;

    /// Re-normalize internally
    /** Flag to refine normalization of each experimental image */
    bool do_norm;
    /** Grey-scale correction values */
    std::vector<double> imgs_scale, imgs_bgmean, refs_avgscale;
    /** Overall average scale (to be forced to one)*/
    double average_scale;

    /// Statistical analysis of the noise distributions
    /** Perform Kolmogorov-Smirnov test on noise distribution */
    bool do_kstest;
    /** Average histogram */
    Matrix1D<double> sumhist;

    /** IN DEVELOPMENT */

    /// Deterministic annealing
    double anneal, anneal_step;

    /// Write out a histogram of the noise distribution for each particle 
    /** Iteration at which to write out histograms */
    int iter_write_histograms;

    /// Each image has its own sigma for the noise
    /** Flag to de per-image noise model */
    bool do_per_image_noise;
    /** Per-image noise models */
    std::vector<double> imgs_noise_sigma;

#ifdef HAVE_FFTW
    /// Use fftw
    /** Use the FFTW libraries instead of Xmipps FFT routines */
    /** FFTW objects */
    xmippFftw forwfftw, backfftw;
#endif

    /** debug flag */
    int debug;

public:
    /// Read arguments from command line
    void read(int argc, char **argv, bool ML3D = false);

    /// Show
    void show(bool ML3D = false);

    /// Usage for ML mode
    void usage();

    /// Extended Usage
    void extended_usage(bool ML3D = false);

    /// Setup lots of stuff
    void produce_Side_info();

    /// Generate initial references from random subset averages
    void generate_initial_references();

    /** Read reference images in memory & set offset vectors
        (This produce_side_info is Selfile-dependent!) */
    void produce_Side_info2(int nr_vols = 1);

    /** Read and write optimal translations to disc
        (not to store them all in memory) */
    void write_offsets(FileName fn, std::vector<double> &data);
    bool read_offsets(FileName fn, std::vector<double> &data);

    /// Calculate probability density distribution for in-plane transformations
    void calculate_pdf_phi();

    /// Fill vector of matrices with all rotations of reference
    void rotate_reference(std::vector< ImageXmippT<double> > &Iref,
                          bool fill_real_space,
                          bool fill_fourier_space,
                          std::vector <std::vector< Matrix2D<double> > > &Mref,
                          std::vector <std::vector< Matrix2D<std::complex<double> > > > &Fref);

    /// Apply reverse rotations to all matrices in vector and fill new matrix with their sum
    void reverse_rotate_reference(std::vector <std::vector< Matrix2D<std::complex<double> > > > &Fnew,
                                  std::vector <std::vector< Matrix2D<double> > > &Mnew, bool real_space,
                                  std::vector<Matrix2D<double> > &Mref);

    /** Calculate which references have projection directions close to
        phi and theta */
    void preselect_directions(float &phi, float &theta,
                              std::vector<double> &pdf_directions);

    /** Pre-calculate which model and phi have significant probabilities
       without taking translations into account! */
    void preselect_significant_model_phi(Matrix2D<double> &Mimg, std::vector<double> &offsets,
                                         std::vector <std::vector< Matrix2D<double > > > &Mref,
                                         Matrix2D<int> &Msignificant,
                                         std::vector<double > &pdf_directions);

    /** Calculate translated matrices for all limited translations
        for each of the flipped variants */
    void calculate_realspace_offsets(Matrix2D<double> &Mimg, std::vector<double > &offsets,
                                     std::vector<double > &pdf_directions,
                                     std::vector<std::vector<Matrix2D<double> > > &Mimg_trans,
                                     Matrix2D<int> &Moffsets, Matrix2D<int> &Moffsets_mirror);

    /** ML-integration over limited translations,
        and with -fast way of selection significant rotations */
    void ML_integrate_locally(Matrix2D<double> &Mimg,
                              std::vector <std::vector< Matrix2D<double> > > &Mref,
                              std::vector <std::vector< Matrix2D<double> > > &Mwsum_imgs,
                              double &wsum_sigma_noise, double &wsum_sigma_offset,
                              std::vector<double> &sumw, std::vector<double> &sumw2, 
                              std::vector<double> &sumwsc, std::vector<double> &sumwsc2, std::vector<double> &sumw_mirror,
                              double &LL, double &fracweight, double &maxweight2, double &opt_scale,
                              double &per_image_mean, double &per_image_sigma,
                              int &opt_refno, double &opt_psi,
                              Matrix1D<double> &opt_offsets,
                              std::vector<double> &opt_offsets_ref,
                              std::vector<double> &pdf_directions);

    /// ML-integration over all (or -fast) translations
    void ML_integrate_complete(Matrix2D<double> &Mimg,
                               std::vector <std::vector< Matrix2D<std::complex<double> > > > &Fref,
                               Matrix2D<int> &Msignificant,
                               std::vector <std::vector< Matrix2D<std::complex<double> > > > &Fwsum_imgs,
                               double &wsum_sigma_noise, double &wsum_sigma_offset,
                               std::vector<double> &sumw, std::vector<double> &sumw2,
                               std::vector<double> &sumwsc, std::vector<double> &sumwsc2, std::vector<double> &sumw_mirror,
                               double &LL, double &fracweight, double &maxweight2, double &opt_scale, 
                               double &per_image_mean, double &per_image_sigma,
                               int &opt_refno, double &opt_psi, int &iopt_psi, int &iopt_flip, 
                               Matrix1D<double> &opt_offsets, std::vector<double> &opt_offsets_ref,
                               std::vector<double > &pdf_directions);

    /// Perform Kolmogorov-Smirnov test
    double performKSTest(Matrix2D<double> &Mimg, FileName &fn_img, bool write_histogram,
                         std::vector <std::vector< Matrix2D<std::complex<double> > > > &Fref,
                         double &per_image_mean, double &per_image_sigma, double &opt_scale,
                         int &opt_refno, int &iopt_psi, int &iopt_flip,
                         Matrix1D<double> &opt_offsets);

    /// Calculate maxCC averages for new model and new model parameters
    void maxCC_search_complete(Matrix2D<double> &Mimg,
                               std::vector <std::vector< Matrix2D<std::complex<double> > > > &Fref,
                               std::vector <std::vector< Matrix2D<double> > > &Mref,
                               double &max_shift,
                               std::vector <std::vector< Matrix2D<double> > > &Msum_imgs,
                               std::vector<double> &sumw, std::vector<double> &sumw_mirror,
                               double &minSQ, int &opt_refno, double &opt_psi,
                               Matrix1D<double> &opt_offsets,
                               std::vector<double> &pdf_directions);

    /// Integrate over all experimental images
    void ML_sum_over_all_images(SelFile &SF, std::vector< ImageXmippT<double> > &Iref, int iter,
                                double &LL, double &sumcorr, DocFile &DFo,
                                std::vector<Matrix2D<double> > &wsum_Mref,
                                double &wsum_sigma_noise, double &wsum_sigma_offset,
                                std::vector<double> &sumw, std::vector<double> &sum2, 
                                std::vector<double> &sumwsc, std::vector<double> &sumwsc2, std::vector<double> &sumw_mirror);

    /// Update all model parameters
    void update_parameters(std::vector<Matrix2D<double> > &wsum_Mref,
                           double &wsum_sigma_noise, double &wsum_sigma_offset, 
			   std::vector<double> &sumw, std::vector<double> &sumw2, 
                           std::vector<double> &sumwsc, std::vector<double> &sumwsc2, std::vector<double> &sumw_mirror,
                           double &sumcorr, double &sumw_allrefs, int refs_per_class=1);

    /// check convergence
    bool check_convergence(std::vector<double> &conv);

    /// Output some parameters to screen
    void output_to_screen(int &iter, double &sumcorr, double &LL);

    /// Write out reference images, selfile and logfile
    void write_output_files(const int iter, DocFile &DFo,
                            double &sumw_allrefs, double &LL, double &avecorr,
                            std::vector<double> &conv);

};
//@}
#endif
