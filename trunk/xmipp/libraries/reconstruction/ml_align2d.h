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

#ifndef _MLALIGN2D_H
#define _MLALIGN2D_H

#include <data/fftw.h>
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
#include <pthread.h>
#include <vector>

#define SIGNIFICANT_WEIGHT_LOW 1e-8
#define SMALLANGLE 1.75
#define DATALINELENGTH 12
class Prog_MLalign2D_prm;

// Thread declaration
void * threadExpectationSingleImage( void * data );

// This structure is needed to pass parameters to threadExpectationSingleImage
typedef struct{
    int thread_id;
    int thread_num;
    Prog_MLalign2D_prm *prm;
    SelFile *SF;
    int *iter;
    double *wsum_sigma_noise;
    double *wsum_sigma_offset;
    double *sumfracweight;
    double *LL;
    std::vector<Matrix2D<std::complex<double> > > *wsumimgs;
    std::vector<Matrix2D<std::complex<double> > > *fref;
    std::vector<Matrix2D<double> > *mref;
    std::vector<Matrix1D<double> > *docfiledata;
    std::vector<double> *sumw;
    std::vector<double> *sumw2;
    std::vector<double> *sumwsc;
    std::vector<double> *sumwsc2;
    std::vector<double> *sumw_mirror;
} structThreadExpectationSingleImage ;

/**@defgroup MLalign2D ml_align2d (Maximum likelihood in 2D)
   @ingroup ReconsLibraryPrograms */
//@{
/** MLalign2D parameters. */
class Prog_MLalign2D_prm
{
public:
    /** Filenames reference selfile/image, fraction docfile & output rootname */
    FileName fn_sel, fn_ref, fn_root, fn_frac, fn_sig, fn_doc, fn_oext, fn_scratch, fn_control;
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
    int oridim, dim, dim2, hdim;
    double ddim2;
    /** Maximum resolution (dig.freq.) */
    double max_resol, scale_factor;
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
    /** Offsets for limited translations */
    std::vector<Matrix1D<double> > Vtrans;
    /** Start all optimal offsets from zero values */
    bool zero_offsets;
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
    /** For initial guess of mindiff */
    double trymindiff_factor;

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
    std::vector<double> imgs_scale, imgs_bgmean, imgs_trymindiff, refs_avgscale;
    /** Optimal refno from previous iteration */
    std::vector<int> imgs_optrefno;
    /** Overall average scale (to be forced to one)*/
    double average_scale;

    /** Threads */
    int threads;

    /** FFTW objects */
    XmippFftw transformer;

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
    void extendedUsage(bool ML3D = false);

    /// Setup lots of stuff
    void produceSideInfo();

    /// Rescale image from dim to oridim and vice versa
    void reScaleImage(Matrix2D<double> &Min, bool down_scale);

    /// Generate initial references from random subset averages
    void generateInitialReferences();

    /** Read reference images in memory & set offset vectors
        (This produce_side_info is Selfile-dependent!) */
    void produceSideInfo2(int nr_vols = 1);

    /** Read and write optimal translations to disc
        (not to store them all in memory) */
    void writeOffsets(FileName fn, std::vector<double> &data);
    bool readOffsets(FileName fn, std::vector<double> &data);

    /// Calculate probability density distribution for in-plane transformations
    void calculatePdfInplane();

    /// Fill vector of matrices with all rotations of reference
    void rotateReference(std::vector< ImageXmippT<double> > &Iref,
                         bool fill_real_space, 
                         std::vector<Matrix2D<double> > &mref, 
                         std::vector<Matrix2D<std::complex<double > > > &fref);

    /// Apply reverse rotations to all matrices in vector and fill new matrix with their sum
    void reverseRotateReference(std::vector<Matrix2D<std::complex<double > > > &fnew, 
                                std::vector<Matrix2D<double> > &Mref);

    /** Calculate which references have projection directions close to
        phi and theta */
    void preselectLimitedDirections(float &phi, float &theta,
                                    std::vector<double> &pdf_directions);

    /** Pre-calculate which model and phi have significant probabilities
       without taking translations into account! */
    void preselectFastSignificant(Matrix2D<double> &Mimg, std::vector<double> &offsets,
                                  std::vector<Matrix2D<double> > &mref,
                                  Matrix2D<int> &Msignificant,
                                  std::vector<double > &pdf_directions);

    /// ML-integration over all (or -fast) translations
    void expectationSingleImage(Matrix2D<double> &Mimg,
                                std::vector<Matrix2D<std::complex<double> > > &fref,
                                std::vector<Matrix2D<std::complex<double> > > &wsumimgs,
                                Matrix2D<int> &Msignificant,
                                double &wsum_sigma_noise, double &wsum_sigma_offset,
                                std::vector<double> &sumw, std::vector<double> &sumw2,
                                std::vector<double> &sumwsc, std::vector<double> &sumwsc2, 
                                std::vector<double> &sumw_mirror,
                                double &LL, double &dLL, double &fracweight, double &sumfracweight, 
                                double &maxweight2, double &opt_scale, double &bgmean, double &trymindiff,
                                int &opt_refno, double &opt_psi, int &iopt_psi, int &iopt_flip, 
                                Matrix1D<double> &opt_offsets, std::vector<double> &opt_offsets_ref,
                                std::vector<double > &pdf_directions);

    /// Integrate over all experimental images
    void expectation(SelFile &SF, std::vector< ImageXmippT<double> > &Iref, int iter,
                     double &LL, double &sumfracweight, DocFile &DFo,
                     std::vector<Matrix2D<double> > &wsum_Mref,
                     double &wsum_sigma_noise, double &wsum_sigma_offset,
                     std::vector<double> &sumw, std::vector<double> &sum2, 
                     std::vector<double> &sumwsc, std::vector<double> &sumwsc2, std::vector<double> &sumw_mirror);

    /// Update all model parameters
    void maximization(std::vector<Matrix2D<double> > &wsum_Mref,
                      double &wsum_sigma_noise, double &wsum_sigma_offset, 
                      std::vector<double> &sumw, std::vector<double> &sumw2, 
                      std::vector<double> &sumwsc, std::vector<double> &sumwsc2, 
                      std::vector<double> &sumw_mirror,
                      double &sumfracweight, double &sumw_allrefs, int refs_per_class=1);

    /// check convergence
    bool checkConvergence(std::vector<double> &conv);

    /// Write out reference images, selfile and logfile
    void writeOutputFiles(const int iter, DocFile &DFo,
                          double &sumw_allrefs, double &LL, double &avefracweight,
                          std::vector<double> &conv);

};
//@}
#endif
