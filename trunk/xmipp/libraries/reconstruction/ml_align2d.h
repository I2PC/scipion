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
#define FOR_ALL_DEFOCUS_GROUPS() for (int ifocus=0; ifocus<nr_focus; ifocus++)
#define SIGNIFICANT_WEIGHT_LOW 1e-8
#define SMALLVALUE 1e-4
#define SMALLANGLE 1.75

/**@name MLalign2D */
//@{
/** MLalign2D parameters. */
class Prog_MLalign2D_prm
{
public:
    /** Filenames reference selfile/image, fraction docfile & output rootname */
    FileName fn_sel, fn_ref, fn_root, fn_frac, fn_sig, fn_doc, fn_ctf, fn_oext;
    /** If true: use maximum cross-correlation instead of maximum likelihood target */
    bool maxCC_rather_than_ML;
    /** Command line */
    string cline;
    /** Sigma value for expected pixel noise */
    double sigma_noise;
    /** sigma-value for origin offsets */
    double sigma_offset;
    /** Vector containing estimated fraction for each model */
    vector<double> alpha_k;
    /** Vector containing estimated fraction for mirror of each model */
    vector<double> mirror_fraction;
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
    vector<double> A2;
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
    vector<Matrix2D<double> > F;
    /** Vector for images to hold references (new & old) */
    vector <ImageXmipp> Iref, Iold, Ictf;
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
    vector<Matrix1D<double> > Vtrans;
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
    /** Vectors to store old phi, theta, xoff and yoff for all images */
    vector<float> imgs_oldphi, imgs_oldtheta, imgs_oldxoff, imgs_oldyoff;
    /** Number of subdirectories to keep for unique offsets filenames */
    int offsets_keepdir;
    /** Flag for using ML3D */
    bool do_ML3D;
    /** Flag for generation of initial models from random subsets */
    bool do_generate_refs;
    /** Flag whether to write offsets to disc */
    bool do_write_offsets;
    /** Vector to store optimal origin offsets (if not written to disc) */
    vector<vector<double> > imgs_offsets;

    /// FOR FOURIER_MODE
    /** If true: use fourier instead of real space maximum likelihood target */
    bool fourier_mode;
    /** Pixel size in Angstroms */
    double sampling;
    /** Vector with number of images per defocuss group (fourier_mode) */
    vector<int> count_defocus;
    /** Flag whether the phases of the experimental images are flipped already */
    bool phase_flipped;
    /** Matrix with resolution value at each Fourier pixel */
    Matrix2D<int> Mresol;
    /** Vectors with sigma2 (for each defocus group) */
    vector<Matrix1D<double> > Vsig, Vctf, Vdec;
    /** pointers for the different ctf-matrices */
    vector<int> pointer_ctf, pointer_i, pointer_j, pointer_sigctf;
    vector<int> pointer_old, pointer_i_old, pointer_j_old;
    /** Number of elements in pointers for the different ctf-matrices */
    int nr_pointer_ctf, nr_pointer_sigctf, nr_pointer_old;
    /** Multiplicative factor for SSNR */
    double reduce_snr;
    /** number of defocus groups */
    int nr_focus;
    /** low and high resolution cutoffs for fourier-mode (in Fourier pixels) */
    int lowres_limit, highres_limit, ini_highres_limit;
    /** Number of Fourier rings to include above SNR<1 */
    int increase_highres_limit;
    /** Do not multiply signal with CTF in the first iteration */
    bool first_iter_noctf;

    /// IN DEVELOPMENT
    /// Deterministic annealing
    double anneal, anneal_step;
    /** Flag to output the log-likelihood target value for MLF mode */
    bool do_output_MLF_LL;
    /** The log-likelihood target function value of the previous iteration */
    double LL_prev_iter;

public:
    /// Read arguments from command line
    void read(int argc, char **argv, bool ML3D = false);

    /// Show
    void show(bool ML3D = false);

    /// Usage for ML mode
    void usage();

    /// Usage for MLF mode
    void MLF_usage();

    /// Extended Usage
    void extended_usage(bool ML3D = false);

    /// Setup lots of stuff
    void produce_Side_info();

    /// Calculate initial sigma2 from average power spectrum of the
    /// experimental images
    void estimate_initial_sigma2();

    /// Calculate Wiener filter for defocus series as defined by Frank
    /// (2nd ed. formula 2.32b on p.60)
    void calculate_wiener_defocus_series(Matrix1D<double> &spectral_signal, int iter);

    /// Generate initial references from random subset averages
    void generate_initial_references();

    /// Read reference images in memory & set offset vectors
    /// (This produce_side_info is Selfile-dependent!)
    void produce_Side_info2(int nr_vols = 1);

    /// Read and write optimal translations to disc
    /// (not to store them all in memory)
    void write_offsets(FileName fn, vector<double> &data);
    bool read_offsets(FileName fn, vector<double> &data);

    /// Calculate probability density distribution for in-plane transformations
    void calculate_pdf_phi();

    /// Fill vector of matrices with all rotations of reference
    void rotate_reference(vector<ImageXmipp> &Iref,
                          bool fill_real_space,
                          bool fill_fourier_space,
                          vector <vector< Matrix2D<double> > > &Mref,
                          vector <vector< Matrix2D<complex<double> > > > &Fref);

    /// Apply reverse rotations to all matrices in vector and fill new matrix with their sum
    void reverse_rotate_reference(vector <vector< Matrix2D<complex<double> > > > &Fnew,
                                  vector <vector< Matrix2D<double> > > &Mnew, bool real_space,
                                  vector<Matrix2D<double> > &Mref);

    /// Calculate which references have projection directions close to
    /// phi and theta
    void preselect_directions(float &phi, float &theta,
                              vector<double> &pdf_directions);

    /// Pre-calculate which model and phi have significant probabilities
    /// without taking translations into account!
    void preselect_significant_model_phi(Matrix2D<double> &Mimg, vector<double> &offsets,
                                         vector <vector< Matrix2D<double > > > &Mref,
                                         Matrix2D<int> &Msignificant,
                                         vector<double > &pdf_directions);

    // Calculate the FT of a translated matrix using a phase shift in
    // Fourier space
    void Fourier_translate2D(const Matrix2D<complex<double> > &Fimg,
                             int focus, Matrix1D<double> &trans,
                             Matrix2D<complex<double> > &Fimg_shift);

    // If not determined yet: search optimal offsets using maxCC
    // Then for all optimal translations, calculate all translated FTs
    // for each of the flipped variants
    void calculate_fourier_offsets(Matrix2D<double> &Mimg, int focus,
                                   vector <vector< Matrix2D<complex<double> > > > &Fref,
                                   Matrix2D<double> &ctf, vector<double> &offsets,
                                   vector<vector<Matrix2D<complex<double> > > > &Fimg_trans,
                                   Matrix2D<int> &Moffsets, Matrix2D<int> &Moffsets_mirror);

    // Calculate translated matrices for all limited translations
    // for each of the flipped variants
    void calculate_realspace_offsets(Matrix2D<double> &Mimg, vector<double > &offsets,
                                     vector<double > &pdf_directions,
                                     vector<vector<Matrix2D<double> > > &Mimg_trans,
                                     Matrix2D<int> &Moffsets, Matrix2D<int> &Moffsets_mirror);

    /// fast MLF integration...
    void MLF_integrate_locally(Matrix2D<double> &Mimg, int focus, bool apply_ctf,
                               vector<vector<Matrix2D<complex<double> > > > &Fref,
                               vector <vector< Matrix2D<complex<double> > > > &Fwsum_imgs,
                               vector <vector< Matrix2D<complex<double> > > > &Fwsum_ctfimgs,
                               Matrix2D<double> &Mwsum_sigma2,
                               double &wsum_sigma_offset, vector<double> &sumw,
                               vector<double> &sumw_mirror,
                               double &LL, double &LL_old, double &fracweight, 
			       int &opt_refno, double &opt_psi, 
			       Matrix1D<double> &opt_offsets,
                               vector<double> &opt_offsets_ref,
                               vector<double > &pdf_directions);

    // ML-integration over limited translations,
    // and with -fast way of selection significant rotations
    void ML_integrate_locally(Matrix2D<double> &Mimg,
                              vector <vector< Matrix2D<double> > > &Mref,
                              vector <vector< Matrix2D<double> > > &Mwsum_imgs,
                              double &wsum_sigma_noise, double &wsum_sigma_offset,
                              vector<double> &sumw, vector<double> &sumw_mirror,
                              double &LL, double &fracweight,
                              int &opt_refno, double &opt_psi,
                              Matrix1D<double> &opt_offsets,
                              vector<double> &opt_offsets_ref,
                              vector<double> &pdf_directions);

    /// ML-integration over all (or -fast) translations
    void ML_integrate_complete(Matrix2D<double> &Mimg,
                               vector <vector< Matrix2D<complex<double> > > > &Fref,
                               Matrix2D<int> &Msignificant,
                               vector <vector< Matrix2D<complex<double> > > > &Fwsum_imgs,
                               double &wsum_sigma_noise, double &wsum_sigma_offset,
                               vector<double> &sumw, vector<double> &sumw_mirror,
                               double &LL, double &fracweight,
                               int &opt_refno, double &opt_psi,
                               Matrix1D<double> &opt_offsets, vector<double> &opt_offsets_ref,
                               vector<double > &pdf_directions);

    /// Calculate maxCC averages for new model and new model parameters
    void maxCC_search_complete(Matrix2D<double> &Mimg,
                               vector <vector< Matrix2D<complex<double> > > > &Fref,
                               vector <vector< Matrix2D<double> > > &Mref,
                               double &max_shift,
                               vector <vector< Matrix2D<double> > > &Msum_imgs,
                               vector<double> &sumw, vector<double> &sumw_mirror,
                               double &minSQ, int &opt_refno, double &opt_psi,
                               Matrix1D<double> &opt_offsets,
                               vector<double> &pdf_directions);

    /// Integrate over all experimental images
    void ML_sum_over_all_images(SelFile &SF, vector<ImageXmipp> &Iref, int iter,
                                double &LL, double &LL_old, double &sumcorr, DocFile &DFo,
                                vector<Matrix2D<double> > &wsum_Mref,
                                vector<Matrix2D<double> > &wsum_ctfMref,
                                double &wsum_sigma_noise, vector<Matrix2D<double> > &Mwsum_sigma2,
                                double &wsum_sigma_offset,
                                vector<double> &sumw, vector<double> &sumw_mirror);

    /// Update all model parameters
    void update_parameters(vector<Matrix2D<double> > &wsum_Mref,
                           vector<Matrix2D<double> > &wsum_ctfMref,
                           double &wsum_sigma_noise, vector<Matrix2D<double> > &Mwsum_sigma2,
                           double &wsum_sigma_offset, vector<double> &sumw,
                           vector<double> &sumw_mirror,
                           double &sumcorr, double &sumw_allrefs,
                           Matrix1D<double> &spectral_signal);

    /// check convergence
    bool check_convergence(vector<double> &conv);

    /// Output some parameters to screen
    void output_to_screen(int &iter, double &sumcorr, double &LL);

    /// Write out reference images, selfile and logfile
    void write_output_files(const int iter, DocFile &DFo,
                            double &sumw_allrefs, double &LL, double &LL_old, double &avecorr,
                            vector<double> &conv);

};
//@}
