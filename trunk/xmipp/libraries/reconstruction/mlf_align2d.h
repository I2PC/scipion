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

/**@name MLFalign2D */
//@{
/** MLFalign2D parameters. */
class Prog_MLFalign2D_prm
{
public:
    /** Filenames reference selfile/image, fraction docfile & output rootname */
    FileName fn_sel, fn_ref, fn_root, fn_frac, fn_sig, fn_doc, fn_ctfdat, fn_oext, fn_scratch;
    /** Command line */
    string cline;
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
    int nr_psi, max_nr_psi;
    /** Number of operations in "flip-array" (depending on do_mirror) */
    int nr_flip;
    /** Sampling rate for in-plane rotation */
    float psi_step;
    /** Total degrees in FOR_ALL_ROTATIONS */
    double psi_max;
    /** Vary psi and translational sampling with resolution */
    bool do_variable_psi, do_variable_trans;
    /** Total number of no-mirror rotations in FOR_ALL_FLIPS */
    int nr_nomirror_flips;
    /** Number of reference images */
    int n_ref;
    /** Total number of experimental images */
    int nr_exp_images;
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
    vector<Matrix2D<double> > F;
    /** Vector for images to hold references (new & old) */
    vector<ImageXmipp> Iref, Iold, Ictf;
    /** Matrices for calculating PDF of (in-plane) translations */
    Matrix2D<double> P_phi, Mr2;
    /** Fast mode */
    bool fast_mode;
    /** Fast mode */
    double C_fast;
    /** Limit translational searches */
    bool limit_trans;
    /** Number of limited translations */
    int nr_trans;
    /** Number for which limited translation is zero */
    int zero_trans;
    /** Offsets for limited translations */
    vector<Matrix1D<double> > Vtrans;
    /** Limited search range for origin offsets */
    int search_shift;
    /** Limit orientational searches */
    bool limit_rot;
    /** Limited search range for projection directions */
    double search_rot;
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

    /** CTFDat file for all images */
    CTFDat ctfdat;
    /** Flag whether to include CTFs in the image formation model */
    bool do_ctf_correction;
    /** Pixel size in Angstroms */
    double sampling;
    /** Vector with number of images per defocuss group (fourier_mode) */
    vector<int> count_defocus;
    /** Flag whether the phases of the experimental images are flipped already */
    bool phase_flipped;
    /** Matrix with resolution value (in Ang) at each Fourier pixel */
    Matrix2D<double> Mresol;
    Matrix2D<int> Mresol_int;
    /** Vectors with sigma2 (for each defocus group) */
    vector<Matrix1D<double> > Vsig, Vctf, Vdec;
    /** Multiplicative factor for SSNR */
    double reduce_snr;
    /** number of defocus groups */
    int nr_focus;
    /** Overall low and high resolution cutoffs for fourier-mode (in Fourier pixels) */
    double lowres_limit, highres_limit, ini_highres_limit;
    /** Do not multiply signal with CTF in the first iteration */
    bool first_iter_noctf;
    /** Divide by CTF (until first zero) instead of wiener filter */
    bool do_divide_ctf;
    /** Include all frequencies in the refinement */
    bool do_include_allfreqs;
    /** Pointers to the 2D matrices (in FourierTransformHalf format) */
    vector<int> pointer_2d, pointer_i, pointer_j;
    int nr_points_prob, nr_points_2d, dnr_points_2d;
    /** Current highest resolution shell */
    int current_highres_shell;

    /// IN DEVELOPMENT
    /** debug flag */
    int debug;

public:
    /// Read arguments from command line
    void read(int argc, char **argv, bool ML3D = false);

    /// Show
    void show(bool ML3D = false);

    /// Usage
    void usage();

    /// Extended Usage
    void extendedUsage(bool ML3D = false);

    /// Setup lots of stuff
    void produceSideInfo();

    /// Calculate initial sigma2 from average power spectrum of the
    /// experimental images
    void estimateInitialNoiseSpectra();

    /// Calculate Wiener filter for defocus series as defined by Frank
    /// (2nd ed. formula 2.32b on p.60)
    void updateWienerFilters(Matrix1D<double> &spectral_signal, int iter);

    /// Vary in-plane and translational sampling rates with resolution
    void setCurrentSamplingRates(double current_probres_limit);

    /// Generate initial references from random subset averages
    void generateInitialReferences();

    /// Read reference images in memory & set offset vectors
    /// (This produce_side_info is Selfile-dependent!)
    void produceSideInfo2(int nr_vols = 1);

    /// Read and write optimal translations to disc
    /// (not to store them all in memory)
    void writeOffsets(FileName fn, vector<double> &data);
    bool readOffsets(FileName fn, vector<double> &data);

    /// Calculate probability density distribution for in-plane transformations
    void calculateInPlanePDF();

    // Append a FourierTransform (in half format!) to a vector
    void appendFTtoVector(const Matrix2D<complex<double> > &Fin,
			  vector<double> &out);

    // get a FourierTransform (in half format!) from a vector
    void getFTfromVector(const vector<double> &in, 
			 const int start_point,
			 Matrix2D<complex<double> > &out, 
			 bool only_real = false);

    /// Fill vector of matrices with all rotations of reference
    void rotateReference(vector<double> &out);

    /// Apply reverse rotations to all matrices in vector and fill new matrix with their sum
    void reverseRotateReference(const vector<double> &in,
				vector<Matrix2D<double > > &out);

    /// Calculate which references have projection directions close to
    /// phi and theta
    void preselectDirections(float &phi, float &theta,
			     vector<double> &pdf_directions);

    /// Pre-calculate which model and phi have significant probabilities
    /// without taking translations into account!
    void preselect_significant_model_phi(Matrix2D<double> &Mimg, vector<double> &offsets,
                                         vector <vector< Matrix2D<double > > > &Mref,
                                         Matrix2D<int> &Msignificant,
                                         vector<double > &pdf_directions);

    // Calculate the FT of a translated matrix using a phase shift in
    // Fourier space
    void fourierTranslate2D(const vector<double> &in,
			    Matrix1D<double> &trans,
			    vector<double> &out,
			    int point_start = 0);

    // If not determined yet: search optimal offsets using maxCC
    // Then for all optimal translations, calculate all translated FTs
    // for each of the flipped variants
    void calculateFourierOffsets(const Matrix2D<double> &Mimg,
				 const vector<double > &offsets,
				 vector<double>  &out,
				 Matrix2D<int> &Moffsets, 
				 Matrix2D<int> &Moffsets_mirror);

    /// Perform expectation step for a single image
    void processOneImage(const Matrix2D<double> &Mimg, 
			 const int focus, bool apply_ctf,
			 const vector<double> &Fref,
			 vector<double> &Fwsum_imgs,
			 vector<double> &Fwsum_ctfimgs,
			 vector<double> &Mwsum_sigma2,
			 double &wsum_sigma_offset, 
			 vector<double> &sumw,
			 vector<double> &sumw_mirror,
			 double &LL, double &fracweight, 
			 int &opt_refno, double &opt_psi,
			 Matrix1D<double> &opt_offsets, 
			 vector<double> &opt_offsets_ref,
			 vector<double > &pdf_directions);

    /// Integrate over all experimental images
    void sumOverAllImages(SelFile &SF, vector<ImageXmipp> &Iref, int iter,
			  double &LL, double &sumcorr, DocFile &DFo,
			  vector<Matrix2D<double> > &wsum_Mref,
			  vector<Matrix2D<double> > &wsum_ctfMref,
			  vector<vector<double> > &Mwsum_sigma2,
			  double &wsum_sigma_offset, 
			  vector<double> &sumw, 
			  vector<double> &sumw_mirror);

    /// Update all model parameters (maximization step)
    void updateParameters(vector<Matrix2D<double> > &wsum_Mref,
			  vector<Matrix2D<double> > &wsum_ctfMref,
			  vector<vector<double> > &Mwsum_sigma2,
			  double &wsum_sigma_offset,
			  vector<double> &sumw, vector<double> &sumw_mirror,
			  double &sumcorr, double &sumw_allrefs,
			  Matrix1D<double> &spectral_signal);

    /// check convergence
    bool checkConvergence(vector<double> &conv);

    /// Write out reference images, selfile and logfile
    void writeOutputFiles(const int iter, DocFile &DFo,
			  double &sumw_allrefs, double &LL, double &avecorr,
			  vector<double> &conv);


};
//@}
