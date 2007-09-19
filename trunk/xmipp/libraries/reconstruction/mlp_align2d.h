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
#include <data/gridding.h>
#include <data/polar.h>
#include <vector>

#define SIGNIFICANT_WEIGHT_LOW 1e-8

/**@name MLPalign2D */
//@{
/** MLPalign2D parameters. */
class Prog_MLPalign2D_prm
{
public:
    /** Filenames reference selfile/image, fraction docfile & output rootname */
    FileName fn_sel, fn_ref, fn_root, fn_frac, fn_doc;
    /** Command line */
    string cline;
    /** Sigma value for expected pixel noise */
    double sigma_noise;
    /** Sigma value for expected origin offsets */
    double sigma_offset;
    /** Vector containing estimated fraction for each model */
    vector<double> alpha_k;
    /** Vector containing estimated fraction for mirror of each model */
    vector<double> mirror_fraction;
    /** Flag for checking mirror images of all references */
    bool do_mirror;
    /** Flag whether to fix estimates for model fractions */
    bool fix_fractions;
    /** Flag whether to fix estimate for sigma of noise */
    bool fix_sigma_noise;
    /** Flag whether to fix estimate for sigma the offsets */
    bool fix_sigma_offset;
    /** Starting iteration */
    int istart;
    /** Number of iterations to be performed */
    int Niter;
    /** dimension of the images */
    int dim;
    /** Number of reference images */
    int nr_ref;
    /** Number of in-plane rotation samples */
    int nr_psi, dnr_psi;
    /** Total number of experimental images */
    int nr_exp_images;
    /** Verbose level:
        1: gives progress bar (=default)
        0: gives no output to screen at all */
    int verb;
    /** Stopping criterium */
    double eps;
    /** SelFile images (working, test and reference set) */
    SelFile SF, SFr;
    /** Vector for images to hold references (new & old) */
    vector <ImageXmipp> Iref, Iold;
    /** Vector for FT rings for all references */
    vector <Polar <complex <double> > > fP_refs;
    /** Vector for sum2 for all references */
    vector < double > sum2_refs;
    /** Empty Polar */
    Polar <complex <double> > fP_zero;
    /** Number of limited translations */
    int nr_trans;
    /** pointers to original x and y translations */
    vector<double> Vxtrans, Vytrans;
    /** Matrices for calculating PDF of (in-plane) translations */
    Matrix2D<double> Mpdf_trans, Mr2;
    /** Limited search range for origin offsets */
    double search_shift;
    /** Limit orientational searches */
    bool limit_rot;
    /** Limited search range for projection directions */
    double search_rot;
    /** Vectors to store old phi, theta, xoff and yoff for all images */
    vector<float> imgs_oldphi, imgs_oldtheta, imgs_oldxoff, imgs_oldyoff;
    /** Flag for using ML3D */
    bool do_ML3D;
    /** Flag for generation of initial models from random subsets */
    bool do_generate_refs;
    /** debug flag */
    bool debug;
    /** Inner and outer radii for polar coordinates */
    int Ri, Ro;
    /** One common kb object for all images! */
    KaiserBessel kb;
    /** Precalculated Voronoi areas for the polar structure*/
    vector<double> voronoi_area;

public:
    /// @defgroup MLPfunctions Functions for MLPalign2D class
    /// @ingroup MLPalign2D
    
    /** Read command line
     * @ingroup MLPfunctions
     *
     * Read arguments from command line
     */
    void read(int argc, char **argv, bool ML3D = false);

    /** Show input information
     * @ingroup MLPfunctions
     *
     * Summarize input information
     */
    void show(bool ML3D = false);


    /** Usage
     * @ingroup MLPfunctions
     *
     * Standard usage
     */
    void usage();

    /** Extended Usage
     * @ingroup MLPfunctions
     *
     * Extended Usage
     */
    void extendedUsage(bool ML3D = false);

    /** Produce side info
     * @ingroup MLPfunctions
     *
     * Setup lots of (selfile-independent) stuff
     */
    void produceSideInfo();

    /** Generate initial reference
     * @ingroup MLPfunctions
     *
     * Generate references from averages of random subsets
     */
    void generateInitialReferences();

    /** Produce side info2
     * @ingroup MLPfunctions
     *
     * Setup lots of (selfile-dependent) stuff
     */
    void produceSideInfo2(int nr_vols = 1);

    /** Preselect directions
     * @ingroup MLPfunctions
     *
     * Calculate which references have projection directions close to
     * phi and theta
     */
    void preselectDirections(float &phi, float &theta,
			     vector<double> &pdf_directions);

    /** Update PDF of the translation
     * @ingroup MLPfunctions
     *
     * Calculate prior probabilities of the translations
     */
    void updatePdfTranslations();

    /** Prepare references
     * @ingroup MLPfunctions
     *
     * Calculate Fourier-transforms of all rings of all references
     * (interpolation based on reverse gridding)
     */
    void calculateFtRingsAllRefs(const vector<ImageXmipp> &Iref,
				 vector< Polar< complex <double> > > &fP_refs,
				 Polar< complex <double> > &fP_zero,
				 vector< double > &sum2_refs,
				 const int &first, const int &last);

    /** Prepare experimental image
     * @ingroup MLPfunctions
     *
     * Calculate Fourier-transforms of all rings of all translated
     * version of the experimental image
     * (interpolation based on reverse gridding)
     */
    void calculateFtRingsAllTransImg(const ImageXmipp &Iexp,
				     vector< Polar< complex <double> > > &fP_trans,
				     vector< Polar< complex <double> > > &fPm_trans,
				     double &Xi2, const int &first, const int &last);


    /** MLP integration over all refs, rots and trans
     * @ingroup MLPfunctions
     *
     * Here the expectation step of the EM-algorithm is performed.
     * and the weighted sums for each image are stored on-the-fly
     *
     */
    void processOneImage(const ImageXmipp &img,
			 const vector < Polar <complex <double> >  > &fP_refs,
			 const vector < double > &sum2_refs,
			 const vector < double > &pdf_directions,
			 vector < Polar <complex <double> > > &fP_wsum_imgs,
			 double &wsum_sigma_noise, double &wsum_sigma_offset, 
			 vector < double > &sumw, vector < double > &sumw_mirror,
			 double &LL, double &fracweight,
			 int &opt_refno, double &opt_psi, double &opt_flip, 
			 double &opt_xoff, double &opt_yoff);

    /** The actual loop over all images
     * @ingroup MLPfunctions
     *
     */
    void sumOverAllImages(SelFile &SF, const vector<ImageXmipp> &Iref,
			  double &LL, double &sumcorr, DocFile &DFo,
			  vector < Polar <complex <double> > > &fP_wsum_imgs,
			  double &wsum_sigma_noise, double &wsum_sigma_offset, 
			  vector <double> &sumw, vector <double> &sumw_mirror);

    /** Update all model parameters
     * @ingroup MLPfunctions
     *
     * Here the maximization step of the EM-algorithm is performed.
     *
     */
    void updateParameters(vector < Polar <complex <double> > > &fP_wsum_imgs,
			  double &wsum_sigma_noise, double &wsum_sigma_offset, 
			  vector <double> &sumw, vector <double> &sumw_mirror,
			  double &sumcorr, double &sumw_allrefs);

    /** Convergence check
     * @ingroup MLPfunctions
     *
     * Convergence is based on signal change in the cartesian-sampled images
     *
     */
    bool checkConvergence(vector<double> &conv);

    /** Write output files
     * @ingroup MLPfunctions
     *
     *  Write out reference images, selfile and logfile
     *
     */
    void writeOutputFiles(const int iter, DocFile &DFo,
                            double &sumw_allrefs, double &LL, double &avecorr,
                            vector<double> &conv);

};
//@}
