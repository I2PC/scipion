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

#ifndef _MLP_ALIGN2D_H
#define _MLP_ALIGN2D_H

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

/**@defgroup MLPalign2D mlp_align2d (Maximum likelihood in polar coordinates for 2D images)
   @ingroup ReconsLibraryPrograms */
#define SIGNIFICANT_WEIGHT_LOW 1e-8

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
    std::vector<double> alpha_k;
    /** Vector containing estimated fraction for mirror of each model */
    std::vector<double> mirror_fraction;
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
    std::vector < ImageXmippT<double> > Iref, Iold;
    /** Vector for FT rings for all references */
    std::vector <Polar <std::complex <double> > > fP_refs;
    /** Vector for sum2 for all references */
    std::vector < double > sum2_refs;
    /** Empty Polar */
    Polar <std::complex <double> > fP_zero;
    /** Number of limited translations */
    int nr_trans;
    /** pointers to original x and y translations */
    std::vector<double> Vxtrans, Vytrans;
    /** Matrices for calculating PDF of (in-plane) translations */
    Matrix2D<double> Mpdf_trans, Mr2;
    /** Limited search range for origin offsets */
    double search_shift;
    /** Limit orientational searches */
    bool limit_rot;
    /** Limited search range for projection directions */
    double search_rot;
    /** Vectors to store old phi, theta, xoff and yoff for all images */
    std::vector<float> imgs_oldphi, imgs_oldtheta, imgs_oldxoff, imgs_oldyoff;
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
    std::vector<double> voronoi_area;

public:
    /** Read command line
     *
     * Read arguments from command line
     */
    void read(int argc, char **argv, bool ML3D = false);

    /** Show input information
     *
     * Summarize input information
     */
    void show(bool ML3D = false);


    /** Usage
     *
     * Standard usage
     */
    void usage();

    /** Extended Usage
     *
     * Extended Usage
     */
    void extendedUsage(bool ML3D = false);

    /** Produce side info
     *
     * Setup lots of (selfile-independent) stuff
     */
    void produceSideInfo();

    /** Generate initial reference
     *
     * Generate references from averages of random subsets
     */
    void generateInitialReferences();

    /** Produce side info2
     *
     * Setup lots of (selfile-dependent) stuff
     */
    void produceSideInfo2(int nr_vols = 1);

    /** Preselect directions
     *
     * Calculate which references have projection directions close to
     * phi and theta
     */
    void preselectDirections(float &phi, float &theta,
                             std::vector<double> &pdf_directions);

    /** Update PDF of the translation
     *
     * Calculate prior probabilities of the translations
     */
    void updatePdfTranslations();

    /** Prepare references
     *
     * Calculate Fourier-transforms of all rings of all references
     * (interpolation based on reverse gridding)
     */
    void calculateFtRingsAllRefs(const std::vector< ImageXmippT<double> > &Iref,
                                 std::vector< Polar< std::complex <double> > > &fP_refs,
                                 Polar< std::complex <double> > &fP_zero,
                                 std::vector< double > &sum2_refs,
                                 const int &first, const int &last);

    /** Prepare experimental image
     *
     * Calculate Fourier-transforms of all rings of all translated
     * version of the experimental image
     * (interpolation based on reverse gridding)
     */
    void calculateFtRingsAllTransImg(const  ImageXmippT<double>  &Iexp,
                                     std::vector< Polar< std::complex <double> > > &fP_trans,
                                     std::vector< Polar< std::complex <double> > > &fPm_trans,
                                     double &Xi2, const int &first, const int &last);


    /** MLP integration over all refs, rots and trans
     *
     * Here the expectation step of the EM-algorithm is performed.
     * and the weighted sums for each image are stored on-the-fly
     *
     */
    void processOneImage(const ImageXmippT<double> &img,
                         const std::vector < Polar <std::complex <double> >  > &fP_refs,
                         const std::vector < double > &sum2_refs,
                         const std::vector < double > &pdf_directions,
                         std::vector < Polar <std::complex <double> > > &fP_wsum_imgs,
                         double &wsum_sigma_noise, double &wsum_sigma_offset,
                         std::vector < double > &sumw, std::vector < double > &sumw_mirror,
                         double &LL, double &fracweight,
                         int &opt_refno, double &opt_psi, double &opt_flip,
                         double &opt_xoff, double &opt_yoff);

    /** The actual loop over all images
     */
    void sumOverAllImages(SelFile &SF, const std::vector< ImageXmippT<double> > &Iref,
                          double &LL, double &sumcorr, DocFile &DFo,
                          std::vector < Polar <std::complex <double> > > &fP_wsum_imgs,
                          double &wsum_sigma_noise, double &wsum_sigma_offset,
                          std::vector <double> &sumw, std::vector <double> &sumw_mirror);

    /** Update all model parameters
     * Here the maximization step of the EM-algorithm is performed.
     *
     */
    void updateParameters(std::vector < Polar <std::complex <double> > > &fP_wsum_imgs,
                          double &wsum_sigma_noise, double &wsum_sigma_offset,
                          std::vector <double> &sumw, std::vector <double> &sumw_mirror,
                          double &sumcorr, double &sumw_allrefs);

    /** Convergence check
     * Convergence is based on signal change in the cartesian-sampled images
     *
     */
    bool checkConvergence(std::vector<double> &conv);

    /** Write output files
     *  Write out reference images, selfile and logfile
     *
     */
    void writeOutputFiles(const int iter, DocFile &DFo,
                          double &sumw_allrefs, double &LL, double &avecorr,
                          std::vector<double> &conv);

};
//@}
#endif
