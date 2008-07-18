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

#ifndef _MLF_TOMO_H
#define _MLF_TOMO_H

#include <data/fftw.h>
#include <data/fft.h>
#include <data/args.h>
#include <data/funcs.h>
#include <data/selfile.h>
#include <data/docfile.h>
#include <data/image.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/geometry.h>
#include "symmetrize.h"
#include <vector>

/**@defgroup mlf_tomo mlf_tomo (Maximum likelihood in Fourier space for tomography)
   @ingroup ReconsLibraryPrograms */
//@{
#define SIGNIFICANT_WEIGHT_LOW 1e-8

/** mlf_tomo parameters. */
class Prog_mlf_tomo_prm
{

public:
    /** Filenames reference selfile/image, fraction docfile & output rootname */
    SelFile SFi, SFr, SFg;
    FileName fn_ref, fn_root, fn_doc, fn_wlist, fn_group, fn_prior;
    /** Flag whether to fix estimates for model fractions */
    bool fix_fractions;
    /** Flag whether to fix estimate for sigma of noise */
    bool fix_sigma_noise;
    /** Calculate FFTW files again if they exist? */
    bool dont_recalc_fftw;
    /** Use TOM Toolbox Euler angle and translation conventions */
    bool use_tom_conventions;
  // For all tomograms: angles, offsets and wedge parameters
    std::vector<double> img_ang1, img_ang2, img_ang3;
    std::vector<double> img_xoff, img_yoff, img_zoff;
    std::vector<double> img_th0, img_thF;

    /** Use FFTW for Fts */
    bool do_fftw;
    /** Sizes of all data vectors (size for complex-double, hsize for double) */
    int size, hsize, fftw_hsize;
    /** Array with booleans whether this point is inside the resolution range */
    bool * is_in_range;
    /** Vector containing estimated fraction for each model */
    std::vector<double> alpha_k;
    /** Matrix1d with number of images per group */
    Matrix1D<double> nr_imgs_per_group;
    /** Vector containing resolution-depedent noise sigma values for each tomogram group */
    std::vector< Matrix1D< double > > sigma_noise;
    /** Starting iteration */
    int istart;
    /** Number of iterations to be performed */
    int Niter;
    /** dimensions of the images  */
    int Xdim, Ydim, Zdim;
    double dim3;
    /** Number of references */
    int nr_ref;
    /** Number of input images */
    double nr_imgs;
    /** Number of groups (e.g. one group for each tomogram) */
    int nr_group;
    /** Regularisation constants */
    double reg0, regF, reg_steps, reg, delta_reg;

    /** Verbose level:
        2: gives debug statements
        1: gives progress bar (=default)
        0: gives no output to screen at all */
    int verb;

    /* High/low and probability-calculation resolution limits */
    double highres, lowres;


    double ang;
    xmippFftw forwfftw, backfftw;

    //// DEBUGGING
    bool debug;

public:

    /// Read arguments from command line
    void read(int argc, char **argv);

    /// Show
    void show();

    /// Usage
    void usage();

    /// Extended Usage
    void extendedUsage();

    /// splitted SF-INdependent side-info calculations
    void produceSideInfo();

    /// splitted SF-dependent side-info calculations
    void produceSideInfo2();

    // Get Transformation matrix (acc to XMIPP or TOM conventions, depending on use_tom_conventions flag)
    // For XMIPP ang1-3 are: rot, tilt, psi
    // For TOM   ang1-3 are: phi, psi, theta
    Matrix2D< double > getTransformationMatrix(double ang1, 
                                               double ang2, 
                                               double ang3, 
                                               double xoff = 0.,
                                               double yoff = 0.,
                                               double zoff = 0.);
    /// Get binary missing wedge (or pyramid) 
    void getMissingWedge(bool * measured,
                         Matrix2D<double> A,
                         const double theta0_alongy, 
                         const double thetaF_alongy,
                         const double theta0_alongx = 0.,
                         const double thetaF_alongx = 0.);

    /// Generate initial references as averages over random subsets
    void generateInitialReferences(double * dataRefs, double * dataWsumWedsPerRef);
 
    /// Calculate FFTWs for images in a selfile and write to disc
    void calculateAllFFTWs();

    /// initial sigma2-estimation for fourier-mode
    void estimateInitialNoiseSpectra(double * dataSigma);

   /// parameters using fourier-space likelihood functions
    void expectationSingleImage(int igroup,
                                double * dataImg,
                                bool   * dataMeasured,
                                double * dataRefs,
                                double * dataSigma,
                                double * dataWsumRefs,
                                double * dataWsumWedsPerRef,
                                double * dataWsumWedsPerGroup,
                                double * dataWsumDist,
                                double * dataSumWRefs,
                                int    & opt_refno, 
                                double & LL, 
                                double & Pmax);

    /// Integrate over all experimental images
    void expectation(double * dataRefs,
                     double * dataSigma,
                     double * dataWsumRefs,
                     double * dataWsumWedsPerRef,
                     double * dataWsumWedsPerGroup,
                     double * dataWsumDist,
                     double * dataSumWRefs,
                     double & LL,
                     double & avePmax,
                     DocFile & DFo);

    /// Update all model parameters
    void maximization(double * dataRefs,
                      double * dataSigma,
                      double * dataWsumRefs,
                      double * dataWsumWedsPerRef,
                      double * dataWsumWedsPerGroup,
                      double * dataWsumDist,
                      double * dataSumWRefs,
                      double & sumw_allrefs,
                      double & avePmax);

    // Regularisation
    void regularise(double * dataRefs,
                    double * dataSigma);

    /// Write out reference images, selfile and logfile
    void writeOutputFiles(int iter, 
                          double * dataRefs,
                          double & sumw_allrefs, 
                          double & LL, 
                          double & avePmax,
                          std::vector<double> &conv);

};
//@}
#endif
