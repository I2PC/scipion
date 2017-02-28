/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
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

#ifndef _PROG_VOLUME_HALVES_RESTORATION
#define _PROG_VOLUME_HALVES_RESTORATION

#include <data/xmipp_program.h>
#include <data/mask.h>
#include <data/xmipp_fftw.h>

/**@defgroup VolumeHalvesRestoration Restoration of a volume given two halves
   @ingroup ReconsLibrary */
//@{


/** Volume restoration Parameters. */
class ProgVolumeHalvesRestoration: public XmippProgram
{
public:
    /** Filename of the two halves and the output root */
    FileName fnV1, fnV2, fnRoot;
    /** Number of iterations of real space denoising */
    int NiterReal;
    /** Number of iterations of Fourier space deconvolution */
    int NiterFourier;
    /** FilterBank */
    double bankStep, bankOverlap;
    /** Initial sigma */
    double sigma0;
    /** Laplacian regularization */
    double lambda;
    /** Weight function */
    int weightFun;
    /** Weight power */
    double weightPower;
public:
    Image<double> V1, V2, V1r, V2r, S, N;
    Mask mask;
    MultidimArray<int> *pMask;
    size_t pMaskSize;
    FourierTransformer transformer, transformer1, transformer2;
    MultidimArray< std::complex<double> > fVol, fV1r, fV2r;
    MultidimArray<double> R2;

    CDF cdfS;
    double sigmaConv1, sigmaConv2; // Sigma for the convolution
public:
    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /// Run
    void run();

    /// Side information
    void produceSideInfo();

    /// Different estimates
    void estimateS();
    void deconvolveS();
    void convolveS();
    void optimizeSigma();
    void significanceRealSpace(const MultidimArray<double> &V1, MultidimArray<double> &V1r);
    void filterBank();
    void filterBand(const MultidimArray< std::complex<double> > &Vin, FourierTransformer &transformer, double w);
};
//@}
#endif
