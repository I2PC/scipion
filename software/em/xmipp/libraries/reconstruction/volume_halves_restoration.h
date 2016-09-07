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

class CDF
{
public:
	MultidimArray<double> x;
	MultidimArray<double> probXLessThanx;
	double minVal, maxVal;
	void calculateCDF(MultidimArray<double> &V, double probStep=0.005);
	double getProbability(double x);
};

/** Volume restoration Parameters. */
class ProgVolumeHalvesRestoration: public XmippProgram
{
public:
    /** Filename of the two halves and the output root */
    FileName fnV1, fnV2, fnRoot;
    /** Apply positivity */
    bool applyPos;
    /** Number of iterations */
    int Niter;
public:
    Image<double> V1, V2, V1r, V2r, S, HS, N;
    Mask mask;
    MultidimArray<int> *pMask;
    size_t pMaskSize;
    FourierTransformer transformer;
    MultidimArray< std::complex<double> > fVol;
    MultidimArray<double> H;
    MultidimArray<int> Ridx;

    CDF cdfHS;
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
    void applyMask(MultidimArray<double> &V);
    void estimateS();
    void estimateHS();
    void estimatePSDn();
    void significanceRealSpace(const MultidimArray<double> &V1, MultidimArray<double> &V1r);
};
//@}
#endif
