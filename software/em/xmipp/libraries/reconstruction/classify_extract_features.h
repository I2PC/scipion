/***************************************************************************
 *
 * Authors:    Tomas Majtner           tmajtner@cnb.csic.es (2017)
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

#ifndef _PROG_EXTRACT_FEATURES
#define _PROG_EXTRACT_FEATURES

#include <data/xmipp_program.h>


class ProgExtractFeatures: public XmippProgram
{
public:
    /** Filename selection file containing the images */
    FileName fnSel;

    /**  Filename output root */
    FileName fnOut;

    /**  Parameter for turning off denoising */
    bool noDenoising;

    /**  Parameter for using LBP */
    bool useLBP;

    /**  Parameter for using entropy features */
    bool useEntropy;

    /**  Parameter for using variance features */
    bool useVariance;

    /**  Parameter for using Zernike moments */
    bool useZernike;

    /** Parameter for ramp */
    bool useRamp;

public:
    ProgExtractFeatures();
    ~ProgExtractFeatures();

    /// Read argument
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /// Function for returning factorial up to n=4
    int facs(int n);

    /// Extracting LBP features
    /// See method at Ojala, Timo, Matti Pietikainen, and Topi Maenpaa.
    /// "Multiresolution gray-scale and rotation invariant texture
    /// classification with local binary patterns." IEEE Transactions on
    /// Pattern Analysis and Machine Intelligence 24.7 (2002): 971-987.
    void extractLBP(const MultidimArray<double> &I, std::vector<double> &fv);

    /// Extracting entropy features
    void extractEntropy(const MultidimArray<double> &I, MultidimArray<double> &Imasked, std::vector<double> &fv);

    /// Extracting variance features
    void extractVariance(const MultidimArray<double> &I, std::vector<double> &fv);

    /// Extracting Zernike moments
    /// See method at A. Tahmasbi, F. Saki, S. B. Shokouhi, Classification
    /// of Benign and Malignant Masses Based on Zernike Moments,
    /// Comput. Biol. Med., vol. 41, no. 8, pp. 726-735, 2011
    /// and also
    /// F. Saki, A. Tahmasbi, H. Soltanian-Zadeh, S. B. Shokouhi,
    /// Fast opposite weight learning rules with application in breast
    /// cancer diagnosis, Comput. Biol. Med., vol. 43, no. 1, pp. 32-41, 2013.
    void extractZernike(const MultidimArray<double> &I, std::vector<double> &fv);

    /// Gray ramp coefficients
    void extractRamp(const MultidimArray<double> &I, std::vector<double> &fv);

    /// Main routine
    void run();

public:
    /**  Masks for entropy features */
    MultidimArray<int> masks[7];

    MultidimArray<int> rampMask;
    FitPoint *fitPoints;
    int NmaskPoints;
};
#endif
