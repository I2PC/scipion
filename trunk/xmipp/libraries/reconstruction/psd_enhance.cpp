/***************************************************************************
 *
 * Authors:
 *
 * Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * and
 *
 * Slavica Jonic (slavica.jonic@impmc.jussieu.fr)
 * UPMC, Paris 6
 *
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

#include "psd_enhance.h"
#include "fourier_filter.h"

#include <data/args.h>
#include <data/filters.h>
#include <data/fft.h>

/* Read parameters --------------------------------------------------------- */
void ProgCTFEnhancePSD::readParams()
{
    XmippMetadataProgram::read(argc, argv);
    center = !checkParam("--dont_center");
    take_log = !checkParam("--dont_log");
    filter_w1 = getDoubleParam("--f1");
    filter_w2 = getDoubleParam("--f2");
    decay_width = getDoubleParam("--decay");
    mask_w1 = getDoubleParam("--m1");
    mask_w2 = getDoubleParam("--m2");
}

/* Usage ------------------------------------------------------------------- */
void ProgCTFEnhancePSD::defineParams()
{
    XmippMetadataProgram::defineParams();
    addParamsLine("  [--dont_center]            : By default, it is assumed that the image");
    addParamsLine("                             : needs to be centered");
    addParamsLine("  [--dont_log]               : Don't take log10 before working");
    addParamsLine("  [--f1 <freq_low=0.05>]     : Low freq. for band pass filtration, max 0.5");
    addParamsLine("  [--f2 <freq_high=0.2>]     : High freq. for band pass filtration, max 0.5");
    addParamsLine("  [--decay <freq_decay=0.02>]: Decay for the transition bands");
    addParamsLine("  [--m1 <freq_low=0.025>]    : Low freq. for frequency mask, max 0.5");
    addParamsLine("  [--m2 <freq_high=0.2>]     : High freq. for frequency mask, max 0.5");
}

/* Show -------------------------------------------------------------------- */
void ProgCTFEnhancePSD::show()
{
    XmippMetadataProgram::show();
    std::cout << "Centering:    " << center      << std::endl
    << "Log10:        " << take_log    << std::endl
    << "Filter w1:    " << filter_w1   << std::endl
    << "Filter w2:    " << filter_w2   << std::endl
    << "Filter decay: " << decay_width << std::endl
    << "Mask w1:      " << mask_w1     << std::endl
    << "Mask w2:      " << mask_w2     << std::endl
    ;
}


/* Apply ------------------------------------------------------------------- */
void ProgCTFEnhancePSD::processImage(const FileName &fnImg, const FileName &fnImgOut,
                                     size_t objId)
{
	Image<double> PSD;
	PSD.read(fnImg);
	if (ZSIZE(PSD())!=1)
		REPORT_ERROR(ERR_MATRIX_DIM,"This program is not intended for volumes");
	apply(PSD());
	PSD.write(fnImgOut);
}

//#define DEBUG
void ProgCTFEnhancePSD::apply(MultidimArray<double> &PSD)
{
    // Take the logarithm
    if (take_log)
        FOR_ALL_ELEMENTS_IN_ARRAY2D(PSD)
        A2D_ELEM(PSD, i, j) = log10(1 + A2D_ELEM(PSD, i, j));

    // Remove single outliers
    if (center)
        CenterFFT(PSD, true);
    MultidimArray<double> aux;
    median_filter3x3(PSD, aux);
    PSD = aux;
    CenterFFT(PSD, false);

    // Reject other outliers
    reject_outliers(PSD, 2);

    // Band pass filter
    CenterFFT(PSD, true);
    ProgFourierFilter Filter;
    Filter.FilterShape = RAISED_COSINE;
    Filter.FilterBand = BANDPASS;
    Filter.w1 = filter_w1;
    Filter.w2 = filter_w2;
    Filter.raised_w = decay_width;
    PSD.setXmippOrigin();
    Filter.generateMask(PSD);
    Filter.applyMaskSpace(PSD);
    STARTINGX(PSD) = STARTINGY(PSD) = 0;
    CenterFFT(PSD, false);

    // Mask the input PSD
    MultidimArray<int> mask;
    mask.resize(PSD);
    Matrix1D<int>    idx(2);  // Indexes for Fourier plane
    Matrix1D<double> freq(2); // Frequencies for Fourier plane
    FOR_ALL_ELEMENTS_IN_ARRAY2D(PSD)
    {
        XX(idx) = j;
        YY(idx) = i;
        FFT_idx2digfreq(PSD, idx, freq);
        if (freq.module() < mask_w1 || freq.module() > mask_w2)
            A2D_ELEM(PSD, i, j) = 0;
        else
            A2D_ELEM(mask, i, j) = 1;
    }

    //Compute the mean and the standard deviation under a tighter mask
    //close to the border and normalize the PSD image
    MultidimArray<int> tighterMask;
    tighterMask.resizeNoCopy(PSD);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(PSD)
    {
        XX(idx) = j;
        YY(idx) = i;
        FFT_idx2digfreq(PSD, idx, freq);
        if (freq.module() > mask_w2*0.9 && freq.module() < mask_w2)
            A2D_ELEM(tighterMask, i, j) = 1;
        else
            A2D_ELEM(tighterMask,i,j) = 0;
    }

    double min_val, max_val, avg, stddev;
    computeStats_within_binary_mask(tighterMask, PSD, min_val, max_val, avg, stddev);
    double istddev=1.0/stddev;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(PSD)
    if (A2D_ELEM(mask, i, j))
        A2D_ELEM(PSD, i, j) = (A2D_ELEM(PSD, i, j) - avg) *istddev;

    // Mask again
    FOR_ALL_ELEMENTS_IN_ARRAY2D(PSD)
    {
        XX(idx) = j;
        YY(idx) = i;
        FFT_idx2digfreq(PSD, idx, freq);
        if (freq.module() < mask_w1 || freq.module() > mask_w2*0.9)
            A2D_ELEM(PSD, i, j) = 0;
    }

    CenterFFT(PSD, true);
}
#undef DEBUG
