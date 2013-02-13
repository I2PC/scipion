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

#include "ctf_enhance_psd.h"
#include "fourier_filter.h"

#include <data/args.h>
#include <data/filters.h>
#include <data/xmipp_fft.h>

/* Read parameters --------------------------------------------------------- */
void ProgCTFEnhancePSD::readParams()
{
    XmippMetadataProgram::readParams();
    method = getParam("--method");
    if (method == "filter")
    {
        filter_w1 = getDoubleParam("--method", 1);
        filter_w2 = getDoubleParam("--method", 2);
        decay_width = getDoubleParam("--method", 3);
    }
    else
    {
        N0 = getIntParam("--method", 1);
        NF = getIntParam("--method", 2);
    }
    mask_w1 = getDoubleParam("--m1");
    mask_w2 = getDoubleParam("--m2");
}

/* Usage ------------------------------------------------------------------- */
void ProgCTFEnhancePSD::defineParams()
{
    each_image_produces_an_output = true;
    defaultComments["-i"].clear();
    defaultComments["-i"].addComment("Metadata with PSDs or a single PSD");
    XmippMetadataProgram::defineParams();
    addUsageLine(
        "Enhances the visibility of the Thon rings in a Power Spectrum Density (PSD).");
    addSeeAlsoLine("ctf_estimate_from_micrograph");
    addParamsLine("[--method <mth=filter>]        : Choose enhancing method");
    addParamsLine("       where <mth>");
    addParamsLine("             filter <freq_low=0.05> <freq_high=0.2> <freq_decay=0.02>: Ad hoc filters. The algorithm is fully described at [[http://www.ncbi.nlm.nih.gov/pubmed/16987671][this article]]");
    addParamsLine("                                                                     : The frequency limits define a raised-cosine bandpass filter, frequencies are normalized to 0.5");
    addParamsLine("             spht   <N0=1> <NF=10> : Spiral phase transform normalization.");
    addParamsLine("                                   : N0 and NF define the minimum and maximum number of fringes in the CTF");
    addParamsLine("==+ Output parameters");
    addParamsLine("  [--dont_center]            : By default, the output is centered");
    addParamsLine("  [--dont_log]               : Don't take log10 before working");
    addParamsLine("==+ Output mask parameters");
    addParamsLine("  [--m1 <freq_low=0.025>]    : Low freq. for frequency mask, max 0.5");
    addParamsLine("  [--m2 <freq_high=0.3>]     : High freq. for frequency mask, max 0.5");
    addExampleLine("xmipp_ctf_enhance_psd -i myPSD.psd -o myPSD_enhanced.psd");
}

/* Show -------------------------------------------------------------------- */
void ProgCTFEnhancePSD::show()
{
    XmippMetadataProgram::show();
    std::cout << "Filter w1:    " << filter_w1 << std::endl << "Filter w2:    "
    << filter_w2 << std::endl << "Filter decay: " << decay_width
    << std::endl << "Mask w1:      " << mask_w1 << std::endl
    << "Mask w2:      " << mask_w2 << std::endl;
}

/* Apply ------------------------------------------------------------------- */
void ProgCTFEnhancePSD::processImage(const FileName &fnImg, const FileName &fnImgOut,
									 const MDRow &rowIn,MDRow &rowOut)
{
    Image<double> PSD;
    PSD.read(fnImg);
    if (ZSIZE(PSD()) != 1)
        REPORT_ERROR(ERR_MATRIX_DIM,"This program is not intended for volumes");
    if (method=="filter")
        applyFilter(PSD());
    else
        applySPHT(PSD());
    PSD.write(fnImgOut);
}

//#define DEBUG
void ProgCTFEnhancePSD::applyFilter(MultidimArray<double> &PSD)
{
    // Take the logarithm
    FOR_ALL_ELEMENTS_IN_ARRAY2D(PSD)
    A2D_ELEM(PSD, i, j) = log10(1 + A2D_ELEM(PSD, i, j));

    // Remove single outliers
    CenterFFT(PSD, true);
    MultidimArray<double> aux;
    medianFilter3x3(PSD, aux);
    PSD = aux;

    // Reject other outliers
    reject_outliers(PSD, 2);

    // Band pass filter
    FourierFilter Filter;
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
    Matrix1D<double> freq(2); // Frequencies for Fourier plane
    double limit0_2 = mask_w1;
    limit0_2 = limit0_2 * limit0_2;
    double limitF_2 = mask_w2;
    limitF_2 = limitF_2 * limitF_2;
    for (int i = STARTINGY(PSD); i <= FINISHINGY(PSD); ++i)
    {
        FFT_IDX2DIGFREQ(i, YSIZE(PSD), YY(freq));
        double freqy2 = YY(freq) * YY(freq);
        for (int j = STARTINGX(PSD); j <= FINISHINGX(PSD); ++j)
        {
            FFT_IDX2DIGFREQ(j, XSIZE(PSD), XX(freq));
            double freq2 = XX(freq) * XX(freq) + freqy2;
            if (freq2 < limit0_2 || freq2 > limitF_2)
                A2D_ELEM(PSD, i, j) = 0;
            else
                A2D_ELEM(mask, i, j) = 1;
        }
    }

    //Compute the mean and the standard deviation under a tighter mask
    //close to the border and normalize the PSD image
    MultidimArray<int> tighterMask;
    tighterMask.resizeNoCopy(PSD);
    limit0_2 = mask_w2 * 0.9;
    limit0_2 = limit0_2 * limit0_2;
    limitF_2 = mask_w2;
    limitF_2 = limitF_2 * limitF_2;
    for (int i = STARTINGY(PSD); i <= FINISHINGY(PSD); ++i)
    {
        FFT_IDX2DIGFREQ(i, YSIZE(PSD), YY(freq));
        double freqy2 = YY(freq) * YY(freq);
        for (int j = STARTINGX(PSD); j <= FINISHINGX(PSD); ++j)
        {
            FFT_IDX2DIGFREQ(j, XSIZE(PSD), XX(freq));
            double freq2 = XX(freq) * XX(freq) + freqy2;
            A2D_ELEM(tighterMask, i, j) = freq2 > limit0_2 && freq2 < limitF_2;
        }
    }

    double avg, stddev;
    computeAvgStdev_within_binary_mask(tighterMask, PSD, avg, stddev);
    double istddev = 1.0 / stddev;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(PSD)
    if (A2D_ELEM(mask, i, j))
        A2D_ELEM(PSD, i, j) = (A2D_ELEM(PSD, i, j) - avg) * istddev;

    // Mask again
    limit0_2 = mask_w1;
    limit0_2 = limit0_2 * limit0_2;
    limitF_2 = mask_w2 * 0.9;
    limitF_2 = limitF_2 * limitF_2;
    for (int i = STARTINGY(PSD); i <= FINISHINGY(PSD); ++i)
    {
        FFT_IDX2DIGFREQ(i, YSIZE(PSD), YY(freq));
        double freqy2 = YY(freq) * YY(freq);
        for (int j = STARTINGX(PSD); j <= FINISHINGX(PSD); ++j)
        {
            FFT_IDX2DIGFREQ(j, XSIZE(PSD), XX(freq));
            double freq2 = XX(freq) * XX(freq) + freqy2;
            if (freq2 < limit0_2 || freq2 > limitF_2)
                A2D_ELEM(PSD, i, j) = 0;
        }
    }

    CenterFFT(PSD, true);
}

void ProgCTFEnhancePSD::applySPHT(MultidimArray<double> &PSD)
{
    FourierTransformer transformer;
    MultidimArray< std::complex<double> > PSDfourier;
    transformer.FourierTransform(PSD, PSDfourier, false);

    transformer.inverseFourierTransform();
}
#undef DEBUG
