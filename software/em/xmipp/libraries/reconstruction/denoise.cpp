/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
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

#include "denoise.h"
#include <data/args.h>
#include <data/wavelet.h>
#include <data/histogram.h>
#include <data/filters.h>

#include <iomanip>

// Empty constructor -------------------------------------------------------
WaveletFilter::WaveletFilter()
{
    DWT_type = "DAUB12";
    denoising_type = REMOVE_SCALE;
    scale = 0;
    output_scale = 0;
    threshold = 50;
    R = -1;
    SNR0 = 1.0 / 10;
    SNRF = 1.0 / 5;
    white_noise = false;
    adjust_range = true;
    verbose = 0;
    dont_denoise = false;
}

// defineParams -------------------------------------------------------------------s
void WaveletFilter::defineParams(XmippProgram *program)
{
    program->addParamsLine("== Wavelet ==");
    program->addParamsLine("  [--wavelet <DWT_type=DAUB12> <mode=remove_scale>]   : Different types of filters using wavelets");
    program->addParamsLine("    where <DWT_type>");
    program->addParamsLine("       DAUB4 DAUB12 DAUB20    : Discrete Wavelet Transform");
    program->addParamsLine("    where <mode>");
    program->addParamsLine("       remove_scale");
    program->addParamsLine("       bayesian <SNR0=0.1> <SNRF=0.2> : Smallest(SNR0) and largest(SNRF) SNR.");
    program->addParamsLine("       soft_thresholding");
    program->addParamsLine("       adaptive_soft");
    program->addParamsLine("       central");
    program->addParamsLine("    alias -w;");
    program->addParamsLine("  [--scale+ <s=0>]             : scale");
    program->addParamsLine("  [--output_scale+ <s=0>]      : output_scale");
    program->addParamsLine("  [--th+ <th=50>]              : threshold of values (%) to remove");
    program->addParamsLine("  [-R+ <r=-1>]                 : Radius to keep, by default half the size");
    program->addParamsLine("  [--white_noise+]             : Select if the noise is white. Used by Bayesian filter.");
}

// Read from command line --------------------------------------------------
void WaveletFilter::readParams(XmippProgram *program)
{

    DWT_type = program->getParam("--wavelet", 0);
    String mode = program->getParam("--wavelet", 1);

    if (mode == "remove_scale")
        denoising_type = REMOVE_SCALE;
    else if (mode == "soft_thresholding")
        denoising_type = SOFT_THRESHOLDING;
    else if (mode == "bayesian")
    {
        SNR0 = program->getDoubleParam("--wavelet", 2);
        SNRF = program->getDoubleParam("--wavelet", 3);
        denoising_type = BAYESIAN;
    }
    else if (mode == "adaptive_soft")
        denoising_type = ADAPTIVE_SOFT;
    else if (mode == "central")
        denoising_type = CENTRAL;
    else
        REPORT_ERROR(ERR_DEBUG_IMPOSIBLE, "Bad argument type, this couldn't happens, check arguments parser!!!");

    scale = program->getIntParam("--scale");
    output_scale = program->getIntParam("--output_scale");
    threshold = program->getDoubleParam("--th");
    R = program->getIntParam("-R");
    white_noise = program->checkParam("--white_noise");
    verbose = program->verbose;
    produceSideInfo();
}

// Produce side info -------------------------------------------------------
void WaveletFilter::produceSideInfo()
{
    if (DWT_type == "DAUB4")
        set_DWT_type(DAUB4);
    else if (DWT_type == "DAUB12")
        set_DWT_type(DAUB12);
    else if (DWT_type == "DAUB20")
        set_DWT_type(DAUB20);
    else
        REPORT_ERROR(ERR_VALUE_INCORRECT, "Unknown DWT type");
}

// Show --------------------------------------------------------------------
void WaveletFilter::show()
{
    if (!verbose)
        return;
    ///Show specific options
    std::cout << "DWT type: " << DWT_type << std::endl;
    std::cout << "Denoising: ";
    switch (denoising_type)
    {
    case REMOVE_SCALE:
        std::cout << " Remove scale " << scale << std::endl;
        break;
    case SOFT_THRESHOLDING:
        std::cout << " Soft thresholding " << threshold << std::endl;
        break;
    case BAYESIAN:
        std::cout << " Bayesian\n";
        std::cout << " SNR between " << SNR0 << " and " << SNRF << std::endl
        << " up to scale " << scale << std::endl;
        if (white_noise)
            std::cout << " Imposing white noise\n";
        break;
    case ADAPTIVE_SOFT:
        std::cout << " Adaptive soft thresholding\n";
        break;
    case CENTRAL:
        std::cout << " Keeping central part " << R << " pixels\n";
        break;
    }
    std::cout << "Output scale: " << output_scale << std::endl;
}

// Denoise volume ----------------------------------------------------------
void WaveletFilter::apply(MultidimArray<double> &img)
{
    if (img.getDim()==2)
    {
        // 2D image denoising
        if (denoising_type == BAYESIAN && adjust_range)
            img.rangeAdjust(0, 1);

        double size2 = log10((double)XSIZE(img)) / log10(2.0);
        if (ABS(size2 - ROUND(size2)) > 1e-6)
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "Input image must be of a size power of 2");
        size2 = log10((double)YSIZE(img)) / log10(2.0);
        if (ABS(size2 - ROUND(size2)) > 1e-6)
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "Input image must be of a size power of 2");
        DWT(img, img);
        Histogram1D hist;
        switch (denoising_type)
        {
        case REMOVE_SCALE:
            clean_quadrant2D(img, scale, "01");
            clean_quadrant2D(img, scale, "10");
            clean_quadrant2D(img, scale, "11");
            break;
        case SOFT_THRESHOLDING:
            compute_hist(img, hist, 100);
            soft_thresholding(img, hist.percentil(threshold));
            break;
        case BAYESIAN:
            estimatedS = bayesian_wiener_filtering2D(img, scale, SNR0, SNRF,
                         white_noise, 0, !dont_denoise);
            break;
        case ADAPTIVE_SOFT:
            adaptive_soft_thresholding2D(img, scale);
            break;
        case CENTRAL:
            DWT_keep_central_part(img, R);
            break;
        }
            if (output_scale != 0)
            {
                int reduction = (int)pow(2.0, output_scale);
                img.resize(YSIZE(img) / reduction, XSIZE(img) / reduction);
            }
            IDWT(img, img);
    }
    else
    {
        // 3D image denoising
        double size2 = log10((double)XSIZE(img)) / log10(2.0);
        if (ABS(size2 - ROUND(size2)) > 1e-6)
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "Input volume must be of a size power of 2");
        size2 = log10((double)YSIZE(img)) / log10(2.0);
        if (ABS(size2 - ROUND(size2)) > 1e-6)
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "Input volume must be of a size power of 2");
        size2 = log10((double)ZSIZE(img)) / log10(2.0);
        if (ABS(size2 - ROUND(size2)) > 1e-6)
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "Input volume must be of a size power of 2");

        DWT(img, img);
        Histogram1D hist;
        switch (denoising_type)
        {
        case REMOVE_SCALE:
            clean_quadrant3D(img, scale, "001");
            clean_quadrant3D(img, scale, "010");
            clean_quadrant3D(img, scale, "011");
            clean_quadrant3D(img, scale, "100");
            clean_quadrant3D(img, scale, "101");
            clean_quadrant3D(img, scale, "110");
            clean_quadrant3D(img, scale, "111");
            break;
        case SOFT_THRESHOLDING:
            compute_hist(img, hist, 100);
            soft_thresholding(img, hist.percentil(threshold));
            break;
        case BAYESIAN:
            estimatedS = bayesian_wiener_filtering3D(img, scale, SNR0, SNRF,
                         white_noise, verbose, !dont_denoise);
            break;
        case ADAPTIVE_SOFT:
            std::cout << "Adaptive soft-thresholding not implemented for imgumes\n";
            break;
        case CENTRAL:
            std::cout << "Keep central part not implemented for volumes\n";
            break;
        }

        if (output_scale != 0)
        {
            int reduction = (int)pow(2.0, output_scale);
            img.resizeNoCopy(ZSIZE(img) / reduction, YSIZE(img) / reduction, XSIZE(img) / reduction);
        }
        IDWT(img, img);

    }
}

void WaveletFilter::denoiseAvgBayesian(MultidimArray<double> &vol)
{
    DWT(vol, vol);
    bayesian_wiener_filtering3D(vol, scale, estimatedS);

    if (output_scale != 0)
    {
        int reduction = (int)pow(2.0, output_scale);
        vol.resizeNoCopy(ZSIZE(vol) / reduction, YSIZE(vol) / reduction, XSIZE(vol) / reduction);
    }
    IDWT(vol, vol);
}
