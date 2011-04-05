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
DenoiseFilter::DenoiseFilter()
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
    Shah_outer = 10;
    Shah_inner = 1;
    Shah_refinement = 1;
    Shah_weight.initZeros(4);
    Shah_weight(1) = Shah_weight(2) = 50;
    Shah_weight(3) = 0.02;
    Shah_edge = false;
    adjust_range = true;
    verbose = 0;
    dont_denoise = false;
}

// defineParams -------------------------------------------------------------------s
void DenoiseFilter::defineParams(XmippProgram *program)
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
    program->addParamsLine("       difussion");
    program->addParamsLine("    alias -w;");
    program->addParamsLine("  [--scale+ <s=0>]             : scale");
    program->addParamsLine("  [--output_scale+ <s=0>]      : output_scale");
    program->addParamsLine("  [--th+ <th=50>]              : threshold of values (%) to remove");
    program->addParamsLine("  [-R+ <r=-1>]                 : Radius to keep, by default half the size");
    program->addParamsLine("  [--white_noise+]             : Select if the noise is white");
    program->addParamsLine("  [--shah_iter+ <outer=10> <inner=1> <refinement=1>]  : Diffusion outer, inner and refinement iterations");
    program->addParamsLine("  [--shah_weight+ <w0=0> <w1=50> <w2=50> <w3=0.02>]:Diffusion weights");
    program->addParamsLine("                             :  w0 = data matching ");
    program->addParamsLine("                             :  w1 = 1st derivative smooth ");
    program->addParamsLine("                             :  w2 = edge strength ");
    program->addParamsLine("                             :  w3 = edge smoothness ");
    program->addParamsLine("  [--shah_only_edge+]         : Produce the edge image of the diffusion");
}

// Read from command line --------------------------------------------------
void DenoiseFilter::readParams(XmippProgram *program)
{

    DWT_type = program->getParam("--denoise", 0);
    String mode = program->getParam("--denoise", 1);

    if (mode == "remove_scale")
        denoising_type = REMOVE_SCALE;
    else if (mode == "soft_thresholding")
        denoising_type = SOFT_THRESHOLDING;
    else if (mode == "bayesian")
    {
        SNR0 = program->getDoubleParam("-d", "bayesian", 0);
        SNRF = program->getDoubleParam("-d", "bayesian", 1);
        denoising_type = BAYESIAN;
    }
    else if (mode == "adaptive_soft")
        denoising_type = ADAPTIVE_SOFT;
    else if (mode == "central")
        denoising_type = CENTRAL;
    else if (mode == "difussion")
    {
        denoising_type = SHAH;
        Shah_outer = program->getIntParam( "--shah_iter", 0);
        Shah_inner = program->getIntParam( "--shah_iter", 1);
        Shah_refinement = program->getIntParam("--shah_iter", 2);
        Shah_weight(0) = program->getDoubleParam("--shah_weight", 0);
        Shah_weight(1) = program->getDoubleParam("--shah_weight", 1);
        Shah_weight(2) = program->getDoubleParam("--shah_weight", 2);
        Shah_weight(3) = program->getDoubleParam("--shah_weight", 3);
        Shah_edge = program->checkParam("--shah_only_edge");
    }
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
void DenoiseFilter::produceSideInfo()
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
void DenoiseFilter::show()
{
  if (!verbose)
    return;
    ///Show specific options
    if (denoising_type != SHAH)
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
    case SHAH:
        std::cout << " Shah difussion\n"
        << " Outer iterations " << Shah_outer << std::endl
        << " Inner iterations " << Shah_inner << std::endl
        << " Refinement interations " << Shah_refinement << std::endl
        << " Weight " << Shah_weight.transpose() << std::endl;
        if (Shah_edge)
            std::cout << " Generating edge image\n";
        break;
    }
    if (denoising_type != SHAH)
        std::cout << "Output scale: " << output_scale << std::endl;
}

// Denoise volume ----------------------------------------------------------
void DenoiseFilter::apply(MultidimArray<double> &img)
{
    if (img.getDim()==2)
    {
        // 2D image denoising
        if (denoising_type == BAYESIAN && adjust_range)
            img.rangeAdjust(0, 1);
        if (denoising_type != SHAH)
        {
            double size2 = log10((double)XSIZE(img)) / log10(2.0);
            if (ABS(size2 - ROUND(size2)) > 1e-6)
                REPORT_ERROR(ERR_MULTIDIM_SIZE, "Input image must be of a size power of 2");
            size2 = log10((double)YSIZE(img)) / log10(2.0);
            if (ABS(size2 - ROUND(size2)) > 1e-6)
                REPORT_ERROR(ERR_MULTIDIM_SIZE, "Input image must be of a size power of 2");
            DWT(img, img);
        }
        double th;
        MultidimArray<double> surface_strength, edge_strength;
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
                         white_noise, verbose, !dont_denoise);
            break;
        case ADAPTIVE_SOFT:
            adaptive_soft_thresholding2D(img, scale);
            break;
        case CENTRAL:
            DWT_keep_central_part(img, R);
            break;
        case SHAH:
            smoothingShah(img, surface_strength, edge_strength,
                           Shah_weight, Shah_outer, Shah_inner,
                           Shah_refinement, adjust_range);
            if (Shah_edge)
                img = edge_strength;
            else
                img = surface_strength;
            break;
        }
        if (denoising_type != SHAH)
        {
            if (output_scale != 0)
            {
                int reduction = (int)pow(2.0, output_scale);
                img.resize(YSIZE(img) / reduction, XSIZE(img) / reduction);
            }
            IDWT(img, img);
        }


    }
    else
    {
        // 3D image denoising
        if (denoising_type == SHAH)
        {
            std::cerr << "Shah denoising is not implemented for volumes\n";
            return;
        }
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
        double th;
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
        case SHAH:
            std::cout << "Shah Difussion not implemented for volumes\n";
            break;
        }

        if (output_scale != 0)
        {
            int reduction = (int)pow(2.0, output_scale);
            img.resize(ZSIZE(img) / reduction, YSIZE(img) / reduction, XSIZE(img) / reduction);
        }
        IDWT(img, img);

    }
}

void DenoiseFilter::denoiseAvgBayesian(MultidimArray<double> &vol)
{
    DWT(vol, vol);
    bayesian_wiener_filtering3D(vol, scale, estimatedS);

    if (output_scale != 0)
    {
        int reduction = (int)pow(2.0, output_scale);
        vol.resize(ZSIZE(vol) / reduction, YSIZE(vol) / reduction, XSIZE(vol) / reduction);
    }
    IDWT(vol, vol);
}
