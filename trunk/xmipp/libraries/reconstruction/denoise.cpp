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
ProgDenoise::ProgDenoise()
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
    tell = 0;
    dont_denoise = false;
}

// defineParams -------------------------------------------------------------------s
void ProgDenoise::defineParams()
{
    XmippMetadataProgram::defineParams();
    addParamsLine("  [-type <DWT_type=DAUB12>]   : Discrete Wavelet Transform");
    addParamsLine("    where <DWT_type> DAUB4 DAUB12 DAUB20");
    addParamsLine("    alias -t;");
    addParamsLine("  [-denoising <mode=remove_scale>]: Denoising action");
    addParamsLine("    where <mode>");
    addParamsLine("       remove_scale");
    addParamsLine("       bayesian <SNR0=0.1> <SNRF=0.2> : Smallest(SNR0) and largest(SNRF) SNR.");
    addParamsLine("       soft_thresholding");
    addParamsLine("       adaptive_soft");
    addParamsLine("       central");
    addParamsLine("       difussion");
    addParamsLine("    alias -d;");
    addParamsLine("==+ Advanced Options ==");
    addParamsLine("  [-scale <s=0>]             : scale");
    addParamsLine("  [-output_scale <s=0>]      : output_scale");
    addParamsLine("  [-th <th=50>]              : threshold of values (%) to remove");
    addParamsLine("  [-R <r=-1>]                : Radius to keep, by default half the size");
    addParamsLine("  [-white_noise]             : Select if the noise is white");
    addParamsLine("  [-outer <it=10>]           : Difussion outer iterations");
    addParamsLine("  [-inner <it=1>]            : Difussion inner iterations");
    addParamsLine("  [-refinement <it=1>]       : Difussion refinement iterations");
    addParamsLine("  [-Shah_weight <w0=0> <w1=50> <w2=50> <w3=0.02>]:Diffusion weights");
    addParamsLine("                             :  w0 = data matching ");
    addParamsLine("                             :  w1 = 1st derivative smooth ");
    addParamsLine("                             :  w2 = edge strength ");
    addParamsLine("                             :  w3 = edge smoothness ");
    addParamsLine("  [-only_edge]               : Produce the edge image of the diffusion");
    addParamsLine("  [-show]                    : Show information about the process");
}

// Read from command line --------------------------------------------------
void ProgDenoise::readParams()
{

    DWT_type = getParam("-type");
    std::string aux = getParam("-denoising");
    if (aux == "remove_scale")
        denoising_type = REMOVE_SCALE;
    else if (aux == "soft_thresholding")
        denoising_type = SOFT_THRESHOLDING;
    else if (aux == "bayesian")
    {
        SNR0 = getDoubleParam("-d", "bayesian", 0);
        SNRF = getDoubleParam("-d", "bayesian", 1);
        denoising_type = BAYESIAN;
    }
    else if (aux == "adaptive_soft")
        denoising_type = ADAPTIVE_SOFT;
    else if (aux == "central")
        denoising_type = CENTRAL;
    else if (aux == "difussion")
        denoising_type = SHAH;
    else
        REPORT_ERROR(ERR_DEBUG_IMPOSIBLE, "Bad argument type, this couldn't happens, check arguments parser!!!");

    scale = getIntParam("-scale");
    output_scale = getIntParam("-output_scale");
    threshold = getDoubleParam("-th");
    R = getIntParam("-R");

    white_noise = checkParam("-white_noise");
    Shah_outer = getIntParam( "-outer");
    Shah_inner = getIntParam( "-inner");
    Shah_refinement = getIntParam("-refinement");
    Shah_weight(0) = getDoubleParam("-Shah_weight", 0);
    Shah_weight(1) = getDoubleParam("-Shah_weight", 1);
    Shah_weight(2) = getDoubleParam("-Shah_weight", 2);
    Shah_weight(3) = getDoubleParam("-Shah_weight", 3);

    Shah_edge = checkParam("-only_edge");
    if (checkParam("-show"))
        tell = 1;

    produceSideInfo();
}

void ProgDenoise::processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
{
	Image<double> img;
    img.readApplyGeo(fnImg,mdIn,objId);
    denoise(img());
}

// Produce side info -------------------------------------------------------
void ProgDenoise::produceSideInfo()
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
void ProgDenoise::show()
{
    XmippMetadataProgram::show();

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
void ProgDenoise::denoise(MultidimArray<double> &img)
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
                         white_noise, tell, !dont_denoise);
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
                         white_noise, tell, !dont_denoise);
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

void ProgDenoise::denoiseAvgBayesian(MultidimArray<double> &vol)
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
