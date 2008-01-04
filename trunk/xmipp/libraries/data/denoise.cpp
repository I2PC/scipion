/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2002)
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

#include "denoise.h"
#include "args.h"
#include "wavelet.h"
#include "histogram.h"
#include "filters.h"

#include <iomanip>

// Empty constructor -------------------------------------------------------
Denoising_parameters::Denoising_parameters(): Prog_parameters()
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

// Read from command line --------------------------------------------------
void Denoising_parameters::read(int argc, char **argv)
{
    Prog_parameters::read(argc, argv);
    DWT_type = getParameter(argc, argv, "-type", "DAUB12");
    std::string aux = getParameter(argc, argv, "-denoising", "remove_scales");
    if (aux == "remove_scales")
        denoising_type = REMOVE_SCALE;
    else if (aux == "soft_thresholding")
        denoising_type = SOFT_THRESHOLDING;
    else if (aux == "bayesian")
        denoising_type = BAYESIAN;
    else if (aux == "adaptive_soft")
        denoising_type = ADAPTIVE_SOFT;
    else if (aux == "central")
        denoising_type = CENTRAL;
    else if (aux == "difussion")
        denoising_type = SHAH;
    else
        denoising_type = REMOVE_SCALE;
    scale = textToInteger(getParameter(argc, argv, "-scale", "0"));
    output_scale = textToInteger(getParameter(argc, argv, "-output_scale", "0"));
    threshold = textToFloat(getParameter(argc, argv, "-th", "50"));
    R = textToInteger(getParameter(argc, argv, "-R", "-1"));
    SNR0 = textToFloat(getParameter(argc, argv, "-SNR0", "0.1"));
    SNRF = textToFloat(getParameter(argc, argv, "-SNRF", "0.2"));
    white_noise = checkParameter(argc, argv, "-white_noise");
    Shah_outer = textToInteger(getParameter(argc, argv, "-outer", "10"));
    Shah_inner = textToInteger(getParameter(argc, argv, "-inner", "1"));
    Shah_refinement = textToInteger(getParameter(argc, argv, "-refinement", "1"));
    if (checkParameter(argc, argv, "-Shah_weight"))
        Shah_weight = getVectorParameter(argc, argv, "-Shah_weight", 4);
    else
    {
        Shah_weight.initZeros(4);
        Shah_weight(1) = Shah_weight(2) = 50;
        Shah_weight(3) = 0.02;
    }
    Shah_edge = checkParameter(argc, argv, "-only_edge");
    if (checkParameter(argc, argv, "-show"))
        tell = 1;

    produce_side_info();
}

// Produce side info -------------------------------------------------------
void Denoising_parameters::produce_side_info()
{
    if (DWT_type == "DAUB4")
        set_DWT_type(DAUB4);
    else if (DWT_type == "DAUB12")
        set_DWT_type(DAUB12);
    else if (DWT_type == "DAUB20")
        set_DWT_type(DAUB20);
    else
        REPORT_ERROR(1, "Denoising::produce_side_info: Unknown DWT type");
}

// Show --------------------------------------------------------------------
void Denoising_parameters::show()
{
    Prog_parameters::show();
    show_specific();
}

void Denoising_parameters::show_specific()
{
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

// Usage -------------------------------------------------------------------s
void Denoising_parameters::usage()
{
    Prog_parameters::usage();
    usage_specific();
}

void Denoising_parameters::usage_specific()
{
    std::cerr << "  [-type <str=DAUB12]        : DWT type. Valid types are:\n"
              << "                               DAUB4, DAUB12, DAUB20\n"
              << "  [-denoising <str=remove_scale] : Denoising action\n"
              << "                               remove_scale\n"
              << "                               bayesian\n"
              << "                               soft_thresholding\n"
              << "                               adaptive_soft\n"
              << "                               central\n"
              << "                               difussion\n"
              << "  [-scale <s=0>]             : scale\n"
              << "  [-output_scale <s=0>]      : output_scale\n"
              << "  [-th <th=50>]              : threshold of values (%) to remove\n"
              << "  [-R <r=-1>]                : Radius to keep, by default half the size\n"
              << "  [-SNR0 <SNR=0.1>]          : Smallest SNR\n"
              << "  [-SNRF <SNR=0.2>]          : Largest SNR\n"
              << "  [-white_noise]             : Select if the noise is white\n"
              << "  [-outer <it=10>]           : Difussion outer iterations\n"
              << "  [-inner <it=1>]            : Difussion inner iterations\n"
              << "  [-refinement <it=1>]       : Difussion refinement iterations\n"
              << "  [-Shah_weight [w0,w1,w2,w3]]:Diffusion weights\n"
              << "                               w0=data matching (=0)\n"
              << "                               w1=1st derivative smooth (=50)\n"
              << "                               w2=edge strength (=50)\n"
              << "                               w3=edge smoothness (=0.02)\n"
              << "  [-only_edge]               : Produce the edge image of the diffusion\n"
              << "  [-show]                    : Show information about the process\n"
     ;
}

// Denoise image -----------------------------------------------------------
void Denoising_parameters::denoise(Matrix2D<double> &img)
{
    if (denoising_type == BAYESIAN && adjust_range)
        img.range_adjust(0, 1);
    if (denoising_type != SHAH)
    {
        double size2 = log10((double)XSIZE(img)) / log10(2.0);
        if (ABS(size2 - ROUND(size2)) > 1e-6)
            REPORT_ERROR(1, "Denoising::denoise: Input image must be of a size power of 2");
        size2 = log10((double)YSIZE(img)) / log10(2.0);
        if (ABS(size2 - ROUND(size2)) > 1e-6)
            REPORT_ERROR(1, "Denoising::denoise: Input image must be of a size power of 2");
        DWT(img, img);
    }
    double th;
    Matrix2D<double> surface_strength, edge_strength;
    histogram1D hist;
    switch (denoising_type)
    {
    case REMOVE_SCALE:
        clean_quadrant(img, scale, "01");
        clean_quadrant(img, scale, "10");
        clean_quadrant(img, scale, "11");
        break;
    case SOFT_THRESHOLDING:
        compute_hist(img, hist, 100);
        soft_thresholding(img, hist.percentil(threshold));
        break;
    case BAYESIAN:
        estimatedS = bayesian_wiener_filtering(img, scale, SNR0, SNRF,
                                               white_noise, tell, !dont_denoise);
        break;
    case ADAPTIVE_SOFT:
        adaptive_soft_thresholding(img, scale);
        break;
    case CENTRAL:
        DWT_keep_central_part(img, R);
        break;
    case SHAH:
        Smoothing_Shah(img, surface_strength, edge_strength,
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

// Denoise volume ----------------------------------------------------------
void Denoising_parameters::denoise(Matrix3D<double> &vol)
{
    if (denoising_type == SHAH)
    {
        std::cerr << "Shah denoising is not implemented for volumes\n";
        return;
    }
    double size2 = log10((double)XSIZE(vol)) / log10(2.0);
    if (ABS(size2 - ROUND(size2)) > 1e-6)
        REPORT_ERROR(1, "Denoising::denoise: Input volume must be of a size power of 2");
    size2 = log10((double)YSIZE(vol)) / log10(2.0);
    if (ABS(size2 - ROUND(size2)) > 1e-6)
        REPORT_ERROR(1, "Denoising::denoise: Input volume must be of a size power of 2");
    size2 = log10((double)ZSIZE(vol)) / log10(2.0);
    if (ABS(size2 - ROUND(size2)) > 1e-6)
        REPORT_ERROR(1, "Denoising::denoise: Input volume must be of a size power of 2");

    DWT(vol, vol);
    double th;
    histogram1D hist;
    switch (denoising_type)
    {
    case REMOVE_SCALE:
        clean_quadrant(vol, scale, "001");
        clean_quadrant(vol, scale, "010");
        clean_quadrant(vol, scale, "011");
        clean_quadrant(vol, scale, "100");
        clean_quadrant(vol, scale, "101");
        clean_quadrant(vol, scale, "110");
        clean_quadrant(vol, scale, "111");
        break;
    case SOFT_THRESHOLDING:
        compute_hist(vol, hist, 100);
        soft_thresholding(vol, hist.percentil(threshold));
        break;
    case BAYESIAN:
        estimatedS = bayesian_wiener_filtering(vol, scale, SNR0, SNRF,
                                               white_noise, tell, !dont_denoise);
        break;
    case ADAPTIVE_SOFT:
        std::cout << "Adaptive soft-thresholding not implemented for volumes\n";
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
        vol.resize(ZSIZE(vol) / reduction, YSIZE(vol) / reduction, XSIZE(vol) / reduction);
    }
    IDWT(vol, vol);
}

void Denoising_parameters::denoise_avg_bayesian(Matrix3D<double> &vol)
{
    DWT(vol, vol);
    bayesian_wiener_filtering(vol, scale, estimatedS);

    if (output_scale != 0)
    {
        int reduction = (int)pow(2.0, output_scale);
        vol.resize(ZSIZE(vol) / reduction, YSIZE(vol) / reduction, XSIZE(vol) / reduction);
    }
    IDWT(vol, vol);
}
