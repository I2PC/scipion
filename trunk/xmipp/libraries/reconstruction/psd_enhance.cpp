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
void Prog_Enhance_PSD_Parameters::read(int argc, char **argv)
{
    Prog_parameters::read(argc, argv);
    center = !checkParameter(argc, argv, "-dont_center");
    take_log = !checkParameter(argc, argv, "-dont_log");
    filter_w1 = textToFloat(getParameter(argc, argv, "-f1", "0.05"));
    filter_w2 = textToFloat(getParameter(argc, argv, "-f2", "0.2"));
    decay_width = textToFloat(getParameter(argc, argv, "-decay", "0.02"));
    mask_w1 = textToFloat(getParameter(argc, argv, "-m1", "0.025"));
    mask_w2 = textToFloat(getParameter(argc, argv, "-m2", "0.2"));
}

/* Usage ------------------------------------------------------------------- */
void Prog_Enhance_PSD_Parameters::usage()
{
    Prog_parameters::usage();
    std::cerr << "  [-dont_center]            : By default, it is assumed that the image\n"
              << "                              needs to be centered\n"
              << "  [-dont_log]               : Don't take log10 before working\n"
              << "  [-f1 <freq_low=0.05>]     : Low freq. for band pass filtration, max 0.5\n"
              << "  [-f2 <freq_high=0.2>]     : High freq. for band pass filtration, max 0.5\n"
              << "  [-decay <freq_decay=0.02>]: Decay for the transition bands\n"
              << "  [-m1 <freq_low=0.025>]    : Low freq. for mask, max 0.5\n"
              << "  [-m2 <freq_high=0.2>      : High freq. for mask, max 0.5\n"
    ;
}

/* Show -------------------------------------------------------------------- */
void Prog_Enhance_PSD_Parameters::show()
{
    Prog_parameters::show();
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
//#define DEBUG
void Prog_Enhance_PSD_Parameters::apply(Matrix2D<double> &PSD)
{
    // Take the logarithm
    if (take_log)
        FOR_ALL_ELEMENTS_IN_MATRIX2D(PSD)
        PSD(i, j) = log10(1 + PSD(i, j));

    // Remove single outliers
    if (center) CenterFFT(PSD, true);
    Matrix2D<double> aux;
    median_filter3x3(PSD, aux);
    PSD = aux;
    CenterFFT(PSD, false);

    // Reject other outliers
    reject_outliers(PSD, 2);

    // Band pass filter
    CenterFFT(PSD, true);
    FourierMask Filter;
    Filter.FilterShape = RAISED_COSINE;
    Filter.FilterBand = BANDPASS;
    Filter.w1 = filter_w1;
    Filter.w2 = filter_w2;
    Filter.raised_w = decay_width;
    PSD.setXmippOrigin();
    Filter.generate_mask(PSD);
    Filter.apply_mask_Space(PSD);
    STARTINGX(PSD) = STARTINGY(PSD) = 0;
    CenterFFT(PSD, false);

    // Mask the input PSD
    Matrix2D<int> mask;
    mask.resize(PSD);
    Matrix1D<int>    idx(2);  // Indexes for Fourier plane
    Matrix1D<double> freq(2); // Frequencies for Fourier plane
    FOR_ALL_ELEMENTS_IN_MATRIX2D(PSD)
    {
        XX(idx) = j;
        YY(idx) = i;
        FFT_idx2digfreq(PSD, idx, freq);
        if (freq.module() < mask_w1 || freq.module() > mask_w2) PSD(i, j) = 0;
        else mask(i, j) = 1;
    }

    //Compute the mean and the standard deviation under the mask
    //and normalize the PSD image
    double min_val, max_val, avg, stddev;
    computeStats_within_binary_mask(mask, PSD, min_val, max_val, avg, stddev);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(PSD)
    if (mask(i, j)) PSD(i, j) = (PSD(i, j) - avg) / stddev;

    // Mask again
    FOR_ALL_ELEMENTS_IN_MATRIX2D(PSD)
    {
        XX(idx) = j;
        YY(idx) = i;
        FFT_idx2digfreq(PSD, idx, freq);
        if (freq.module() < mask_w1 || freq.module() > mask_w2*0.9) PSD(i, j) = 0;
        else mask(i, j) = 1;
    }

    CenterFFT(PSD, true);
}
#undef DEBUG
