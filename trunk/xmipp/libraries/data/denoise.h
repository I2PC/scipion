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

#ifndef DENOISE_H
#define DENOISE_H

#include <string>

#include "progs.h"

/// @defgroup Denoise Image denoising
/// @ingroup DataLibraryPrograms

/// Parameters for denoise program
/// @ingroup Denoise
class Denoising_parameters: public Prog_parameters
{
public:
    typedef enum
    {
        REMOVE_SCALE,
        SOFT_THRESHOLDING,
        BAYESIAN,
        ADAPTIVE_SOFT,
        CENTRAL,
        SHAH
    } Denoising_type;

    /** Wavelet type.
     *
     * Valid types DAUB4, DAUB12, DAUB20
     */
    std::string DWT_type;

    /** Denoising type.
     *
     * Valid types are REMOVE_SCALE, SOFT_THRESHOLDING, BAYESIAN, ADAPTIVE_SOFT,
     * CENTRAL, SHAH.
     */
    Denoising_type denoising_type;

    /** Scale to which the denoising is applied.
     *
     * It is used by remove scale, adaptive soft
     */
    int scale;

    /** Output scale.
     *
     * All wavelet coefficients up to this output scale will be removed. Thus,
     * the output image is reduced by a factor that is related to this output
     * scale. For instance, if the output scale is 0, then no reduction is done.
     * If the output scale is 1, then the output image is reduced by 1/2, if the
     * output scale is 2, then it is reduced by 1/4, ...
     */
    int output_scale;

    /** Threshold for soft thresholding.
     */
    double threshold;

    /** Radius for central
     */
    int R;

    /** Smallest SNR
     */
    double SNR0;

    /** Largest SNR
     */
    double SNRF;

    /** White noise
     */
    bool white_noise;

    /** Shah number of outer iterations
     */
    int Shah_outer;

    /** Shah number of inner iterations
     */
    int Shah_inner;

    /** Shah number of refinement iterations
     */
    int Shah_refinement;

    /** Shah weight.
     *
     * w0 = data matching (=0)
     * w1 = 1st derivative smooth (=50)
     * w2 = edge strength (=50)
     * w3 = edge smoothness (=0.02)
     */
    Matrix1D< double > Shah_weight;

    /** Produce Shah edge instead of Shah smooth.
     */
    bool Shah_edge;

    /** Adjust range in Shah or Bayesian denoising.
     */
    bool adjust_range;

    /** Debugging level.
     *
     * 0, no dubug
     */
    int tell;

    /** Don't denoise.
     *
     * This is used in the Bayesian method, where the estimatedS parameters can
     * be averaged.
     */
    bool dont_denoise;

public:
    /** EstimatedS of the Bayesian method.
     */
    Matrix1D< double > estimatedS;

    /** Empty constructor
     */
    Denoising_parameters();

    /** Read parameters from command line
     */
    void read(int argc, char** argv);

    /** Produce side info.
     *
     * The DWT type is translated and set
     */
    void produce_side_info();

    /** Show parameters. This function calls show_specific.
     */
    void show();

    /** Show specific
     */
    void show_specific();

    /** Usage. This function calls usage_specific.
     */
    void usage();

    /** Show specific parameters.
     */
    void usage_specific();

    /** Denoise an image.
     */
    void denoise(Matrix2D< double >& img);

    /** Denoise a volume.
     */
    void denoise(Matrix3D< double >& vol);

    /** Denoise a volume using a precalculated estimate of the bayesian
     * parameters.
     */
    void denoise_avg_bayesian(Matrix3D< double >& vol);
};

#endif
