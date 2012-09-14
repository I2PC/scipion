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

#ifndef DENOISE_H
#define DENOISE_H

#include <string>
#include <data/filters.h>

/// @defgroup Denoise Image denoising
/// @ingroup ReconsLibrary

/// Parameters for denoise program
/// @ingroup Denoise
class WaveletFilter: public XmippFilter
{
public:
    typedef enum
    {
        REMOVE_SCALE,
        SOFT_THRESHOLDING,
        BAYESIAN,
        ADAPTIVE_SOFT,
        CENTRAL
    } DenoisingType;

    /** Wavelet type.
     *
     * Valid types DAUB4, DAUB12, DAUB20
     */
    String DWT_type;

    /** Denoising type.
     *
     * Valid types are REMOVE_SCALE, SOFT_THRESHOLDING, BAYESIAN, ADAPTIVE_SOFT,
     * CENTRAL.
     */
    DenoisingType denoising_type;

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

    /** Adjust range in Bayesian denoising.
     */
    bool adjust_range;

    /** Verbosity used by denoiser.
     * By default will be 0, programs using denoiser should set this
     */
    int verbose;

    /** Don't denoise.
     *
     * This is used in the Bayesian method, where the estimatedS parameters can
     * be averaged.
     */
    bool dont_denoise;

    static void defineParams(XmippProgram *program);
    void readParams(XmippProgram *program);

public:
    /** EstimatedS of the Bayesian method.
     */
    Matrix1D< double > estimatedS;

    /** Empty constructor
     */
    WaveletFilter();

    /** Produce side info.
     *
     * The DWT type is translated and set
     */
    void produceSideInfo();

    void show();

    /** Denoise an image.
     */
    void apply(MultidimArray< double >& img);

    /** Denoise a volume using a precalculated estimate of the bayesian
     * parameters.
     */
    void denoiseAvgBayesian(MultidimArray< double >& vol);
};

#endif
