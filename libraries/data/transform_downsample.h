/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#ifndef _DOWNSAMPLE
#define _DOWNSAMPLE

#include "xmipp_program.h"

///@defgroup MicrographDownsample Micrograph Downsample
/// @ingroup ReconsLibrary
//@{
/** Downsample parameters. */
class ProgTransformDownsample: public XmippMetadataProgram
{
public:
    /// Downsampling factor
    double step;

    typedef enum {KER_RECTANGLE, FOURIER, SMOOTH} TDownsamplingMethod;

    /** Downsampling method*/
    TDownsamplingMethod method;

    /// Number of Threads used in the Fourier Transform
    int nThreads;

public:
    // Kernel
    MultidimArray<double> kernel;
public:
    /** Read input parameters.*/
    void readParams();

    /// Define params
    void defineParams();

    /** Really downsample.*/
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);
};

/** Downsample a micrograph using a rectangular kernel.
 * The input and output micrographs must be already open.
 */
void downsampleKernel(const ImageGeneric &M, double step, ImageGeneric &Mp);

/** Downsample a micrograph in Fourier space.
 * The input and output micrographs must be already open.
 */
void downsampleFourier(const ImageGeneric &M, double step, ImageGeneric &Mp, int nThreads);

/** Downsample a micrograph using smooth and color dithering.
 * The input and output micrographs must be already open.
 */
void downsampleSmooth(const ImageGeneric &M, ImageGeneric &Mp);
//@}
#endif
