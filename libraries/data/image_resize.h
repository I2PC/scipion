/***************************************************************************
 * Authors:     Joaquin Oton    (joton@cnb.csic.es)
 *              J.M. De la Rosa (jmdelarosa@cnb.csic.es)
 *
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

#ifndef IMAGE_RESIZE_H
#define IMAGE_RESIZE_H

#include "transformations.h"
#include "xmipp_image_generic.h"
#include "metadata_extension.h"
#include "xmipp_fftw.h"
#include "xmipp_program.h"
#include "matrix2d.h"


typedef enum { RESIZE_NONE, RESIZE_FACTOR, RESIZE_FOURIER, RESIZE_PYRAMID_EXPAND, RESIZE_PYRAMID_REDUCE } ScaleType;
#define INTERP_FOURIER -1

class ProgImageResize: public XmippMetadataProgram
{
public:
  /** Constructor and destructor, just to avoid vtable undefined references errors */
  ProgImageResize();
  ~ProgImageResize();

protected:
    ScaleType scale_type;

    int             splineDegree, dim, pyramid_level, fourier_threads;
    bool            isVol, apply_geo, temporaryOutput;
    //Matrix2D<double> R, T, S, A, B;
    Matrix1D<double>          shiftV, rotV, scaleV;
    ImageGeneric img, imgOut;

    void defineParams();
    void readParams();
    void preProcess();
    void postProcess();
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

};
#endif //IMAGE_RESIZE_H
