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

#ifndef TRANSFORMGEOMETRY_H
#define TRANSFORMGEOMETRY_H

#include "transformations.h"
#include "xmipp_image_generic.h"
#include "metadata_extension.h"
#include "xmipp_fftw.h"
#include "xmipp_program.h"
#include "matrix2d.h"


class ProgTransformGeometry: public XmippMetadataProgram
{
public:
    /** Constructor and destructor, just to avoid vtable undefined references errors */
    ProgTransformGeometry();
    ~ProgTransformGeometry();

protected:
    int             splineDegree, dim;
    bool            applyTransform, inverse, wrap, isVol, flip, mdVol;
    Matrix2D<double> R, A, B, T;
    ImageGeneric img, imgOut;
    String matrixStr; // To read directly the matrix

    void defineParams();
    void readParams();

    /** Calculate the rotation matrix according to the parameters
     */
    void calculateRotationMatrix();
    void preProcess();
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

};
#endif //TRANSFORMGEOMETRY_H
