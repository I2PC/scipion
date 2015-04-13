/***************************************************************************
 *
 * Authors: J.M. de la Rosa Trevin   (jmdelarosa@cnb.csic.es) x 0.95
 *          Joaquin Oton             (joton@cnb.csic.es)      x 0.05
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

#ifndef IMAGE_OPERATE_H
#define IMAGE_OPERATE_H

#include "metadata.h"
#include "xmipp_image.h"
#include "xmipp_strings.h"
#include "xmipp_program.h"

/// @defgroup ImageOperate Mathematical operations of images and volumes
/// @ingroup DataLibrary
//@{
/**@name Image operation procedures
   */
//@{
/**This define the prototype of binary operations on images
   the result will be left in op1 */
typedef void ImageBinaryOperator(Image<double> &op1, const Image<double> &op2);

/**This define the prototype of unary operations on images
   the result will be left in op */
typedef void ImageUnaryOperator(Image<double> &op);

/** Operate program */
class ProgOperate: public XmippMetadataProgram
{
private:
    //Functions pointer to binary operation
    ImageBinaryOperator * binaryOperator;
    //Functions pointer to unary operation
    ImageUnaryOperator * unaryOperator;
    FileName fn2;
    MetaData md2;
    MDIterator md2Iterator;
    Image<double> img2;
    bool isValue;
    double value;
    FileName file_or_value;
protected:
    /// Define parameters
    void defineParams();

    /// Read input parameters
    void readParams();

    /// Process one image
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);
}
;
//@}
#endif
