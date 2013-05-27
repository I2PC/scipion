/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
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

#ifndef XMIPP_IMAGE_EXTENSION_H_
#define XMIPP_IMAGE_EXTENSION_H_

#include "xmipp_image.h"

/** @defgroup ImageExtension Extension to Image class
 *  @ingroup Images
 *  @{
 */

/** What is the size of an image */
void getImageSize(const FileName &filename, size_t &Xdim, size_t &Ydim, size_t &Zdim,
                  size_t &Ndim);
/** Get basic information from image file */
void getImageInfo(const FileName &filename, size_t &Xdim, size_t &Ydim, size_t &Zdim,
                  size_t &Ndim, DataType &datatype);
void getImageInfo(const FileName &name, ImageInfo &imgInfo);

/** Get datatype information from image file*/
void getImageDatatype(const FileName &name, DataType &datatype);
DataType getImageDatatype(const FileName &name);

/** A filename is an image? */
bool isImage(const FileName &name);

/** Check if the file is at least as large as needed.
 *  If error = true, then report an error in case of check fail. */
bool checkImageFileSize(const FileName &name, const ImageInfo &imgInfo, bool error = false);
bool checkImageFileSize(const FileName &name, bool error = false);

/** Check if image corners have the same variance as the central part of the image.
 * The test of equality of variances F-Snedecor is applied with 99.99% of confidence
 */
bool checkImageCorners(const FileName &name);
//@}

#endif /* XMIPP_IMAGE_EXTENSION_H_ */
