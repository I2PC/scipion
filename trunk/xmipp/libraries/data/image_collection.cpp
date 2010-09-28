/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#include "image_collection.h"

ImageCollection::ImageCollection(const MetaData &md)
{
  copyMetadata(md);
}

ImageCollection::ImageCollection(const FileName &fnImage)
{
  Image<double> image;
  image.read(fnImage, false);
  if (image().ndim == 1)
    addImage(fnImage);
  else
  {
    FileName fnTemp;
    for (int i = 0; i < image().ndim; ++i)
    {
      fnTemp = integerToString(i, 5);
      fnTemp += "@" + fnImage;
    }
  }

}

void ImageCollection::addImage(const FileName &fnImage)
{

}
