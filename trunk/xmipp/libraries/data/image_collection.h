/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosam@cnb.csic.es)
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

#ifndef IMAGE_COLLECTION_H_
#define IMAGE_COLLECTION_H_

#include "metadata.h"
#include "image.h"

/** Class to represent a collection of images.
 * This class is an abstraction layer to handle(iterate, read, write) a group
 * of images. An image collection could be:
 * - A single image, the simple case.
 * - Metadata containing images(not in stack).
 * - Images on a stack.
 * - Images on a Metadata from one(or several) stacks.
 * - Slices of a volume.
 */

class ImageCollection: public MetaData
{
private:
  FileName lastInStack, lastOutStack;

  ///Open stack to start reading images
  void openStack(const FileName & fnStack);
  ///Close stack
  void closeStack(const FileName & fnStack);
  ///Add image to metadata
  void addImage(const FileName &fnImage);

public:
  ///Constructors and destructor
  ImageCollection(const MetaData &md);
  ImageCollection(const FileName &fnImage);
  ~ImageCollection();

};//end of class ImageCollection

#endif /* IMAGE_COLLECTION_H_ */
