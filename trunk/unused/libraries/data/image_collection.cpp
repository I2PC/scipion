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

#include "xmipp_image_collection.h"

ImageCollection::ImageCollection(int mode)
{
    this->mode = mode;
}

ImageCollection::ImageCollection(const MetaData &md, int mode)
{
    copyMetadata(md);
    this->mode = mode;
}

ImageCollection::ImageCollection(const FileName &fnImage, int mode)
{
    Image<char> image;
    image.read(fnImage, false);
    if (image().ndim == 1)
    {
        addObject();
        setValue(MDL_IMAGE, fnImage);
        setValue(MDL_ENABLED, 1);
    }
    else
    {
        FileName fnTemp;
        for (size_t i = 0; i < image().ndim; ++i)
        {
            fnTemp.compose(i, fnImage);
            addObject();
            setValue(MDL_IMAGE, fnImage);
            setValue(MDL_ENABLED, 1);
        }
    }
    this->mode = mode;
}

ImageCollection::~ImageCollection()
{
    std::map<FileName, ImageFHandler*>::iterator it;
    Image<char> image;
    for (it = openedStacks.begin(); it != openedStacks.end(); ++it)
        image.closeFile(it->second);
}

ImageFHandler* ImageCollection::getStackHandle(const FileName & fnStack)
{
    std::map<FileName, ImageFHandler*>::iterator it;
    it = openedStacks.find(fnStack);
    if (it != openedStacks.end())
        return it->second;
    Image<char> image;
    return (openedStacks[fnStack] = image.openFile(fnStack, mode));
}

/** This is a wrap of Image::read */
int ImageCollection::readImage(Image<double> &image, const FileName &name, bool readdata, int select_img,
                               bool apply_geo, bool only_apply_shifts, MDRow * row, bool mapData)
{
    if (name.isInStack())
    {
        FileName stackName;
        int imgno;
        name.decompose(imgno, stackName);
        ImageFHandler * fIH = getStackHandle(stackName);
        image._read(stackName, fIH, readdata, imgno, apply_geo, only_apply_shifts, row, mapData);
    }
    else
        image.read(name, readdata, select_img, apply_geo, only_apply_shifts, row, mapData);
}

/** This is a wrap of Image::write */
void ImageCollection::writeImage(Image<double> &image, const FileName &name, int select_img, bool isStack)
{
    if (name.isInStack() || select_img != -1)
    {
        FileName stackName;
        int imgno;
        name.decompose(imgno, stackName);
        if (select_img == -1)
            select_img = imgno;
        ImageFHandler * fIH = getStackHandle(stackName);
        image._write(stackName, fIH, select_img, true, mode);
    }
    else
        image.write(name, select_img, isStack, mode);
}

