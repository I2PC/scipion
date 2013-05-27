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

#include "xmipp_image_extension.h"
#include "xmipp_error.h"


void getImageSize(const FileName &filename, size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim)
{
    Image<char> img;
    img.read(filename, HEADER);
    img.getDimensions(Xdim, Ydim, Zdim, Ndim);
}

void getImageInfo(const FileName &filename, size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim, DataType &datatype)
{
    Image<char> img;
    img.read(filename, HEADER);
    img.getDimensions(Xdim, Ydim, Zdim, Ndim);
    datatype = img.datatype();
}

void getImageInfo(const FileName &name, ImageInfo &imgInfo)
{
    Image<char> image;
    image.getInfo(name, imgInfo);
}

/* Return the datatype of the image file.
 */
void getImageDatatype(const FileName &name, DataType &datatype)
{
    Image<char> image;
    image.read(name, HEADER);
    datatype = image.datatype();
}
DataType getImageDatatype(const FileName &name)
{
    Image<char> image;
    image.read(name, HEADER);
    return image.datatype();
}

bool isImage(const FileName &name)
{
    Image<char> I;
    return I.isImage(name);
}

bool checkImageFileSize(const FileName &name, const ImageInfo &imgInfo, bool error)
{
    FileName ext = name.getExtension();
    FileName dataFname;

    if ( ext.contains("hed"))
        dataFname = name.removeLastExtension().addExtension("img");
    else if (ext.contains("inf"))
        dataFname = name.removeLastExtension();
    else if (ext.contains("tif") || ext.contains("jpg"))
        return true;
    else
        dataFname = name;

    size_t expectedSize = imgInfo.adim.nzyxdim*gettypesize(imgInfo.datatype) + imgInfo.offset;
    size_t actualSize = dataFname.removeBlockNameOrSliceNumber().removeFileFormat().getFileSize();
    bool result = (actualSize >= expectedSize);

    if (error && !result)
        REPORT_ERROR(ERR_IO_SIZE, formatString("Image Extension: File %s has wrong size.\n"
                                               "Expected size (at least) %u bytes. Actual size %u bytes.", name.c_str(), expectedSize, actualSize));

    return result;
}

bool checkImageFileSize(const FileName &name, bool error)
{
    ImageInfo imgInfo;
    getImageInfo(name, imgInfo);

    return checkImageFileSize(name, imgInfo, error);
}

bool checkImageCorners(const FileName &name)
{
    const int windowSize=11;
    const int windowSize_2=windowSize/2;
    const double Flower=0.4871; // MATLAB: N=11*11-1; finv(0.00005,N,N)
    // const double Fupper=2.0530; // MATLAB: N=11*11-1; finv(0.99995,N,N)

    ImageGeneric I;
    I.readMapped(name);
    size_t Xdim, Ydim, Zdim;
    I.getDimensions(Xdim,Ydim,Zdim);
    if (Zdim>1)
        return true;
    if (Xdim<=2*windowSize || Ydim<=2*windowSize)
        return true;

    MultidimArray<double> window;
    size_t i=Ydim/2;
    size_t j=Xdim/2;
    I().window(window,0,0,i-windowSize_2,j-windowSize_2,0,0,i+windowSize_2,j+windowSize_2);
    double stddev0=window.computeStddev();
    double var0=stddev0*stddev0;

    i=windowSize_2;
    j=windowSize_2;
    I().window(window,0,0,i-windowSize_2,j-windowSize_2,0,0,i+windowSize_2,j+windowSize_2);
    double stddev1=window.computeStddev();
    double var1=stddev1*stddev1;
    double F=var1/var0;
    if (F<Flower)//|| F>Fupper)
        return false;

    i=Ydim-1-windowSize_2;
    j=windowSize_2;
    I().window(window,0,0,i-windowSize_2,j-windowSize_2,0,0,i+windowSize_2,j+windowSize_2);
    stddev1=window.computeStddev();
    var1=stddev1*stddev1;
    F=var1/var0;
    if (F<Flower)//|| F>Fupper)
        return false;

    i=windowSize_2;
    j=Xdim-1-windowSize_2;
    I().window(window,0,0,i-windowSize_2,j-windowSize_2,0,0,i+windowSize_2,j+windowSize_2);
    stddev1=window.computeStddev();
    var1=stddev1*stddev1;
    F=var1/var0;
    if (F<Flower)//|| F>Fupper)
        return false;

    i=Ydim-1-windowSize_2;
    j=Xdim-1-windowSize_2;
    I().window(window,0,0,i-windowSize_2,j-windowSize_2,0,0,i+windowSize_2,j+windowSize_2);
    stddev1=window.computeStddev();
    var1=stddev1*stddev1;
    F=var1/var0;
    if (F<Flower)//|| F>Fupper)
        return false;

    return true;
}
