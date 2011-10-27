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

#include "xmipp_image.h"
#include "xmipp_image_generic.h"
#include "xmipp_error.h"


ImageGeneric::ImageGeneric(DataType _datatype)
{
    init();
    setDatatype(_datatype);
}

ImageGeneric::ImageGeneric(const FileName &filename)
{
    init();
    read(filename);
}

ImageGeneric::ImageGeneric(const ImageGeneric &img)
{
    init();
    copy(img);
}

ImageGeneric::~ImageGeneric()
{
    delete image;
    delete data;
}

void ImageGeneric::init()
{
    image = NULL;
    data = NULL;
    datatype = Unknown_Type;
}

void ImageGeneric::clear()
{
    if (image != NULL)
    {
        image->clear();
        delete image;
        delete data;
        init();
    }
}

void  ImageGeneric::copy(const ImageGeneric &img)
{
    if (img.datatype != Unknown_Type)
    {
        setDatatype(img.datatype);
#define COPY(type) (*(Image<type>*)image) = (*(Image<type>*)img.image);

        SWITCHDATATYPE(datatype, COPY);
#undef COPY

    }

}

void ImageGeneric::getInfo(const FileName &name, ImageInfo &imgInfo)
{
    Image<char> image;
    image.getInfo(name, imgInfo);
}

void ImageGeneric::getInfo(ImageInfo &imgInfo) const
{
	image->getInfo(imgInfo);
}

void ImageGeneric::getImageType(const FileName &imgName, DataType &datatype)
{
    Image<char> Im;
    Im.read(imgName, HEADER);
    datatype = Im.datatype();
}

void ImageGeneric::getImageType(const FileName &imgName, DataType &datatype, bool &swap)
{
    Image<char> Im;
    Im.read(imgName, HEADER);
    datatype = Im.datatype();
    swap = (Im.getSwap() > 0);
}

void ImageGeneric::setDatatype(DataType imgType)
{
    clear();
    datatype = imgType;
    switch (datatype)
    {
    case Float:
        {
            Image<float> *imT = new Image<float>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case Double:
        {
            Image<double> *imT = new Image<double>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case UInt:
        {
            Image<unsigned int> *imT = new Image<unsigned int>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case Int:
        {
            Image<int> *imT = new Image<int>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case UShort:
        {
            Image<unsigned short> *imT = new Image<unsigned short>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case Short:
        {
            Image<short> *imT = new Image<short>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case UChar:
        {
            Image<unsigned char> *imT = new Image<unsigned char>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case SChar:
        {
            Image<char> *imT = new Image<char>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case Unknown_Type:
        REPORT_ERROR(ERR_IMG_UNKNOWN,"Datatype of the image file is unknown.");
        break;
    default:
        REPORT_ERROR(ERR_NOT_IMPLEMENTED, "Datatype not implemented.");
        break;
    }
}

int ImageGeneric::read(const FileName &name, DataMode datamode, size_t select_img,
                       bool mapData)
{
    bool swap;
    DataType datatype;
    getImageType(name, datatype, swap);
    setDatatype(datatype);
    return image->read(name, datamode, select_img, mapData && !swap);
}

int ImageGeneric::readMapped(const FileName &name, size_t select_img, int mode)
{
    bool swap;
    DataType datatype;
    getImageType(name, datatype, swap);
    setDatatype(datatype);

    return image->read(name, DATA, select_img, !swap, mode);
}

int ImageGeneric::readOrReadMapped(const FileName &name, size_t select_img, int mode)
{
    try
    {
        return read(name, DATA, select_img, false);
    }
    catch (XmippError &xe)
    {
        if (xe.__errno == ERR_MEM_NOTENOUGH)
        {
            reportWarning("ImageBase::readOrReadMapped: Not enough memory to allocate. \n"
                          " Proceeding to map image from file.");
            return readMapped(name, select_img, mode);
        }
        else
            throw xe;
    }
}

int ImageGeneric::readPreview(const FileName &name, int Xdim, int Ydim, int select_slice, size_t select_img)
{
    DataType datatype;
    getImageType(name, datatype);
    setDatatype(datatype);

    return image->readPreview(name, Xdim, Ydim, select_slice, select_img);
}

int ImageGeneric::readOrReadPreview(const FileName &name, int Xdim, int Ydim, int select_slice, size_t select_img, bool mapData)
{
    ImageInfo imInfo;
    this->getInfo(name, imInfo);
    getInfo(name, imInfo);
    setDatatype(imInfo.datatype);

    image->readOrReadPreview(name, Xdim, Ydim, select_slice, select_img, mapData);
}


void  ImageGeneric::mapFile2Write(int Xdim, int Ydim, int Zdim, const FileName &_filename,
                                  bool createTempFile, size_t select_img, bool isStack,int mode, int swapWrite)
{
    image->setDataMode(HEADER); // Use this to ask rw* which datatype to use
    if (swapWrite > 0)
        image->swapOnWrite();
    image->mapFile2Write(Xdim,Ydim,Zdim,_filename,createTempFile, select_img, isStack, mode);

    DataType writeDT = image->datatype();
    if ( writeDT != datatype)
    {
        setDatatype(writeDT);
        image->mapFile2Write(Xdim,Ydim,Zdim,_filename,createTempFile, select_img, isStack, mode);
    }
}

int ImageGeneric::readApplyGeo(const FileName &name, const MDRow &row, bool only_apply_shifts,
                               DataMode datamode, size_t select_img, bool wrap)
{
    DataType datatype;
    getImageType(name, datatype);
    setDatatype(datatype);
    return image->readApplyGeo(name, row, only_apply_shifts, datamode, select_img, wrap);
}

int ImageGeneric::readApplyGeo(const FileName &name, const MetaData &md, size_t objId, bool only_apply_shifts, DataMode datamode,
                               size_t select_img, bool wrap)
{
    DataType datatype;
    getImageType(name, datatype);
    setDatatype(datatype);
    return image->readApplyGeo(name, md, objId, only_apply_shifts, datamode, select_img, wrap);
}

/** Read an image from metadata, filename is taken from MDL_IMAGE */
int ImageGeneric::readApplyGeo(const MetaData &md, size_t objId, bool only_apply_shifts,
                               DataMode datamode, size_t select_img, bool wrap)
{
    FileName name;
    md.getValue(MDL_IMAGE, name, objId/*md.firstObject()*/);
    DataType datatype;
    getImageType(name, datatype);
    setDatatype(datatype);
    return image->readApplyGeo(name, md, objId, only_apply_shifts, datamode, select_img, wrap);
}

/** Apply geometry in refering metadata to the image */
void ImageGeneric::applyGeo(const MetaData &md, size_t objId, bool only_apply_shifts, bool wrap)
{
    image->applyGeo(md, objId, only_apply_shifts, wrap);
}

ImageGeneric& ImageGeneric::operator=(const ImageGeneric &img)
{
    copy(img);
    return *this;
}

bool ImageGeneric::operator==(const ImageGeneric &i1) const
{
    return(*(this->data) == *(i1.data));
}

void ImageGeneric::print() const
{
    std::cout << *image;
}

void ImageGeneric::toString(String &s) const
{
    std::stringstream ss;
    ss << *image;
    s = ss.str();
}

void ImageGeneric::add(const ImageGeneric &img)
{
    if (datatype != img.datatype)
        REPORT_ERROR(ERR_TYPE_INCORRECT, "Images have different datatypes");

#define ADD(type) MultidimArray<type> & kk = *((MultidimArray<type>*) data->im);\
                     MultidimArray<type> & pp = *((MultidimArray<type>*) img.data->im);\
                     kk += pp;

    SWITCHDATATYPE(datatype, ADD);
#undef ADD

}

void ImageGeneric::subtract(const ImageGeneric &img)
{
    if (datatype != img.datatype)
        REPORT_ERROR(ERR_TYPE_INCORRECT, "Images have different datatypes");

#define MINUS(type) MultidimArray<type> & kk = *((MultidimArray<type>*) data->im);\
                     MultidimArray<type> & pp = *((MultidimArray<type>*) img.data->im);\
                     kk -= pp;

    SWITCHDATATYPE(datatype, MINUS);
#undef MINUS

}


void createEmptyFile(const FileName &filename, int Xdim, int Ydim, int Zdim,
                     size_t select_img, bool isStack, int mode, int _swapWrite)
{
    ImageGeneric image;
    size_t found = filename.find_first_of("%");
    String strType = "";

    if (found == String::npos)
        image.setDatatype(Float);
    else
    {
        strType = filename.substr(found+1).c_str();
        image.setDatatype(datatypeString2Int(strType));
    }

    image.mapFile2Write(Xdim, Ydim, Zdim, filename, false, select_img, isStack, mode, _swapWrite);
}
