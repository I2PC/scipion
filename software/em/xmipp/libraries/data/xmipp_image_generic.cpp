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
#include "xmipp_image_extension.h"
#include "xmipp_error.h"
#include "xmipp_fftw.h"
#include "xmipp_types.h"
#include "xvsmooth.h"


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
    datatype = DT_Unknown;
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
    if (img.datatype != DT_Unknown)
    {
        setDatatype(img.datatype);
#define COPY(type) (*(Image<type>*)image) = (*(Image<type>*)img.image);

        SWITCHDATATYPE(datatype, COPY);
#undef COPY

    }

}

void ImageGeneric::getDimensions(size_t & xdim, size_t & ydim, size_t & Zdim) const
{
    size_t Ndim;
    image->getDimensions(xdim, ydim, Zdim, Ndim);
}

void ImageGeneric::getDimensions(size_t & xdim, size_t & ydim, size_t & Zdim, size_t & Ndim) const
{
    image->getDimensions(xdim, ydim, Zdim, Ndim);
}

void ImageGeneric::getDimensions(ArrayDim &aDim) const
{
    image->getDimensions(aDim);
}

void ImageGeneric::getInfo(ImageInfo &imgInfo) const
{
    image->getInfo(imgInfo);
}

void ImageGeneric::setDatatype(DataType imgType)
{
    if (imgType == datatype)
        return;

    clear();
    datatype = imgType;
    switch (datatype)
    {
    case DT_Float:
        {
            Image<float> *imT = new Image<float>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case DT_Double:
        {
            Image<double> *imT = new Image<double>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case DT_UInt:
        {
            Image<unsigned int> *imT = new Image<unsigned int>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case DT_Int:
        {
            Image<int> *imT = new Image<int>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case DT_UShort:
        {
            Image<unsigned short> *imT = new Image<unsigned short>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case DT_Short:
        {
            Image<short> *imT = new Image<short>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case DT_UChar:
        {
            Image<unsigned char> *imT = new Image<unsigned char>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case DT_SChar:
        {
            Image<char> *imT = new Image<char>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case DT_ULong:
        {
            Image<unsigned long> *imT = new Image<unsigned long>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case DT_Long:
        {
            Image<long> *imT = new Image<long>;
            image = imT;
            data = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), datatype);
        }
        break;
    case DT_Unknown:
        REPORT_ERROR(ERR_IMG_UNKNOWN,"Datatype of the image file is unknown.");
        break;
    default:
        REPORT_ERROR(ERR_NOT_IMPLEMENTED, "Datatype not implemented.");
        break;
    }
}

void ImageGeneric::setDatatype(const FileName &name, ImageInfo &imgInf)
{
    getImageInfo(name, imgInf);
    checkImageFileSize(name, imgInf, true);
    setDatatype(imgInf.datatype);

}

#define SET_DATATYPE(name) ImageInfo imgInf; setDatatype(name, imgInf);

void ImageGeneric::setDatatype(const FileName &name)
{
    SET_DATATYPE(name);
}

int ImageGeneric::read(const FileName &name, DataMode datamode, size_t select_img,
                       bool mapData)
{
    SET_DATATYPE(name);
    return image->read(name, datamode, select_img, mapData && !imgInf.swap);
}

int ImageGeneric::readMapped(const FileName &name, size_t select_img, int mode)
{
    SET_DATATYPE(name);
    return image->read(name, DATA, select_img, !imgInf.swap, mode);
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

double getScale(ImageInfo imgInf, size_t &xdim, size_t &ydim)
{
    double scale = 0;
    // If only xdim is passed, it is the higher allowable size, for any dimension
    if (ydim == 0 && imgInf.adim.xdim < imgInf.adim.ydim)
    {
        ydim = xdim;
        scale = ((double) ydim)/((double) imgInf.adim.ydim);
        xdim = (int)(scale*imgInf.adim.xdim);
    }
    else
    {
        scale = ((double) xdim)/((double) imgInf.adim.xdim);
        if (ydim == 0)
            ydim = (int)(scale*imgInf.adim.ydim);
    }
    return scale;
}

int ImageGeneric::readPreview(const FileName &name, size_t xdim, size_t ydim, int select_slice, size_t select_img)
{
    SET_DATATYPE(name);
    double scale = getScale(imgInf, xdim, ydim);
    //std::cerr << "DEBUG_JM: readPreview, scale: " << scale << std::endl;

    if (scale < 0.1)
      return readPreviewSmooth(name, xdim, ydim, select_slice, select_img);
    return image->readPreview(name, xdim, ydim, select_slice, select_img);
}

int ImageGeneric::readOrReadPreview(const FileName &name, size_t xdim, size_t ydim, int select_slice, size_t select_img,
                                    bool mapData, bool wrap)
{
    SET_DATATYPE(name);
    double scale = getScale(imgInf, xdim, ydim);
    //std::cerr << "DEBUG_JM: readOrReadPreview, scale: " << scale << std::endl;

    if (scale < 0.1)
      return readPreviewSmooth(name, xdim, ydim, select_slice, select_img);
    return image->readOrReadPreview(name, xdim, ydim, select_slice, select_img, !imgInf.swap && mapData);
}

void ImageGeneric::getPreview(ImageGeneric &imgOut, int xdim, int ydim, int select_slice, size_t select_img)
{
    imgOut.setDatatype(getDatatype());
    image->getPreview(imgOut.image, xdim, ydim, select_slice, select_img);
}


int ImageGeneric::readPreviewFourier(const FileName &name, size_t xdim, size_t ydim, int select_slice, size_t select_img)
{
    ImageGeneric ig;
    int result = ig.read(name, DATA, select_img);
    ImageInfo ii;
    ig.getInfo(ii);
    setDatatype(ii.datatype);
    switch (select_slice)
    {
    case CENTRAL_SLICE:
    case ALL_SLICES:
        select_slice = (int)((float)ii.adim.zdim / 2.0);
        break;
    default:
        select_slice--;
        break;
    }

    getScale(ii, xdim, ydim);

    ig().getSlice(select_slice, data);


    data->setXmippOrigin();
    selfScaleToSizeFourier((int)ydim, (int)xdim, *data, 1);
    return result;
}

int ImageGeneric::readPreviewSmooth(const FileName &name, size_t xdim, size_t ydim, int select_slice, size_t select_img)
{
  //std::cerr << "DEBUG_JM: readPreviewSmooth" << std::endl;
    ImageGeneric ig;
    int result = ig.read(name, DATA, select_img);
    ig.convert2Datatype(DT_UChar);
    ImageInfo ii;
    ig.getInfo(ii);

    getScale(ii, xdim, ydim);

    setDatatype(DT_UChar);
    resize(xdim, ydim, 1, 1, false);

    byte *inputArray = (byte*) MULTIDIM_ARRAY_GENERIC(ig).getArrayPointer();
    byte *outputArray = (byte*) MULTIDIM_ARRAY_GENERIC(*this).getArrayPointer();

    SmoothResize(inputArray, outputArray, ii.adim.xdim, ii.adim.ydim, xdim, ydim);
    return result;
}


void  ImageGeneric::mapFile2Write(int xdim, int ydim, int Zdim, const FileName &_filename,
                                  bool createTempFile, size_t select_img, bool isStack,int mode, int swapWrite)
{
    image->setDataMode(HEADER); // Use this to ask rw* which datatype to use
    if (swapWrite > 0)
        image->swapOnWrite();
    image->mapFile2Write(xdim,ydim,Zdim,_filename,createTempFile, select_img, isStack, mode);

    DataType writeDT = image->datatype();

    if ( writeDT != datatype)
    {
        setDatatype(writeDT);

        image->mapFile2Write(xdim,ydim,Zdim,_filename,createTempFile, select_img, isStack, mode);
    }
}

int ImageGeneric::readApplyGeo(const FileName &name, const MDRow &row,
                               const ApplyGeoParams &params)
{
    setDatatype(name);
    return image->readApplyGeo(name, row, params);
}

int ImageGeneric::readApplyGeo(const FileName &name, const MetaData &md, size_t objId,
                               const ApplyGeoParams &params)
{
    setDatatype(name);
    return image->readApplyGeo(name, md, objId, params);
}

void ImageGeneric::mirrorY()
{
    return image->mirrorY();
}
/** Read an image from metadata, filename is taken from MDL_IMAGE */
int ImageGeneric::readApplyGeo(const MetaData &md, size_t objId,
                               const ApplyGeoParams &params)
{
    FileName name;
    md.getValue(MDL_IMAGE, name, objId/*md.firstObject()*/);
    setDatatype(name);
    return image->readApplyGeo(name, md, objId, params);
}

/** Apply geometry in referring metadata to the image */
void ImageGeneric::applyGeo(const MetaData &md, size_t objId,
                            const ApplyGeoParams &params)
{
    image->applyGeo(md, objId, params);
}


void ImageGeneric::convert2Datatype(DataType _datatype, CastWriteMode castMode)
{
    if (_datatype == datatype || _datatype == DT_Unknown)
        return;

    ArrayDim aDim;
    data->getDimensions(aDim);

    ImageBase * newImage;
    MultidimArrayGeneric * newMAG;

#define CONVERTTYPE(type) Image<type> *imT = new Image<type>; \
        newImage = imT;\
        newMAG = new MultidimArrayGeneric((MultidimArrayBase*) &(imT->data), _datatype);\
        MultidimArray<type>* pMAG;\
        newMAG->getMultidimArrayPointer(pMAG);\
        if (castMode == CW_CAST)\
         data->getImage(*pMAG);\
        else\
        {\
         pMAG->resize(aDim);\
         double min, max;\
         data->computeDoubleMinMax(min, max);\
         ((Image<double>*) image)->getCastConvertPageFromT(0, (char*)pMAG->data, _datatype, aDim.nzyxdim, min, max, castMode);\
        }\

    SWITCHDATATYPE(_datatype, CONVERTTYPE)

#undef CONVERTTYPE

    /* aDimFile must be set in order to movePointerTo can be used.
     * If movePointerTo has been used before convert2Datatype, then
     * only the pointed image/slice is converted, and the images/slices left
     * are lost. This is why we set the new dimensions to the new ImageGeneric Object*/
    newImage->setADimFile(aDim);

    clear();
    datatype = _datatype;
    image = newImage;
    data = newMAG;
}


void ImageGeneric::reslice(AxisView face, ImageGeneric &imgOut)
{
    ArrayDim aDim, aDimOut;
    data->getDimensions(aDim);

    char axis;
    bool reverse;

    aDimOut = aDim;

    if (face == VIEW_Y_NEG || face == VIEW_Y_POS)
    {
        axis = 'Y';
        aDimOut.ydim = aDim.zdim;
        aDimOut.zdim = aDim.ydim;
        reverse = (face == VIEW_Y_NEG);
    }
    else if (face == VIEW_X_NEG || face == VIEW_X_POS)
    {
        axis = 'X';
        aDimOut.xdim = aDim.zdim;
        aDimOut.zdim = aDim.xdim;
        reverse = (face == VIEW_X_NEG);
    }
    else
      REPORT_ERROR(ERR_VALUE_INCORRECT, formatString("reslice: not supported view %d", (int)face));

    DataType dtype = getDatatype();
    imgOut.setDatatype(dtype);
    imgOut().resize(aDimOut, false);
    imgOut.image->setADimFile(aDimOut);

    MultidimArrayGeneric imTemp;

    int index;

    for (size_t k = 0; k < aDimOut.zdim; k++)
    {
        imTemp.aliasSlice(MULTIDIM_ARRAY_GENERIC(imgOut), k);
        index = k + (aDimOut.zdim - 1 - 2*k) * (int)reverse;
        MULTIDIM_ARRAY_GENERIC(*this).getSlice(index, &imTemp, axis, !reverse);
    }

}

void ImageGeneric::reslice(AxisView face)
{
    ImageGeneric imgOut;
    reslice(face, imgOut);
    clear();
    datatype = imgOut.getDatatype();
    image = imgOut.image;
    data = imgOut.data;
    imgOut.image = NULL;
    imgOut.data = NULL;
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
    String s;
    toString(s);
    std::cout << s <<std::endl;
}

void ImageGeneric::toString(String &s) const
{
    std::stringstream ss;
    if (image == NULL)
        ss << "Xmipp::ImageGeneric: Uninitialized image";
    else
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

void ImageGeneric::multiply(const ImageGeneric &img)
{

#define MULTIPLYIMG(type) MultidimArray<type> & kk = *((MultidimArray<type>*) data->im);\
                     MultidimArray<type> & pp = *((MultidimArray<type>*) img.data->im);\
                     kk *= pp;

    SWITCHDATATYPE(datatype, MULTIPLYIMG);
#undef MULTIPLYIMG

}

void ImageGeneric::multiply(const double value)
{

#define MULTIPLY(type) MultidimArray<type> & kk = *((MultidimArray<type>*) data->im);\
                     kk *= value;

    SWITCHDATATYPE(datatype, MULTIPLY);
#undef MULTIPLY

}

void ImageGeneric::divide(const double value)
{

#define DIVIDE(type) MultidimArray<type> & kk = *((MultidimArray<type>*) data->im);\
                     kk /= value;

    SWITCHDATATYPE(datatype, DIVIDE);
#undef DIVIDE

}

void createEmptyFile(const FileName &filename, int xdim, int ydim, int Zdim,
                     size_t select_img, bool isStack, int mode, int _swapWrite)
{
    ImageGeneric image;
    size_t found = filename.find_first_of("%");
    String strType = "";

    if (found == String::npos)
        image.setDatatype(DT_Float);
    else
    {
        strType = filename.substr(found+1).c_str();
        image.setDatatype(str2Datatype(strType));
    }

    image.mapFile2Write(xdim, ydim, Zdim, filename, false, select_img, isStack, mode, _swapWrite);
}
