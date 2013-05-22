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

#include "xmipp_image_base.h"


/** TIFF Data Header
*/
struct TIFFDirHead
{                                   // Header for each Directory in TIFF
    unsigned short  bitsPerSample;
    unsigned short  samplesPerPixel;
    unsigned int   imageWidth;
    unsigned int   imageLength;
    uint16           imageSampleFormat;
    unsigned short  resUnit;
    float            xTiffRes,yTiffRes;
    unsigned int subFileType;
    uint16 pNumber, pTotal; // pagenumber and total number of pages of current directory
};

/**
 * castTiffTile2T
 *
 * write the content of a tile from a TIFF file to an image array
 */
void ImageBase::castTiffTile2T(
    size_t offset ,
    char* tif_buf,
    unsigned int x, unsigned int y,
    unsigned int imageWidth, unsigned int imageLength,
    unsigned int tileWidth, unsigned int tileLength,
    unsigned short samplesPerPixel,
    DataType datatype)
{
    int typeSize = gettypesize(datatype);
    unsigned int i, j;
    unsigned int x_max = x + tileWidth,
                         y_max = y + tileLength;

    if (x_max > imageWidth)
        x_max = imageWidth;
    if (y_max > imageLength)
        y_max = imageLength;


    for (j = y; j < y_max; j++)
        for (i = x; i < x_max; i++)
            setPage2T(offset+(j*imageWidth + i), (char*) tif_buf+((j-y)*samplesPerPixel*typeSize*tileWidth+(i-x)*samplesPerPixel*typeSize), datatype, (size_t) 1);
}

/** castTiffLine2T
 *
 * write the content of a line from a TIFF file to an image array
*/
void ImageBase::castTiffLine2T(
    size_t offset,
    char* tif_buf,
    unsigned int y,
    unsigned int imageWidth, unsigned int imageLength,
    unsigned short samplesPerPixel,
    DataType datatype)
{
    unsigned int x;
    int typeSize = gettypesize(datatype);

    for (x = 0; x < imageWidth; x++)
        setPage2T(offset+(y*imageWidth + x), (char*) tif_buf+(samplesPerPixel*typeSize * x), datatype, (size_t) 1);
}

/** Determine datatype of the TIFF format file.
  * @ingroup TIFF
*/
DataType ImageBase::datatypeTIFF(TIFFDirHead dHead)
{
    DataType datatype;

    switch (dHead.bitsPerSample)
    {
    case 8:
        if (dHead.imageSampleFormat == SAMPLEFORMAT_INT)
            datatype = DT_SChar;
        else
            datatype = DT_UChar;
        break;
    case 16:
        if (dHead.imageSampleFormat == SAMPLEFORMAT_INT)
            datatype = DT_Short;
        else
            datatype = DT_UShort;

        //        else if (dHead.imageSampleFormat == SAMPLEFORMAT_UINT ||
        //                 dHead.imageSampleFormat == SAMPLEFORMAT_IEEEFP ) //Don't know why
        //            datatype = DT_UShort;
        //        //        else if (dHead.imageSampleFormat == 0     ||
        //        //                 dHead.imageSampleFormat == 32767 ) // Format 0 and 32767 are not declared in TIFF 6.0 specifications ¿?
        //        else
        //        datatype = DT_UShort;
        break;
    case 32:
        if (dHead.imageSampleFormat == SAMPLEFORMAT_INT)
            datatype = DT_Int;
        else if (dHead.imageSampleFormat == SAMPLEFORMAT_UINT )
            datatype = DT_UInt;
        else if (dHead.imageSampleFormat == SAMPLEFORMAT_IEEEFP )
            datatype = DT_Float;
        else
            datatype = DT_Unknown;
        break;
    default:
        datatype = DT_Unknown;
        //        REPORT_ERROR(ERR_TYPE_INCORRECT,"rwTIFF: Unsupported TIFF sample format.");
        break;
    }
    return datatype;
}

/**
 *  Read TIFF format files.
*/
int ImageBase::readTIFF(size_t select_img, bool isStack)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readTIFF: Reading TIFF file\n");
#endif

    //    TIFFSetWarningHandler(NULL); // Switch off warning messages

    char*  tif_buf = NULL;
    std::vector<TIFFDirHead> dirHead;
    TIFFDirHead dhRef;

    /* Get TIFF image properties */
    do
    {
        dhRef.imageSampleFormat = SAMPLEFORMAT_VOID;
        if (TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE,  &dhRef.bitsPerSample) == 0)
            REPORT_ERROR(ERR_IO_NOREAD,"rwTIFF: Error reading TIFFTAG_BITSPERSAMPLE");
        if (TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL,&dhRef.samplesPerPixel) == 0)
            dhRef.samplesPerPixel = 1;

        if (TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,     &dhRef.imageWidth) == 0)
            REPORT_ERROR(ERR_IO_NOREAD,"rwTIFF: Error reading TIFFTAG_IMAGEWIDTH");
        if (TIFFGetField(tif, TIFFTAG_IMAGELENGTH,    &dhRef.imageLength) == 0)
            REPORT_ERROR(ERR_IO_NOREAD,"rwTIFF: Error reading TIFFTAG_IMAGELENGTH");
        if (TIFFGetField(tif, TIFFTAG_SUBFILETYPE,    &dhRef.subFileType) == 0)
            dhRef.subFileType = 0; // Some scanners does not provide this label. So, we set this to zero
        //            REPORT_ERROR(ERR_IO_NOREAD,"rwTIFF: Error reading TIFFTAG_SUBFILETYPE");
        TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT,   &dhRef.imageSampleFormat);
        TIFFGetField(tif, TIFFTAG_RESOLUTIONUNIT, &dhRef.resUnit);
        TIFFGetField(tif, TIFFTAG_XRESOLUTION,    &dhRef.xTiffRes);
        TIFFGetField(tif, TIFFTAG_YRESOLUTION,    &dhRef.yTiffRes);
        TIFFGetField(tif, TIFFTAG_PAGENUMBER,     &dhRef.pNumber, &dhRef.pTotal);

        if ((dhRef.subFileType & 0x00000001) != 0x00000001) //add image if not a thumbnail
            dirHead.push_back(dhRef);
    }
    while(TIFFReadDirectory(tif));

    swap = TIFFIsByteSwapped(tif);

    //Check select_img is lower than stack size
    if (select_img > dirHead.size())
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, formatString("readTIFF (%s): Image number %lu exceeds stack size %lu", filename.c_str(), select_img, dirHead.size()));
    else if (select_img == ALL_IMAGES)// Check images dimensions. Need to be the same
    {
        for (size_t i = 1; i < dirHead.size(); i++)
        {
            if (dirHead[0].imageLength != dirHead[i].imageLength || \
                dirHead[0].imageWidth != dirHead[i].imageWidth)
                dirHead.resize(i);
            //REPORT_ERROR(ERR_IMG_NOREAD, formatString("readTIFF: %s file contains %lu images with, at least,"\
            //             " two of them with different dimensions. Try to read them individually.",filename.c_str(), dirHead.size()));
        }
    }

    // Calculate x,y space dimension resolution
    double xRes, yRes;

    switch (dirHead[0].resUnit)
    {
    case RESUNIT_NONE:
        {
            xRes = yRes = -1;
            break;
        }
    case RESUNIT_INCH:
        {
            xRes = 2.54e8/dirHead[0].xTiffRes;
            yRes = 2.54e8/dirHead[0].yTiffRes;
            break;
        }
    case RESUNIT_CENTIMETER:
        {
            xRes = 1e8/dirHead[0].xTiffRes;
            yRes = 1e8/dirHead[0].yTiffRes;
            break;
        }
    }

    // cast image data to image class datatypes
    ArrayDim aDim;

    size_t   imgStart = IMG_INDEX(select_img);

    aDim.xdim = (int) dirHead[imgStart].imageWidth;
    aDim.ydim = (int) dirHead[imgStart].imageLength;
    aDim.zdim = 1;
    aDim.ndim = replaceNsize = (select_img == ALL_IMAGES)? dirHead.size() : 1;
    setDimensions(aDim);

    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : aDim.ndim;

    DataType datatype = datatypeTIFF(dirHead[0]);

    //Set main header
    MDMainHeader.clear();
    MDMainHeader.setValue(MDL_SAMPLINGRATE_X, xRes);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_Y, yRes);
    MDMainHeader.setValue(MDL_DATATYPE,(int) datatype);

    //Read header only
    if( dataMode < DATA ||( dataMode == _HEADER_ALL && aDim.ndim > 1) )
        return 0;

    /* As we cannot mmap a TIFF File, when this option is passed we are going to mmap
     * the multidimarray of Image
     */
    if (mmapOnRead)
    {
        mmapOnRead = false;
        if (aDim.nzyxdim*gettypesize(datatype) > tiff_map_min_size)
            mdaBase->setMmap(true);
    }

    // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
    //if memory already allocated use it (no resize allowed)
    mdaBase->coreAllocateReuse();

    size_t pad = aDim.yxdim;
    int imReaded = 0;

    MD.clear();
    MD.resize(aDim.ndim,MDL::emptyHeader);

    uint32 rowsperstrip;
    tsize_t scanline;

    unsigned int x, y;
    // Dimensions of tiles
    unsigned int tileWidth, tileLength;

    for (size_t i = imgStart; i < imgEnd; ++i)
    {
        TIFFSetDirectory(tif,(tdir_t) i);

        // If samplesPerPixel is higher than 3 it means there are extra samples, as associated alpha data
        // Greyscale images are usually samplesPerPixel=1
        // RGB images are usually samplesPerPixel=3 (this is only implemented for untiled 8-bit tiffs)
        if (dirHead[i].samplesPerPixel > 3)
            dirHead[i].samplesPerPixel = 1;

        if (TIFFIsTiled(tif))
        {
            TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth);
            TIFFGetField(tif, TIFFTAG_TILELENGTH,&tileLength);
            tif_buf = (char*)_TIFFmalloc(TIFFTileSize(tif));
        }
        else
        {
            TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);
            scanline = TIFFScanlineSize(tif);
            tif_buf = (char*)_TIFFmalloc(scanline);
        }
        if (tif_buf == 0)
        {
            TIFFError(TIFFFileName(tif), "No space for strip buffer");
            exit(-1);
        }

        /* Start to convert the TIFF image to type T */

        datatype = datatypeTIFF(dirHead[i]);

        if (TIFFIsTiled(tif))
        {
            for (y = 0; y < dirHead[0].imageLength; y += tileLength)
                for (x = 0; x < dirHead[0].imageWidth; x += tileWidth)
                {
                    TIFFReadTile(tif, tif_buf, x, y, 0, 0);
                    if (swap)
                        swapPage((char*)tif_buf, TIFFTileSize(tif)*sizeof(unsigned char), datatype);

                    castTiffTile2T((pad*imReaded), tif_buf, x, y,
                                   dirHead[i].imageWidth, dirHead[i].imageLength,
                                   tileWidth, tileLength,
                                   dirHead[i].samplesPerPixel,
                                   datatype);
                }
        }
        else
        {
            for (y = 0; y < dirHead[i].imageLength; y++)
            {
                TIFFReadScanline(tif, tif_buf, y);
                castTiffLine2T((pad*imReaded), tif_buf, y,
                               dirHead[i].imageWidth, dirHead[i].imageLength,
                               dirHead[i].samplesPerPixel,
                               datatype);
            }
        }

        ++imReaded;
        _TIFFfree(tif_buf);
    }
    return 0;
}

/**
 * Write TIFF format files.
*/
int ImageBase::writeTIFF(size_t select_img, bool isStack, int mode, String bitDepth, CastWriteMode castMode)
{
#undef DEBUG

    if (isComplexT())
    {
        REPORT_ERROR(ERR_TYPE_INCORRECT,"rwTIFF: Complex images are not supported by TIFF format.");
        return 0;
    }

    ArrayDim aDim;
    mdaBase->getDimensions(aDim);

    // Volumes are not supported
    if (aDim.zdim > 1)
        REPORT_ERROR(ERR_MULTIDIM_DIM, "writeTIFF does not support volumes.");
    // TIFF cannot open a file to read/write at the same time, so the file must be overwritten
    //    if (mode != WRITE_OVERWRITE)
    //        REPORT_ERROR(ERR_VALUE_INCORRECT, "writeTIFF: LIBTIFF cannot modify an existing file, only overwrite it.");

    TIFFDirHead dhMain; // Main header

    //Selection of output datatype

    DataType wDType,myTypeID = myT();

    if (bitDepth == "")
    {
        castMode = CW_CAST;
        switch(myTypeID)
        {
        case DT_Double:
        case DT_Float:
            wDType = DT_Float;
            dhMain.imageSampleFormat = SAMPLEFORMAT_IEEEFP;
            break;
        case DT_Int:
            castMode = CW_CONVERT;
        case DT_UInt:
            wDType = DT_UInt;
            dhMain.imageSampleFormat = SAMPLEFORMAT_UINT;
            break;
        case DT_Short:
            castMode = CW_CONVERT;
        case DT_UShort:
            wDType = DT_UShort;
            dhMain.imageSampleFormat = SAMPLEFORMAT_UINT;
            break;
        case DT_SChar:
            castMode = CW_CONVERT;
        case DT_UChar:
            wDType = DT_UChar;
            dhMain.imageSampleFormat = SAMPLEFORMAT_UINT;
            break;
        default:
            wDType = DT_Unknown;
            REPORT_ERROR(ERR_TYPE_INCORRECT,formatString("rwTIFF: cannot write TIFF format from %s\
                         datatype.",datatype2Str(myTypeID).c_str()));
        }
    }
    else
    {
        // Default Value
        wDType = (bitDepth == "default") ? DT_UChar : datatypeRAW(bitDepth);

        switch(wDType)
        {
        case DT_Float:
            dhMain.imageSampleFormat = SAMPLEFORMAT_IEEEFP;
            break;
        case DT_UInt:
        case DT_UShort:
        case DT_UChar:
            dhMain.imageSampleFormat = SAMPLEFORMAT_UINT;
            break;
        default:
            wDType = DT_Unknown;
            REPORT_ERROR(ERR_TYPE_INCORRECT,formatString("rwTIFF: TIFF format does not support %s " \
                         "datatype.",datatype2Str(myTypeID).c_str()));
        }
    }

    if (mmapOnWrite)
    {
        /* As we cannot mmap a TIFF File, when this option is passed we are going to mmap
         * the multidimarray of Image
         */
        mmapOnWrite = false;
        dataMode = DATA;
        MDMainHeader.setValue(MDL_DATATYPE,(int) myTypeID);

        if (mdaBase->nzyxdim*gettypesize(myTypeID) > tiff_map_min_size)
            mdaBase->setMmap(true);

        // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
        //if memory already allocated use it (no resize allowed)

        mdaBase->coreAllocateReuse();

        return 0;
    }

    int nBytes = gettypesize(wDType);

    /* Set TIFF image properties */
    dhMain.bitsPerSample = (unsigned short int) nBytes*8;
    dhMain.samplesPerPixel = 1;
    dhMain.imageWidth  = aDim.xdim;
    dhMain.imageLength = aDim.ydim;
    dhMain.resUnit = RESUNIT_CENTIMETER;

    double aux;

    if (!MDMainHeader.empty())
    {
        dhMain.xTiffRes = (MDMainHeader.getValue(MDL_SAMPLINGRATE_X, aux)) ? (float) 1e8/aux : 0. ;
        dhMain.yTiffRes = (MDMainHeader.getValue(MDL_SAMPLINGRATE_X, aux)) ? (float) 1e8/aux : 0. ;
    }

    size_t imgStart = (mode == WRITE_APPEND)? replaceNsize : IMG_INDEX(select_img);

    size_t bufferSize, datasize_n;
    bufferSize = aDim.xdim*nBytes;
    datasize_n = aDim.zyxdim;

    char*  tif_buf;

    if ((tif_buf = (char*)_TIFFmalloc(bufferSize)) == 0)
    {
        TIFFError(TIFFFileName(tif), "No space for strip buffer");
        exit(-1);
    }

    //Write each image in a directory
    for (size_t i = 0; i < aDim.ndim; i++ )
    {
        TIFFSetDirectory(tif,(tdir_t) i + imgStart);

        // Image header
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE,  dhMain.bitsPerSample);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL,dhMain.samplesPerPixel);
        TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT,   dhMain.imageSampleFormat);
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH,     dhMain.imageWidth);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH,    dhMain.imageLength);
        TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, dhMain.resUnit);
        TIFFSetField(tif, TIFFTAG_XRESOLUTION,    dhMain.xTiffRes);
        TIFFSetField(tif, TIFFTAG_YRESOLUTION,    dhMain.yTiffRes);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC,    PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_COMPRESSION,    COMPRESSION_NONE);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP,   (uint32) dhMain.imageLength);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG,   PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_SOFTWARE,       "Xmipp 3.0");

        if (aDim.ndim == 1 && isStack == false)
        {
            TIFFSetField(tif, TIFFTAG_SUBFILETYPE, (unsigned int) 0x0);
            TIFFSetField(tif, TIFFTAG_PAGENUMBER, (uint16) 0, (uint16) 0);
        }
        else
        {
            TIFFSetField(tif, TIFFTAG_SUBFILETYPE, (unsigned int) 0x2);
            TIFFSetField(tif, TIFFTAG_PAGENUMBER, (uint16) i, (uint16) aDim.ndim);
        }

        double min0, max0;

        if (castMode != CW_CAST)
            mdaBase->computeDoubleMinMaxRange(min0, max0, i*datasize_n, datasize_n);

        for (uint32 y = 0; y < aDim.ydim; y++)
        {
            if (castMode == CW_CAST)
                getPageFromT(i*datasize_n + y*aDim.xdim, (char *)tif_buf, wDType, (size_t) aDim.xdim);
            else
                getCastConvertPageFromT(i*datasize_n + y*aDim.xdim,
                                        (char *)tif_buf, wDType, (size_t) aDim.xdim, min0, max0, castMode);

            TIFFWriteScanline(tif, tif_buf,y,0);
        }

        if (aDim.ndim >1)
            TIFFWriteDirectory(tif);
    }

    _TIFFfree(tif_buf);
    return(0);
}
//@}
