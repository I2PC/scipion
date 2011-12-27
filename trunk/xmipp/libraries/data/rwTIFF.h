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

#ifndef RWTIFF_H_
#define RWTIFF_H_

///@defgroup TIFF TIFF File format
///@ingroup ImageFormats
//@{

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
void castTiffTile2T(
    T * ptrDest ,
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
            castPage2T((char*) tif_buf+((j-y)*samplesPerPixel*typeSize*tileWidth+(i-x)*samplesPerPixel*typeSize), ptrDest+(j*imageWidth + i), datatype, (size_t) 1);
}

/** castTiffLine2T
 *
 * write the content of a line from a TIFF file to an image array
*/
void castTiffLine2T(
    T * ptrDest,
    char* tif_buf,
    unsigned int y,
    unsigned int imageWidth, unsigned int imageLength,
    unsigned short samplesPerPixel,
    DataType datatype)
{
    unsigned int x;
    int typeSize = gettypesize(datatype);

    for (x = 0; x < imageWidth; x++)
        castPage2T((char*) tif_buf+(samplesPerPixel*typeSize * x), ptrDest+(y*imageWidth + x), datatype, (size_t) 1);
}

/** Determine datatype of the TIFF format file.
  * @ingroup TIFF
*/
DataType datatypeTIFF(TIFFDirHead dHead)
{
    DataType datatype;

    switch (dHead.bitsPerSample)
    {
    case 8:
        if (dHead.imageSampleFormat == SAMPLEFORMAT_INT)
            datatype = SChar;
        else
            datatype = UChar;
        break;
    case 16:
        if (dHead.imageSampleFormat == SAMPLEFORMAT_INT)
            datatype = Short;
        else
            datatype = UShort;

        //        else if (dHead.imageSampleFormat == SAMPLEFORMAT_UINT ||
        //                 dHead.imageSampleFormat == SAMPLEFORMAT_IEEEFP ) //Don't know why
        //            datatype = UShort;
        //        //        else if (dHead.imageSampleFormat == 0     ||
        //        //                 dHead.imageSampleFormat == 32767 ) // Format 0 and 32767 are not declared in TIFF 6.0 specifications ¿?
        //        else
        //        datatype = UShort;
        break;
    case 32:
        if (dHead.imageSampleFormat == SAMPLEFORMAT_INT)
            datatype = Int;
        else if (dHead.imageSampleFormat == SAMPLEFORMAT_UINT )
            datatype = UInt;
        else if (dHead.imageSampleFormat == SAMPLEFORMAT_IEEEFP )
            datatype = Float;
        else
            datatype = Unknown_Type;
        break;
    default:
        datatype = Unknown_Type;
        //        REPORT_ERROR(ERR_TYPE_INCORRECT,"rwTIFF: Unsupported TIFF sample format.");
        break;
    }
    return datatype;
}

/**
 *  Read TIFF format files.
*/
int readTIFF(size_t select_img, bool isStack=false)
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

    uint32 rowsperstrip;
    tsize_t scanline;

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
                REPORT_ERROR(ERR_IMG_NOREAD, formatString("readTIFF: %s file contains %lu images with, at least,"\
                             " two of them with different dimensions. Try to read them individually.",filename.c_str(), dirHead.size()));
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
    int _xDim,_yDim,_zDim;
    size_t _nDim;

    size_t   imgStart = IMG_INDEX(select_img);

    _xDim = (int) dirHead[imgStart].imageWidth;
    _yDim = (int) dirHead[imgStart].imageLength;
    _zDim = 1;
    _nDim = (select_img == ALL_IMAGES)? dirHead.size() : 1;
    replaceNsize = _nDim;
    setDimensions(_xDim, _yDim, 1, _nDim);

    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : _nDim;

    DataType datatype = datatypeTIFF(dirHead[0]);

    //Set main header
    MDMainHeader.clear();
    MDMainHeader.setValue(MDL_SAMPLINGRATEX, xRes);
    MDMainHeader.setValue(MDL_SAMPLINGRATEY, yRes);
    MDMainHeader.setValue(MDL_DATATYPE,(int) datatype);

    //Read header only
    if( dataMode < DATA ||( dataMode == _HEADER_ALL && _nDim > 1) )
        return 0;

    /* As we cannot mmap a TIFF File, when this option is passed we are going to mmap
     * the multidimarray of Image
     */
    if (mmapOnRead)
    {
        mmapOnRead = false;
        if (data.nzyxdim*gettypesize(datatype) > tiff_map_min_size)
            data.setMmap(true);
    }

    // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
    //if memory already allocated use it (no resize allowed)
    data.coreAllocateReuse();

    size_t pad = _xDim * _yDim;
    int imReaded = 0;

    MD.clear();
    MD.resize(_nDim,MDL::emptyHeader);

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

                    castTiffTile2T(data.data+(pad*imReaded), tif_buf, x, y,
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
                castTiffLine2T(data.data+(pad*imReaded), tif_buf, y,
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
int writeTIFF(size_t select_img, bool isStack=false, int mode=WRITE_OVERWRITE, String bitDepth="", CastWriteMode castMode = CW_CAST)
{
#undef DEBUG

    //    if (mode == WRITE_REPLACE)
    //        REPORT_ERROR(ERR_TYPE_INCORRECT,"rwTIFF: Images cannot be replaced in TIFF file.");
    if (typeid(T) == typeid(std::complex<double>))
    {
        REPORT_ERROR(ERR_TYPE_INCORRECT,"rwTIFF: Complex images are not supported by TIFF format.");
        return 0;
    }

    int Xdim, Ydim, Zdim;
    size_t Ndim;
    getDimensions(Xdim, Ydim, Zdim, Ndim);

    // Volumes are not supported
    if (Zdim > 1)
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
        case Double:
        case Float:
            wDType = Float;
            dhMain.imageSampleFormat = SAMPLEFORMAT_IEEEFP;
            break;
        case Int:
            castMode = CW_CONVERT;
        case UInt:
            wDType = UInt;
            dhMain.imageSampleFormat = SAMPLEFORMAT_UINT;
            break;
        case Short:
            castMode = CW_CONVERT;
        case UShort:
            wDType = UShort;
            dhMain.imageSampleFormat = SAMPLEFORMAT_UINT;
            break;
        case SChar:
            castMode = CW_CONVERT;
        case UChar:
            wDType = UChar;
            dhMain.imageSampleFormat = SAMPLEFORMAT_UINT;
            break;
        default:
            wDType = Unknown_Type;
            REPORT_ERROR(ERR_TYPE_INCORRECT,(std::string)"ERROR: IMAGIC format does not write from " \
                         + typeid(T).name() + "type.");
        }
    }
    else
    {
        // Default Value
        wDType = (bitDepth == "default") ? UChar : datatypeRAW(bitDepth);

        switch(wDType)
        {
        case Float:
            dhMain.imageSampleFormat = SAMPLEFORMAT_IEEEFP;
            break;
        case UInt:
        case UShort:
        case UChar:
            dhMain.imageSampleFormat = SAMPLEFORMAT_UINT;
            break;
        default:
            wDType = Unknown_Type;
            REPORT_ERROR(ERR_TYPE_INCORRECT,"ERROR: incorrect TIFF bits depth value.");
        }
    }

    if (mmapOnWrite)
    {
        MDMainHeader.setValue(MDL_DATATYPE,(int) wDType);
        if (!checkMmapT(wDType))
        {
            if (dataMode < DATA ) // This means ImageGeneric wants to know which DataType must use in mapFile2Write
                return 0;
            else
                REPORT_ERROR(ERR_MMAP, "File datatype and image declaration not compatible with mmap.");
        }
        else
            dataMode = DATA;

        /* As we cannot mmap a TIFF File, when this option is passed we are going to mmap
         * the multidimarray of Image
         */
        mmapOnWrite = false;
        if (data.nzyxdim*gettypesize(wDType) > tiff_map_min_size)
            data.setMmap(true);

        // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
        //if memory already allocated use it (no resize allowed)

        data.coreAllocateReuse();

        return 0;
    }

    int nBytes = gettypesize(wDType);

    /* Set TIFF image properties */
    dhMain.bitsPerSample = (unsigned short int) nBytes*8;
    dhMain.samplesPerPixel = 1;
    dhMain.imageWidth  = Xdim;
    dhMain.imageLength = Ydim;
    dhMain.resUnit = RESUNIT_CENTIMETER;

    double aux;

    if (!MDMainHeader.empty())
    {
        dhMain.xTiffRes = (MDMainHeader.getValue(MDL_SAMPLINGRATEX, aux)) ? (float) 1e8/aux : 0. ;
        dhMain.yTiffRes = (MDMainHeader.getValue(MDL_SAMPLINGRATEX, aux)) ? (float) 1e8/aux : 0. ;
    }

    //    size_t imgStart = 0;
    //    size_t imgEnd = Ndim;

    size_t imgStart = (mode == WRITE_APPEND)? replaceNsize : IMG_INDEX(select_img);



    size_t bufferSize, datasize_n;
    bufferSize = Xdim*nBytes;
    datasize_n = Xdim*Ydim*Zdim;

    double min0, max0;

    char*  tif_buf;

    if ((tif_buf = (char*)_TIFFmalloc(bufferSize)) == 0)
    {
        TIFFError(TIFFFileName(tif), "No space for strip buffer");
        exit(-1);
    }

    //    if (mode == WRITE_REPLACE)
    //        tiffsetdirectory(tif,select_img);

    //Write each image in a directory
    for (size_t i = 0; i < Ndim; i++ )
    {
        TIFFSetDirectory(tif,(tdir_t) i + imgStart);
        //        TIFFSetDirectory(tif, select_img);


        //        if (mode == WRITE_REPLACE)
        //          TIFFReadDirectory(tif);

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

        if (Ndim == 1 && isStack == false)
        {
            TIFFSetField(tif, TIFFTAG_SUBFILETYPE, (unsigned int) 0x0);
            TIFFSetField(tif, TIFFTAG_PAGENUMBER, (uint16) 0, (uint16) 0);
        }
        else
        {
            TIFFSetField(tif, TIFFTAG_SUBFILETYPE, (unsigned int) 0x2);
            TIFFSetField(tif, TIFFTAG_PAGENUMBER, (uint16) i, (uint16) Ndim);
        }

        data.computeDoubleMinMaxRange(min0, max0, i*datasize_n, datasize_n);

        for (uint32 y = 0; y < Ydim; y++)
        {
            if (castMode == CW_CAST)
                castPage2Datatype(MULTIDIM_ARRAY(data)+i*datasize_n + y*Xdim, (char *)tif_buf, wDType, (size_t) Xdim);
            else
                castConvertPage2Datatype(MULTIDIM_ARRAY(data)+i*datasize_n + y*Xdim,
                                         (char *)tif_buf, wDType, (size_t) Xdim,min0 ,max0, castMode);

            TIFFWriteScanline(tif, tif_buf,y,0);
        }

        if (Ndim >1)
            TIFFWriteDirectory(tif);
    }

    _TIFFfree(tif_buf);
    return(0);
}
//@}
#endif /* RWTIFF_H_ */
