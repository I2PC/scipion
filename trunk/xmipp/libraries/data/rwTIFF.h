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

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <iostream>

///@defgroup TIFF TIFF File format
///@ingroup ImageFormats

/*
//
// castTiffTile2T
//
// write the content of a tile from a TIFF file to an image array
//
*/
void castTiffTile2T(
    T * ptrDest ,
    unsigned char* tif_buf,
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

/*
//
// castTiffLine2T
//
// write the content of a line from a TIFF file to an image array
//
*/
void castTiffLine2T(
    T * ptrDest,
    unsigned char* tif_buf,
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


// I/O prototypes
/** TIFF Reader
  * @ingroup TIFF
*/
int readTIFF(int img_select, bool isStack=false)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readTIFF: Reading TIFF file\n");
#endif

    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "r") ) == NULL )
        return(-1);


    TIFF*           tif     = NULL;
    unsigned char*  tif_buf = NULL;
    int             raw_fd;
    unsigned char*  raw_buf = NULL;

    unsigned short  bitsPerSample;
    unsigned short  samplesPerPixel;
    unsigned int   imageWidth;
    unsigned int   imageLength;
    unsigned short  imageSampleFormat;
    unsigned short  imageBitsPerSample;
    unsigned int    tileWidth;
    unsigned int    tileLength;
    unsigned short  resUnit;
    float            xTiffRes,yTiffRes;
    unsigned int subFileType;
    uint16 pNumber, pTotal;
    uint32 rowsperstrip;
    tsize_t scanline;

    unsigned int    x, y;

    /* Open TIFF image */
    if ((tif = TIFFOpen(filename.c_str(), "r")) == NULL)
        REPORT_ERROR(1501,"rwTIFF: There is a problem opening the TIFF file.");

    do
    {
        std::cout << TIFFCurrentDirectory(tif) <<std::endl;

        /* Get TIFF image properties */
        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE,  &bitsPerSample);
        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL,&samplesPerPixel);
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,     &imageWidth);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH,    &imageLength);
        TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT,   &imageSampleFormat);
        TIFFGetField(tif, TIFFTAG_RESOLUTIONUNIT, &resUnit);
        TIFFGetField(tif, TIFFTAG_XRESOLUTION,    &xTiffRes);
        TIFFGetField(tif, TIFFTAG_YRESOLUTION,    &yTiffRes);
        TIFFGetField(tif, TIFFTAG_SUBFILETYPE, &subFileType);
        TIFFGetField(tif, TIFFTAG_PAGENUMBER, &pNumber, &pTotal);


        swap = TIFFIsByteSwapped(tif);
    }
    while(TIFFReadDirectory(tif));



    // Calculate x,y space dimension resolution
    double xRes, yRes;

    switch (resUnit)
    {
    case RESUNIT_NONE:
        {
            xRes = yRes = zeroD;
            break;
        }
    case RESUNIT_INCH:
        {
            xRes = 2.54e8/xTiffRes;
            yRes = 2.54e8/yTiffRes;
            break;
        }
    case RESUNIT_CENTIMETER:
        {
            xRes = 1e8/xTiffRes;
            yRes = 1e8/yTiffRes;
            break;
        }
    }


    // Set datatype
    DataType datatype;

    switch (bitsPerSample)
    {
    case 8:
        datatype = UChar;
        break;
    case 16:
        if (imageSampleFormat == SAMPLEFORMAT_INT)
            datatype = Short;
        else if (imageSampleFormat == SAMPLEFORMAT_UINT )
            datatype = UShort;
        else if (imageSampleFormat == 0)
            datatype = UShort;
        break;
    default:
        REPORT_ERROR(1000,"rwTIFF: Unsupported TIFF sample format.");
        break;
    }


    // cast image date to image class datatypes
    int _xDim,_yDim,_zDim;
    unsigned long int _nDim;

    _xDim = (int) imageWidth;
    _yDim = (int) imageLength;
    _zDim = (int) 1;
    _nDim = (int) 1;

    // Map the parameters
    data.setDimensions(_xDim, _yDim, 1, _nDim);

    //Set main header
    MDMainHeader.removeObjects();
    MDMainHeader.setColumnFormat(false);
    MDMainHeader.addObject();
    MDMainHeader.setValue(MDL_SAMPLINGRATEX, xRes);
    MDMainHeader.setValue(MDL_SAMPLINGRATEY, yRes);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    if( dataflag < 0 )
    {
        fclose(fimg);
        return 0;
    }


    // If samplesPerPixel is higher than 3 it means there are extra samples, as associated alpha data
    // Greyscale images are usually samplesPerPixel=1
    // RGB images are usually samplesPerPixel=3 (this is only implemented for untiled 8-bit tiffs)
    if (samplesPerPixel > 3)
        samplesPerPixel = 1;

    if (TIFFIsTiled(tif))
    {
        TIFFGetField(tif, TIFFTAG_TILEWIDTH,       &tileWidth);
        TIFFGetField(tif, TIFFTAG_TILELENGTH,      &tileLength);
        tif_buf = (unsigned char*)_TIFFmalloc(TIFFTileSize(tif));
        if (samplesPerPixel != 1)
        {
            std::cerr<<"rwTIFF ERROR: samplePerPixel is not 1: not yet implemented for RGB images";
            exit(1);
        }
    }
    else
    {
        TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip);
        scanline = TIFFScanlineSize(tif);
        tif_buf = (unsigned char*)_TIFFmalloc(scanline);
    }
    if (tif_buf == 0)
    {
        TIFFError(TIFFFileName(tif), "No space for strip buffer");
        exit(-1);
    }


    // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
    //if memory already allocated use it (no resize allowed)
    data.coreAllocateReuse();

    /* Start to convert the TIFF image to type T */

    if (TIFFIsTiled(tif))
    {
        for (y = 0; y < imageLength; y += tileLength)
            for (x = 0; x < imageWidth; x += tileWidth)
            {
                TIFFReadTile(tif, tif_buf, x, y, 0, 0);
                if (swap)
                    swapPage((char*)tif_buf, TIFFTileSize(tif)*sizeof(unsigned char), datatype);

                castTiffTile2T(data.data, tif_buf, x, y,
                               imageWidth, imageLength,
                               tileWidth, tileLength,
                               samplesPerPixel,
                               datatype);
            }
    }
    else
    {
        for (y = 0; y < imageLength; y++)
        {
            TIFFReadScanline(tif, tif_buf, y);
            castTiffLine2T(data.data, tif_buf, y,
                           imageWidth, imageLength,
                           samplesPerPixel,
                           datatype);
        }
    }
    _TIFFfree(tif_buf);

    TIFFClose(tif);

    //     return 0;


    unsigned long   imgStart=0;
    unsigned long   imgEnd =_nDim;
    if (img_select != -1)
    {
        imgStart=img_select;
        imgEnd=img_select+1;
    }

    //    MD.removeObjects();
    //    for ( i=imgStart; i<imgEnd; i++ )
    //for(int i=0;i< Ndim;i++)
    {
        MD.addObject();

        MD.setValue(MDL_ORIGINX, zeroD);
        MD.setValue(MDL_ORIGINY, zeroD);
        MD.setValue(MDL_ORIGINZ,  zeroD);
        MD.setValue(MDL_ANGLEROT, zeroD);
        MD.setValue(MDL_ANGLETILT,zeroD);
        MD.setValue(MDL_ANGLEPSI, zeroD);
        MD.setValue(MDL_WEIGHT,   oneD);
        MD.setValue(MDL_FLIP,     falseb);
    }

    //
    //


    //
    //    DataType datatype;
    //    //    dataHeaders[0].isSigned = false;
    //    int TIA_DT;
    //    if (img_select==-1)
    //    {
    //        TIA_DT = dataHeaders[0].DATA_TYPE;
    //        offset = header->pDATA_OFFSET[0] + TIAdataSIZE;
    //    }
    //    else
    //    {
    //        TIA_DT = dataHeaders[img_select].DATA_TYPE;
    //        offset = header->pDATA_OFFSET[img_select] + TIAdataSIZE;
    //    }
    //
    //
    //    switch ( TIA_DT )
    //    {
    //    case 1:
    //        datatype = UChar;
    //        break;
    //    case 2:
    //        datatype = UShort;
    //        //        datatype = Short;
    //        break;
    //    case 3:
    //        datatype = UInt;
    //        break;
    //    case 4:
    //        datatype = SChar;
    //        break;
    //    case 5:
    //        datatype = Short;
    //        //        dataHeaders[0].isSigned = true;
    //        break;
    //    case 6:
    //        datatype = Int;
    //        break;
    //    case 7:
    //        datatype = Float;
    //        break;
    //    case 8:
    //        datatype = Double;
    //        break;
    //    case 9:
    //        datatype = ComplexFloat;
    //        break;
    //    case 10:
    //        datatype = ComplexDouble;
    //        break;
    //    default:
    //        datatype = Unknown_Type;
    //        break;
    //    }
    //
    //    MDMainHeader.removeObjects();
    //    MDMainHeader.setColumnFormat(false);
    //    MDMainHeader.addObject();
    //    MDMainHeader.setValue(MDL_SAMPLINGRATEX,(double)dataHeaders[0].PIXEL_WIDTH);
    //    MDMainHeader.setValue(MDL_SAMPLINGRATEY,(double)dataHeaders[0].PIXEL_HEIGHT);
    //    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);
    //
    //    if( dataflag == -2 )
    //    {
    //        fclose(fimg);
    //        return 0;
    //    }
    //

    //
    //    //#define DEBUG
    //#ifdef DEBUG
    //
    //    MDMainHeader.write(std::cerr);
    //    MD.write(std::cerr);
    //#endif
    //
    //    delete header;
    //    size_t pad = TIAdataSIZE;
    //
    //    readData(fimg, img_select, datatype, pad);
    //
    //
    //    if ( !mmapOn )
    //        fclose(fimg);

    return 0;
}



/** TIFF Writer
  * @ingroup TIFF
*/
int writeTIFF(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE)
{
    REPORT_ERROR(6001, "ERROR: writeTIFF is not implemented.");
    return(-1);
}
#endif /* RWTIA_H_ */
