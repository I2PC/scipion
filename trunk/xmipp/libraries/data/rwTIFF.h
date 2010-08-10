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

///@definegroup TIFF TIFF File format
///@ingroup ImageFormats

// I/O prototypes
/** TIFF Reader
  * @ingroup TIFF
*/
int readTIFF(int img_select, bool isStack=false)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readTIA: Reading TIA file\n");
#endif

    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "r") ) == NULL )
        return(-1);


    TIFF*           tif     = NULL;
    unsigned char*  tif_buf = NULL;
    int             raw_fd;
    unsigned char*  raw_buf = NULL;
    FILE*           inf     = NULL;

    unsigned short  bitsPerSample;
    unsigned short  samplesPerPixel;
    unsigned int   imageWidth;
    unsigned int   imageLength;
    unsigned short  imageSampleFormat;
    unsigned short  imageBitsPerSample;
    unsigned int    tileWidth;
    unsigned int    tileLength;
    int             byte_swapped;
    size_t          len;
    uint32 rowsperstrip;
    tsize_t scanline;

    unsigned int    x, y;
    char     zsMinval[64], zsMaxval[64];
    zsMinval[0] = '\0';
    zsMaxval[0] = '\0';

    /* Open TIFF image */
    if ((tif = TIFFOpen(filename.c_str(), "r")) == NULL)
        REPORT_ERROR(1501,"rwTIFF: There is a problem opening the TIFF file.");

    /* Get TIFF image properties */
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE,   &bitsPerSample);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,      &imageWidth);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH,     &imageLength);
    TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT,    &imageSampleFormat);

    swap = TIFFIsByteSwapped(tif);


    // Sjors 6nov07: some scanners set samplesPerPixel to a very high value
    // If this happens, set to 1 and write a warning message...
    // Greyscale images are usually samplesPerPixel=1
    // RGB images are usually samplesPerPixel=3 (this is only implemented for untiled 8-bit tiffs)
    if (samplesPerPixel != 1 && samplesPerPixel != 3)
    {
        std::cerr <<"WARNING! This tif has a strange value for samplesPerPixel (i.e. "<<samplesPerPixel<<")"<<std::endl;
        std::cerr <<"         Some scanners do not set this value correctly"<<std::endl;
        std::cerr <<"         Setting samplePerPixel to 1 and continuing execution ... "<<std::endl;
        samplesPerPixel = 1;
    }

    if (TIFFIsTiled(tif))
    {
        TIFFGetField(tif, TIFFTAG_TILEWIDTH,       &tileWidth);
        TIFFGetField(tif, TIFFTAG_TILELENGTH,      &tileLength);
        tif_buf = (unsigned char*)_TIFFmalloc(TIFFTileSize(tif));
        if (samplesPerPixel != 1)
        {
            std::cerr<<"ERROR: samplePerPixel is not 1: not implemented for tiled images";
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
















    /* Calculates the TIFF image size */
    len = (imageLength * imageWidth * bitsPerSample * samplesPerPixel) / 8;
    initFile(argv[2], len);

    /* Mapping from file to memory (it will save primary memory) */
    raw_fd  = open(argv[2], O_RDWR);
    raw_buf = (unsigned char*)mmap(NULL, len, PROT_READ | PROT_WRITE, MAP_SHARED,
                                   raw_fd, 0);
    /* Writes the TIFF image properties to the .inf file */
    inf     = fopen(strcat(argv[2], ".inf"), "w");
    fprintf(inf, "# Bits per sample\n");
    fprintf(inf, "bitspersample= %d\n", bitsPerSample);
    fprintf(inf, "# Samples per pixel\n");
    fprintf(inf, "samplesperpixel= %d\n", samplesPerPixel);
    fprintf(inf, "# Image width\n");
    fprintf(inf, "Xdim= %d\n", imageWidth);
    fprintf(inf, "# Image length\n");
    fprintf(inf, "Ydim= %d\n", imageLength);

    /* Start to convert the TIFF image to the RAW file */

    if (TIFFIsTiled(tif))
    {
        for (y = 0; y < imageLength; y += tileLength)
            for (x = 0; x < imageWidth; x += tileWidth)
            {
                TIFFReadTile(tif, tif_buf, x, y, 0, 0);
                rawWriteTile(raw_buf, tif_buf, x, y,
                             imageWidth, imageLength,
                             tileWidth, tileLength,
                             bitsPerSample, imageSampleFormat,
                             byte_swapped,
                             zsMinval, zsMaxval);
            }
    }
    else
    {
        for (y = 0; y < imageLength; y++)
        {
            TIFFReadScanline(tif, tif_buf, y);
            rawWriteLine(raw_buf, tif_buf, y,
                         imageWidth, imageLength,
                         bitsPerSample, samplesPerPixel, imageSampleFormat,
                         byte_swapped,
                         zsMinval, zsMaxval);
        }
    }
    _TIFFfree(tif_buf);
    munmap((char *) raw_buf, len);

    TIFFClose(tif);
    close(raw_fd);

    fprintf(inf, "# Minimum value=%s\n", zsMinval);
    fprintf(inf, "# Maximum value=%s\n", zsMaxval);
    if (imageSampleFormat == SAMPLEFORMAT_INT)
        fprintf(inf, "# Signed short?\nis_signed=TRUE\n");
    else
        fprintf(inf, "# Signed short?\nis_signed=FALSE\n");

    fclose(inf);
    return 0;






    int _xDim,_yDim,_zDim;
    unsigned long int _nDim;

    _xDim = (int) imageWidth;
    _yDim = (int) imageLength;
    _zDim = (int) 1;
    _nDim = (int) 1;

    // Map the parameters
    data.setDimensions(_xDim, _yDim, 1, _nDim);

    unsigned long   imgStart=0;
    unsigned long   imgEnd =_nDim;
    if (img_select != -1)
    {
        imgStart=img_select;
        imgEnd=img_select+1;
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
    //    MD.removeObjects();
    //    for ( i=imgStart; i<imgEnd; i++ )
    //        //for(int i=0;i< Ndim;i++)
    //    {
    //        MD.addObject();
    //        double aux;
    //        if(MDMainHeader.getValue(MDL_SAMPLINGRATEX,aux))
    //        {
    //            aux = ROUND(dataHeaders[i].CalibrationElementX - \
    //                        dataHeaders[i].CalibrationOffsetX/aux - data.xdim/2);
    //            MD.setValue(MDL_ORIGINX, aux);
    //        }
    //        if(MDMainHeader.getValue(MDL_SAMPLINGRATEY,aux))
    //        {
    //            aux = ROUND(dataHeaders[i].CalibrationElementY - \
    //                        dataHeaders[i].CalibrationOffsetY/aux -data.ydim/2);
    //            MD.setValue(MDL_ORIGINY, aux);
    //        }
    //        MD.setValue(MDL_ORIGINZ,  zeroD);
    //
    //        MD.setValue(MDL_ANGLEROT, zeroD);
    //        MD.setValue(MDL_ANGLETILT,zeroD);
    //        MD.setValue(MDL_ANGLEPSI, zeroD);
    //        MD.setValue(MDL_WEIGHT,   oneD);
    //        MD.setValue(MDL_FLIP,     falseb);
    //    }
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
