/*******************************************************************************/
/* tiff2raw.c                             08/08/2000 */
/*******************************************************************************/
/* TIFF to raw converter.                            */
/* The TIFF file should be in the tile oriented format                         */
/* It will make two files: the raw file and another one (filename.inf) with    */
/* the TIFF image size, the bits per sample and the samples per pixel          */
/*******************************************************************************/
/* Centro Nacional de Biotecnología                                            */
/* Carlos Manzanares Sancho                                                    */
/*******************************************************************************/
/*
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <iostream>

#include "tiffio.h"

/*
// initFile
//
// Creates a file with the specified size
// (to speed up the process it makes a 100000 bytes buffer)
//
*/
int initFile(char *filename, size_t len)
{
    unsigned char*  buffer = NULL;
    size_t          i;
    FILE *                                  fd;

    buffer = (unsigned char*)malloc(sizeof(unsigned char) * 100000);
    if (buffer == NULL) return(-1);
    for (int i = 0; i < 100000; i++) buffer[i] = 0;

    fd = fopen(filename, "w");
    if (fd == NULL)
    {
        fprintf(stderr, "Error opening file %s", filename);
        exit(1);
    }
    for (i = 0; i < len / 100000; i++)
        fwrite(buffer, sizeof(unsigned char), 100000, fd);
    fwrite(buffer, sizeof(unsigned char), len % 100000 , fd);

    fclose(fd);

    return(0);
}

/*
//
// rawWriteTile
//
// write the content of a tile from a TIFF file to an image array
//
*/
void rawWriteTile(
    unsigned char* raw_buf, unsigned char* tif_buf,
    unsigned int x, unsigned int y,
    unsigned int imageWidth, unsigned int imageLength,
    unsigned int tileWidth, unsigned int tileLength,
    unsigned short bitsPerSample, unsigned short imageSampleFormat,
    int byte_swapped,
    char * zsMinval, char * zsMaxval)
{
    unsigned int i, j;
    unsigned int x_max = x + tileWidth,
                         y_max = y + tileLength;
    static double minval = 10e10;
    static double maxval = -10e10;
    unsigned short int  uiVal;
    short int            iVal;
    unsigned char ucVal;

    unsigned char * aux_pointer;

    uiVal = 0;
    iVal = 0;
    ucVal = 0;


    if (x_max > imageWidth)  x_max = imageWidth;
    if (y_max > imageLength) y_max = imageLength;
    if (bitsPerSample == 16)
    {
        if (imageSampleFormat == SAMPLEFORMAT_INT) aux_pointer = (unsigned char *) & iVal;
        else aux_pointer = (unsigned char *) & uiVal;
    }

    for (j = y; j < y_max; j++)
        for (i = x; i < x_max; i++)
        {
            switch (bitsPerSample)
            {
            case  8:
                ucVal = tif_buf[(j-y)*tileWidth+(i-x)];
                if (((double)ucVal) < minval) minval = (double)ucVal;
                else if (((double)ucVal) > maxval) maxval = (double)ucVal;
                raw_buf[j*imageWidth + i] = ucVal;
                break;
            case 16:
                if (imageSampleFormat == SAMPLEFORMAT_INT)
                {

                    if (!byte_swapped)
                    {
                        aux_pointer[0] = tif_buf[((j-y)*tileWidth+(i-x))*2];
                        aux_pointer[1] = tif_buf[((j-y)*tileWidth+(i-x))*2+1];
                    }
                    else
                    {
                        aux_pointer[0] = tif_buf[((j-y)*tileWidth+(i-x))*2+1];
                        aux_pointer[1] = tif_buf[((j-y)*tileWidth+(i-x))*2];
                    }
                    if ((double)iVal < minval) minval = (double)iVal;
                    else if ((double)iVal > maxval) maxval = (double)iVal;
                }
                else
                {

                    if (!byte_swapped)
                    {
                        aux_pointer[0] = tif_buf[((j-y)*tileWidth+(i-x))*2];
                        aux_pointer[1] = tif_buf[((j-y)*tileWidth+(i-x))*2+1];
                    }
                    else
                    {
                        aux_pointer[0] = tif_buf[((j-y)*tileWidth+(i-x))*2+1];
                        aux_pointer[1] = tif_buf[((j-y)*tileWidth+(i-x))*2];
                    }
                    if ((double)uiVal < minval) minval = (double)uiVal;
                    else if ((double)uiVal > maxval) maxval = (double)uiVal;
                }
                raw_buf[(j*imageWidth + i)*2]   = aux_pointer[0];
                raw_buf[(j*imageWidth + i)*2+1] = aux_pointer[1];
                break;
            }//switch
        }//for
    sprintf(zsMinval, " %f", minval);
    sprintf(zsMaxval, " %f", maxval);
}

/*
//
// rawWriteLine
//
// write the content of a line from a TIFF file to an image array
//
*/
void rawWriteLine(
    unsigned char* raw_buf, unsigned char* tif_buf,
    unsigned int y,
    unsigned int imageWidth, unsigned int imageLength,
    unsigned short bitsPerSample, unsigned short samplesPerPixel, unsigned short imageSampleFormat,
    int byte_swapped,
    char * zsMinval, char * zsMaxval)
{
    unsigned int x;
    static double minval = 10e10;
    static double maxval = -10e10;
    unsigned short int  uiVal;
    short int            iVal;
    unsigned char ucVal;

    unsigned char * aux_pointer;

    uiVal = 0;
    iVal = 0;
    ucVal = 0;

    if (bitsPerSample == 16)
    {
        if (imageSampleFormat == SAMPLEFORMAT_INT) aux_pointer = (unsigned char *) & iVal;
        else aux_pointer = (unsigned char *) & uiVal;
        if (samplesPerPixel!=1)
        {
            std::cerr<<"ERROR: samplePerPixel is not 1: not implemented for 16-bit images";
            exit(1);
        }
    }
    for (x = 0; x < imageWidth; x++)
    {
        switch (bitsPerSample)
        {
        case  8:
            ucVal = tif_buf[samplesPerPixel * x];
            if (((double)ucVal) < minval) minval = (double)ucVal;
            else if (((double)ucVal) > maxval) maxval = (double)ucVal;
            raw_buf[y*imageWidth + x] = ucVal;
            break;
        case 16:
            if (imageSampleFormat == SAMPLEFORMAT_INT)
            {
                if (!byte_swapped)
                {
                    aux_pointer[0] = tif_buf[x*2];
                    aux_pointer[1] = tif_buf[x*2+1];
                }
                else
                {
                    aux_pointer[0] = tif_buf[x*2+1];
                    aux_pointer[1] = tif_buf[x*2];
                }
                if ((double)iVal < minval) minval = (double)iVal;
                else if ((double)iVal > maxval) maxval = (double)iVal;
            }
            else
            {
                if (!byte_swapped)
                {
                    aux_pointer[0] = tif_buf[x*2];
                    aux_pointer[1] = tif_buf[x*2+1];
                }
                else
                {
                    aux_pointer[0] = tif_buf[x*2+1];
                    aux_pointer[1] = tif_buf[x*2];
                }
                if ((double)uiVal < minval) minval = (double)uiVal;
                else if ((double)uiVal > maxval) maxval = (double)uiVal;
            }
            raw_buf[(y*imageWidth+x)*2]   = aux_pointer[0];
            raw_buf[(y*imageWidth+x)*2+1] = aux_pointer[1];
            break;
        }//switch
    }//for
    sprintf(zsMinval, " %f", minval);
    sprintf(zsMaxval, " %f", maxval);
}

/*
//
// main
//
// usage: tiff2raw inputfile outputfile
// inputfile is the TIFF filename
// outputfile is the RAW filename
// (outputfile.inf will be the properties file)
//
*/
int main(int argc, char *argv[])
{
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
    if (argc < 3)
    {
        fprintf(stderr, "Usage: tiff2raw inputfile outputfile\n");
        exit(-1);
    }


    /* Open TIFF image */
    tif = TIFFOpen(argv[1], "r");
    if (tif == NULL) exit(-1);

    /* Get TIFF image properties */
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE,   &bitsPerSample);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,      &imageWidth);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH,     &imageLength);
    TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT,    &imageSampleFormat);

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
    byte_swapped = TIFFIsByteSwapped(tif);

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
}

/* Colimate menu =========================================================== */
/*Colimate:
   PROGRAM Tiff2Raw {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Tiff2Raw/Help/tiff2raw.html";
      help="Convert TIFF micrographs into raw images suitable for marking";
      OPEN MENU Tiff2Raw;
      COMMAND LINES {
         + usual: xmipp_tiff2raw $FILEIN $FILEOUT
      }
      PARAMETER DEFINITIONS {
         $FILEIN {
            label="Input TIFF micrograph";
            type=file existing;
         }
         $FILEOUT {
            label="Output RAW micrograph";
            type=file;
         }
      }
   }
   MENU Tiff2Raw {
      "I/O parameters"
      $FILEIN
      $FILEOUT
   }
*/
