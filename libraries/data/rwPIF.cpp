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

#define PIFHEADERSIZE sizeof(PIFDataHeader) // size of EM file header

int  ImageBase::readPIF(size_t select_img)
{
    PIFMainHeader  mainHeader;
    if ( fread( &mainHeader, PIFHEADERSIZE, 1, fimg ) != 1 )
        REPORT_ERROR(ERR_IO_NOREAD, formatString("rwPIF: cannot read Spider main header from file %s"
                     ". Error message: %s", filename.c_str() ,strerror(errno)));

    // Check Machine endianess
    bool isLE = IsLittleEndian();

    // Determine byte order and swap bytes if from different-endian machine
    swap = (isLE == mainHeader.endianNess);

    if (swap)
        swapPage((char *) &mainHeader, PIFHEADERSIZE, DT_Int);


    DataType datatype;
    switch (mainHeader.mode)
    {
    case 0:
        datatype = DT_SChar;
        break;
    case 1:
    case 7: // Actually it's floatint*2, but it isn't moved to float
        datatype = DT_Short;
        break;
    case 2:
    case 46: // This case is not documented, but I think this is floatInt*4
        datatype = DT_Int; // We aren't moving to float actually
        break;
    case 3:
    case 8: // Actually it's complex floatint*2, but it isn't moved to float
        datatype = DT_CShort;
        break;
    case 4:
        datatype = DT_CInt; // We aren't moving to float actually
        break;
    case 9:
        datatype = DT_Float;
        break;
    case 10:
        datatype = DT_CFloat;
        break;
    case 97:
        REPORT_ERROR(ERR_NOT_IMPLEMENTED, "readPIF: Reading from colormap datatype file not implemented.\n"
                     "If you're interested in this feature, please contact us at xmipp@cnb.csic.es \n"
                     "and send us an example image file to help implementing.");
        break;
    default:
        REPORT_ERROR(ERR_NOT_IMPLEMENTED, formatString("readPIF: Reading from datatype code %d file not implemented.\n"
                     "If you're interested in this feature, please contact us at xmipp@cnb.csic.es \n"
                     "and send us an example image file to help implementing.", mainHeader.mode));
        break;
    }
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);


    // Check selected image
    if (select_img > (size_t)mainHeader.numImages)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, formatString("readPIF (%s): Image number %lu exceeds stack size %lu" ,filename.c_str(),select_img, mainHeader.numImages));

    // Setting image dimensions
    ArrayDim aDim;

    aDim.ndim = (select_img == ALL_IMAGES)? mainHeader.numImages : 1;
    aDim.zdim = mainHeader.nz;
    aDim.ydim = mainHeader.ny;
    aDim.xdim = mainHeader.nx;

    replaceNsize = aDim.ndim;

    setDimensions(aDim);


    size_t imgStart = IMG_INDEX(select_img);
    size_t datatypeSize = gettypesize(datatype);
    offset = PIFHEADERSIZE*2;// + (aDimFile.zyxdim*datatypeSize + PIFHEADERSIZE)*imgStart;

    if (dataMode == HEADER || (dataMode == _HEADER_ALL && aDim.ndim > 1)) // Stop reading if not necessary
        return 0;

    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : aDim.ndim;

    size_t imageSize = PIFHEADERSIZE + aDimFile.zyxdim*datatypeSize;
    size_t pad       = PIFHEADERSIZE;
    size_t headerOffset  = PIFHEADERSIZE;

    MD.clear();
    MD.resize(imgEnd - imgStart,MDL::emptyHeader);

    PIFDataHeader dataHeader;

    for (size_t n = 0, i = imgStart; i < imgEnd; ++i, ++n )
    {
        if (fseek( fimg, headerOffset + i*imageSize, SEEK_SET ) != 0)//fseek return 0 on success
            REPORT_ERROR(ERR_IO, formatString("rwPIF: error seeking %lu to read image %lu", headerOffset, i));

        if ( fread( &dataHeader, PIFHEADERSIZE, 1, fimg ) != 1 )
            REPORT_ERROR(ERR_IO_NOREAD, formatString("rwPIF: cannot read PIF image %lu header", i));

        if ( swap )
            swapPage((char *) &dataHeader, PIFHEADERSIZE, DT_Int);

        if (dataMode == _HEADER_ALL || dataMode == _DATA_ALL)
        {
            MD[n].setValue(MDL_ORIGIN_X, (double)dataHeader.nxstart);
            MD[n].setValue(MDL_ORIGIN_Y, (double)dataHeader.nystart);
            MD[n].setValue(MDL_ORIGIN_Z, (double)dataHeader.nzstart);
            MD[n].setValue(MDL_ANGLE_ROT, (double)dataHeader.alpha);
            MD[n].setValue(MDL_ANGLE_TILT, (double)dataHeader.beta);
            MD[n].setValue(MDL_ANGLE_PSI, (double)dataHeader.gamma);
            MD[n].setValue(MDL_SAMPLINGRATE_X, (double)dataHeader.xlength/aDim.xdim);
            MD[n].setValue(MDL_SAMPLINGRATE_Y, (double)dataHeader.ylength/aDim.ydim);
            MD[n].setValue(MDL_SAMPLINGRATE_Z, (double)dataHeader.zlength/aDim.zdim);
        }
    }

    if (dataMode < DATA)   // Don't read  data if not necessary but read the header
        return 0;

    //offset should point to the begin of the data
    readData(fimg, select_img, datatype, pad );

    return(0);
}

int ImageBase::writePIF(size_t select_img, bool isStack, int mode)
{
    REPORT_ERROR(ERR_IO_NOWRITE, "ERROR: writePIF is not implemented.");
    return(-1);
}
