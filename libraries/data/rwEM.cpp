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

#define EMHEADERSIZE 512 // size of EM file header

int ImageBase::readEM(size_t select_img)
{
    // EM File formats does not support stacks
    if (select_img > FIRST_IMAGE)
        REPORT_ERROR(ERR_ARG_INCORRECT, "readEM: EM file format does not support stacks.");

    EMHead header;

    if ( fread( &header, EMHEADERSIZE, 1, fimg ) != 1 )
        REPORT_ERROR(ERR_IO_NOREAD, formatString("rwEM: cannot read EM main header from file %s"
                     ". Error message: %s", filename.c_str() ,strerror(errno)));

    // endian: If machine is SGI, OS-9 or MAC: Big Endian, otherwise Litle Endian
    // Check Machine endianess
    bool isLE = IsLittleEndian();

    if (header.machine == 0 || header.machine == 3 || header.machine == 5)
        swap = isLE;
    else if (header.machine == 1 || header.machine == 2 || header.machine == 4 || header.machine == 6)
        swap = !isLE;
    else
        REPORT_ERROR(ERR_IMG_UNKNOWN, "rwEM: Unknown source machine to determine Endianness");

    if (swap)
        swapPage((char *) &header, EMHEADERSIZE - 256, DT_UInt); // EMHEADERSIZE - 256 is to exclude userdata from swapping

    // Setting image dimensions
    replaceNsize = 1;

    ArrayDim aDim;
    aDim.ndim = 1;
    aDim.xdim = header.xdim;
    aDim.ydim= header.ydim;
    aDim.zdim = header.zdim;

    setDimensions(aDim);


    DataType datatype;
    switch (header.datatype)
    {
    case 1:
        datatype = DT_SChar;
        break;
    case 2:
        datatype = DT_Short;
        break;
    case 4:
        datatype = DT_Long;
        break;
    case 3:
    case 5:
        datatype = DT_Float;
        break;
    case 8:
        datatype = DT_CFloat;
        break;
    case 9:
        datatype = DT_CDouble;
        break;
    default:
        REPORT_ERROR(ERR_ARG_INCORRECT, formatString("readEM: Unknown datatype value: %c", header.datatype));
        break;
    }
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    offset = EMHEADERSIZE;

    // If only header is read: return
    if (dataMode==HEADER || (dataMode == _HEADER_ALL && aDim.ndim > 1)) // Stop reading if not necessary
        return 0;

    size_t   imgStart = IMG_INDEX(select_img);
    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : aDim.ndim;

    MD.clear();
    MD.resize(imgEnd - imgStart,MDL::emptyHeader);

    /* As MRC does not support stacks, we use the geometry stored in the header
    for any image when we simulate the file is a stack.*/
    if (dataMode == _HEADER_ALL || dataMode == _DATA_ALL)
    {
        for ( size_t i = 0; i < imgEnd - imgStart; ++i )
        {
            MD[i].setValue(MDL_SHIFT_X, 0.);
            MD[i].setValue(MDL_SHIFT_Y, 0.);
            MD[i].setValue(MDL_SHIFT_Z, 0.);
            MD[i].setValue(MDL_ORIGIN_X, 0.);
            MD[i].setValue(MDL_ORIGIN_Y, 0.);
            MD[i].setValue(MDL_ORIGIN_Z, 0.);
        }
    }

    if ( dataMode < DATA )   // Don't read the individual header and the data if not necessary
        return 0;

    readData(fimg, select_img, datatype, 0);

    return(0);
}


int ImageBase::writeEM(size_t select_img, bool isStack, int mode)
{
    REPORT_ERROR(ERR_IO_NOWRITE, "ERROR: writeEM is not implemented.");
    return(-1);
}
