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


DataType ImageBase::datatypeRAW(String strDT)
{
    DataType datatype;

    if(strDT=="uint8")
        datatype = DT_UChar;
    else if (strDT=="int8")
        datatype = DT_SChar;
    else if (strDT=="uint16")
        datatype = DT_UShort;
    else if (strDT=="int16")
        datatype = DT_Short;
    else if (strDT=="uint32")
        datatype = DT_UInt;
    else if (strDT=="int32")
        datatype = DT_Int;
    else if (strDT=="long")
        datatype = DT_Long;
    else if (strDT=="float")
        datatype = DT_Float;
    else if (strDT=="double")
        datatype = DT_Double;
    else if (strDT=="cint16")
        datatype = DT_CShort;
    else if (strDT=="cint32")
        datatype = DT_CInt;
    else if (strDT=="cfloat")
        datatype = DT_CFloat;
    else if (strDT=="cdouble")
        datatype = DT_CDouble;
    else if (strDT=="bool")
        datatype = DT_Bool;
    else
        datatype = DT_Unknown;

    return datatype;
}


int ImageBase::readRAW(size_t select_img, bool isStack)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readRAW: Reading RAW file\n");
#endif

    int _xDim,_yDim,_zDim;
    size_t _nDim;
    int rPos;

    size_t found;
    FileName infolist;
    StringVector info;
    DataType datatype;

    found = filename.find_first_of("#");
    infolist = filename.substr(found+1);
    filename = filename.substr(0,found);
    infolist.toLowercase();
    splitString(infolist,",",info, false);

    if (info.size() < 4)
        REPORT_ERROR(ERR_ARG_MISSING, (String) " Cannot open file " + filename +
                     ". Not enough header arguments.");

    _xDim = textToInteger(info[0]);
    _yDim = textToInteger(info[1]);

    if (atoi(info[3].c_str()) == 0 && info[3]!="0") // Check if zdim is not included
    {
        rPos = 2;
        _zDim = 1;
    }
    else
    {
        rPos = 3;
        _zDim = textToInteger(info[2]);
    }

    offset = (size_t)textToInteger(info[rPos]);
    datatype = datatypeRAW(info[rPos+1]);
    _nDim = 1;

    // Check the reverse argument
    swap = (info.back() == "r");

    // Map the parameters
    setDimensions(_xDim, _yDim, _zDim, _nDim);

    MDMainHeader.setValue(MDL_SAMPLINGRATE_X,(double) -1);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_Y,(double) -1);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    if (dataMode==HEADER || (dataMode == _HEADER_ALL && _nDim > 1)) // Stop reading if not necessary
        return 0;

    size_t   imgStart = IMG_INDEX(select_img);
    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : _nDim;

    MD.clear();
    MD.resize(imgEnd - imgStart,MDL::emptyHeader);

    if( dataMode < DATA )
        return 0;

    size_t pad = 0;
    readData(fimg, select_img, datatype, pad);

    return(0);
}
