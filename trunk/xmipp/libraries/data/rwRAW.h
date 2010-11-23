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

#ifndef RWRAW_H_
#define RWRAW_H_

///@defgroup RAW RAW Data type
///@ingroup ImageFormats

// I/O prototypes
/** RAW Reader
  * @ingroup RAW
*/

DataType datatypeRAW(std::string strDT)
{
    DataType datatype;

    if(strDT=="uint8")
        datatype = UChar;
    else if (strDT=="int8")
        datatype = SChar;
    else if (strDT=="uint16")
        datatype = UShort;
    else if (strDT=="int16")
        datatype = Short;
    else if (strDT=="uint32")
        datatype = UInt;
    else if (strDT=="int32")
        datatype = Int;
    else if (strDT=="long")
        datatype = Long;
    else if (strDT=="float")
        datatype = Float;
    else if (strDT=="double")
        datatype = Double;
    else if (strDT=="cint16")
        datatype = ComplexShort;
    else if (strDT=="cint32")
        datatype = ComplexInt;
    else if (strDT=="cfloat")
        datatype = ComplexFloat;
    else if (strDT=="cdouble")
        datatype = ComplexDouble;
    else if (strDT=="bool")
        datatype = Bool;
    else
        datatype = Unknown_Type;

    return datatype;
}

///@defgroup RAW RAW File format
///@ingroup ImageFormats

// I/O prototypes
/** RAW Reader
  * @ingroup RAW
*/

int readRAW(int img_select,bool isStack=false)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readRAW: Reading RAW file\n");
#endif

    int _xDim,_yDim,_zDim;
    unsigned long int _nDim;
    int rPos;

    size_t found;
    FileName infolist;
    std::vector<std::string> info;
    DataType datatype;

    found = filename.find_first_of("#");
    infolist = filename.substr(found+1);
    filename = filename.substr(0,found);
    infolist.toLowercase();
    splitString(infolist,",",info, false);

    if (info.size() < 4)
        REPORT_ERROR(ERR_ARG_MISSING, (std::string) " Cannot open file " + filename +
                     ". Not enough header arguments.");

    _xDim = textToInteger(info[0]);
    _yDim = textToInteger(info[1]);

    if (atoi(info[3].c_str()) == 0) // Check if zdim is not included
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
    _nDim = (int) 1;

    // Check the reverse argument
    if (info.back() == "r")
            swap = true;
        else
            swap = false;


    // Map the parameters
    	data.setDimensions(_xDim, _yDim, _zDim, _nDim);

    unsigned long   imgStart=0;
    unsigned long   imgEnd =_nDim;
    if (img_select != -1)
    {
        imgStart=img_select;
        imgEnd=img_select+1;
    }

    MDMainHeader.setValue(MDL_SAMPLINGRATEX,(double) -1);
    MDMainHeader.setValue(MDL_SAMPLINGRATEY,(double) -1);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    if( dataflag < 0 )
        return 0;

    MD.clear();
    MD.resize(imgEnd - imgStart);
    for ( i = imgStart; i<imgEnd; ++i )
    {
        MD[i-imgStart].setValue(MDL_ORIGINX, zeroD);
        MD[i-imgStart].setValue(MDL_ORIGINY, zeroD);
        MD[i-imgStart].setValue(MDL_ORIGINZ,  zeroD);
        MD[i-imgStart].setValue(MDL_ANGLEROT, zeroD);
        MD[i-imgStart].setValue(MDL_ANGLETILT,zeroD);
        MD[i-imgStart].setValue(MDL_ANGLEPSI, zeroD);
        MD[i-imgStart].setValue(MDL_WEIGHT,   oneD);
        MD[i-imgStart].setValue(MDL_FLIP,     falseb);
        MD[i-imgStart].setValue(MDL_SCALE,    oneD);
    }

    size_t pad = 0;
    readData(fimg, img_select, datatype, pad);

    return(0);
}


#endif /* RWINF_H_ */
