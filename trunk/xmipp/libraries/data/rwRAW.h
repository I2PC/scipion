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



// I/O prototypes
int readRAW(int img_select,bool isStack=false)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readTIA: Reading TIA file\n");
#endif

    int _xDim,_yDim,_zDim, __depth;
    unsigned long int _nDim;
    bool __is_signed;

    FileName fn_inf;

    fn_inf = filename.add_extension("inf");
    FILE *fh_inf = fopen(fn_inf.c_str(), "r");
    if (!fh_inf)
        REPORT_ERROR(1, (std::string)"Micrograph::open_micrograph: Cannot find " +
                     fn_inf);
    _xDim = textToInteger(getParameter(fh_inf, "Xdim"));
    _yDim = textToInteger(getParameter(fh_inf, "Ydim"));
    __depth = textToInteger(getParameter(fh_inf, "bitspersample"));
    if (checkParameter(fh_inf, "offset"))
        offset = textToInteger(getParameter(fh_inf, "offset"));
    else
        offset = 0;
    if (checkParameter(fh_inf, "is_signed"))
        __is_signed = (getParameter(fh_inf, "is_signed") == "true" ||
                       getParameter(fh_inf, "is_signed") == "TRUE");
    else
        __is_signed = false;
    if (checkParameter(fh_inf, "reversed"))
        swap = (getParameter(fh_inf, "reversed") == "true" ||
                getParameter(fh_inf, "reversed") == "TRUE");
    else
    	swap = false;
    fclose(fh_inf);

    _zDim = (int) 1;
    _nDim = (int) 1;

    // Map the parameters
    data.setDimensions(_xDim, _yDim, _zDim, _nDim);

    unsigned long   imgStart=0;
    unsigned long   imgEnd =_nDim;
    if (img_select != -1)
    {
        imgStart=img_select;
        imgEnd=img_select+1;
    }

    DataType datatype;

    switch ( __depth )
    {
    case 8:
        datatype = UChar;
        break;
    case 16:
        if (__is_signed)
            datatype = Short;
        else
            datatype = UShort;
        break;
    case 32:
        datatype = Float;
        break;
    default:
        REPORT_ERROR(1, "Micrograph::open_micrograph: depth is not 8, 16 nor 32");
    }

    MDMainHeader.removeObjects();
    MDMainHeader.setColumnFormat(false);
    MDMainHeader.addObject();
    MDMainHeader.setValue(MDL_SAMPLINGRATEX,(double) -1);
    MDMainHeader.setValue(MDL_SAMPLINGRATEY,(double) -1);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    if( dataflag == -2 )
    {
//        fclose(fimg);
        return 0;
    }

    MD.removeObjects();
    for ( i=imgStart; i<imgEnd; i++ )
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

    //#define DEBUG
#ifdef DEBUG

    MDMainHeader.write(std::cerr);
    MD.write(std::cerr);
#endif

    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "r") ) == NULL )
        return(-1);

    size_t pad = 0;

    readData(fimg, img_select, datatype, pad);

    if ( !mmapOn )
        fclose(fimg);

    return(0);
}

int writeRAW(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE)
{
    REPORT_ERROR(6001, "ERROR: writeRAW is not implemented.");
    return(-1);
}


#endif /* RWRAW_H_ */
