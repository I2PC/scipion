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

#ifndef RWINF_H_
#define RWINF_H_

///@defgroup INF INF File format
///@ingroup ImageFormats

// I/O prototypes
/** INF Reader
  * @ingroup INF
*/
int readINF(int img_select,bool isStack=false)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readINF: Reading INF file\n");
#endif

    int _xDim,_yDim,_zDim, __depth;
    unsigned long int _nDim;
    bool __is_signed;

    _xDim = textToInteger(getParameter(fhed, "Xdim"));
    _yDim = textToInteger(getParameter(fhed, "Ydim"));
    __depth = textToInteger(getParameter(fhed, "bitspersample"));
    if (checkParameter(fhed, "offset"))
        offset = textToInteger(getParameter(fhed, "offset"));
    else
        offset = 0;
    if (checkParameter(fhed, "is_signed"))
        __is_signed = (getParameter(fhed, "is_signed") == "true" ||
                       getParameter(fhed, "is_signed") == "TRUE");
    else
        __is_signed = false;
    if (checkParameter(fhed, "endianess") &&
        (getParameter(fhed, "endianess") == "big" || getParameter(fhed, "endianess") == "BIG"))
        swap = true;
    else
        swap = false;

    if (IsBigEndian())
        swap = !swap;


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
        REPORT_ERROR(ERR_TYPE_INCORRECT, "rwINF::read: depth is not 8, 16 nor 32");
    }

    MDMainHeader.setValue(MDL_SAMPLINGRATEX,(double) -1);
    MDMainHeader.setValue(MDL_SAMPLINGRATEY,(double) -1);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    if( dataflag < 0 )
        return 0;

    MD.clear();
    MD.resize(imgEnd - imgStart);
    for ( i = imgStart; i < imgEnd; ++i )
    {
        MD[i-imgStart].setValue(MDL_ORIGINX, zeroD);
        MD[i-imgStart].setValue(MDL_ORIGINY, zeroD);
        MD[i-imgStart].setValue(MDL_ORIGINZ,  zeroD);
        MD[i-imgStart].setValue(MDL_ANGLEROT, zeroD);
        MD[i-imgStart].setValue(MDL_ANGLETILT,zeroD);
        MD[i-imgStart].setValue(MDL_ANGLEPSI, zeroD);
        MD[i-imgStart].setValue(MDL_WEIGHT,   oneD);
        MD[i-imgStart].setValue(MDL_FLIP,     falseb);
    }

    //#define DEBUG
#ifdef DEBUG
    MDMainHeader.write(std::cerr);
    MD.write(std::cerr);
#endif


    size_t pad = 0;
    readData(fimg, img_select, datatype, pad);
    return(0);
}

/** INF Writer
  * @ingroup INF
*/
int writeINF(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE)
{
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG writeINF: Writing INF file\n");
    printf("DEBUG writeINF: File %s\n", filename.c_str());
#endif
#undef DEBUG

    int Xdim = XSIZE(data);
    int Ydim = YSIZE(data);
    int Zdim = ZSIZE(data);
    int Ndim = NSIZE(data);

    int _depth;
    bool _is_signed;

    // Volumes and stacks are not supported
    if (Zdim > 1 || Ndim > 1)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "rwINF::write does not support neither volumes nor stacks.");

    if (mode != WRITE_OVERWRITE)
        REPORT_ERROR(ERR_ARG_INCORRECT, "rwINF::write only can overwrite image files,"
                     "neither append nor replace.");

    DataType wDType;

    if (typeid(T)==typeid(double)||typeid(T)==typeid(float))
    {
        wDType = Float;
        _is_signed = true;
    }
    else if (typeid(T)==typeid(int))
    {
        wDType = Int;
        _is_signed = true;
    }
    else if (typeid(T)==typeid(unsigned int))
    {
        wDType = UInt;
        _is_signed = false;

    }
    else if (typeid(T)==typeid(short))
    {
        wDType = Short;
        _is_signed = true;
    }
    else if (typeid(T)==typeid(unsigned short))
    {
        wDType = UShort;
        _is_signed = false;
    }
    else if (typeid(T)==typeid(char))
    {
        wDType = SChar;
        _is_signed = true;
    }
    else if (typeid(T)==typeid(unsigned char))
    {
        wDType = UChar;
        _is_signed = false;
    }
    else
        REPORT_ERROR(ERR_TYPE_INCORRECT,(std::string)"ERROR: rwINF does not write from " + typeid(T).name() + "type.");

    _depth = gettypesize(wDType);

    //locking
    struct flock fl;

    fl.l_type   = F_WRLCK;  /* F_RDLCK, F_WRLCK, F_UNLCK    */
    fl.l_whence = SEEK_SET; /* SEEK_SET, SEEK_CUR, SEEK_END */
    fl.l_start  = 0;        /* Offset from l_whence         */
    fl.l_len    = 0;        /* length, 0 = to EOF           */
    fl.l_pid    = getpid(); /* our PID                      */

    // Lock Header file
    fl.l_type   = F_WRLCK;
    fcntl(fileno(fhed), F_SETLKW, &fl); /* locked */

    /* Write INF file ==================================*/
    fprintf(fhed,"# Bits per sample\n");
    fprintf(fhed,"bitspersample= %d\n",_depth*8);
    fprintf(fhed,"# Samples per pixel\n");
    fprintf(fhed,"samplesperpixel= 1\n");
    fprintf(fhed,"# Image width\n");
    fprintf(fhed,"Xdim= %d\n", Xdim);
    fprintf(fhed,"# Image length\n");
    fprintf(fhed,"Ydim= %d\n",Ydim);
    fprintf(fhed,"# offset in bytes (zero by default)\n");
    fprintf(fhed,"offset= 0\n");
    fprintf(fhed,"# Is a signed or Unsigned int (by default true)\n");
    if (_is_signed)
        fprintf(fhed,"is_signed = true\n");
    else
        fprintf(fhed,"is_signed = false\n");
    fprintf(fhed,"# Byte order\n");
    if (IsBigEndian())
        fprintf(fhed,"endianess = big\n");
    else
        fprintf(fhed,"endianess = little\n");

    //Unlock Header file
    fl.l_type = F_UNLCK;
    fcntl(fileno(fhed), F_SETLK, &fl);

    // Lock Image file
    fl.l_type   = F_WRLCK;
    fcntl(fileno(fimg), F_SETLKW, &fl);

    /* Write Image file ==================================*/
    size_t datasize_n;
    datasize_n = Xdim*Ydim*Zdim;

    writePageAsDatatype(fimg, wDType, datasize_n);

    // Unlock Image file
    fl.l_type   = F_UNLCK;
    fcntl(fileno(fimg), F_SETLK, &fl);

    return(0);
}

#endif /* RWINF_H_ */
