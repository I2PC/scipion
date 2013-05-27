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
#include "args.h"

///@defgroup INF INF File format
///@ingroup ImageFormats

// I/O prototypes
/** INF Reader
  * @ingroup INF
*/
int ImageBase::readINF(size_t select_img,bool isStack)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readINF: Reading INF file\n");
#endif

    int _xDim,_yDim,_zDim, __depth;
    size_t _nDim;
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

    _zDim = 1;
    _nDim = 1;

    DataType datatype;
    switch ( __depth )
    {
    case 8:
        if (__is_signed)
            datatype = DT_SChar;
        else
            datatype = DT_UChar;
        break;
    case 16:
        if (__is_signed)
            datatype = DT_Short;
        else
            datatype = DT_UShort;
        break;
    case 32:
        datatype = DT_Float;
        break;
    default:
        REPORT_ERROR(ERR_TYPE_INCORRECT, "rwINF::read: depth is not 8, 16 nor 32");
    }

    MDMainHeader.setValue(MDL_SAMPLINGRATE_X,(double) -1);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_Y,(double) -1);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    // Map the parameters
    setDimensions(_xDim, _yDim, _zDim, _nDim);

    size_t   imgStart = IMG_INDEX(select_img);
    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : _nDim;

    if (dataMode == HEADER || (dataMode == _HEADER_ALL && _nDim > 1)) // Stop reading if not necessary
        return 0;

    MD.clear();
    MD.resize(imgEnd - imgStart,MDL::emptyHeader);

    //#define DEBUG
#ifdef DEBUG

    MDMainHeader.write(std::cerr);
    MD.write(std::cerr);
#endif

    if( dataMode < DATA )
        return 0;

    size_t pad = 0;
    readData(fimg, select_img, datatype, pad);
    return(0);
}

/** INF Writer
  * @ingroup INF
*/
int ImageBase::writeINF(size_t select_img, bool isStack, int mode, String bitDepth, bool adjust)
{
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG writeINF: Writing INF file\n");
    printf("DEBUG writeINF: File %s\n", filename.c_str());
#endif
#undef DEBUG

    size_t Xdim, Ydim, Zdim, Ndim;
    getDimensions(Xdim, Ydim, Zdim, Ndim);

    int _depth;
    bool _is_signed = false;

    // Volumes and stacks are not supported
    if (Zdim > 1 || Ndim > 1)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "rwINF::write does not support neither volumes nor stacks.");

    if (mode == WRITE_APPEND)
        REPORT_ERROR(ERR_ARG_INCORRECT, "rwINF::write only can overwrite or replace image files,"
                     "not append.");

    DataType wDType,myTypeID = myT();
    CastWriteMode castMode = CW_CAST;

    if (bitDepth != "")
    {
        myTypeID = (bitDepth == "default") ? DT_Float : datatypeRAW(bitDepth);
        castMode = (adjust)? CW_ADJUST : CW_CONVERT;
    }

    switch(myTypeID)
    {
    case DT_Double:
    case DT_UInt:
    case DT_Int:
        if (bitDepth != "")
            REPORT_ERROR(ERR_TYPE_INCORRECT,"ERROR: incorrect RAW bits depth value.");
    case DT_Float:
        wDType = DT_Float;
        break;
    case DT_UShort:
        wDType = DT_UShort;
        _is_signed = false;
        break;
    case DT_Short:
        wDType = DT_Short;
        _is_signed = true;
        break;
    case DT_UChar:
        wDType = DT_UChar;
        _is_signed = false;
        break;
    case DT_SChar:
        wDType = DT_SChar;
        _is_signed = true;
        break;
    default:
        wDType = DT_Unknown;
        REPORT_ERROR(ERR_TYPE_INCORRECT,(std::string)"ERROR: Unsupported data type by RAW format.");
    }

    if (mmapOnWrite)
    {
        MDMainHeader.setValue(MDL_DATATYPE,(int) wDType);
        if (!checkMmapT(wDType))
        {
            if (dataMode < DATA && castMode == CW_CAST) // This means ImageGeneric wants to know which DataType must use in mapFile2Write
                return 0;
            else //Mapping is an extra. When not available, go on and do not report an error.
            {
                /* In this case we cannot map the file because required and feasible datatypes are
                 * not compatible. Then we denote to MapFile2Write the same incoming datatype to
                 * keep using this Image object as usual, without mapping on write.
                 */
                mmapOnWrite = false;
                dataMode = DATA;
                MDMainHeader.setValue(MDL_DATATYPE,(int) myTypeID);

                // In case Image size great then, at least, map the multidimarray
                if (mdaBase->nzyxdim*gettypesize(myTypeID) > tiff_map_min_size)
                    mdaBase->setMmap(true);

                // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
                //if memory already allocated use it (no resize allowed)
                mdaBase->coreAllocateReuse();

                return 0;
            }
        }
        else
            dataMode = DATA;
    }

    _depth = gettypesize(wDType);

    //locking

    // Lock Header file
    FileLock flock;
    flock.lock(fhed);

    /* Write INF file ==================================*/
    fprintf(fhed,"# Bits per sample\n");
    fprintf(fhed,"bitspersample= %d\n",_depth*8);
    fprintf(fhed,"# Samples per pixel\n");
    fprintf(fhed,"samplesperpixel= 1\n");
    fprintf(fhed,"# Image width\n");
    fprintf(fhed,"Xdim= %d\n", (int)Xdim);
    fprintf(fhed,"# Image length\n");
    fprintf(fhed,"Ydim= %d\n",(int)Ydim);
    fprintf(fhed,"# offset in bytes (zero by default)\n");
    fprintf(fhed,"offset= 0\n");
    fprintf(fhed,"# Is a signed or Unsigned int (by default true)\n");
    if (_is_signed)
        fprintf(fhed,"is_signed= true\n");
    else
        fprintf(fhed,"is_signed= false\n");
    fprintf(fhed,"# Byte order\n");
    if ( swapWrite^IsBigEndian() )
        fprintf(fhed,"endianess= big\n");
    else
        fprintf(fhed,"endianess= little\n");

    //Unlock Header file
    flock.unlock();

    /* Write Image file ==================================*/
    size_t datasize_n, datasize;
    datasize_n = Xdim*Ydim*Zdim;
    datasize = datasize_n * gettypesize(wDType);

    // Lock Image file
    flock.lock(fimg);

    if (mmapOnWrite)
    {
        mappedOffset = 0;
        mappedSize = mappedOffset + datasize;
        fseek(fimg, datasize-1, SEEK_SET);
        fputc(0, fimg);
        mmapFile();
    }
    else
        writeData(fimg, 0, wDType, datasize_n, castMode);

    // Unlock Image file
    flock.unlock();

    return(0);
}
