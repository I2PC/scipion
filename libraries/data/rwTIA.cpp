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

#define TIASIZE    30 // Size of the TIA header without pDATA_OFFSET

///@defgroup TIA TIA File format
///@ingroup ImageFormats

/** TIA Header
  * @ingroup TIA
*/
struct TIAhead
{
    short int endianess;
    short int SeriesID;
    short int SeriesVersion;
    int DATA_TYPE_ID;
    int TagTypeID;
    int TotalNumberElements;
    int NUMBER_IMAGES;
    int OFFSET_ARRAY_OFFSET;
    int numberdimensions;
    int * pDATA_OFFSET;
};

#define TIAdataSIZE    50 // Size of the TIA data header to be read

/** TIA Data Header
  * @ingroup TIA
*/
struct TIAdataHead
{
    double      CalibrationOffsetX;   //CalibrationOffsetX
    double      PIXEL_WIDTH;          //CalibrationDeltaX
    int         CalibrationElementX;  //CalibrationElementX
    double      CalibrationOffsetY;   //CalibrationOffsetY
    double      PIXEL_HEIGHT;         //CalibrationDeltaY
    int          CalibrationElementY; //CalibrationElementY
    short int   DATA_TYPE;            //DataType
    int         IMAGE_WIDTH;          //ArraySizeX
    int         IMAGE_HEIGHT;         //ArraySizeY
    short int   DATA_TYPE_SIZE;
    std::string  DATA_TYPE_SIZE_STRING;
    bool        isSigned;
};

// I/O prototypes
/** TIA Reader
  * @ingroup TIA
*/
int ImageBase::readTIA(int select_img,bool isStack)
{
    TIAhead * header = new TIAhead;

    xmippFREAD(&header->endianess, sizeof(short int), 1, fimg, swap );

    // Set Endianess
    if (header->endianess == 18761)
        swap = 0;
    else
        swap = 1;
    if (IsBigEndian())
        swap = !swap;

    xmippFREAD(&header->SeriesID, sizeof(short int), 1, fimg, swap );
    xmippFREAD(&header->SeriesVersion, sizeof(short int), 1, fimg, swap);
    xmippFREAD(&header->DATA_TYPE_ID, sizeof(int), 1, fimg, swap);
    xmippFREAD(&header->TagTypeID, sizeof(int), 1, fimg, swap );
    xmippFREAD(&header->TotalNumberElements, sizeof(int), 1, fimg, swap );
    xmippFREAD(&header->NUMBER_IMAGES, sizeof(int), 1, fimg, swap );
    xmippFREAD(&header->OFFSET_ARRAY_OFFSET, sizeof(int), 1, fimg, swap );
    xmippFREAD(&header->numberdimensions, sizeof(int), 1, fimg, swap );

    // Check data type
    if (header->DATA_TYPE_ID != 16674)
        REPORT_ERROR(ERR_TYPE_INCORRECT, "ERROR: readTIA only processes images in real space");

    fseek(fimg, header->OFFSET_ARRAY_OFFSET, SEEK_SET);
    header->pDATA_OFFSET = (int *) askMemory(header->NUMBER_IMAGES * sizeof(int));
    xmippFREAD(header->pDATA_OFFSET, sizeof(int), header->NUMBER_IMAGES, fimg, swap);

    TIAdataHead* dataHeaders = new TIAdataHead [header->NUMBER_IMAGES];

    // Read all the image headers
    for (int i = 0; i < header->NUMBER_IMAGES; i++)
    {
        fseek(fimg, header->pDATA_OFFSET[i], SEEK_SET);
        xmippFREAD(&(dataHeaders[i].CalibrationOffsetX), sizeof(double), 1, fimg, swap);
        xmippFREAD(&dataHeaders[i].PIXEL_WIDTH, sizeof(double), 1, fimg, swap);
        xmippFREAD(&dataHeaders[i].CalibrationElementX, sizeof(int), 1, fimg, swap);
        xmippFREAD(&dataHeaders[i].CalibrationOffsetY, sizeof(double), 1, fimg, swap);
        xmippFREAD(&dataHeaders[i].PIXEL_HEIGHT, sizeof(double), 1, fimg, swap);
        xmippFREAD(&dataHeaders[i].CalibrationElementY, sizeof(int), 1, fimg, swap);
        xmippFREAD(&dataHeaders[i].DATA_TYPE, sizeof(short int), 1, fimg, swap);
        xmippFREAD(&dataHeaders[i].IMAGE_WIDTH, sizeof(int), 1, fimg, swap);
        xmippFREAD(&dataHeaders[i].IMAGE_HEIGHT, sizeof(int), 1, fimg, swap);
    }

    int _xDim,_yDim;
    size_t _nDim;

    size_t   imgStart = IMG_INDEX(select_img);
    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : header->NUMBER_IMAGES;

    if (select_img >  header->NUMBER_IMAGES)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, formatString("readTIA: Image number %lu exceeds stack size %lu", select_img, header->NUMBER_IMAGES));
    else if (select_img == ALL_IMAGES)
    {
        for (int i = 1; i < header->NUMBER_IMAGES; i++)    // Check images dimensions. Need to be the same
        {
            if (dataHeaders[0].IMAGE_HEIGHT != dataHeaders[i].IMAGE_HEIGHT || \
                dataHeaders[0].IMAGE_WIDTH != dataHeaders[i].IMAGE_WIDTH  || \
                dataHeaders[0].DATA_TYPE != dataHeaders[i].DATA_TYPE)
                REPORT_ERROR(ERR_IMG_NOREAD, "readTIA: images in TIA file with different dimensions and data types are not supported");
        }
        _xDim = dataHeaders[0].IMAGE_WIDTH;
        _yDim = dataHeaders[0].IMAGE_HEIGHT;
        _nDim = (size_t) header->NUMBER_IMAGES;
    }
    else
    {
        _xDim = dataHeaders[imgStart].IMAGE_WIDTH;
        _yDim = dataHeaders[imgStart].IMAGE_HEIGHT;
        _nDim = 1;
    }

    // Map the parameters
    setDimensions(_xDim, _yDim, 1, _nDim);

    DataType datatype;
    //    dataHeaders[0].isSigned = false;
    int tiaDT;
    tiaDT = dataHeaders[imgStart].DATA_TYPE;
    offset = header->pDATA_OFFSET[imgStart] + TIAdataSIZE;

    switch ( tiaDT )
    {
    case 1:
        datatype = DT_UChar;
        break;
    case 2:
        datatype = DT_UShort;
        //        datatype = DT_Short;
        break;
    case 3:
        datatype = DT_UInt;
        break;
    case 4:
        datatype = DT_SChar;
        break;
    case 5:
        datatype = DT_Short;
        //        dataHeaders[0].isSigned = true;
        break;
    case 6:
        datatype = DT_Int;
        break;
    case 7:
        datatype = DT_Float;
        break;
    case 8:
        datatype = DT_Double;
        break;
    case 9:
        datatype = DT_CFloat;
        break;
    case 10:
        datatype = DT_CDouble;
        break;
    default:
        datatype = DT_Unknown;
        break;
    }

    MDMainHeader.setValue(MDL_SAMPLINGRATE_X,(double)dataHeaders[0].PIXEL_WIDTH);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_Y,(double)dataHeaders[0].PIXEL_HEIGHT);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    if (dataMode == HEADER || (dataMode == _HEADER_ALL && _nDim > 1)) // Stop reading if not necessary
    {
        delete header;
        return 0;
    }

    MD.clear();
    MD.resize(imgEnd - imgStart,MDL::emptyHeader);
    double aux;
    for ( size_t i = 0; i < imgEnd - imgStart; ++i )
    {
        if (dataMode == _HEADER_ALL || dataMode == _DATA_ALL)
        {
            if(MDMainHeader.getValue(MDL_SAMPLINGRATE_X,aux))
            {
                MD[i].setValue(MDL_SHIFT_X, dataHeaders[i].CalibrationOffsetX/aux);
                aux = ROUND(dataHeaders[i].CalibrationElementX - \
                            dataHeaders[i].CalibrationOffsetX/aux - _xDim/2);
                MD[i].setValue(MDL_ORIGIN_X, aux);
            }
            if(MDMainHeader.getValue(MDL_SAMPLINGRATE_Y,aux))
            {
                MD[i].setValue(MDL_SHIFT_Y, dataHeaders[i].CalibrationOffsetY/aux);
                aux = ROUND(dataHeaders[i].CalibrationElementY - \
                            dataHeaders[i].CalibrationOffsetY/aux -_yDim/2);
                MD[i].setValue(MDL_ORIGIN_Y, aux);
            }
        }
    }

    delete header;

    if( dataMode < DATA )
        return 0;

    size_t pad = TIAdataSIZE;
    readData(fimg, select_img, datatype, pad);

    return(0);
}
