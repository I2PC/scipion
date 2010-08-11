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

#ifndef RWTIA_H_
#define RWTIA_H_

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
int readTIA(int img_select,bool isStack=false, double dStddev=5)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readTIA: Reading TIA file\n");
#endif

    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "r") ) == NULL )
        return(-1);

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
        REPORT_ERROR(6001, "ERROR: readTIA only processes images in real space");

    fseek(fimg, header->OFFSET_ARRAY_OFFSET, SEEK_SET);
    header->pDATA_OFFSET = (int *) askMemory(header->NUMBER_IMAGES * sizeof(int));
    xmippFREAD(header->pDATA_OFFSET, sizeof(int), header->NUMBER_IMAGES, fimg, swap);

    TIAdataHead* dataHeaders = new TIAdataHead [header->NUMBER_IMAGES];

    // Read all the image headers
    for (i = 0; i < header->NUMBER_IMAGES; i++)
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

    // Check images dimensions. Need to be the same
    int _xDim,_yDim,_zDim;
    unsigned long int _nDim;

    if (img_select==-1)
    {
        for (i = 1; i < header->NUMBER_IMAGES; i++)
        {
            if (dataHeaders[0].IMAGE_HEIGHT != dataHeaders[i].IMAGE_HEIGHT || \
                dataHeaders[0].IMAGE_WIDTH != dataHeaders[i].IMAGE_WIDTH  || \
                dataHeaders[0].DATA_TYPE != dataHeaders[i].DATA_TYPE)
                REPORT_ERROR(6001, "readTIA: images in TIA file with different dimensions and data types are not supported");
        }
        _xDim = (int) dataHeaders[0].IMAGE_WIDTH;
        _yDim = (int) dataHeaders[0].IMAGE_HEIGHT;
        _zDim = (int) 1;
        _nDim = (int) header->NUMBER_IMAGES;
    }
    else
    {
        _xDim = (int) dataHeaders[img_select].IMAGE_WIDTH;
        _yDim = (int) dataHeaders[img_select].IMAGE_HEIGHT;
        _zDim = (int) 1;
        _nDim = (int) 1;
    }

    // Map the parameters
    data.setDimensions(_xDim, _yDim, 1, _nDim);

    unsigned long   imgStart=0;
    unsigned long   imgEnd =_nDim;
    if (img_select != -1)
    {
        imgStart=img_select;
        imgEnd=img_select+1;
    }

    DataType datatype;
    //    dataHeaders[0].isSigned = false;
    int TIA_DT;
    if (img_select==-1)
    {
        TIA_DT = dataHeaders[0].DATA_TYPE;
        offset = header->pDATA_OFFSET[0] + TIAdataSIZE;
    }
    else
    {
        TIA_DT = dataHeaders[img_select].DATA_TYPE;
        offset = header->pDATA_OFFSET[img_select] + TIAdataSIZE;
    }

    switch ( TIA_DT )
    {
    case 1:
        datatype = UChar;
        break;
    case 2:
        datatype = UShort;
        //        datatype = Short;
        break;
    case 3:
        datatype = UInt;
        break;
    case 4:
        datatype = SChar;
        break;
    case 5:
        datatype = Short;
        //        dataHeaders[0].isSigned = true;
        break;
    case 6:
        datatype = Int;
        break;
    case 7:
        datatype = Float;
        break;
    case 8:
        datatype = Double;
        break;
    case 9:
        datatype = ComplexFloat;
        break;
    case 10:
        datatype = ComplexDouble;
        break;
    default:
        datatype = Unknown_Type;
        break;
    }

    MDMainHeader.removeObjects();
    MDMainHeader.setColumnFormat(false);
    MDMainHeader.addObject();
    MDMainHeader.setValue(MDL_SAMPLINGRATEX,(double)dataHeaders[0].PIXEL_WIDTH);
    MDMainHeader.setValue(MDL_SAMPLINGRATEY,(double)dataHeaders[0].PIXEL_HEIGHT);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    if( dataflag == -2 )
    {
        fclose(fimg);
        return 0;
    }

    MD.removeObjects();
    for ( i=imgStart; i<imgEnd; i++ )
        //for(int i=0;i< Ndim;i++)
    {
        MD.addObject();
        double aux;
        if(MDMainHeader.getValue(MDL_SAMPLINGRATEX,aux))
        {
            aux = ROUND(dataHeaders[i].CalibrationElementX - \
                        dataHeaders[i].CalibrationOffsetX/aux - data.xdim/2);
            MD.setValue(MDL_ORIGINX, aux);
        }
        if(MDMainHeader.getValue(MDL_SAMPLINGRATEY,aux))
        {
            aux = ROUND(dataHeaders[i].CalibrationElementY - \
                        dataHeaders[i].CalibrationOffsetY/aux -data.ydim/2);
            MD.setValue(MDL_ORIGINY, aux);
        }
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

    delete header;
    size_t pad = TIAdataSIZE;

    readData(fimg, img_select, datatype, pad);

    if (dataflag == 1)
    {
        if (dStddev == NULL)
            dStddev = 5;

        double temp, avg, stddev;
        double size = YXSIZE(data);

        avg = 0;
        stddev = 0;

        for ( int n=imgStart; n<imgEnd; n++ )
        {
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(data)
            {
                temp = abs(DIRECT_NZYX_ELEM(data,n,0,i,j));
                avg += temp;
                stddev += temp * temp;
            }
            avg /= size;
            stddev = stddev/size - avg * avg;
            stddev *= size/(size -1);
            stddev = sqrt(stddev);

            double low  = (avg - dStddev * stddev);
            double high = (avg + dStddev * stddev);

            FOR_ALL_ELEMENTS_IN_ARRAY3D(data)
            {
                if (abs(DIRECT_NZYX_ELEM(data,n,0,i,j)) < low)
                    DIRECT_NZYX_ELEM(data,n,0,i,j) = (T) low;
                else if (abs(DIRECT_NZYX_ELEM(data,n,0,i,j)) > high)
                    DIRECT_NZYX_ELEM(data,n,0,i,j) = (T) high;
            }
        }
    }

    if ( !mmapOn )
    	fclose(fimg);

    return(0);
}

/** TIA Writer
  * @ingroup TIA
*/
int writeTIA(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE)
{
    REPORT_ERROR(6001, "ERROR: writeTIA is not implemented.");
    return(-1);
}


#endif /* RWTIA_H_ */
