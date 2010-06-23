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

#ifndef RWDM3_H_
#define RWDM3_H_


struct DM3head
{
    int fileVersion;
    int fileLength;
    int byteOrder;
    char sorted;
    char open;
    int nTags;
    MetaData tags;
};


int readDM3(int img_select,bool isStack=false)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readDM3: Reading DM3 file\n");
#endif

    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "r") ) == NULL )
        return(-1);

    DM3head * header = new DM3head;
    int dummy;

    // Check Machine endianess
    bool isLE = IsLittleEndian();

    xmippFREAD(&header->fileVersion, sizeof(int), 1, fimg, isLE);
    xmippFREAD(&dummy, sizeof(int), 1, fimg, isLE);
    xmippFREAD(&header->byteOrder, sizeof(int), 1, fimg, isLE);

    // Set Endianess
    swap = (isLE^header->byteOrder);

    if ( header->fileVersion!=3 )
        REPORT_ERROR(6001, "readDM3: Input file is not Digital Micrograph 3 format.");

    xmippFREAD(&header->sorted, sizeof(char), 1, fimg, false);
    xmippFREAD(&header->open, sizeof(char), 1, fimg, false);
    xmippFREAD(&header->nTags, sizeof(int), 1, fimg, isLE);


    //FIXME: I am here

    (header->tags).clear();

    int depLevel=0, imCount=0, imCountF=0;

    for (int n=1;n<=header->nTags;n++)
        readTagDM3(fimg, header, depLevel, isLE);

    header->tags.write("images.txt");



    /*

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
        for (i = 1; i < header->NUMBER_IMAGES; i++)
        {
            if (dataHeaders[0].IMAGE_HEIGHT != dataHeaders[i].IMAGE_HEIGHT || \
                dataHeaders[0].IMAGE_WIDTH != dataHeaders[i].IMAGE_WIDTH  || \
                dataHeaders[0].DATA_TYPE != dataHeaders[i].DATA_TYPE)
                REPORT_ERROR(6001, "readTIA: images in TIA file with different dimensions and data types are not supported");
        }


        int _xDim,_yDim,_zDim;
        unsigned long int _nDim;
        _xDim = (int) dataHeaders[0].IMAGE_WIDTH;
        _yDim = (int) dataHeaders[0].IMAGE_HEIGHT;
        _zDim = (int) 1;
        _nDim = (int) header->NUMBER_IMAGES;



        // Map the parameters
        if (img_select==-1)
            data.setDimensions(_xDim, _yDim, 1, _nDim);
        else
            data.setDimensions(_xDim, _yDim, 1, 1);

        unsigned long   imgStart=0;
        unsigned long   imgEnd =_nDim;
        if (img_select != -1)
        {
            imgStart=img_select;
            imgEnd=img_select+1;
        }

        DataType datatype;
        dataHeaders[0].isSigned = false;
        switch ( dataHeaders[0].DATA_TYPE )
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
            dataHeaders[0].isSigned = true;
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

        offset = header->pDATA_OFFSET[0] + TIAdataSIZE;
        size_t pad = TIAdataSIZE;


        //#define DEBUG
    #ifdef DEBUG

        MDMainHeader.write(std::cerr);
        MD.write(std::cerr);
    #endif

        delete header;
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

        fclose(fimg);

        return(0);*/
}


double readTagDM3(FILE *fimg,
                  DM3head *header,
                  int depLevel,
                  bool isLE)
{
    depLevel++;

    /* Header Tag ============================================================== */

    unsigned char cdTag;
    unsigned int  idTag;
    unsigned short int ltName;
    xmippFREAD(&cdTag,sizeof (unsigned char),1,fimg,false); // Identification tag: 20 = tag dir,  21 = tag
    xmippFREAD(&ltName,sizeof(unsigned short int), 1,fimg,isLE); // Length of the tag name
    idTag = int(cdTag);

    char * tagName;
    tagName =  new char[ltName+1];
    xmippFREAD(tagName,ltName,1,fimg,false); // Tag name
    tagName[ltName] = '\0';



    //    for (int n=1;n<=depLevel;n++)
    //  printf("%d.",index[n-1]);

    /* Reading tags ===================================================================*/
    if (idTag == 20)  // Tag directory
    {
        //  printf("- Dir: %s\n",tagName);
        unsigned char dummy;
        unsigned int nTags;
        xmippFREAD(&dummy,sizeof(unsigned char),1,fimg,false); // 1 = sorted (normally = 1)
        xmippFREAD(&dummy,sizeof(unsigned char),1,fimg,false); //  0 = closed, 1 = open (normally = 0)
        xmippFREAD(&nTags,sizeof(int),1,fimg,isLE);             //  number of tags in tag directory


        if (strcmp(tagName,"ImageList")==0)    // Number of images
        {
            header->tags.addObject();
            //                header->tags.setValue(MDL_DM3_LEVEL, depLevel);
            //                header->tags.setValue(MDL_DM3_IDTAG, idTag);
            //                header->tags.setValue(MDL_DM3_TAGNAME, tagName);
        }
        else if (strcmp(tagName,"Dimensions")==0)
        {
            //            newIndex[depLevel] = 1;
            header->tags.setValue(MDL_DIMY,(int) readTagDM3(fimg, header, depLevel, isLE));

            //            newIndex[depLevel] = 2;
            header->tags.setValue(MDL_DIMX,(int) readTagDM3(fimg, header, depLevel, isLE));

            header->tags.addObject();
            return 0;
        }

        for (int n=1;n<=nTags;n++) // Rest of directories
        {
            readTagDM3(fimg, header, depLevel, isLE);
        }
        return 0;
    }
    else if (idTag == 21)    // Tag
    {
        //  printf("- Tag: %s ",tagName);

        unsigned int nnum;
        char buf[4]; // to read %%%% symbols

        xmippFREAD(&buf,1,4,fimg,false); // To read %%%% symbols
        xmippFREAD(&nnum,sizeof(unsigned int),1,fimg,isLE); // Size of info array

        unsigned int * info;
        info = new unsigned int[nnum];
        xmippFREAD(info,sizeof(unsigned int),nnum,fimg,isLE); // Reading of Info


        /* Tag classification  =======================================*/

        if (nnum == 1)   // Single entry tag
        {
            double tagValue = 0;

            FREADTagValueDM3(&tagValue,info[0],1,fimg);
            //   printf(" = %s\n", sprintfTagValue(&tagValue,info[0]));

            if (strcmp(tagName,"DataType")==0)
            {
                header->tags.setValue(MDL_DATATYPE,(int) tagValue);
            }

            return tagValue;
        }
        else if(nnum == 3 && info[0]==20)   // Tag array
        {             /*nnum = 3
                                                                                                                                                                                                                                                                                 info(0) = 20
                                                                                                                                                                                                                                                                                 info(1) = number type for all values
                                                                                                                                                                                                                                                                                 info(2) = info(nnum) = size of array*/

            if (strcmp(tagName,"Data")==0)    // Locating the image data
            {
                header->tags.setValue(MDL_HEADERSIZE, (int) ftell(fimg));
            }

            // Jump the array values
            int k;
            if(info[1] == 2 || info[1] == 4)
                k = 2;
            else if(info[1] == 3 || info[1] == 5)
                k = 4;
            else if(info[1] == 10 )
                k = 1;

            //   printf("Position before array = %d \n ", ftell(fh_in));
            //    printf("(Array size = %d) \n", info[nnum-1]*k);

            fseek( fimg, ftell(fimg)+(info[nnum-1])*k , SEEK_SET );

            //   printf("Position = %d \n ", ftell(fh_in));
            //  }else{
            //   int arrayValue[info[2]];
            //   memset(arrayValue,0,sizeof(arrayValue));
            //   FREADTagValue(&arrayValue,info[1],info[2],fh_in,bigEndian);
            //   for (int n=0;info[2];n++){
            //    printf("%d  ", arrayValue[n]);
            //   } std::cout << "\n\n";

            return 0;

        }
        else if (info[0]==20 && info[1] == 15)    // Tag Group array
        {
            /*nnum = size of array
                     info(0) = 20 (array)
                     info(1) = 15 (group)
                     info(2) = 0 (always 0)
                     info(3) = number of values in group
                     info(2*i+3) = number type for value i
                     info(nnum) = size of info array*/

            //   printf(" = (");

            int nBytes=0, k;

            double fieldValue;


            for (int n=1;n<=info[3];n++)
            {
                fieldValue=0;

                FREADTagValueDM3(&fieldValue,info[3+2*n],1,fimg);

                if(info[3+2*n] == 2 || info[3+2*n] == 4)
                    k = 2;
                else if(info[3+2*n] == 3 || info[3+2*n] == 5)
                    k = 4;
                else if(info[3+2*n] == 10 )
                    k = 1;
                nBytes+=k;
                //    printf( "%s", sprintfTagValue(&fieldValue,info[3+2*n]));
                //    if(n<info[3]) printf(","); else printf(")\n");
            }

            // Jump the array values

            //   printf("Reading of array values\n");
            //  printf("Position = %d \n ", ftell(fh_in));

            fseek( fimg, ftell(fimg)+(info[nnum-1]-1)*nBytes , SEEK_SET );
            return 0;
        }
        else if (info[0] == 15)    // Tag Group  (struct)
        {
            /* info[1] = length group name (always =0)
                    info[2] = number of entries in group */
            /*nnum = size of info array
            info(1) = 0fh
            info(3) = number of values in group
            info(2*i+3) = number type for value i
            Other info entries are always zero*/

            //   printf(" = [");

            for (int n=1;n<=info[2];n++)
            {
                double fieldValue=0;

                FREADTagValueDM3(&fieldValue,info[2+2*n],1,fimg);

                //    printf( "%s", sprintfTagValue(&fieldValue,info[2+2*n]));
                //    if(n<info[2]) printf(","); else printf("]\n");
            }
            //     printf("Position = %d \n ", ftell(fh_in));

            return 0;
        }
    }
}

void FREADTagValueDM3(double *fieldValue,
                      int numberType,
                      int n,
                      FILE* fimg)
{

    DataType datatype;


    switch(numberType)
    {
    case 2:      // (02h =  2  i2* signed    (short)
        datatype = Short;
        break;
    case 3:          // 03h =  3  i4* signed    (long)
        datatype = Int;
        break;
    case 4:       //  04h =  4  i2* unsigned  (ushort) or unicode string
        datatype = UShort;
        break;
    case 5:        //  05h =  5  i4* unsigned  (ulong)
        datatype = UInt;
        break;
    case 6:        //  06h =  6  f4*           (float)
        datatype = Float;
        break;
    case 7:        //  07h =  7  f8*           (double)
        datatype = Double;
        break;
    case 8:        //  08h =  8  i1            (boolean)
        datatype = Bool;
        break;
    case 9:        //  0ah = 10  i1
        datatype = SChar;
        break;
    case 10:        //  0ah = 10  i1
        datatype = SChar;
        break;
    default:
        datatype = Unknown_Type;
        break;
    }

    size_t datatypesize=gettypesize(datatype);

//    char* page = new char [datatypesize];

    char* page = (char*) fieldValue;

    xmippFREAD(page, datatypesize*n, 1, fimg, swap);
    //    memcpy(fieldValue, page, datatypesize*n);
//    fieldValue = (double *) page;
//    for(int i=0; i<n;i++)
//        fieldValue[i]=(double) ptr[i];
    //    castPage2T(page, fieldValue, datatype, datatypesize);
}



#endif /* RWDM3_H_ */
