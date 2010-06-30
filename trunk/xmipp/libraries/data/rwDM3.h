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
    int nIm;
};


struct DM3dataHead
{
    double      CalibrationOffsetX;  //CalibrationOffsetX
    double      pixelWidth;          //CalibrationDeltaX
    int         CalibrationElementX;    //CalibrationElementX
    double      CalibrationOffsetY;   //CalibrationOffsetY
    double      pixelHeight;           //CalibrationDeltaY
    int          CalibrationElementY;    //CalibrationElementY
    short int   dataType;     //DataType
    int         imageWidth;            //ArraySizeX
    int         imageHeight;           //ArraySizeY
    short int   dataTypeSize;
    int          headerSize;
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

    header->tags.addLabel(MDL_DM3_NODEID);
    header->tags.addLabel(MDL_DM3_PARENTID);
    header->tags.addLabel(MDL_DM3_IDTAG);
    header->tags.addLabel(MDL_DM3_TAGNAME);
    header->tags.addLabel(MDL_DM3_TAGCLASS);
    header->tags.addLabel(MDL_DM3_SIZE);
    header->tags.addLabel(MDL_DM3_NUMBER_TYPE);
    header->tags.addLabel(MDL_DM3_VALUE);


    int nodeID=0, parentID=0, imCount=0, imCountF=0;

    for (int n=1;n<=header->nTags;n++)
        readTagDM3(fimg, header, parentID, nodeID, isLE);

    header->tags.write("images.txt");


    std::vector<long int> vIm;
    header->tags.findObjects(vIm, MDValueEqual(MDL_DM3_TAGNAME,(std::string)"DataType"));

    //    int parentId =parentDM3(header->tags, vIm[1], 2);


    header->nIm = 0;
    std::vector<DM3dataHead> dataHeaders;

    int value;

    for (int n = 0; n < vIm.size(); n++)
    {
        header->tags.goToObject(vIm[n]);
        header->tags.getValue(MDL_DM3_VALUE, value);

        if (value != 23)
        {
            parentID = parentDM3(header->tags, vIm[n], 2);

            nodeID = gotoTagDM3(header->tags, parentID, "ImageData Data");
            header->tags.getValue(MDL_DM3_NUMBER_TYPE, value);
            dataHeaders[header->nIm].dataType = (short int) value;

            header->tags.getValue(MDL_DM3_VALUE, value);
            dataHeaders[header->nIm].headerSize = (int) value;

            nodeID = gotoTagDM3(header->tags, parentID, "ImageData Dimensions");
            header->tags.nextObject();
            header->tags.getValue(MDL_DM3_VALUE, value);
            dataHeaders[header->nIm].imageHeight = (int) value;

            header->tags.nextObject();
            header->tags.getValue(MDL_DM3_VALUE, value);
            dataHeaders[header->nIm].imageWidth = (int) value;

            header->tags.nextObject();
            header->tags.getValue(MDL_DM3_VALUE, value);
            dataHeaders[header->nIm].dataTypeSize = (short int) value;

            nodeID = gotoTagDM3(header->tags, parentID, "ImageTags Acquisition Frame CCD Dimensions");
            header->tags.nextObject();
            header->tags.getValue(MDL_DM3_VALUE, value);
            dataHeaders[header->nIm].imageHeight = (int) value;














//            xmippFREAD(&(dataHeaders[i].CalibrationOffsetX), sizeof(double), 1, fimg, swap);
//            xmippFREAD(&dataHeaders[i].pixelWidth, sizeof(double), 1, fimg, swap);
//            xmippFREAD(&dataHeaders[i].CalibrationElementX, sizeof(int), 1, fimg, swap);
//            xmippFREAD(&dataHeaders[i].CalibrationOffsetY, sizeof(double), 1, fimg, swap);
//            xmippFREAD(&dataHeaders[i].pixelHeight, sizeof(double), 1, fimg, swap);
//            xmippFREAD(&dataHeaders[i].CalibrationElementY, sizeof(int), 1, fimg, swap);
//            xmippFREAD(&dataHeaders[i].dataType, sizeof(short int), 1, fimg, swap);
//            xmippFREAD(&dataHeaders[i].imageWidth, sizeof(int), 1, fimg, swap);
//            xmippFREAD(&dataHeaders[i].imageHeight, sizeof(int), 1, fimg, swap);
        }

    }



    //    nodeDim = gotoTagDM3(MD,parentId, "ImageData Dimensions");





    //    printDM3(header->tags);


    std::string names("perro gato raton");

    std::istringstream iss(names);
    std::vector<std::string> vSnames;
    //    splitString(names, " ", vSnames, false);

    std::string subs;
    while (iss >> subs)
    {
        std::cout << subs <<std::endl;

    }

    //    for (int i = 0; i < vSnames.size(); i++)
    // {
    //     std::cout << vSnames[i] <<std::endl;
    // }


    //    header->NUMBER_IMAGES = (int) dtemp;






    /*

        // Check data type
        if (header->DATA_TYPE_ID != 16674)
            REPORT_ERROR(6001, "ERROR: readTIA only processes images in real space");






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
                  int parentId,
                  int &nodeId,
                  bool isLE)
{


    /* Header Tag ============================================================== */

    unsigned char cdTag;
    int  idTag;
    unsigned short int ltName;
    xmippFREAD(&cdTag,sizeof (unsigned char),1,fimg,false); // Identification tag: 20 = tag dir,  21 = tag
    xmippFREAD(&ltName,sizeof(unsigned short int), 1,fimg,isLE); // Length of the tag name
    idTag = int(cdTag);

    char * tagName;
    std::string stagName;

    tagName =  new char[ltName+1];
    xmippFREAD(tagName,ltName,1,fimg,false); // Tag name
    tagName[ltName] = '\0';

    stagName = tagName;

    std::cout << tagName <<std::endl;

    header->tags.addObject();

    nodeId++;
    header->tags.setValue(MDL_DM3_NODEID, nodeId);
    header->tags.setValue(MDL_DM3_PARENTID, parentId);
    header->tags.setValue(MDL_DM3_IDTAG, idTag);
    header->tags.setValue(MDL_DM3_TAGNAME, stagName);



    //    for (int n=1;n<=depLevel;n++)
    //  printf("%d.",index[n-1]);

    /* Reading tags ===================================================================*/
    if (idTag == 20)  // Tag directory
    {
        //  printf("- Dir: %s\n",tagName);
        unsigned char dummy;
        int nTags;
        xmippFREAD(&dummy,sizeof(unsigned char),1,fimg,false); // 1 = sorted (normally = 1)
        xmippFREAD(&dummy,sizeof(unsigned char),1,fimg,false); //  0 = closed, 1 = open (normally = 0)
        xmippFREAD(&nTags,sizeof(int),1,fimg,isLE);             //  number of tags in tag directory


        header->tags.setValue(MDL_DM3_TAGCLASS,(std::string) "Dir");
        header->tags.setValue(MDL_DM3_SIZE, nTags);



        if (strcmp(tagName,"ImageList")==0)    // Number of images
        {

        }
        else if (strcmp(tagName,"Dimensions")==0)
        {
            //            newIndex[depLevel] = 1;
            //            header->tags.setValue(MDL_DIMY,(int) readTagDM3(fimg, header, depLevel, isLE));

            //            newIndex[depLevel] = 2;
            //            header->tags.setValue(MDL_DIMX,(int) readTagDM3(fimg, header, depLevel, isLE));

            //            header->tags.addObject();
            //            return 0;
        }

        parentId = nodeId;
        for (int n=1;n<=nTags;n++) // Rest of directories
        {
            readTagDM3(fimg, header, parentId, nodeId, isLE);
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

        int * info;
        info = new int[nnum];
        xmippFREAD(info,sizeof(unsigned int),nnum,fimg,isLE); // Reading of Info


        /* Tag classification  =======================================*/

        if (nnum == 1)   // Single entry tag
        {
            header->tags.setValue(MDL_DM3_TAGCLASS,(std::string) "Single");


            double tagValue = 0;

            FREADTagValueDM3(&tagValue,info[0],1,fimg);
            //   printf(" = %s\n", sprintfTagValue(&tagValue,info[0]));
            std::cout << "double " << tagValue <<std::endl;
            float fcaca;
            memcpy(&fcaca, &tagValue, sizeof(short));
            std::cout << "float " << fcaca <<std::endl;
            double caca2 = (double ) fcaca;
            std::cout << "double float " << caca2 <<std::endl;
            short scaca;
            memcpy(&scaca, &tagValue, sizeof(short));
            std::cout << scaca <<std::endl;
            caca2 = (double ) scaca;
            std::cout << "double short " << caca2 <<std::endl;
            int caca;
            memcpy(&caca, &tagValue, sizeof(int));
            std::cout << "int " << caca <<std::endl;
            caca2 = (double) caca;
            std::cout << "double int " << caca2 <<std::endl;

            std::vector<double> vtagValue(1);

            vtagValue.assign(1,tagValue);

            header->tags.setValue(MDL_DM3_NUMBER_TYPE, info[0]);
            header->tags.setValue(MDL_DM3_VALUE, vtagValue);


            if (strcmp(tagName,"DataType")==0)
            {
                //                header->tags.setValue(MDL_DATATYPE,(int) tagValue);
            }

            return tagValue;
        }
        else if(nnum == 3 && info[0]==20)   // Tag array
        {             /*nnum = 3
                                                                                                                                                                     info(0) = 20
                                                                                                                                                                     info(1) = number type for all values
                                                                                                                                                                     info(2) = info(nnum) = size of array*/

            header->tags.setValue(MDL_DM3_TAGCLASS,(std::string) "Array");

            header->tags.setValue(MDL_DM3_NUMBER_TYPE, info[1]);
            header->tags.setValue(MDL_DM3_SIZE, info[nnum-1]);
            std::vector<double> vtagValue(1);
            vtagValue.assign(1,(double) ftell(fimg));

            header->tags.setValue(MDL_DM3_VALUE, vtagValue);

            if (strcmp(tagName,"Data")==0)    // Locating the image data
            {

            }

            // Jump the array values
            int k;
            if(info[1] == 2 || info[1] == 4)
                k = 2;
            else if(info[1] == 3 || info[1] == 5)
                k = 4;
            else if(info[1] == 10 )
                k = 1;

            fseek( fimg, ftell(fimg)+(info[nnum-1])*k , SEEK_SET );

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


            header->tags.setValue(MDL_DM3_TAGCLASS, (std::string) "GroupArray");
            header->tags.setValue(MDL_DM3_SIZE, info[3]);


            int nBytes=0, k;

            double fieldValue;


            for (int n=1;n<=info[3];n++)
            {
                fieldValue=0;

                FREADTagValueDM3(&fieldValue,info[3+2*n],1,fimg);
                std::cout << fieldValue <<std::endl;
                float fcaca;
                memcpy(&fcaca, &fieldValue, sizeof(short));
                std::cout << fcaca <<std::endl;
                short scaca;
                memcpy(&scaca, &fieldValue, sizeof(short));
                std::cout << scaca <<std::endl;
                int caca;
                memcpy(&caca, &fieldValue, sizeof(int));
                std::cout << caca <<std::endl;

                if(info[3+2*n] == 2 || info[3+2*n] == 4)
                    k = 2;
                else if(info[3+2*n] == 3 || info[3+2*n] == 5)
                    k = 4;
                else if(info[3+2*n] == 10 )
                    k = 1;
                nBytes+=k;
            }

            // Jump the array values

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


            header->tags.setValue(MDL_DM3_TAGCLASS, (std::string) "Group");
            header->tags.setValue(MDL_DM3_SIZE, info[2]);

            std::vector<double> vtagValue(info[2]);

            for (int n=1;n<=info[2];n++)
            {
                double fieldValue=0;

                FREADTagValueDM3(&fieldValue,info[2+2*n],1,fimg);

                vtagValue.assign(n,fieldValue);

                std::cout << fieldValue <<std::endl;
                float fcaca;
                memcpy(&fcaca, &fieldValue, sizeof(short));
                std::cout << fcaca <<std::endl;
                short scaca;
                memcpy(&scaca, &fieldValue, sizeof(short));
                std::cout << scaca <<std::endl;
                int caca;
                memcpy(&caca, &fieldValue, sizeof(int));
                std::cout << caca <<std::endl;
            }

            header->tags.setValue(MDL_DM3_VALUE, vtagValue);

            return 0;
        }
    }
}

void FREADTagValueDM3(double *fieldValue,
                      int numberType,
                      int n,
                      FILE* fimg)
{

    DataType datatype = datatypeDM3(numberType);
    size_t datatypesize=gettypesize(datatype);

    xmippFREAD(fieldValue, datatypesize, n, fimg, swap);


    switch(numberType)
    {
    case 2:      // (02h =  2  i2* signed    (short)
        {
            datatype = Short;
            short* sValue = (short*) fieldValue;
            *fieldValue = (double) *sValue;
            break;
        }
    case 3:          // 03h =  3  i4* signed    (long)
        {
            datatype = Int;
            int* iValue = (int*) fieldValue;
            *fieldValue = (double) *iValue;
            break;
        }
    case 4:       //  04h =  4  i2* unsigned  (ushort) or unicode string
        {
            datatype = UShort;
            unsigned short* usValue = (unsigned short*) fieldValue;
            *fieldValue = (double) *usValue;
            break;
        }
    case 5:        //  05h =  5  i4* unsigned  (ulong)
        {
            datatype = UInt;
            unsigned int* uiValue = (unsigned int*) fieldValue;
            *fieldValue = (double) *uiValue;
            break;
        }
    case 6:        //  06h =  6  f4*           (float)
        {
            datatype = Float;
            float* fValue = (float*) fieldValue;
            *fieldValue = (double) *fValue;
            break;
        }
    case 7:        //  07h =  7  f8*           (double)
        {
            datatype = Double;
            //            double* caca = (double*) fieldValue;
            break;
        }
    case 8:        //  08h =  8  i1            (boolean)
        {
            datatype = Bool;
            bool* bValue = (bool*) fieldValue;
            *fieldValue = (double) *bValue;
            break;
        }
    case 9:        //  0ah = 10  i1
    case 10:        //  0ah = 10  i1
        {
            datatype = SChar;
            char* cValue = (char*) fieldValue;
            *fieldValue = (double) *cValue;
            break;
        }
    default:
        {
            datatype = Unknown_Type;
            break;
        }
    }
}


DataType datatypeDM3(int nType)
{
    DataType datatype;

    switch(nType)
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
    return datatype;
}

int parentDM3(MetaData &MD, int nodeId, int depth = 1)
{
    for (int n = 0; n < depth; n++)
    {
        MD.gotoFirstObject(MDValueEqual(MDL_DM3_NODEID,nodeId));
        MD.getValue(MDL_DM3_PARENTID,nodeId);

        if (nodeId == 0)
            break;
    }
    return nodeId;
}

double gotoTagDM3(MetaData &MD, int nodeId, std::string tagsList)
{
    std::istringstream iss(tagsList);
    std::string tag;

    MDMultiQuery queries;

    while (iss >> tag)
    {
        queries.clear();
        std::cout << queries.queryString <<std::endl;
        queries.addAndQuery(MDValueEqual(MDL_DM3_PARENTID,nodeId));
        std::cout << queries.queryString <<std::endl;
        queries.addAndQuery(MDValueEqual(MDL_DM3_TAGNAME,tag));
        std::cout << queries.queryString <<std::endl;

        MD.gotoFirstObject(queries);

        MD.getValue(MDL_DM3_NODEID,nodeId);
    }
    return nodeId;
}

int space;

void printDM3node(MetaData &MD, long int objId)
{
    MD.goToObject(objId);

    std::string tag;
    MD.getValue(MDL_DM3_TAGNAME,tag);

    int nodeId;
    MD.getValue(MDL_DM3_NODEID,nodeId);

    for (int i = 0; i < space; i++)
        std::cout << " ";

    std::cout << tag << std::endl;
    std::vector<long int> vObjs;

    space += 3;

    MD.findObjects(vObjs,MDValueEqual(MDL_DM3_PARENTID, nodeId));
    for (int i = 0; i < vObjs.size(); i++)
    {
        printDM3node(MD, vObjs[i]);
    }

    space -= 3;

}
void printDM3(MetaData MD)
{

    std::vector<long int> vObjs;

    space = 0;

    MD.findObjects(vObjs,MDValueEqual(MDL_DM3_PARENTID, 0));

    for (int i = 0; i < vObjs.size(); i++)
    {
        printDM3node(MD, vObjs[i]);
    }
}

#endif /* RWDM3_H_ */
