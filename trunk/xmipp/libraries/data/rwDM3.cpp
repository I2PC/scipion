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

#include "image.h"

///@defgroup DM3 DM3 File format
///@ingroup ImageFormats

/** DM3 Header
  * @ingroup DM3
*/
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

/** DM3 Data Header
  * @ingroup DM3
*/
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
    bool   flip;
};

/** DM3 Determine DM3 datatype
  * @ingroup DM3
*/
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

/** DM3 Low level reader
  * @ingroup DM3
*/
void FREADTagValueDM3(double *fieldValue,int numberType,int n, FILE* fimg, int swap)
{
    DataType datatype = datatypeDM3(numberType);
    size_t datatypesize=gettypesize(datatype);

    xmippFREAD(fieldValue, datatypesize, n, fimg, swap);

    switch(numberType)
    {
    case 2:      // (02h =  2  i2* signed    (short)
        {
            short* sValue = (short*) fieldValue;
            *fieldValue = (double) *sValue;
            break;
        }
    case 3:          // 03h =  3  i4* signed    (long)
        {
            int* iValue = (int*) fieldValue;
            *fieldValue = (double) *iValue;
            break;
        }
    case 4:       //  04h =  4  i2* unsigned  (ushort) or unicode string
        {
            unsigned short* usValue = (unsigned short*) fieldValue;
            *fieldValue = (double) *usValue;
            break;
        }
    case 5:        //  05h =  5  i4* unsigned  (ulong)
        {
            unsigned int* uiValue = (unsigned int*) fieldValue;
            *fieldValue = (double) *uiValue;
            break;
        }
    case 6:        //  06h =  6  f4*           (float)
        {
            float* fValue = (float*) fieldValue;
            *fieldValue = (double) *fValue;
            break;
        }
    case 7:        //  07h =  7  f8*           (double)
        {
            //            double* caca = (double*) fieldValue;
            break;
        }
    case 8:        //  08h =  8  i1            (boolean)
        {
            bool* bValue = (bool*) fieldValue;
            *fieldValue = (double) *bValue;
            break;
        }
    case 9:        //  0ah = 10  i1
    case 10:        //  0ah = 10  i1
        {
            char* cValue = (char*) fieldValue;
            *fieldValue = (double) *cValue;
            break;
        }
    default:
        {
            break;
        }
    }
}

/** DM3 Get DM3 parent
  * @ingroup DM3
*/
int parentDM3(MetaData &MD, int nodeId, int depth = 1)
{
    for (int n = 0; n < depth; n++)
    {
        MD.gotoFirstObject(MDValueEQ(MDL_DM3_NODEID,nodeId));
        MD.getValue(MDL_DM3_PARENTID,nodeId);

        if (nodeId == 0)
            break;
    }
    return nodeId;
}

/** DM3 Go to tag
  * @ingroup DM3
*/
double gotoTagDM3(MetaData &MD, int nodeId, std::string tagsList)
{
    std::string tag;
    std::vector<std::string> vTags;
    splitString(tagsList,",",vTags, false);

    MDValueEQ queryParentId(MDL_DM3_PARENTID,-1), queryTagname(MDL_DM3_TAGNAME,tag);
    MDMultiQuery queries;

    queries.addAndQuery(queryParentId);
    queries.addAndQuery(queryTagname);

    for (int n = 0; n < vTags.size(); n++)
    {
        tag = vTags[n];

        queryParentId.setValue(nodeId);
        queryTagname.setValue(tag);

        MD.gotoFirstObject(queries);
        MD.getValue(MDL_DM3_NODEID,nodeId);
    }
    return nodeId;
}

int space;
/** DM3 Print DM3 node
  * @ingroup DM3
*/
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
    MD.findObjects(vObjs,MDValueEQ(MDL_DM3_PARENTID, nodeId));

    space += 3;

    for (int i = 0; i < vObjs.size(); i++)
        printDM3node(MD, vObjs[i]);

    space -= 3;

}

/** DM3 Print DM3 header
  * @ingroup DM3
*/
void printDM3(MetaData MD)
{
    std::vector<long int> vObjs;
    space = 0;
    MD.findObjects(vObjs,MDValueEQ(MDL_DM3_PARENTID, 0));

    for (int i = 0; i < vObjs.size(); i++)
        printDM3node(MD, vObjs[i]);
}

/** DM3 Tag reader
  * @ingroup DM3
*/
double readTagDM3(FILE *fimg, DM3head *header, int parentId, int &nodeId, bool isLE, int swap)
{
    /* Header Tag ============================================================== */
    int  idTag;
    unsigned char cdTag;
    unsigned short int ltName;

    xmippFREAD(&cdTag,sizeof (unsigned char),1,fimg,isLE); // Identification tag: 20 = tag dir,  21 = tag
    xmippFREAD(&ltName,sizeof(unsigned short int), 1,fimg,isLE); // Length of the tag name
    idTag = int(cdTag);

    char * tagName;
    std::string stagName;

    tagName =  new char[ltName+1];
    xmippFREAD(tagName,ltName,1,fimg,isLE); // Tag name
    tagName[ltName] = '\0';

    stagName = tagName;

    header->tags.addObject();

    nodeId++;
    header->tags.setValue(MDL_DM3_NODEID, nodeId);
    header->tags.setValue(MDL_DM3_PARENTID, parentId);
    header->tags.setValue(MDL_DM3_IDTAG, idTag);
    header->tags.setValue(MDL_DM3_TAGNAME, stagName);

    /* Reading tags ===================================================================*/
    if (idTag == 20)  // Tag directory
    {
        //  printf("- Dir: %s\n",tagName);
        unsigned char dummy;
        int nTags;
        xmippFREAD(&dummy,sizeof(unsigned char),1,fimg,isLE); // 1 = sorted (normally = 1)
        xmippFREAD(&dummy,sizeof(unsigned char),1,fimg,isLE); //  0 = closed, 1 = open (normally = 0)
        xmippFREAD(&nTags,sizeof(int),1,fimg,isLE);             //  number of tags in tag directory


        header->tags.setValue(MDL_DM3_TAGCLASS,(std::string) "Dir");
        header->tags.setValue(MDL_DM3_SIZE, nTags);

        parentId = nodeId;
        for (int n=1;n<=nTags;n++) // Rest of directories
        {
            readTagDM3(fimg, header, parentId, nodeId, isLE, swap);
        }
        return 0;
    }
    else if (idTag == 21)    // Tag
    {
        unsigned int nnum;
        char buf[4]; // to read %%%% symbols

        xmippFREAD(&buf,1,4,fimg,isLE); // To read %%%% symbols
        xmippFREAD(&nnum,sizeof(unsigned int),1,fimg,isLE); // Size of info array

        int * info;
        info = new int[nnum];
        xmippFREAD(info,sizeof(unsigned int),nnum,fimg,isLE); // Reading of Info

        /* Tag classification  =======================================*/

        if (nnum == 1)   // Single entry tag
        {
            header->tags.setValue(MDL_DM3_TAGCLASS,(std::string) "Single");

            double tagValue = 0;

            FREADTagValueDM3(&tagValue,info[0],1,fimg, swap);

            std::vector<double> vtagValue(1);
            vtagValue.assign(1,tagValue);

            header->tags.setValue(MDL_DM3_NUMBER_TYPE, info[0]);
            header->tags.setValue(MDL_DM3_VALUE, vtagValue);

            return tagValue;
        }
        else if(nnum == 3 && info[0]==20)   // Tag array
        {
            /*nnum = 3
            info(0) = 20
            info(1) = number type for all values
            info(2) = info(nnum) = size of array*/

            header->tags.setValue(MDL_DM3_TAGCLASS,(std::string) "Array");

            header->tags.setValue(MDL_DM3_NUMBER_TYPE, info[1]);
            header->tags.setValue(MDL_DM3_SIZE, info[nnum-1]);
            std::vector<double> vtagValue(1);
            vtagValue.assign(1,(double) ftell(fimg));

            header->tags.setValue(MDL_DM3_VALUE, vtagValue);

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
                FREADTagValueDM3(&fieldValue,info[3+2*n],1,fimg, swap);

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
                FREADTagValueDM3(&fieldValue,info[2+2*n],1,fimg, swap);
                vtagValue.assign(n,fieldValue);
            }
            header->tags.setValue(MDL_DM3_VALUE, vtagValue);
            return 0;
        }
    }
}

/** DM3 Reader
  * @ingroup DM3
*/
template<typename T>
int Image<T>::readDM3(int img_select,bool isStack)
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
        REPORT_ERROR(ERR_IO_NOREAD, "readDM3: Input file is not Digital Micrograph 3 format.");

    xmippFREAD(&header->sorted, sizeof(char), 1, fimg, isLE);
    xmippFREAD(&header->open, sizeof(char), 1, fimg, isLE);
    xmippFREAD(&header->nTags, sizeof(int), 1, fimg, isLE);

    header->tags.addLabel(MDL_DM3_NODEID);
    header->tags.addLabel(MDL_DM3_PARENTID);
    header->tags.addLabel(MDL_DM3_IDTAG);
    header->tags.addLabel(MDL_DM3_TAGNAME);
    header->tags.addLabel(MDL_DM3_TAGCLASS);
    header->tags.addLabel(MDL_DM3_SIZE);
    header->tags.addLabel(MDL_DM3_NUMBER_TYPE);
    header->tags.addLabel(MDL_DM3_VALUE);

    int nodeID=0, parentID=0;

    for (int n=1;n<=header->nTags;n++)
        readTagDM3(fimg, header, parentID, nodeID, isLE, swap);

    //#define DEBUG
#ifdef DEBUG

    header->tags.write("images.txt");
    printDM3(header->tags);
#endif


    std::vector<long int> vIm;
    header->tags.findObjects(vIm, MDValueEQ(MDL_DM3_TAGNAME,(std::string)"DataType"));

    header->nIm = 0;
    std::vector<DM3dataHead> dataHeaders;
    DM3dataHead dhRef;

    std::vector<double> vValue;
    int iValue;

    // Read all the image headers
    for (int n = 0; n < vIm.size(); n++)
    {
        header->tags.goToObject(vIm[n]);
        header->tags.getValue(MDL_DM3_VALUE, vValue);

        if (vValue[0] != 23) //avoid thumb images
        {
            dataHeaders.push_back(dhRef);

            parentID = parentDM3(header->tags, vIm[n], 2);

            nodeID = gotoTagDM3(header->tags, parentID, "ImageData,Data");
            header->tags.getValue(MDL_DM3_NUMBER_TYPE, iValue);
            dataHeaders[header->nIm].dataType = (short int) iValue;

            header->tags.getValue(MDL_DM3_VALUE, vValue);
            dataHeaders[header->nIm].headerSize = (int) vValue[0];

            nodeID = gotoTagDM3(header->tags, parentID, "ImageData,Dimensions");
            header->tags.nextObject();
            header->tags.getValue(MDL_DM3_VALUE, vValue);
            dataHeaders[header->nIm].imageHeight = (int) vValue[0];

            header->tags.nextObject();
            header->tags.getValue(MDL_DM3_VALUE, vValue);
            dataHeaders[header->nIm].imageWidth = (int) vValue[0];

            header->tags.nextObject();
            header->tags.getValue(MDL_DM3_VALUE, vValue);
            dataHeaders[header->nIm].dataTypeSize = (short int) vValue[0];

            nodeID = gotoTagDM3(header->tags, parentID, "ImageTags,Acquisition,Frame,CCD,Pixel Size (um)");
            //            header->tags.nextObject();
            header->tags.getValue(MDL_DM3_VALUE, vValue);
            dataHeaders[header->nIm].pixelHeight = vValue[0]*1e4;
            dataHeaders[header->nIm].pixelWidth  = vValue[1]*1e4;

            //TODO: Do I have to include FLIP?!?!? which? vertical or horizontal?

            header->nIm++;
        }

    }

    // Check images dimensions. Need to be the same
    for (int i = 1; i < header->nIm; i++)
    {
        if (dataHeaders[0].imageHeight != dataHeaders[i].imageHeight || \
            dataHeaders[0].imageWidth != dataHeaders[i].imageWidth  || \
            dataHeaders[0].dataType != dataHeaders[i].dataType)
            REPORT_ERROR(ERR_IMG_NOREAD, "readDM3: images in DM3 file with different \
                         dimensions and data types are not currently supported. Try to read them individually.");
    }


    int _xDim,_yDim;
    unsigned long int _nDim;
    _xDim = (int) dataHeaders[0].imageWidth;
    _yDim = (int) dataHeaders[0].imageHeight;
    _nDim = (int) header->nIm;

    //FIXME: Code is not totally implemented to load automatically multiple images if they have same size.

    // Map the parameters
    if (img_select==-1)
    {
        if (_nDim>1)
            REPORT_ERROR(ERR_IO_NOREAD, "readDM3: Reading multiple \
                         images at once in DM3 file are not currently supported. Try to read them individually.");
    }
    else
        _nDim=1;

    data.setDimensions(_xDim, _yDim, 1, _nDim);


    unsigned long   imgStart=0;
    unsigned long   imgEnd =_nDim;

    if (img_select != -1)
    {
        imgStart=img_select;
        imgEnd=img_select+1;
    }

    DataType datatype = datatypeDM3(dataHeaders[0].dataType);


    MDMainHeader.setValue(MDL_SAMPLINGRATEX,(double)dataHeaders[0].pixelWidth);
    MDMainHeader.setValue(MDL_SAMPLINGRATEY,(double)dataHeaders[0].pixelHeight);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    MD.clear();
    MD.resize(imgEnd - imgStart);
    for ( i = imgStart; i < imgEnd; i++ )
    {
        MD[i-imgStart].setValue(MDL_ORIGINX,  zeroD);
        MD[i-imgStart].setValue(MDL_ORIGINY,  zeroD);
        MD[i-imgStart].setValue(MDL_ORIGINZ,  zeroD);

        MD[i-imgStart].setValue(MDL_ANGLEROT, zeroD);
        MD[i-imgStart].setValue(MDL_ANGLETILT,zeroD);
        MD[i-imgStart].setValue(MDL_ANGLEPSI, zeroD);
        MD[i-imgStart].setValue(MDL_WEIGHT,   oneD);
        MD[i-imgStart].setValue(MDL_FLIP,     falseb);
    }

    offset = (unsigned long) dataHeaders[imgStart].headerSize;
    size_t pad = 0;
    delete header;

    if( dataflag < 0 )
    {
        fclose(fimg);
        return 0;
    }


    //#define DEBUG
#ifdef DEBUG

    MDMainHeader.write(std::cerr);
    MD.write(std::cerr);
#endif


    readData(fimg, img_select, datatype, pad);

    if ( !mmapOn )
        fclose(fimg);

    return(0);
}

/** DM3 Writer
  * @ingroup DM3
*/
template<typename T>
int Image<T>::writeDM3(int img_select, bool isStack, int mode)
{
    REPORT_ERROR(ERR_IO_NOWRITE, "ERROR: writeDM3 is not implemented.");
    return(-1);
}

