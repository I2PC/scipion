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
    size_t       headerSize;
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
        datatype = DT_Short;
        break;
    case 3:          // 03h =  3  i4* signed    (long)
        datatype = DT_Int;
        break;
    case 4:       //  04h =  4  i2* unsigned  (ushort) or unicode string
        datatype = DT_UShort;
        break;
    case 5:        //  05h =  5  i4* unsigned  (ulong)
        datatype = DT_UInt;
        break;
    case 6:        //  06h =  6  f4*           (float)
        datatype = DT_Float;
        break;
    case 7:        //  07h =  7  f8*           (double)
        datatype = DT_Double;
        break;
    case 8:        //  08h =  8  i1            (boolean)
        datatype = DT_Bool;
        break;
    case 9:        //  0ah = 10  i1
        datatype = DT_SChar;
        break;
    case 10:        //  0ah = 10  i1
        datatype = DT_SChar;
        break;
    default:
        datatype = DT_Unknown;
        break;
    }
    return datatype;
}

/** DM3 Low level reader
  * @ingroup DM3
*/
void FREADTagValueDM3(double *fieldValue,int numberType,int n,FILE* fimg, bool swap)
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

/** DM3 Tag reader
  * @ingroup DM3
*/
double readTagDM3(FILE *fimg, DM3head *header, int parentId, int &nodeId, bool isLE, bool swap)
{
    /* Header Tag ============================================================== */
    unsigned char cdTag;
    unsigned short int ltName;

    xmippFREAD(&cdTag,sizeof (unsigned char),1,fimg,false); // Identification tag: 20 = tag dir,  21 = tag
    xmippFREAD(&ltName,sizeof(unsigned short int), 1,fimg,isLE); // Length of the tag name
    int idTag = int(cdTag);

    std::string stagName;

    char * tagName =  new char[ltName+1];
    xmippFREAD(tagName,1,ltName,fimg,false); // Tag name
    tagName[ltName] = '\0';
    stagName = tagName;
    delete [] tagName;

    size_t id = header->tags.addObject();

    nodeId++;
    header->tags.setValue(MDL_DM3_NODEID, nodeId, id);
    header->tags.setValue(MDL_DM3_PARENTID, parentId, id);
    header->tags.setValue(MDL_DM3_IDTAG, idTag, id);
    header->tags.setValue(MDL_DM3_TAGNAME, stagName, id);

    /* Reading tags ===================================================================*/
    if (idTag == 20)  // Tag directory
    {
        unsigned char dummy;
        int nTags;
        xmippFREAD(&dummy,sizeof(unsigned char),1,fimg,false); // 1 = sorted (normally = 1)
        xmippFREAD(&dummy,sizeof(unsigned char),1,fimg,false); //  0 = closed, 1 = open (normally = 0)
        xmippFREAD(&nTags,sizeof(int),1,fimg,isLE);             //  number of tags in tag directory


        header->tags.setValue(MDL_DM3_TAGCLASS,(std::string) "Dir", id);
        header->tags.setValue(MDL_DM3_SIZE, nTags, id);

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

        xmippFREAD(&buf,1,4,fimg,false); // To read %%%% symbols
        xmippFREAD(&nnum,sizeof(unsigned int),1,fimg,isLE); // Size of info array

        int * info;
        info = new int[nnum];
        xmippFREAD(info,sizeof(unsigned int),nnum,fimg,isLE); // Reading of Info

        /* Tag classification  =======================================*/

        if (nnum == 1)   // Single entry tag
        {
            header->tags.setValue(MDL_DM3_TAGCLASS,(std::string) "Single", id);

            double tagValue = 0;

            FREADTagValueDM3(&tagValue,info[0],1,fimg,swap);

            std::vector<double> vtagValue(1);
            vtagValue.assign(1,tagValue);

            header->tags.setValue(MDL_DM3_NUMBER_TYPE, info[0], id);
            header->tags.setValue(MDL_DM3_VALUE, vtagValue, id);
            delete []info;
            return tagValue;
        }
        else if(nnum == 3 && info[0]==20)   // Tag array
        {
            /*nnum = 3
            info(0) = 20
            info(1) = number type for all values
            info(2) = info(nnum) = size of array*/

            header->tags.setValue(MDL_DM3_TAGCLASS,(std::string) "Array", id);

            header->tags.setValue(MDL_DM3_NUMBER_TYPE, info[1], id);
            header->tags.setValue(MDL_DM3_SIZE, info[nnum-1], id);
            std::vector<double> vtagValue(1);
            vtagValue.assign(1,(double) ftell(fimg));

            header->tags.setValue(MDL_DM3_VALUE, vtagValue, id);

            // Jump the array values
            int k;
            if(info[1] == 2 || info[1] == 4)
                k = 2;
            else if(info[1] == 3 || info[1] == 5)
                k = 4;
            else if(info[1] == 10 )
                k = 1;
            else
                REPORT_ERROR(ERR_ARG_INCORRECT,"");

            fseek( fimg, ftell(fimg)+(info[nnum-1])*k , SEEK_SET );
            delete []info;
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

            header->tags.setValue(MDL_DM3_TAGCLASS, (std::string) "GroupArray", id);
            header->tags.setValue(MDL_DM3_SIZE, info[3], id);
            int nBytes=0, k;
            double fieldValue;
            for (int n=1;n<=info[3];n++)
            {
                fieldValue=0;
                FREADTagValueDM3(&fieldValue,info[3+2*n],1,fimg,swap);

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
            delete []info;
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

            header->tags.setValue(MDL_DM3_TAGCLASS, (std::string) "Group", id);
            header->tags.setValue(MDL_DM3_SIZE, info[2], id);
            std::vector<double> vtagValue(info[2]);
            for (int n=1;n<=info[2];n++)
            {
                double fieldValue=0;
                FREADTagValueDM3(&fieldValue,info[2+2*n],1,fimg,swap);
                vtagValue.assign(n,fieldValue);
            }
            header->tags.setValue(MDL_DM3_VALUE, vtagValue, id);
            delete []info;
            return 0;
        }
        delete []info;
    }
    return 0;
}

/** DM3 Get DM3 parent
  * @ingroup DM3
*/
int parentDM3(MetaData &MD, int nodeId, int depth = 1)
{
    for (int n = 0; n < depth; n++)
    {
        MD.getValue(MDL_DM3_PARENTID, nodeId, MD.firstObject(MDValueEQ(MDL_DM3_NODEID, nodeId)));

        if (nodeId == 0)
            break;
    }
    return nodeId;
}

/** DM3 Go to tag
  * @ingroup DM3
*/
size_t gotoTagDM3(MetaData &MD, int &nodeId, const std::string &tagsList)
{
    std::string tag;
    std::vector<std::string> vTags;
    splitString(tagsList,",",vTags, false);
    size_t id=0;

    MDValueEQ queryParentId(MDL_DM3_PARENTID,-1), queryTagname(MDL_DM3_TAGNAME,tag);
    MDMultiQuery queries;

    queries.addAndQuery(queryParentId);
    queries.addAndQuery(queryTagname);

    for (size_t n = 0; n < vTags.size(); n++)
    {
        tag = vTags[n];

        queryParentId.setValue(nodeId);
        queryTagname.setValue(tag);
        id = MD.firstObject(queries);
        MD.getValue(MDL_DM3_NODEID, nodeId, id);
    }
    return id;
}

int space;
/** DM3 Print DM3 node
  * @ingroup DM3
*/
void printDM3node(MetaData &MD, size_t id)
{
    std::string tag;
    MD.getValue(MDL_DM3_TAGNAME, tag, id);

    int nodeId;
    MD.getValue(MDL_DM3_NODEID, nodeId, id);

    for (int i = 0; i < space; i++)
        std::cout << " ";

    std::cout << tag << std::endl;

    std::vector<size_t> vObjs;
    MD.findObjects(vObjs, MDValueEQ(MDL_DM3_PARENTID, nodeId));

    space += 3;

    for (size_t i = 0; i < vObjs.size(); i++)
        printDM3node(MD, vObjs[i]);

    space -= 3;

}

/** DM3 Print DM3 header
  * @ingroup DM3
*/
void printDM3(MetaData MD)
{
    std::vector<size_t> vObjs;
    space = 0;
    MD.findObjects(vObjs,MDValueEQ(MDL_DM3_PARENTID, 0));

    for (size_t i = 0; i < vObjs.size(); i++)
        printDM3node(MD, vObjs[i]);
}

/** DM3 Reader
  * @ingroup DM3
*/
int ImageBase::readDM3(size_t select_img,bool isStack)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readDM3: Reading DM3 file\n");
#endif

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

    xmippFREAD(&header->sorted, sizeof(char), 1, fimg, false);
    xmippFREAD(&header->open, sizeof(char), 1, fimg, false);
    xmippFREAD(&header->nTags, sizeof(int), 1, fimg, isLE);

    header->tags.addLabel(MDL_DM3_NODEID);
    header->tags.addLabel(MDL_DM3_PARENTID);
    header->tags.addLabel(MDL_DM3_IDTAG);
    header->tags.addLabel(MDL_DM3_TAGNAME);
    header->tags.addLabel(MDL_DM3_TAGCLASS);
    header->tags.addLabel(MDL_DM3_SIZE);
    header->tags.addLabel(MDL_DM3_NUMBER_TYPE);
    header->tags.addLabel(MDL_DM3_VALUE);

    int nodeID = 0, parentID = 0;

    for (int n=1;n<=header->nTags;n++)
        readTagDM3(fimg, header, parentID, nodeID, isLE, swap);

    //#define DEBUG
#ifdef DEBUG

    header->tags.write("images.txt");
    printDM3(header->tags);
#endif


    //std::vector<size_t> vIm;
    //header->tags.findObjects(vIm, );

    header->nIm = 0;
    std::vector<DM3dataHead> dataHeaders;
    DM3dataHead dhRef;

    std::vector<double> vValue;
    int iValue;
    size_t id;
    //Initialize query for later use
    MDValueEQ queryNodeId(MDL_DM3_NODEID, -1);

    for (MDIterator iter(header->tags, MDValueEQ(MDL_DM3_TAGNAME,(String)"DataType")); iter.hasNext(); iter.moveNext())
        // Read all the image headers
        //for (int n = 0; n < vIm.size(); n++)
    {
        //header->tags.goToObject(vIm[n]);
        header->tags.getValue(MDL_DM3_VALUE, vValue, iter.objId);

        if (vValue[0] != 23) //avoid thumb images
        {
            dataHeaders.push_back(dhRef);

            parentID = parentDM3(header->tags, iter.objId, 2);

            nodeID = parentID;
            id = gotoTagDM3(header->tags, nodeID, "ImageData,Data");
            header->tags.getValue(MDL_DM3_NUMBER_TYPE, iValue, id);
            dataHeaders[header->nIm].dataType = (short int) iValue;

            header->tags.getValue(MDL_DM3_VALUE, vValue, id);
            dataHeaders[header->nIm].headerSize = (size_t) vValue[0];

            nodeID = parentID;
            id = gotoTagDM3(header->tags, nodeID, "ImageData,Dimensions");
            nodeID++;
            queryNodeId.setValue(nodeID);
            id = header->tags.firstObject(queryNodeId);
            header->tags.getValue(MDL_DM3_VALUE, vValue, id);
            dataHeaders[header->nIm].imageWidth = (int) vValue[0];

            nodeID++;
            queryNodeId.setValue(nodeID);
            id = header->tags.firstObject(queryNodeId);
            header->tags.getValue(MDL_DM3_VALUE, vValue, id);
            dataHeaders[header->nIm].imageHeight = (int) vValue[0];

            nodeID++;
            queryNodeId.setValue(nodeID);
            id = header->tags.firstObject(queryNodeId);
            header->tags.getValue(MDL_DM3_VALUE, vValue, id);
            dataHeaders[header->nIm].dataTypeSize = (short int) vValue[0];

            nodeID = parentID;
            id = gotoTagDM3(header->tags, nodeID, "ImageTags,Acquisition,Frame,CCD,Pixel Size (um)");
            //            header->tags.nextObject();
            header->tags.getValue(MDL_DM3_VALUE, vValue, id);
            dataHeaders[header->nIm].pixelHeight = vValue[0]*1e4;
            dataHeaders[header->nIm].pixelWidth  = vValue[1]*1e4;

            //TODO: Do I have to include FLIP?!?!? which? vertical or horizontal?
            header->nIm++;
        }

    }

    if (dataHeaders.size() == 0)
        REPORT_ERROR(ERR_IMG_NOREAD,formatString("readDM3: Image information not found in file %s",filename.c_str()));

    int _xDim,_yDim;
    size_t _nDim;
    _xDim = dataHeaders[0].imageWidth;
    _yDim = dataHeaders[0].imageHeight;
    _nDim = header->nIm;

    // Map the parameters
    if (select_img >  _nDim)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, formatString("readDM3: Image number %lu exceeds stack size %lu", select_img, _nDim));
    else if (select_img == ALL_IMAGES)
    {
        // Check images dimensions. Need to be the same
        for (size_t i = 1; i < _nDim ; i++)
        {
            if (dataHeaders[0].imageHeight != dataHeaders[i].imageHeight || \
                dataHeaders[0].imageWidth != dataHeaders[i].imageWidth  || \
                dataHeaders[0].dataType != dataHeaders[i].dataType)
                REPORT_ERROR(ERR_IMG_NOREAD, "readDM3: images in DM3 file with different \
                             dimensions and data types are not currently supported. Try to read them individually.");
        }
        //FIXME: Code is not totally implemented to load automatically multiple images if they have same size.
        if (_nDim > 1)
            REPORT_ERROR(ERR_IO_NOREAD, "readDM3: Reading multiple \
                         images at once in DM3 file are not currently supported. Try to read them individually.");
    }
    else
        _nDim = 1;

    setDimensions(_xDim, _yDim, 1, _nDim);

    size_t   imgStart = IMG_INDEX(select_img);
    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : _nDim;

    DataType datatype = datatypeDM3(dataHeaders[0].dataType);

    MDMainHeader.setValue(MDL_SAMPLINGRATE_X,(double)dataHeaders[0].pixelWidth);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_Y,(double)dataHeaders[0].pixelHeight);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    if (dataMode == HEADER) // Stop reading if not necessary
    {
        delete header;
        return 0;
    }

    MD.clear();
    MD.resize(imgEnd - imgStart,MDL::emptyHeader);
    offset = dataHeaders[imgStart].headerSize;
    delete header;

    if( dataMode < DATA )
        return 0;

#undef DEBUG
#ifdef DEBUG

    MDMainHeader.write(std::cerr);
    MD.write(std::cerr);
#endif

    size_t pad = 0;
    readData(fimg, select_img, datatype, pad);

    return(0);
}
