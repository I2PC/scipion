/*
 * metadata_extension.h
 *
 *  Created on: May 12, 2010
 *      Author: roberto
 */
#include "metadata_extension.h"

#ifndef __linux__
#define MAXDOUBLE __DBL_MAX__
#endif


/*----------   Statistics --------------------------------------- */
void getStatistics(MetaData &MT_in, Image<double> & _ave, Image<double> & _sd, double& _min,
                   double& _max, bool apply_geo)
{
    MetaData MT(MT_in); //copy constructor so original MT is not changed
    _min = MAXDOUBLE;
    _max = -MAXDOUBLE;
    bool first = true;
    int n = 0;
    //Remove disabled images if present
    if (MT.containsLabel(MDL_ENABLED))
        MT.removeObjects(MDValueEQ(MDL_ENABLED, -1));
    // Calculate Mean
    if (MT.isEmpty())
        REPORT_ERROR(ERR_MD_OBJECTNUMBER, "There is no selected images in Metadata.");

    Image<double> image, tmpImg;
    double min, max, avg, stddev;
    FileName fnImg;
    FOR_ALL_OBJECTS_IN_METADATA(MT)
    {
        if (apply_geo)
            image.readApplyGeo(MT,__iter.objId, HEADER);
        else
        {
            MT.getValue(MDL_IMAGE,fnImg,__iter.objId);
            image.read(fnImg);
        }
        image().computeStats(avg, stddev, min, max);
        if (min < _min)
            _min = min;
        if (max > _max)
            _max = max;
        if (first)
        {
            _ave = image;
            first = false;
        }
        else
        {
            _ave() += image();
        }
        n++;
    }

    if (n > 0)
        _ave() /= n;
    _sd = _ave;
    _sd().initZeros();
    // Calculate SD
    FOR_ALL_OBJECTS_IN_METADATA(MT)
    {
        if (apply_geo)
            image.readApplyGeo(MT,__iter.objId, HEADER);
        else
        {
            MT.getValue(MDL_IMAGE,fnImg,__iter.objId);
            image.read(fnImg);
        }
        tmpImg() = ((image() - _ave()));
        tmpImg() *= tmpImg();
        _sd() += tmpImg();
    }
    _sd() /= (n - 1);
    _sd().selfSQRT();
}

/*----------   Statistics --------------------------------------- */
void getStatistics(MetaData &MT_in, double& _ave, double& _sd, double& _min,
                   double& _max, bool apply_geo)
{
    MetaData MT(MT_in); //copy constructor so original MT is not changed
    _min = MAXDOUBLE;
    _max = -MAXDOUBLE;
    _ave = _sd = 0;
    int n = 0;
    //Remove disabled images if present
    if (MT.containsLabel(MDL_ENABLED))
        MT.removeObjects(MDValueEQ(MDL_ENABLED, -1));
    // Calculate Mean
    if (MT.isEmpty())
        REPORT_ERROR(ERR_MD_OBJECTNUMBER, "There is no selected images in Metadata.");

    ImageGeneric image;
    double min, max, avg, stddev;
    FileName fnImg;
    FOR_ALL_OBJECTS_IN_METADATA(MT)
    {
        if (apply_geo)
            image.readApplyGeo(MT,__iter.objId);
        else
        {
            MT.getValue(MDL_IMAGE,fnImg,__iter.objId);
            image.read(fnImg, DATA, ALL_IMAGES, true);
        }
        image().computeStats(avg, stddev, min, max);

        if (min < _min)
            _min = min;
        if (max > _max)
            _max = max;

        _ave += avg;
        _sd += stddev;

        n++;
    }

    _ave /= n;
    _sd /= n;
}

void ImgSize(const MetaData &MD, int &Xdim, int &Ydim, int &Zdim, size_t &Ndim)
{
    if (!MD.isEmpty())
    {
        FileName fn_img;
        MD.getValue(MDL_IMAGE, fn_img, MD.firstObject());
        SingleImgSize(fn_img, Xdim, Ydim, Zdim, Ndim);

    }
    else
        REPORT_ERROR(ERR_MD_NOOBJ, "Can not read image size from empty metadata");
}

void ImgSize(const MetaData &MD, int &Xdim, int &Ydim, int &Zdim, size_t &Ndim, DataType &datatype)
{
    if (!MD.isEmpty())
    {
        FileName fn_img;
        MD.getValue(MDL_IMAGE, fn_img, MD.firstObject());
        SingleImgSize(fn_img, Xdim, Ydim, Zdim, Ndim, datatype);

    }
    else
        REPORT_ERROR(ERR_MD_NOOBJ, "Can not read image size from empty metadata");
}

void ImgSize(const FileName &filename, int &Xdim, int &Ydim, int &Zdim, size_t &Ndim)
{
    ImgSize(MetaData(filename), Xdim, Ydim, Zdim, Ndim);
}

bool ImgCompare(const FileName &filename1, const FileName &filename2)
{
    ImageGeneric Im1,Im2;
    Im1.read(filename1);
    Im2.read(filename2);
    bool result =  *(Im1.data)  == *(Im2.data);
    return( result);
}


void getBlocksInMetaDataFile(const FileName &inFile, StringVector& blockList)
{
    blockList.clear();
    if (!inFile.isMetaData())
        return;

    std::ifstream is(inFile.data(), std::ios_base::in);
    int state=0;
    String candidateBlock, line;
    while (!is.eof())
    {
        getline(is,line);
        trim(line);
        switch (state)
        {
        case 0:
            if (line.find("data_")==0)
            {
                state=1;
                candidateBlock=line.substr(5,line.size()-5);
            }
            break;
        case 1:
            if (line.find("loop_")==0)
            {
                state=0;
                blockList.push_back(candidateBlock);
            }
            break;
        }
    }
}

int MaxFileNameLength(MetaData &MD)
{
    int maxLength=0;
    FOR_ALL_OBJECTS_IN_METADATA(MD)
    {
        FileName fnImg;
        MD.getValue(MDL_IMAGE, fnImg, __iter.objId);
        int length=fnImg.length();
        maxLength=XMIPP_MAX(length,maxLength);
    }
    return maxLength;
}

void mpiSelectPart(MetaData &md, int rank, int size, int &num_img_tot)
{
    num_img_tot = md.size();
    MetaData aux(md);
    md.selectSplitPart(aux, size, rank);
}

void readMetaDataWithTwoPossibleImages(const FileName &fn, MetaData &MD)
{
    if (fn.isStar1(true))
        MD.read(fn);
    else
    {
        // Try to read a one or two column file
        std::ifstream fhIn;
        fhIn.open(fn.c_str());
        if (!fhIn)
            REPORT_ERROR(ERR_IO_NOTEXIST,fn);
        MD.clear();
        String line;
        size_t id;
        while (!fhIn.eof())
        {
            getline(fhIn,line);
            std::vector<std::string> tokens;
            tokenize(line, tokens, " \t");
            switch (tokens.size())
            {
            case 0:
                break;
            case 1:
                id = MD.addObject();
                MD.setValue(MDL_IMAGE,tokens[0], id);
                break;
            case 2:
                id = MD.addObject();
                MD.setValue(MDL_IMAGE,tokens[0], id);
                MD.setValue(MDL_ASSOCIATED_IMAGE1,tokens[1], id);
                break;
            default:
                REPORT_ERROR(ERR_MD_OBJECTNUMBER,
                             (String)"Invalid number of objects in line:"+line);
            }
        }
        fhIn.close();
    }
}

/* Substitute ------------------------------------------------------------- */
void substituteOriginalImages(const FileName &fn, const FileName &fnOrig, const FileName &fnOut,
                              MDLabel label, bool skipFirstBlock)
{
    // Read the original files
    MetaData MDorig(fnOrig);
    if (MDorig.containsLabel(MDL_ENABLED))
        MDorig.removeObjects(MDValueEQ(MDL_ENABLED, -1));
    StringVector filesOrig;
    MDorig.getColumnValues(MDL_IMAGE,filesOrig);
    MDorig.clear(); // Save memory

    // Read the blocks available
    StringVector blocks;
    getBlocksInMetaDataFile(fn, blocks);

    // Delete the output file if it exists
    if (exists(fnOut))
        unlink(fnOut.c_str());

    // Process each block
    for (int b=0; b<blocks.size(); b++)
    {
        MetaData MD;
        MD._read(fn,NULL,blocks[b]);
        if (MD.containsLabel(label) && (!skipFirstBlock || b!=0))
        {
            FileName fnImg;
            size_t stkNo;
            String stkName;
            FOR_ALL_OBJECTS_IN_METADATA(MD)
            {
                MD.getValue(label, fnImg, __iter.objId);
                fnImg.decompose(stkNo,stkName);
                MD.setValue(label, filesOrig[stkNo], __iter.objId);
            }
        }
        MD._write(fnOut,blocks[b],MD_APPEND);
    }
}
