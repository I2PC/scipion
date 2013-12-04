/*
 * metadata_extension.h
 *
 *  Created on: May 12, 2010
 *      Author: roberto
 */
#include "metadata_extension.h"
#include "xmipp_image_extension.h"
#include "xmipp_fftw.h"
#include "xmipp_image_convert.h"

#ifndef __linux__
#define MAXDOUBLE __DBL_MAX__
#endif


/*----------   Statistics --------------------------------------- */
//Copy of the Metadata is required to remove disabled objects before computing stats
void getStatistics(MetaData MD, Image<double> & _ave, Image<double> & _sd, double& _min,
                   double& _max, bool apply_geo, MDLabel image_label)
{
    _min = MAXDOUBLE;
    _max = -MAXDOUBLE;
    bool first = true;
    int n = 0;
    //Remove disabled images if present
    MD.removeDisabled();
    // Calculate Mean
    if (MD.isEmpty())
        REPORT_ERROR(ERR_MD_OBJECTNUMBER, "There is no selected images in Metadata.");

    Image<double> image, tmpImg;
    double min=0, max=0, avg, stddev;
    FileName fnImg;
    FOR_ALL_OBJECTS_IN_METADATA(MD)
    {
        MD.getValue(image_label,fnImg,__iter.objId);
        if (apply_geo)
            image.readApplyGeo(fnImg, MD,__iter.objId);
        else
            image.read(fnImg);
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
    FOR_ALL_OBJECTS_IN_METADATA(MD)
    {
        MD.getValue(image_label,fnImg,__iter.objId);
        if (apply_geo)
            image.readApplyGeo(fnImg, MD,__iter.objId);
        else
            image.read(fnImg);
        tmpImg() = ((image() - _ave()));
        tmpImg() *= tmpImg();
        _sd() += tmpImg();
    }
    _sd() /= (n - 1);
    _sd().selfSQRT();
}

//Copy of the Metadata is required to remove disabled objects before computing stats
void getAverageApplyGeo(MetaData MD, MultidimArray<double> & _ave, MDLabel image_label)
{
    bool first = true;
    int n = 0;
    MD.removeDisabled();
    // Calculate Mean
    if (MD.isEmpty())
        return;

    Image<double> image;
    FileName fnImg;

    FOR_ALL_OBJECTS_IN_METADATA(MD)
    {
        MD.getValue(image_label,fnImg,__iter.objId);
        image.readApplyGeo(fnImg, MD,__iter.objId);

        if (first)
        {
            _ave = image();
            first = false;
        }
        _ave += image();
        n++;
    }

    if (n > 0)
        _ave /= n;
}

/*----------   Statistics --------------------------------------- */
//Copy of the Metadata is required to remove disabled objects before computing stats
void getStatistics(MetaData MD, double& _ave, double& _sd, double& _min,
                   double& _max, bool apply_geo, MDLabel image_label)
{
    _min = MAXDOUBLE;
    _max = -MAXDOUBLE;
    _ave = _sd = 0;
    int n = 0;
    //Remove disabled images if present
    MD.removeDisabled();
    // Calculate Mean
    if (MD.isEmpty())
        REPORT_ERROR(ERR_MD_OBJECTNUMBER, "There is no selected images in Metadata.");

    ImageGeneric image;
    double min, max, avg, stddev;
    FileName fnImg;
    FOR_ALL_OBJECTS_IN_METADATA(MD)
    {
        MD.getValue(image_label,fnImg,__iter.objId);
        if (apply_geo)
            image.readApplyGeo(fnImg, MD,__iter.objId);
        else
            image.read(fnImg, DATA, ALL_IMAGES, true);
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

/* Get Fourier statistics ------------------------------------------------- */
void getFourierStatistics(MetaData &MDin, double sam, MetaData &MDout,
                          bool do_dpr, double max_sam, MDLabel image_label)
{
    MetaData MDaux;
    std::vector<MetaData> vMD;
    MDaux.randomize(MDin);
    MDaux.split(2,vMD,image_label);
    MetaData &MD1 = vMD.at(0);
    MetaData &MD2 = vMD.at(1);

    Image<double> I1, I2, Id;
    double dummy;
    getStatistics(MD1,I1,Id,dummy,dummy,true, image_label);
    getStatistics(MD2,I2,Id,dummy,dummy,true, image_label);
    I1().setXmippOrigin();
    I2().setXmippOrigin();

    MultidimArray<double> freq, frc, dpr, frc_noise, ssnr, error_l2;
    frc_dpr(I1(), I2(), sam, freq, frc, frc_noise, dpr,error_l2,do_dpr);

    MDout.clear();
    FOR_ALL_ELEMENTS_IN_ARRAY1D(freq)
    {
        if (i>0)
        {
            size_t id=MDout.addObject();
            if(max_sam >=0 && ((1./dAi(freq, i))<max_sam) )
            {
                if(do_dpr)
                    dAi(dpr, i)=0.;
                dAi(frc, i)=0.;
            }
            MDout.setValue(MDL_RESOLUTION_FREQ,dAi(freq, i),id);
            MDout.setValue(MDL_RESOLUTION_FRC,dAi(frc, i),id);
            if(do_dpr)
                MDout.setValue(MDL_RESOLUTION_DPR,dAi(dpr, i),id);
            MDout.setValue(MDL_RESOLUTION_ERRORL2,dAi(error_l2, i),id);
            MDout.setValue(MDL_RESOLUTION_FRCRANDOMNOISE,dAi(frc_noise, i),id);
            MDout.setValue(MDL_RESOLUTION_FREQREAL,1./dAi(freq, i),id);
        }
    }
}

void getImageSize(const MetaData &MD, size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim, MDLabel image_label)
{
    if (!MD.isEmpty())
    {
        FileName fn_img;
        MD.getValue(image_label, fn_img, MD.firstObject());
        getImageSize(fn_img, Xdim, Ydim, Zdim, Ndim);

    }
    else
        REPORT_ERROR(ERR_MD_NOOBJ, "Can not read image size from empty metadata");
}

void getImageInfo(const MetaData &md, size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim, DataType &datatype, MDLabel image_label)
{
    if (!md.isEmpty())
    {
        FileName fn_img;
        md.getValue(image_label, fn_img, md.firstObject());
        getImageInfo(fn_img, Xdim, Ydim, Zdim, Ndim, datatype);

    }
    else
        REPORT_ERROR(ERR_MD_NOOBJ, "Can not read image size from empty metadata");
}

void getImageInfo(const MetaData &md, ImageInfo &imgInfo, MDLabel image_label)
{
    if (!md.isEmpty())
    {
        FileName fn_img;
        md.getValue(image_label, fn_img, md.firstObject());
        getImageInfo(fn_img, imgInfo);
    }
    else
        REPORT_ERROR(ERR_MD_NOOBJ, "Can not read image size from empty metadata");
}

void getImageSizeFromFilename(const FileName &filename, size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim, MDLabel image_label)
{
    if (filename.hasImageExtension())
        getImageSize(filename, Xdim, Ydim, Zdim, Ndim);
    else
    {
        MetaData mdi;
        mdi.setMaxRows(1);
        mdi.read(filename);
        getImageSize(mdi, Xdim, Ydim, Zdim, Ndim, image_label);
    }
}

bool compareImage(const FileName &filename1, const FileName &filename2)
{
    ImageGeneric Im1,Im2;
    Im1.read(filename1);
    Im2.read(filename2);
    bool result =  *(Im1.data)  == *(Im2.data);
    return( result);
}

bool compareImageSize(const FileName &filename1, const FileName &filename2)
{
    size_t x,y,z, X, Y, Z, n, N;
    getImageSizeFromFilename(filename1,x,y,z,n);
    getImageSizeFromFilename(filename2,X,Y,Z,N);
    return (x==X && y == Y && z == Z && n == N);
}

bool compareTwoMetadataFiles(const FileName &fn1, const FileName &fn2)
{
    StringVector blockList;
    getBlocksInMetaDataFile(fn1,blockList);
    FileName fn_1aux, fn_2aux;
    MetaData md1, md2;

    for (StringVector::iterator it= blockList.begin();
         it!=blockList.end(); it++)
    {
        fn_1aux.compose(*it, fn1);
        fn_2aux.compose(*it, fn2);
        md1.read(fn_1aux);
        md2.read(fn_2aux);
        if (!(md1 == md2))
            return false;
    }

    return true;
}

void copyImages(const MetaData &md, const char * output, bool independent, MDLabel image_label)
{
    ProgConvImg conv;
    FileName out, oroot;
    if (independent)
        oroot = output;
    else
        out = output;
    //Image convert dont use applyGeo
    conv.setup(new MetaData(md), out, oroot, false, image_label);
    conv.tryRun();
}

int maxFileNameLength(const MetaData &md, MDLabel image_label)
{
    int maxLength=0;
    FOR_ALL_OBJECTS_IN_METADATA(md)
    {
        FileName fnImg;
        md.getValue(image_label, fnImg, __iter.objId);
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

void readMetaDataWithTwoPossibleImages(const FileName &fn, MetaData &md)
{
    if (fn.isStar1(true))
        md.read(fn);
    else
    {
        // Try to read a one or two column file
        std::ifstream fhIn;
        fhIn.open(fn.c_str());
        if (!fhIn)
            REPORT_ERROR(ERR_IO_NOTEXIST,fn);
        md.clear();
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
                id = md.addObject();
                md.setValue(MDL_IMAGE,tokens[0], id);
                break;
            case 2:
                id = md.addObject();
                md.setValue(MDL_IMAGE,tokens[0], id);
                md.setValue(MDL_IMAGE1,tokens[1], id);
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
    MetaData mdorig(fnOrig);
    if (mdorig.containsLabel(MDL_ENABLED))
        mdorig.removeObjects(MDValueEQ(MDL_ENABLED, -1));
    StringVector filesOrig;
    mdorig.getColumnValues(MDL_IMAGE, filesOrig);
    mdorig.clear(); // Save memory
    FileName auxFn;

    // Read the blocks available
    StringVector blocks;
    getBlocksInMetaDataFile(fn, blocks);

    // Delete the output file if it exists
    fnOut.deleteFile();

    // Process each block
    for (size_t b=0; b<blocks.size(); b++)
    {
        MetaData md;
        md.read(blocks[b]+"@"+fn);
        if (md.containsLabel(label) && (!skipFirstBlock || b!=0))
        {
            FileName fnImg;
            size_t stkNo;
            String stkName;
            FOR_ALL_OBJECTS_IN_METADATA(md)
            {
                md.getValue(label, fnImg, __iter.objId);
                fnImg.decompose(stkNo,stkName);
                md.setValue(label, filesOrig[stkNo], __iter.objId);
            }
        }
        auxFn.compose(blocks[b],fnOut);
        md.write(auxFn, MD_APPEND);
        //MD._write(fnOut,blocks[b],MD_APPEND);
    }
}
