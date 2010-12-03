/*
 * metadata_extension.h
 *
 *  Created on: May 12, 2010
 *      Author: roberto
 */
#include "metadata_extension.h"

/*----------   Statistics --------------------------------------- */
void getStatistics(MetaData &MT_in, Image<double> & _ave, Image<double> & _sd, double& _min,
                   double& _max, bool apply_geo)
{
    MetaData MT(MT_in); //copy constructor so original MT is not changed
    _min = MAXFLOAT;
    _max = 0.;
    bool first = true;
    int n = 0;
    // Calculate Mean
    if (MT.isEmpty())
    {
        std::cerr << "Empty inputFile File\n";
        exit(1);
    }

    FileName image_name;
    int _enabled;
    FOR_ALL_OBJECTS_IN_METADATA(MT)
    {
        MT.getValue(MDL_IMAGE, image_name);
        MT.getValue(MDL_ENABLED, _enabled);
        if (_enabled == (-1) || image_name == "")
            continue;
        Image<double> image;
        image.read(image_name,true,-1,apply_geo, false);
        double min, max, avg, stddev;
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
        MT.getValue(MDL_IMAGE, image_name);
        MT.getValue(MDL_ENABLED, _enabled);
        if (_enabled == (-1) || image_name == "")
            continue;

        Image<double> image, tmpImg;
        image.read(image_name);
        tmpImg() = ((image() - _ave()));
        tmpImg() *= tmpImg();
        _sd() += tmpImg();
    }
    _sd() /= (n - 1);
    _sd().selfSQRT();
}

void SingleImgSize(const FileName &filename, int &Xdim, int &Ydim, int &Zdim, unsigned long &Ndim)
{
    Image<char> img;
    img.read(filename, false);
    img.getDimensions(Xdim, Ydim, Zdim, Ndim);
}

void ImgSize(const MetaData &MD, int &Xdim, int &Ydim, int &Zdim, unsigned long &Ndim)
{
    if (!MD.isEmpty())
    {
        FileName fn_img;
        MD.getValue(MDL_IMAGE, fn_img);
        SingleImgSize(fn_img, Xdim, Ydim, Zdim, Ndim);

    }
    else
        REPORT_ERROR(ERR_MD_NOOBJ, "Can not read image size from empty metadata");
}

void ImgSize(const FileName &filename, int &Xdim, int &Ydim, int &Zdim, unsigned long  &Ndim)
{
    ImgSize(MetaData(filename), Xdim, Ydim, Zdim, Ndim);
}

int MaxFileNameLength(MetaData &MD)
{
    int maxLength=0;
    FOR_ALL_OBJECTS_IN_METADATA(MD)
    {
        FileName fnImg;
        MD.getValue(MDL_IMAGE,fnImg);
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
    	std::string line;
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
    			MD.addObject();
    			MD.setValue(MDL_IMAGE,tokens[0]);
    			break;
    		case 2:
    			MD.addObject();
    			MD.setValue(MDL_IMAGE,tokens[0]);
    			MD.setValue(MDL_ASSOCIATED_IMAGE1,tokens[1]);
    			break;
    		default:
    			REPORT_ERROR(ERR_MD_OBJECTNUMBER,
    				(std::string)"Invalid number of objects in line:"+line);
    		}
    	}
    	fhIn.close();
    }
}
