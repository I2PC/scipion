/***************************************************************************
 * 
 * Authors:      J.R. Bilbao-Castro (jrbcast@ace.ual.es)
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

#include "metadata.h"

//-----Constructors and related functions ------------
void MetaData::clear()
{
    path.clear();
    comment.clear();
    fastStringSearch.clear();
    fastStringSearchLabel = MDL_UNDEFINED;
    activeLabels.clear();
    ignoreLabels.clear();
    isColumnFormat = true;
    inFile = FileName::FileName();
    //TODO: DROP THE TABLE FROM DB IF EXISTS
}//close clear

void MetaData::init(const std::vector<MDLabel> *labelsVector)
{
    clear();
    //TODO: ASK FOR UNIQUE TABLE NAME TO DB
    //TODO: CREATE NEW TABLE ON ON DB
    if (labelsVector != NULL)
    {
        this->activeLabels = *labelsVector;
        //TODO: CREATE COLUMNS ON DB
    }
}//close init

void MetaData::copyInfo(const MetaData &md)
{
    if (this == &md) //not sense to copy same metadata
        return;
    this->setComment(md.getComment());
    this->setPath(md.getPath());
    this->isColumnFormat = md.isColumnFormat;
    this->inFile = md.inFile;
    this->fastStringSearchLabel = md.fastStringSearchLabel;
    this->activeLabels = md.activeLabels;
    this->ignoreLabels = md.ignoreLabels;

}//close copyInfo

void MetaData::copyMetadata(const MetaData &md)
{
    if (this == &md) //not sense to copy same metadata
        return;
    init(&md.activeLabels);
    copyInfo(md);
    //TODO COPY DATA FROM MD TO THIS METADATA
}

MetaData::MetaData()
{
    init(NULL);
}//close MetaData default Constructor

MetaData::MetaData(const std::vector<MDLabel> *labelsVector)
{
    init(labelsVector);
}//close MetaData default Constructor

MetaData::MetaData(const FileName &fileName, const std::vector<MDLabel> *labelsVector)
{
    init(labelsVector);
    //TODO: Read values from File and store in MD
}//close MetaData from file Constructor

MetaData::MetaData(const MetaData &md)
{
    copyMetadata(md);
}//close MetaData copy Constructor

MetaData& MetaData::operator =(const MetaData &md)
{
    copyMetadata(md);
    return *this;
}//close metadata operator =

MetaData::~MetaData()
{
    clear();
}//close MetaData Destructor

//-------- Getters and Setters ----------

MetaDataContainer * MetaData::getObject(long int objectID) const
{
    REPORT_ERROR(-55, "getObject not implemented");
    return NULL;
}

bool MetaData::getColumnFormat() const
{
    return isColumnFormat;
}
/** Set to false for row format (parameter files)
 *  @ingroup GettersAndSetters
 *  set to true  for column format (this is the default) (docfiles)
 *
 */
void MetaData::setColumnFormat(bool column)
{
    isColumnFormat = column;
}
std::string MetaData::getPath()   const
{
    return path;
}

void MetaData::setPath(std::string newPath)
{
    path = (newPath == "") ? std::string(getcwd(NULL, 0)) : newPath;
}

std::string MetaData::getComment() const
{
    return  comment;
}

void MetaData::setComment(const std::string newComment )
{
    comment = newComment;
}

FileName MetaData::getFilename()
{
    return inFile;
}

std::vector<MDLabel> MetaData::getActiveLabels() const
{
    return activeLabels;
}

int MetaData::MaxStringLength( const MDLabel thisLabel) const
{
    REPORT_ERROR(-55, "getMaxStringsetValue, not implemented yet");
}

bool MetaData::setValue(const std::string &name, const std::string &value, long int objectID)
{
    REPORT_ERROR(-55, "setValue, not implemented yet");
}

bool MetaData::setValue(const std::string &name, const float &value, long int objectID)
{
    std::cerr << "Do not use setValue with floats, use double"<< std::endl;
    std::cerr << "Floats are banned from metadata class"<< std::endl;
    exit(1);
}
bool MetaData::setValue(const std::string &name, const char &value, long int objectID)
{
    std::cerr << "Do not use setValue with char, use string"<< std::endl;
    std::cerr << "chars are banned from metadata class"<< std::endl;
    exit(1);
}

bool MetaData::isEmpty() const
{
    REPORT_ERROR(-55, "isEmpty, not implemented yet");
    return true;

}

size_t MetaData::size() const
{
    REPORT_ERROR(-55, "isEmpty, not implemented yet");
    return 0;
}

long int MetaData::addObject(long int objectID)
{
    REPORT_ERROR(-55, "addObject not yet implemented");
    return -1;
}

void MetaData::importObjects(const MetaData &md, const std::vector<long int> &objectsToAdd)
{
    REPORT_ERROR(-55, "importObjects not yet implemented");
}

bool MetaData::removeObject(long int objectID)
{
    REPORT_ERROR(-55, "removeObject not yet implemented");
    return false;
}

void MetaData::removeObjects(const std::vector<long int> &toRemove)
{
    REPORT_ERROR(-55, "removeObjects not yet implemented");
}

//----------Iteration functions -------------------
long int MetaData::firstObject()
{
    REPORT_ERROR(-55, "firstObject not yet implemented");
}
long int MetaData::nextObject()
{
    REPORT_ERROR(-55, "nextObject not yet implemented");
}
long int MetaData::lastObject()
{
    REPORT_ERROR(-55, "lastObject not yet implemented");
}
long int MetaData::goToObject(long int objectID)
{
    REPORT_ERROR(-55, "goToObject not yet implemented");
}

//--------------IO functions -----------------------
void MetaData::write(const std::string &fileName)
{
    REPORT_ERROR(-55, "write not yet implemented");
}
void MetaData::read(FileName infile, std::vector<MDLabel> * labelsVector)
{
    REPORT_ERROR(-55, "read not yet implemented");
}

//-------------Set Operations ----------------------
void MetaData::union_(const MetaData &MD, const MDLabel thisLabel){}
void MetaData::unionAll(const MetaData &MD){}
void MetaData::aggregate(MetaData MDIn,
               MDLabel aggregateLabel,
               MDLabel entryLabel,
               MDLabel operationLabel){}

void MetaData::merge(const FileName &fn){}
void MetaData::intersection(MetaData & minuend, MetaData & , MDLabel thisLabel){}
void MetaData::substraction(MetaData & minuend, MetaData & subtrahend,
                  MDLabel thisLabel){}


/*----------   Statistics --------------------------------------- */
#include "metadata_extension.h"
void get_statistics(MetaData MT_in, Image<double> & _ave, Image<double> & _sd, double& _min,
                    double& _max, bool apply_geo)
{
    MetaData MT(MT_in); //copy constructor so original MT is not changed
    _min = MAXFLOAT;
    _max = 0;
    bool first = true;
    int n = 0;
    // Calculate Mean
    long int ret = MT.firstObject();
    if (ret == MetaData::NO_OBJECTS_STORED)
    {
        std::cerr << "Empty inputFile File\n";
        exit(1);
    }
    FileName image_name;
    int _enabled;
    do
    {
        MT.getValue(MDL_IMAGE, image_name);
        MT.getValue(MDL_ENABLED, _enabled);
        if (_enabled == (-1) || image_name == "")
            continue;
        Image<double> image;
        MetaDataContainer * aux = MT.getObject();
        image.read(image_name,true,
                   -1,apply_geo,
                   false);
        double min, max, avg, stddev;
        image().computeStats(avg, stddev, min, max);
        if (_min > min)
            _min = min;
        if (_max < max)
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
    while (MT.nextObject() != MetaData::NO_MORE_OBJECTS);

    if (n > 0)
        _ave() /= n;
    _sd = _ave;
    _sd().initZeros();
    // Calculate SD
    MT.firstObject();
    do
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
    while (MT.nextObject() != MetaData::NO_MORE_OBJECTS);
    _sd() /= (n - 1);
    _sd().selfSQRT();
}

void ImgSize(MetaData &MD, int &Xdim, int &Ydim, int &Zdim, int &Ndim)
{
    if (!MD.isEmpty())
    {
        FileName fn_img;
        Image<double> img;
        MD.getValue(MDL_IMAGE, fn_img);
        img.read(fn_img, false);
        img.getDimensions(Xdim, Ydim, Zdim, Ndim);
    }
    else
        REPORT_ERROR(-1, "Can not read image size from empty metadata");
}
