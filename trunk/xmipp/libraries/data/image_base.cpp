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

#include "image_base.h"


void ImageBase::init()
{
    clearHeader();
    dataflag = -1;
    if (isComplexT())
        transform = Standard;
    else
        transform = NoTransform;
    i = 0;
    filename = "";
    offset = 0;
    swap = 0;
    replaceNsize=0;
    mmapOnRead = mmapOnWrite = false;
    mappedSize = 0;
    mFd    = NULL;
}

void ImageBase::clearHeader()
{
    MDMainHeader.clear();
    MD.clear();
    //Just to ensure there is an empty MDRow
    MD.push_back(MDMainHeader);
}

/** General read function
 */
int ImageBase::read(const FileName &name, bool readdata, int select_img, bool mapData)
{
    ImageFHandler* hFile = openFile(name);
    int err = _read(name, hFile, readdata, select_img, false, false, NULL, mapData);
    closeFile(hFile);

    return err;
}

/** New mapped file */
void ImageBase::newMappedFile(int Xdim, int Ydim, int Zdim, int Ndim, const FileName &_filename,
                              bool createTempFile)
{
    clear();
    mmapOnWrite = true;
    setDimensions(Xdim, Ydim, Zdim, Ndim);
    MD.resize(Ndim);
    filename = _filename;
    FileName fnToOpen;
    if (createTempFile)
    {
        tempFilename.initUniqueName("temp_XXXXXX");
        fnToOpen = tempFilename + ":" + _filename.getExtension();
    }
    else
    	fnToOpen=_filename;

    ImageFHandler *hFile = openFile(fnToOpen, WRITE_OVERWRITE);
    _write(fnToOpen, hFile, -1, false, WRITE_OVERWRITE);
    closeFile(hFile);
}

/** General read function
 */
int ImageBase::readApplyGeo(const FileName &name, bool readdata, int select_img,
                            bool only_apply_shifts, MDRow * row)
{
    ImageFHandler* hFile = openFile(name);
    int err = _read(name, hFile, readdata, select_img, only_apply_shifts, row, false);
    closeFile(hFile);

    return err;
}

/** Read an image from metadata, filename is provided
*/
int ImageBase::readApplyGeo(const FileName &name, const MetaData &md, size_t objId, bool readdata,
                            int select_img, bool only_apply_shifts)
{
    ImageFHandler* hFile = openFile(name);
    MDRow row;
    md.getRow(row, objId);
    int err = _read(name, hFile, readdata, select_img, true, only_apply_shifts, &row, false);
    closeFile(hFile);

    return err;
}

/** Read an image from metadata, filename is taken from MDL_IMAGE
 */
int ImageBase::readApplyGeo(const MetaData &md, size_t objId, bool readdata, int select_img,
                            bool only_apply_shifts)
{
    MDRow row;
    md.getRow(row, objId);
    FileName name;
    row.getValue(MDL_IMAGE, name);
    ImageFHandler* hFile = openFile(name);
    int err = _read(name, hFile, readdata, select_img, true, only_apply_shifts, &row, false);
    closeFile(hFile);

    return err;
}

void ImageBase::write(const FileName &name, int select_img, bool isStack,
                      int mode, bool adjust)
{
    // If image is already mapped to file then close the file and clear.
    if (mmapOnWrite && mappedSize > 0)
    {
        munmapFile();
        if (tempFilename!="")
        {
            if (std::rename(tempFilename.c_str(),name.c_str())!=0)
                REPORT_ERROR(ERR_IO, "Error renaming the file.");
        }
        return;
    }

    const FileName &fname = (name == "") ? filename : name;

    /* If the filename is in stack we will suppose you want to write this,
     * even if you have not set the flags to.
     */
    if ( fname.isInStack() && !isStack && mode == WRITE_OVERWRITE)
    {
        isStack = true;
        mode = WRITE_APPEND;
    }
    else if (!isStack && mode != WRITE_OVERWRITE)
        mode = WRITE_OVERWRITE;

    ImageFHandler* hFile = openFile(fname, mode);
    _write(fname, hFile, select_img, isStack, mode, adjust);
    closeFile(hFile);
}

void ImageBase::swapPage(char * page, size_t pageNrElements, DataType datatype)
{
    unsigned long datatypesize = gettypesize(datatype);
#ifdef DEBUG

    std::cerr<<"DEBUG swapPage: Swapping image data with swap= "
    << swap<<" datatypesize= "<<datatypesize
    << " pageNrElements " << pageNrElements
    << " datatype " << datatype
    <<std::endl;
    ;
#endif

    // Swap bytes if required
    if ( swap == 1 )
    {
        if ( datatype >= ComplexShort )
            datatypesize /= 2;
        for ( unsigned long i=0; i<pageNrElements; i+=datatypesize )
            swapbytes(page+i, datatypesize);
    }
    else if ( swap > 1 )
    {
        for ( unsigned long i=0; i<pageNrElements; i+=swap )
            swapbytes(page+i, swap);
    }
}

/** Get Rot angle
    *
    * @code
    * std::cout << "First Euler angle " << I.rot() << std::endl;
    * @endcode
    */
double ImageBase::rot(const long int n) const
{
    double dummy = 0;
    MD[n].getValue(MDL_ANGLEROT, dummy);
    return dummy;
}

/** Get Tilt angle
 *
 * @code
 * std::cout << "Second Euler angle " << I.tilt() << std::endl;
 * @endcode
 */
double ImageBase::tilt(const long int n) const
{
    double dummy = 0;
    MD[n].getValue(MDL_ANGLETILT, dummy);
    return dummy;
}

/** Get Psi angle
 *
 * @code
 * std::cout << "Third Euler angle " << I.psi() << std::endl;
 * @endcode
 */
double ImageBase::psi(const long int n) const
{
    double dummy = 0;
    MD[n].getValue(MDL_ANGLEPSI, dummy);
    return dummy;
}

/** Get Xoff
 *
 * @code
 * std::cout << "Origin offset in X " << I.Xoff() << std::endl;
 * @endcode
 */
double ImageBase::Xoff(const long int n) const
{
    double dummy = 0;
    MD[n].getValue(MDL_ORIGINX, dummy);
    return dummy;
}

/** Get Yoff
 *
 * @code
 * std::cout << "Origin offset in Y " << I.Yoff() << std::endl;
 * @endcode
 */
double ImageBase::Yoff(const long int n) const
{
    double dummy = 0;
    MD[n].getValue(MDL_ORIGINY, dummy);
    return dummy;
}

/** Get Zoff
 *
 * @code
 * std::cout << "Origin offset in Z " << I.Zoff() << std::endl;
 * @endcode
 */
double ImageBase::Zoff(const long int n) const
{
    double dummy = 0;
    MD[n].getValue(MDL_ORIGINZ, dummy);
    return dummy;
}

/** Get Weight
*
* @code
* std::cout << "weight= " << I.weight() << std::endl;
* @endcode
*/
double ImageBase::weight(const long int n) const
{
    double dummy = 1;
    MD[n].getValue(MDL_WEIGHT, dummy);
    return dummy;
}

/** Get Scale factor
*
* @code
* std::cout << "scale= " << I.scale() << std::endl;
* @endcode
*/
double ImageBase::scale(const long int n) const
{
    double dummy = 1;
    MD[n].getValue(MDL_SCALE, dummy);
    return dummy;
}


/** Get Flip
*
* @code
* std::cout << "flip= " << flip() << std::endl;
* @endcode
*/
bool ImageBase::flip(const long int n) const
{
    bool dummy = false;
    MD[n].getValue(MDL_FLIP, dummy);
    return dummy;
}

/** Data type
    *
    * @code
    * std::cout << "datatype= " << dataType() << std::endl;
    * @endcode
    */
DataType ImageBase::dataType() const
{
    int dummy;
    MDMainHeader.getValue(MDL_DATATYPE, dummy);
    return (DataType)dummy;
}

/** Sampling RateX
*
* @code
* std::cout << "sampling= " << samplingRateX() << std::endl;
* @endcode
*/
double ImageBase::samplingRateX() const
{
    double dummy = 1.;
    MDMainHeader.getValue(MDL_SAMPLINGRATEX, dummy);
    return dummy;
}

/** Set Euler angles in image header
 */
void ImageBase::setEulerAngles(double rot, double tilt, double psi,
                               long int n)
{
    MD[n].setValue(MDL_ANGLEROT, rot);
    MD[n].setValue(MDL_ANGLETILT, tilt);
    MD[n].setValue(MDL_ANGLEPSI, psi);
}

/** Get Euler angles from image header
 */
void ImageBase::getEulerAngles(double &rot, double &tilt, double &psi,
                               long int n)
{
    MD[n].getValue(MDL_ANGLEROT, rot);
    MD[n].getValue(MDL_ANGLETILT, tilt);
    MD[n].getValue(MDL_ANGLEPSI, psi);
}

/** Set origin offsets in image header
     */
void ImageBase::setShifts(double xoff, double yoff, double zoff, long int n)
{
    MD[n].setValue(MDL_ORIGINX, xoff);
    MD[n].setValue(MDL_ORIGINY, yoff);
    MD[n].setValue(MDL_ORIGINZ, zoff);
}
/** Get origin offsets from image header
  */
void ImageBase::getShifts(double &xoff, double &yoff, double &zoff, long int n) const
{
    MD[n].getValue(MDL_ORIGINX, xoff);
    MD[n].getValue(MDL_ORIGINY, yoff);
    MD[n].getValue(MDL_ORIGINZ, zoff);
}

/** Open file function
  * Open the image file and returns its file hander.
  */
ImageFHandler* ImageBase::openFile(const FileName &name, int mode) const
{
    ImageFHandler* hFile = new ImageFHandler;
    FileName fileName, headName = "";
    FileName ext_name = name.getFileFormat();

    int dump;
    name.decompose(dump, fileName);

    fileName = fileName.removeFileFormat();

    size_t found = fileName.find_first_of("%");
    if (found!=std::string::npos)
        fileName = fileName.substr(0, found) ;

    hFile->exist = exists(fileName);

    std::string wmChar;

    switch (mode)
    {
    case WRITE_READONLY:
        if (!hFile->exist)
            REPORT_ERROR(ERR_IO_NOTEXIST,(std::string) "Cannot read file "
                         + fileName + ". It does not exist" );
        wmChar = "r";
        break;
    case WRITE_OVERWRITE:
        wmChar = "w";
        break;
    case WRITE_APPEND:
    case WRITE_REPLACE:
        if (hFile->exist)
            wmChar = "r+";
        else
            wmChar = "w+";
        break;

    }

    if (ext_name.contains("tif"))
    {
        TIFFSetWarningHandler(NULL); // Switch off warning messages
        if ((hFile->tif = TIFFOpen(fileName.c_str(), wmChar.c_str())) == NULL)
            REPORT_ERROR(ERR_IO_NOTOPEN,"rwTIFF: There is a problem opening the TIFF file.");
        hFile->fimg = NULL;
        hFile->fhed = NULL;
    }
    else
    {
        hFile->tif = NULL;

        if (ext_name.contains("img") || ext_name.contains("hed"))
        {
            fileName = fileName.withoutExtension();
            headName = fileName.addExtension("hed");
            fileName = fileName.addExtension("img");
        }
        else if (ext_name.contains("raw"))
        {
            if (mode != WRITE_READONLY || exists(fileName.addExtension("inf")) )
            {
                headName = fileName.addExtension("inf");
                ext_name = "inf";
            }
            else
                ext_name = "raw";
        }
        else if (ext_name.contains("inf"))
        {
            headName = fileName;
            fileName = fileName.withoutExtension();
            ext_name = "inf";
        }

        // Open image file
        if ( ( hFile->fimg = fopen(fileName.c_str(), wmChar.c_str()) ) == NULL )
            REPORT_ERROR(ERR_IO_NOTOPEN,(std::string)"Image::openFile cannot open: " + fileName);

        if (headName != "")
        {
            if ( ( hFile->fhed = fopen(headName.c_str(), wmChar.c_str()) ) == NULL )
                REPORT_ERROR(ERR_IO_NOTOPEN,(std::string)"Image::openFile cannot open: " + headName);
        }
        else
            hFile->fhed = NULL;

    }
    hFile->fileName = fileName;
    hFile->headName = headName;
    hFile->ext_name = ext_name;

    return hFile;
}

/** Close file function.
  * Close the image file according to its name and file handler.
  */
void ImageBase::closeFile(ImageFHandler* hFile)
{
    FileName ext_name;
    FILE* fimg, *fhed;
    TIFF* tif;

    if (hFile != NULL)
    {
        ext_name = hFile->ext_name;
        fimg = hFile->fimg;
        fhed = hFile->fhed;
        tif  = hFile->tif;
    }
    else
    {
        ext_name = filename.getFileFormat();
        fimg = this->fimg;
        fhed = this->fhed;
        tif  = this->tif;
    }

    if (ext_name.contains("tif"))
        TIFFClose(tif);
    else
    {
        if (fclose(fimg) != 0 )
            REPORT_ERROR(ERR_IO_NOCLOSED,(std::string)"Can not close image file "+ filename);

        if (fhed != NULL &&  fclose(fhed) != 0 )
            REPORT_ERROR(ERR_IO_NOCLOSED,(std::string)"Can not close header file of "
                         + filename);
    }
    delete hFile;
}
