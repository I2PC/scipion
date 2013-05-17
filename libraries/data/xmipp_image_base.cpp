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
#include "xmipp_image.h"
#include "xmipp_error.h"

//This is needed for static memory allocation

void ImageBase::init()
{
    clearHeader();

    filename = tempFilename = dataFName = "";
    fimg = fhed = NULL;
    hFile = NULL;
    tif = NULL;
    dataMode = DATA;
    stayOpen = false;
    transform = isComplexT() ? Standard : NoTransform;
    filename.clear();
    offset = 0;
    swap = swapWrite = 0;
    replaceNsize = 0;
    _exists = mmapOnRead = mmapOnWrite = false;
    mFd        = 0;
    mappedSize = mappedOffset = virtualOffset = 0;
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
int ImageBase::read(const FileName &name, DataMode datamode, size_t select_img,
                    bool mapData, int mode)
{
    if (!mapData)
        mode = WRITE_READONLY; //TODO: Check if openfile other than readonly is necessary

    hFile = openFile(name, mode);
    int err = _read(name, hFile, datamode, select_img, mapData);
    closeFile(hFile);

    return err;
}

int ImageBase::readMapped(const FileName &name, size_t select_img, int mode)
{
    read(name, HEADER);
    bool swap = this->swap > 0;
    return read(name, DATA, select_img, !swap, mode);
}

int ImageBase::readOrReadMapped(const FileName &name, size_t select_img, int mode)
{
    try
    {
        return read(name, DATA, select_img, false, mode);
    }
    catch (XmippError &xe)
    {
        if (xe.__errno == ERR_MEM_NOTENOUGH)
        {
            reportWarning("ImageBase::readOrReadMapped: Not enough memory to allocate. \n"
                          " Proceeding to map image from file.");
            return readMapped(name, select_img, mode);
        }
        else
            throw xe;
    }
}

int ImageBase::readOrReadPreview(const FileName &name, size_t Xdim, size_t Ydim, int select_slice, size_t select_img,
                                 bool mapData)
{
    read(name, HEADER);
    size_t imXdim, imYdim, imZdim, imNdim;
    getDimensions(imXdim, imYdim, imZdim, imNdim);

    if (imXdim != Xdim || imYdim != Ydim)
        return readPreview(name, Xdim, Ydim, select_slice, select_img);
    else
    {
        int ret = read(name, DATA, select_img, mapData);
        if (select_slice != ALL_SLICES)
            movePointerTo(select_slice);
        return ret;
    }

}

/** New mapped file */
void ImageBase::mapFile2Write(size_t Xdim, size_t Ydim, size_t Zdim, const FileName &_filename,
                              bool createTempFile, size_t select_img, bool isStack, int mode)
{
    /** If XMIPP_MMAP is not defined this function is supposed to create
     *  the empty file only */
#ifdef XMIPP_MMAP
    mmapOnWrite = true;
#endif

    setDimensions(Xdim, Ydim, Zdim, 1); // Images with Ndim >1 cannot be mapped to image file
    MD.resize(1);
    filename = _filename;
    FileName fnToOpen;
    if (createTempFile)
    {
        tempFilename.initUniqueName("temp_XXXXXX");
        fnToOpen = tempFilename + ":" + _filename.getExtension();
    }
    else
        fnToOpen=_filename;

    /* If the filename is in stack or an image is selected, we will suppose
     * you want to write this, even if you have not set the flags to.
     */
    if ( (filename.isInStack() || select_img > ALL_IMAGES) && mode == WRITE_OVERWRITE)
    {
        isStack = true;
        mode = WRITE_REPLACE;
    }

    hFile = openFile(fnToOpen, mode);
    _write(filename, hFile, select_img, isStack, mode);
    closeFile(hFile);
}

/** General read function
 */
/** Macros for dont type */
#define GET_ROW()               MDRow row; md.getRow(row, objId)

#define READ_AND_RETURN()        ImageFHandler* hFile = openFile(name);\
                                  int err = _read(name, hFile, params.datamode, params.select_img); \
                                  applyGeo(row, params.only_apply_shifts, params.wrap); \
                                  closeFile(hFile); \
                                  return err

#define APPLY_GEO()        MDRow row; md.getRow(row, objId); \
                           applyGeo(row, params.only_apply_shifts, params.wrap) \

void ImageBase::applyGeo(const MetaData &md, size_t objId,
const ApplyGeoParams &params)
{
    APPLY_GEO();
}

int ImageBase::readApplyGeo(const FileName &name, const MDRow &row,
                            const ApplyGeoParams &params)
{
    READ_AND_RETURN();
}

/** Read an image from metadata, filename is provided
*/
int ImageBase::readApplyGeo(const FileName &name, const MetaData &md, size_t objId,
                            const ApplyGeoParams &params)
{
    GET_ROW();
    READ_AND_RETURN();
}

/** Read an image from metadata, filename is taken from MDL_IMAGE
 */
int ImageBase::readApplyGeo(const MetaData &md, size_t objId,
                            const ApplyGeoParams &params)
{
    GET_ROW();
    FileName name;
    row.getValue(MDL_IMAGE, name);
    READ_AND_RETURN();
}


void ImageBase::write(const FileName &name, size_t select_img, bool isStack,
                      int mode, CastWriteMode castMode, int _swapWrite)
{
    const FileName &fname = (name.empty()) ? filename : name;

    if (mmapOnWrite && mappedSize > 0)
    {
        bool hasTempFile = !tempFilename.empty();
        if (hasTempFile && fname.isInStack())
            mmapOnWrite = !(mmapOnRead = true); // We change mmap mode from write to read to allow writing the image into a stack.
        else
        {
            if (_swapWrite > 0)
                REPORT_ERROR(ERR_ARG_INCORRECT, "Cannot swap endianess on writing if file is already mapped.");
            munmapFile();
            if (hasTempFile && std::rename(tempFilename.c_str(), fname.c_str()) != 0)
                REPORT_ERROR(ERR_IO, formatString("Error renaming file '%s' to '%s'.", tempFilename.c_str(), fname.c_str()));
            return;
        }
    }
    // Swap the endianess of the image file when writing
    swapWrite = _swapWrite;

    /* If the filename is in stack we will suppose you want to write this,
     * even if you have not set the flags to.
     */
    if ( fname.isInStack() && mode == WRITE_OVERWRITE)
    {
        isStack = true;
        mode = WRITE_REPLACE;
    }
    //    else if (!isStack && mode != WRITE_OVERWRITE)
    //        mode = WRITE_OVERWRITE;

    hFile = openFile(fname, mode);
    _write(fname, hFile, select_img, isStack, mode, castMode);
    closeFile(hFile);
}

void ImageBase::swapPage(char * page, size_t pageNrElements, DataType datatype, int swap)
{
    size_t datatypesize = gettypesize(datatype);
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
        if ( datatype >= DT_CShort )
            datatypesize /= 2;
        for ( size_t i = 0; i < pageNrElements; i += datatypesize )
            swapbytes(page+i, datatypesize);
    }
    else if ( swap > 1 )
    {
        for (size_t i=0; i<pageNrElements; i+=swap )
            swapbytes(page+i, swap);
    }
}

/** Get Rot angle
    *
    * @code
    * std::cout << "First Euler angle " << I.rot() << std::endl;
    * @endcode
    */
double ImageBase::rot(const size_t n) const
{
    double dummy = 0;
    MD[n].getValue(MDL_ANGLE_ROT, dummy);
    return dummy;
}

/** Get Tilt angle
 *
 * @code
 * std::cout << "Second Euler angle " << I.tilt() << std::endl;
 * @endcode
 */
double ImageBase::tilt(const size_t n) const
{
    double dummy = 0;
    MD[n].getValue(MDL_ANGLE_TILT, dummy);
    return dummy;
}

/** Get Psi angle
 *
 * @code
 * std::cout << "Third Euler angle " << I.psi() << std::endl;
 * @endcode
 */
double ImageBase::psi(const size_t n) const
{
    double dummy = 0;
    MD[n].getValue(MDL_ANGLE_PSI, dummy);
    return dummy;
}

/** Get Xoff
 *
 * @code
 * std::cout << "Origin offset in X " << I.Xoff() << std::endl;
 * @endcode
 */
double ImageBase::Xoff(const size_t n) const
{
    double dummy = 0;
    MD[n].getValue(MDL_SHIFT_X, dummy);
    return dummy;
}

/** Get Yoff
 *
 * @code
 * std::cout << "Origin offset in Y " << I.Yoff() << std::endl;
 * @endcode
 */
double ImageBase::Yoff(const size_t n) const
{
    double dummy = 0;
    MD[n].getValue(MDL_SHIFT_Y, dummy);
    return dummy;
}

/** Get Zoff
 *
 * @code
 * std::cout << "Origin offset in Z " << I.Zoff() << std::endl;
 * @endcode
 */
double ImageBase::Zoff(const size_t n) const
{
    double dummy = 0;
    MD[n].getValue(MDL_SHIFT_Z, dummy);
    return dummy;
}

/** Get Weight
*
* @code
* std::cout << "weight= " << I.weight() << std::endl;
* @endcode
*/
double ImageBase::weight(const size_t n) const
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
double ImageBase::scale(const size_t n) const
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
bool ImageBase::flip(const size_t n) const
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
DataType ImageBase::datatype() const
{
    int dummy;
    MDMainHeader.getValue(MDL_DATATYPE, dummy);
    return (DataType)dummy;
}

void ImageBase::setDatatype(DataType datatype)
{
    MDMainHeader.setValue(MDL_DATATYPE, (int) datatype);
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
    MDMainHeader.getValue(MDL_SAMPLINGRATE_X, dummy);
    return dummy;
}

/** Set Euler angles in image header
 */
void ImageBase::setEulerAngles(double rot, double tilt, double psi,
                               const size_t n)
{
    MD[n].setValue(MDL_ANGLE_ROT, rot);
    MD[n].setValue(MDL_ANGLE_TILT, tilt);
    MD[n].setValue(MDL_ANGLE_PSI, psi);
}

/** Get Euler angles from image header
 */
void ImageBase::getEulerAngles(double &rot, double &tilt, double &psi,
                               const size_t n) const
{
    MD[n].getValue(MDL_ANGLE_ROT, rot);
    MD[n].getValue(MDL_ANGLE_TILT, tilt);
    MD[n].getValue(MDL_ANGLE_PSI, psi);
}

/** Set origin offsets in image header
     */
void ImageBase::setShifts(double xoff, double yoff, double zoff, const size_t n)
{
    MD[n].setValue(MDL_SHIFT_X, xoff);
    MD[n].setValue(MDL_SHIFT_Y, yoff);
    MD[n].setValue(MDL_SHIFT_Z, zoff);
}
/** Get origin offsets from image header
  */
void ImageBase::getShifts(double &xoff, double &yoff, double &zoff, const size_t n) const
{
    MD[n].getValue(MDL_SHIFT_X, xoff);
    MD[n].getValue(MDL_SHIFT_Y, yoff);
    MD[n].getValue(MDL_SHIFT_Z, zoff);
}

void ImageBase::getDimensions(size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim) const
{
    Xdim = XSIZE(*mdaBase);
    Ydim = YSIZE(*mdaBase);
    Zdim = ZSIZE(*mdaBase);
    Ndim = NSIZE(*mdaBase);
}

/** Get Image dimensions
 */
void ImageBase::getInfo(ImageInfo &imgInfo) const
{
    imgInfo.offset = offset;
    imgInfo.datatype = datatype();
    imgInfo.swap = getSwap() > 0;
    imgInfo.adim = aDimFile ;
}

void ImageBase::getInfo(const FileName &name, ImageInfo &imgInfo)
{
    read(name, HEADER);
    getInfo(imgInfo);
}

/** Open file function
  * Open the image file and returns its file hander.
  */
ImageFHandler* ImageBase::openFile(const FileName &name, int mode) const
{
    if (name.empty())
        REPORT_ERROR(ERR_PARAM_INCORRECT, "ImageBase::openFile Cannot open an empty Filename.");

    ImageFHandler* hFile = new ImageFHandler;
    FileName fileName, headName = "";
    FileName ext_name = name.getFileFormat();

    // Remove image number
    size_t dump;
    name.decompose(dump, fileName);

    fileName = fileName.removeFileFormat();

    size_t found = fileName.find_first_of("%");
    if (found!=String::npos)
        fileName = fileName.substr(0, found) ;

    hFile->exist = fileName.exists() && fileName.getFileSize() > 0;
    hFile->mode = mode;

    String wmChar;

    switch (mode)
    {
    case WRITE_READONLY:
        if (!hFile->exist)
            REPORT_ERROR(ERR_IO_NOTEXIST, formatString("Cannot read file '%s'. It doesn't exist", name.c_str()));
        wmChar = "r";
        break;
    case WRITE_OVERWRITE:
        wmChar = "w";
        break;
    case WRITE_APPEND:
    case WRITE_REPLACE:
        wmChar = (hFile->exist) ? "r+" : "w+";
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
    else if (ext_name.contains("hdf5"))
    {
        if ((hFile->fhdf5 = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)) == -1 )
            REPORT_ERROR(ERR_IO_NOTOPEN,"ImageBase::openFile: There is a problem opening the HDF5 file.");
        hFile->fimg = NULL;
        hFile->fhed = NULL;
        hFile->tif = NULL;
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
            if (mode != WRITE_READONLY || fileName.addExtension("inf").exists() )
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
        if ( (hFile->fimg = fopen(fileName.c_str(), wmChar.c_str())) == NULL )
        {
            if (errno == EACCES)
                REPORT_ERROR(ERR_IO_NOPERM,formatString("Image::openFile: permission denied when opening %s",fileName.c_str()));
            else
                REPORT_ERROR(ERR_IO_NOTOPEN,formatString("Image::openFile cannot open: %s", fileName.c_str()));
        }


        if (headName != "")
        {
            if ((hFile->fhed = fopen(headName.c_str(), wmChar.c_str()))  == NULL )
            {
                if (errno == EACCES)
                    REPORT_ERROR(ERR_IO_NOPERM,formatString("Image::openFile: permission denied when opening %s",headName.c_str()));
                else
                    REPORT_ERROR(ERR_IO_NOTOPEN,formatString("Image::openFile cannot open: %s",headName.c_str()));
            }

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
void ImageBase::closeFile(ImageFHandler* hFile) const
{
    FileName ext_name, fileName;
    FILE* fimg, *fhed;
    TIFF* tif;
    hid_t fhdf5;

    if (hFile != NULL)
    {
        fileName = hFile->fileName;
        ext_name = hFile->ext_name;
        fimg = hFile->fimg;
        fhed = hFile->fhed;
        tif  = hFile->tif;
        fhdf5 = hFile->fhdf5;
    }
    else
    {
        fileName = filename;
        ext_name = filename.getFileFormat();
        fimg = this->fimg;
        fhed = this->fhed;
        tif  = this->tif;
        fhdf5 = this->fhdf5;

    }

    if (ext_name.contains("tif"))
    {
        TIFFClose(tif);
        /* Since when creating a TIFF file without adding an image the file is 8 bytes
         * and this same file returns an error when trying to open again, we are going
         * to suppose that under 8 bytes this is empty.
        */
        if (fileName.getFileSize() < 9)
            filename.deleteFile();
    }
    else if (ext_name.contains("hdf5"))
    {
             H5Fclose(fhdf5);
    }
    else
    {
        if (fclose(fimg) != 0 )
            REPORT_ERROR(ERR_IO_NOCLOSED,(String)"Can not close image file "+ filename);

        if (fhed != NULL &&  fclose(fhed) != 0 )
            REPORT_ERROR(ERR_IO_NOCLOSED,(String)"Can not close header file of "
                         + filename);
    }
    delete hFile;
}

/* Internal read image file method.
 */
int ImageBase::_read(const FileName &name, ImageFHandler* hFile, DataMode datamode, size_t select_img,
                     bool mapData)
{
    // Temporary Error to find old select_img == -1
    if (select_img == (size_t) -1)
        REPORT_ERROR(ERR_DEBUG_TEST, "To select all images use ALL_IMAGES macro, or FIRST_IMAGE macro.");

    int err = 0;
    dataMode = datamode;

    // If MultidimArray pointer has been moved to a slice/image different from zero, then reset it.
    // This check must be done prior to mappedSize check, since mappedSlice is a trick over data pointer
    if ( virtualOffset != 0)
        movePointerTo(ALL_SLICES);
    // If Image has been previously used with mmap, then close the previous file
    if (mappedSize != 0)
        munmapFile();

    // Check whether to map the data or not
#ifdef XMIPP_MMAP

    mmapOnRead = mapData;
#endif

    FileName ext_name = hFile->ext_name;
    fimg = hFile->fimg;
    fhed = hFile->fhed;
    tif  = hFile->tif;
    fhdf5 = hFile->fhdf5;

    size_t image_num;
    name.decompose(image_num, filename);
    filename = name;
    dataFName = hFile->fileName;

    if (image_num != ALL_IMAGES)
        select_img = image_num;

#undef DEBUG
    //#define DEBUG
#ifdef DEBUG

    std::cerr << "READ\n" <<
    "name="<<name <<std::endl;
    std::cerr << "ext= "<<ext_name <<std::endl;
    std::cerr << " now reading: "<< filename <<" dataflag= "<<dataMode
    << " select_img "  << select_img << std::endl;
#endif
#undef DEBUG

    //Just clear the header before reading
    MDMainHeader.clear();
    //Set the file pointer at beginning
    if (fimg != NULL)
        fseek(fimg, 0, SEEK_SET);
    if (fhed != NULL)
        fseek(fhed, 0, SEEK_SET);

    if (ext_name.contains("spi") || ext_name.contains("xmp")  ||
        ext_name.contains("stk") || ext_name.contains("vol"))
        err = readSPIDER(select_img);
    else if (ext_name.contains("mrcs")||ext_name.contains("st"))//mrc stack MUST go BEFORE plain MRC
        err = readMRC(select_img,true);
    else if (ext_name.contains("mrc")||ext_name.contains("map"))//mrc
        err = readMRC(select_img,false);
    else if (ext_name.contains("img") || ext_name.contains("hed"))//
        err = readIMAGIC(select_img);//imagic is always an stack
    else if (ext_name.contains("ser"))//TIA
        err = readTIA(select_img,false);
    else if (ext_name.contains("dm3"))//DM3
        err = readDM3(select_img,false);
    else if (ext_name.contains("em"))//EM
        err = readEM(select_img);
    else if (ext_name.contains("pif"))//PIF
        err = readPIF(select_img);
    else if (ext_name.contains("inf"))//RAW with INF file
        err = readINF(select_img,false);
    else if (ext_name.contains("raw"))//RAW without INF file
        err = readRAW(select_img,false);
    else if (ext_name.contains("tif"))//TIFF
        err = readTIFF(select_img,false);
    else if (ext_name.contains("spe"))//SPE
        err = readSPE(select_img,false);
    else if (ext_name.contains("jpg"))//SPE
        err = readJPEG(select_img);
    else if (ext_name.contains("hdf5"))//SPE
        err = readHDF5(select_img);
    else
        err = readSPIDER(select_img);

    // Negative errors are bad.
    return err;
}

/* Internal write image file method.
 */
void ImageBase::_write(const FileName &name, ImageFHandler* hFile, size_t select_img,
                       bool isStack, int mode, CastWriteMode castMode)
{

    // Temporary Error to find old select_img == -1
    if (select_img == (size_t) -1)
        REPORT_ERROR(ERR_DEBUG_TEST, "To select all images use ALL_IMAGES macro, or FIRST_IMAGE macro.");

    int err = 0;

    // if image is mapped to file then close the file and clear
    if (mmapOnWrite && mappedSize > 0)
    {
        munmapFile();
        return;
    }

    filename = name;
    dataFName = hFile->fileName;
    _exists = hFile->exist;
    fimg = hFile->fimg;
    fhed = hFile->fhed;
    tif  = hFile->tif;

    FileName ext_name = hFile->ext_name;

    size_t aux;
    FileName filNamePlusExt;
    name.decompose(aux, filNamePlusExt);

    if (select_img == ALL_IMAGES)
        select_img = aux;

    /// Datatype info must be from filename after "%" symbol
    size_t found = filNamePlusExt.find_first_of("%");
    String imParam = "";
    if (found!=String::npos)
    {
        imParam =  filNamePlusExt.substr(found+1).c_str();
        filNamePlusExt = filNamePlusExt.substr(0, found) ;
    }

    //#define DEBUG
#ifdef DEBUG

    std::cerr << "write" <<std::endl;
    std::cerr<<"extension for write= "<<ext_name<<std::endl;
    std::cerr<<"filename= "<<filename<<std::endl;
    std::cerr<<"mode= "<<mode<<std::endl;
    std::cerr<<"isStack= "<<isStack<<std::endl;
    std::cerr<<"select_img= "<<select_img<<std::endl;
#endif
#undef DEBUG
    // Check that image is not empty
    if (getSize() < 1)
        REPORT_ERROR(ERR_MULTIDIM_EMPTY,"write Image ERROR: image is empty!");

    replaceNsize = 0;//reset replaceNsize in case image is reused
    if(isStack && select_img == ALL_IMAGES && mode == WRITE_REPLACE)
        REPORT_ERROR(ERR_VALUE_INCORRECT,"Please specify object to be replaced");
    else if (_exists && (mode == WRITE_REPLACE || mode == WRITE_APPEND))
    {
        // CHECK FOR INCONSISTENCIES BETWEEN data.xdim and x, etc???
        size_t Xdim, Ydim, Zdim, _Xdim, _Ydim, _Zdim, Ndim, _Ndim;
        Xdim = Ydim = Zdim = _Xdim = _Ydim = _Zdim = Ndim = _Ndim = 0;
        Image<char> auxI;
        auxI._read(filNamePlusExt, hFile, HEADER, ALL_IMAGES);

        this->getDimensions(Xdim, Ydim, Zdim, Ndim);
        auxI.getDimensions(_Xdim, _Ydim, _Zdim, _Ndim);

        if(auxI.getSize()>1)
        {
            replaceNsize = _Ndim;

            /** If we are going to changes all images, then swap of the file may be changed,
             *  otherwise, original swap remains. */
            if (select_img > ALL_IMAGES || Ndim < replaceNsize)
                swapWrite = auxI.swap;

            if(Xdim != _Xdim ||
               Ydim != _Ydim ||
               Zdim != _Zdim)
            {
                REPORT_ERROR(ERR_MULTIDIM_SIZE,formatString(
                                 "ImageBase::Write: images source and target have different sizes:\n"
                                 "Image source to be written (x,y,z,n) = %d %d %d %lu\n"
                                 "Image file target %s (x,y,z,n) = %d %d %d %lu",
                                 Xdim,Ydim,Zdim,Ndim,dataFName.c_str(),_Xdim,_Ydim,_Zdim,_Ndim));
            }
        }
    }
    else if(!_exists && mode == WRITE_APPEND)
    {
        ;
    }
    else if (mode == WRITE_READONLY)//If new file we are in the WRITE_OVERWRITE mode
    {
        REPORT_ERROR(ERR_ARG_INCORRECT, formatString("File %s  opened in read-only mode. Cannot write.", name.c_str()));
    }
    /*
     * SELECT FORMAT
     */
    //Set the file pointer at beginning
    if (fimg != NULL)
        fseek(fimg, 0, SEEK_SET);
    if (fhed != NULL)
        fseek(fhed, 0, SEEK_SET);

    if(ext_name.contains("spi") || ext_name.contains("xmp") ||
       ext_name.contains("vol"))
        err = writeSPIDER(select_img,isStack,mode);
    else if (ext_name.contains("stk"))
        err = writeSPIDER(select_img,true,mode);
    //    else if (ext_name.contains("mrcs"))
    //        writeMRC(select_img,true,mode,imParam,castMode);
    else if (ext_name.contains("mrc")||ext_name.contains("map")
             ||ext_name.contains("mrcs")||ext_name.contains("st"))
        writeMRC(select_img,isStack,mode,imParam,castMode);
    else if (ext_name.contains("img") || ext_name.contains("hed"))
        writeIMAGIC(select_img,mode,imParam,castMode);
    else if (ext_name.contains("dm3"))
        writeDM3(select_img,false,mode);
    else if (ext_name.contains("em"))
        writeEM(select_img,false,mode);
    else if (ext_name.contains("pif"))
        writePIF(select_img,false,mode);
    else if (ext_name.contains("ser"))
        writeTIA(select_img,false,mode);
    else if (ext_name.contains("raw") || ext_name.contains("inf"))
        writeINF(select_img,false,mode,imParam,castMode);
    else if (ext_name.contains("tif"))
        writeTIFF(select_img,isStack,mode,imParam,castMode);
    else if (ext_name.contains("spe"))
        writeSPE(select_img,isStack,mode);
    else if (ext_name.contains("jpg"))
        writeJPEG(select_img);
    else if (ext_name.contains("hdf5"))
        writeHDF5(select_img);
    else
        err = writeSPIDER(select_img,isStack,mode);

    if ( err < 0 )
    {
        std::cerr << " Filename = " << filename << " Extension= " << ext_name << std::endl;
        REPORT_ERROR(ERR_IO_NOWRITE, "Error writing file");
    }

    /* If initially the file did not existed, once the first image is written,
     * then the file exists
     */
    if (!_exists)
        hFile->exist = _exists = true;
}


/** Show image properties
      */
std::ostream& operator<<(std::ostream& o, const ImageBase& I)
{
    o << std::endl;
    DataType * fileDT = NULL;
    if (!I.filename.empty())
    {
        o << "--- File information ---" << std::endl;
        o << "Filename       : " << I.filename << std::endl;
        o << "Endianess      : ";
        if (I.swap^IsLittleEndian())
            o << "Little"  << std::endl;
        else
            o << "Big" << std::endl;

        o << "Reversed       : ";
        if (I.swap)
            o << "True"  << std::endl;
        else
            o << "False" << std::endl;
        fileDT = new DataType;
        *fileDT = I.datatype();
        o << "Data type      : " << datatype2StrLong(*fileDT) << std::endl;
        o << "Data offset    : " << I.offset << std::endl;
    }

    o << "--- Image information ---" << std::endl;

    DataType myDT = I.myT();
    if ((fileDT == NULL || myDT != *fileDT) && I.dataMode >= DATA )
        o << "Memory datatype: " << datatype2StrLong(I.myT()) << std::endl;
    o << "Image type     : ";
    if (I.isComplex())
        o << "Fourier-space image" << std::endl;
    else
        o << "Real-space image" << std::endl;

    size_t xdim, ydim, zdim, ndim;
    I.getDimensions(xdim, ydim, zdim, ndim);
    o << "Dimensions     : " << ndim << " x " << zdim << " x " << ydim << " x " << xdim;
    o << "  ((N)Objects x (Z)Slices x (Y)Rows x (X)Columns)" << std::endl;

    double sampling;
    I.MDMainHeader.getValue(MDL_SAMPLINGRATE_X, sampling);
    if (sampling > 1e-300)
    {
        o << "Sampling rate  : " << std::endl;
        o << "                 X-rate (Angstrom/pixel) = " << sampling << std::endl;
        I.MDMainHeader.getValue(MDL_SAMPLINGRATE_Y, sampling);
        o << "                 Y-rate (Angstrom/pixel) = " << sampling << std::endl;
        I.MDMainHeader.getValue(MDL_SAMPLINGRATE_Z, sampling);
        o << "                 Z-rate (Angstrom/pixel) = " << sampling << std::endl;
    }

    std::stringstream oGeo;

    if (I.individualContainsLabel(MDL_ANGLE_ROT))
    {
        oGeo << "Euler angles   : " << std::endl;
        oGeo << "                 Phi   (rotation around Z axis)                  = " << I.rot() << std::endl;
        oGeo << "                 Theta (tilt, second rotation around new Y axis) = " << I.tilt() << std::endl;
        oGeo << "                 Psi   (third rotation around new Z axis)        = " << I.psi() << std::endl;
    }
    if (I.individualContainsLabel(MDL_SHIFT_X))
    {
        oGeo << "Origin Offsets : " << std::endl;
        oGeo << "                 Xoff  (origin offset in X-direction) = " << I.Xoff() << std::endl;
        oGeo << "                 Yoff  (origin offset in Y-direction) = " << I.Yoff() << std::endl;
        oGeo << "                 Zoff  (origin offset in Z-direction) = " << I.Zoff() << std::endl;
    }
    if (I.individualContainsLabel(MDL_SCALE))
        oGeo << "Scale          : " <<I.scale() << std::endl;
    if (I.individualContainsLabel(MDL_WEIGHT))
        oGeo << "Weight         : " << I.weight() << std::endl;
    if (I.individualContainsLabel(MDL_FLIP))
        oGeo << "Flip           : " << I.flip() << std::endl;

    if (!oGeo.str().empty())
        o << "--- Geometry ---" << std::endl << oGeo.str();

    return o;
}
