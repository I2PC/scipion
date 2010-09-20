#include "image.h"

template<typename T>
Image<T>::Image()
{
    mmapOn = false;
    clear();
}

template<typename T>
Image<T>::Image(int Xdim, int Ydim, int Zdim, int Ndim)
{
    mmapOn = false;
    clear();
    data.resize(Ndim, Zdim, Ydim, Xdim);
    MD.resize(Ndim);
}

/** Clear.
  * Initialize everything to 0
  */
template<typename T>
void Image<T>::clear()
{
    if (mmapOn)
    {
        munmap(data.data-offset,mappedSize);
        close(mFd);
        data.data = NULL;
    }
    else
        data.clear();

    dataflag = -1;
    if (isComplexT())
        transform = Standard;
    else
        transform = NoTransform;
    i = 0;
    filename = "";
    offset = 0;
    swap = 0;
    clearHeader();
    replaceNsize=0;
    mmapOn = false;
}

/** Clear the header of the image
 */
template<typename T>
void Image<T>::clearHeader()
{
    MDMainHeader.clear();
    MD.clear();
}

/** Check whether image is complex based on T
 */
template<typename T>
bool Image<T>::isComplexT() const
{
    return ( typeid(T) == typeid(std::complex<double>) ||
             typeid(T) == typeid(std::complex<float>) );
}

/** Check whether image is complex based on transform
  */
template<typename T>
bool Image<T>::isComplex() const
{
    return !(transform==NoTransform);
}

/** Destructor.
 */
template<typename T>
Image<T>::~Image()
{
    clear();
}

/** Specific read functions for different file formats
  */

/** Is this file an image
 *
 *  Check whether a real-space image can be read
 *
 */
template<typename T>
bool Image<T>::isImage(const FileName &name)
{
    return !read(name, false);
}

/** Is this file a real-valued image
 *
 *  Check whether a real-space image can be read
 *
 */
template<typename T>
bool Image<T>::isRealImage(const FileName &name)
{
    return (isImage(name) && !isComplex());
}

/** Is this file a complex image
 *
 *  Check whether a fourier-space (complex) image can be read
 *
 */
template<typename T>
bool Image<T>::isComplexImage(const FileName &name)
{
    return (isImage(name) && isComplex());
}

/** Rename the image
  */
template<typename T>
void Image<T>::rename (const FileName &name)
{
    filename = name;
}

/** General read function
 * you can read a single image from a single image file
 * or a single image file from an stack, in the second case
 * the select slide may come in the image name or in the select_img parameter
 * file name takes precedence over select_img
 * If -1 is given the whole object is read
 *
 */
template<typename T> int Image<T>::read(const FileName &name, bool readdata, int select_img,
         bool apply_geo, bool only_apply_shifts, MDRow * row, bool mapData)
{
    //const MetaData &docFile = *docFilePtr;
    //std::vector<MDLabel> &activeLabels = *activeLabelsPtr;

    int err = 0;
    // Check whether to read the data or only the header
    dataflag = ( readdata ) ? 1 : -1;

    // Check whether to map the data or not
    mmapOn = mapData;

    FileName ext_name = name.get_file_format();
    size_t found;
    filename = name;
    found = filename.find_first_of("@");
    if (found != std::string::npos)
    {
        select_img =  atoi(filename.substr(0, found).c_str());
        filename = filename.substr(found+1) ;
    }

    double imParam = NULL;
    found = filename.find_first_of("%");
    if (found != std::string::npos)
    {
        imParam = atof(filename.substr(found+1).c_str());
        filename = filename.substr(0, found) ;
    }

    filename = filename.remove_file_format();

    if(ext_name.contains("inf"))
        filename = filename.without_extension();
    else if (ext_name.contains("raw") && exists(filename.add_extension("inf")))
        ext_name = "inf";


#undef DEBUG
    //#define DEBUG
#ifdef DEBUG

    std::cerr << "READ\n" <<
    "name="<<name <<std::endl;
    std::cerr << "ext= "<<ext_name <<std::endl;
    std::cerr << " now reading: "<< filename <<" dataflag= "<<dataflag
    << " select_img "  << select_img << std::endl;
#endif
#undef DEBUG

    //Just clear the header before reading
    MDMainHeader.clear();

    if (ext_name.contains("spi") || ext_name.contains("xmp") )//mrc stack MUST go BEFORE plain MRC
        err = readSPIDER(select_img,true);
    else if (ext_name.contains("mrcs"))//mrc stack MUST go BEFORE plain MRC
        err = readMRC(select_img,true);
    else if (ext_name.contains("mrc"))//mrc
        err = readMRC(select_img,false);
    else if (ext_name.contains("img") || ext_name.contains("hed"))//
        err = readIMAGIC(select_img);//imagic is always an stack
    else if (ext_name.contains("ser"))//TIA
        err = readTIA(select_img,false, imParam);
    else if (ext_name.contains("dm3"))//DM3
        err = readDM3(select_img,false);
    else if (ext_name.contains("inf"))//RAW with INF file
        err = readINF(select_img,false);
    else if (ext_name.contains("raw"))//RAW without INF file
        err = readRAW(select_img,false);
    else if (ext_name.contains("tif") || ext_name.contains("tiff"))//TIFF
        err = readTIFF(select_img,false);
    else if (ext_name.contains("spe"))//SPE
        err = readSPE(select_img,false);
    else
        err = readSPIDER(select_img,true);

    //This implementation does not handle stacks,
    //read in a block
    if (row != NULL)
    {
        if (data.ndim != 1)
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "Header overwriting not available for stacks!!!");
        MDLabel label;

        for (MDRow::const_iterator it = row->begin(); it != row->end(); ++it)
        {
            label = (*it)->label;
            if (MD[select_img].containsLabel(label))
                *(MD[select_img].getObject(label)) = *(*it);
            else
                MD[select_img].push_back(new MDObject(*(*it)));
        }
    }

    //apply geo has not been defined for volumes
    if(this->data.getDim()>2)
        apply_geo=false;

    if (readdata && (apply_geo || only_apply_shifts))
    {
        Matrix2D< double > A = getTransformationMatrix(only_apply_shifts);
        if (!A.isIdentity())
        {
            MultidimArray<T> tmp = (*this)();
            applyGeometry(BSPLINE3, (*this)(), tmp, A, IS_INV, WRAP);
        }
    }

    // Negative errors are bad.
    return err;
}

/** General write function
 * select_img= which slice should I replace
 * overwrite = 0, append slice
 * overwrite = 1 overwrite slice
 */
template<typename T>
void Image<T>::write(FileName name, int select_img, bool isStack, int mode)
{
    int err = 0;

    if (name == "")
        name = filename;

    FileName ext_name = name.get_file_format();
    size_t found;
    filename = name;
    found=filename.find_first_of("@");
    FileName filNamePlusExt;
    if (found!=std::string::npos)
    {
        //select_img = atoi(filename.substr(0, found).c_str());
        filename   =      filename.substr(found+1) ;
    }
    int imParam = NULL;
    found=filename.find_first_of("%");
    if (found!=std::string::npos)
    {
        imParam =  atoi(filename.substr(found+1).c_str());
        filename = filename.substr(0, found) ;
    }
    filNamePlusExt = filename;
    found=filename.find_first_of(":");
    if ( found!=std::string::npos)
        filename   = filename.substr(0, found);

    if(ext_name.contains("inf"))
        filename = filename.without_extension();

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

    // CHECK FOR INCONSISTENCIES BETWEEN data.xdim and x, etc???
    int Xdim, Ydim, Zdim, Ndim;
    this->getDimensions(Xdim,Ydim, Zdim, Ndim);

    _exists = exists(filename);
    Image<T> auxI;
    replaceNsize=0;//reset replaceNsize in case image is reused
    if(select_img==-1 && mode==WRITE_REPLACE)
        REPORT_ERROR(ERR_VALUE_INCORRECT,"writeSPIDER: Please specify object to be replaced");
    else if(!_exists && mode==WRITE_REPLACE)
    {
        std:: stringstream replace_number;
        replace_number << select_img;
        REPORT_ERROR(ERR_IO_NOTEXIST,(std::string)"Cannot replace object number: "
                     + replace_number.str()
                     + " in file " +filename
                     + ". It does not exist");
    }
    else if (_exists && (mode==WRITE_REPLACE || mode==WRITE_APPEND))
    {
        auxI.dataflag = -2;
        auxI.read(filNamePlusExt,false);
        int _Xdim, _Ydim, _Zdim, _Ndim;
        auxI.getDimensions(_Xdim,_Ydim, _Zdim, _Ndim);
        replaceNsize=_Ndim;
        if(Xdim!=_Xdim ||
           Ydim!=_Ydim ||
           Zdim!=_Zdim
          )
            REPORT_ERROR(ERR_MULTIDIM_SIZE,"write: target and source objects have different size");
        if(mode==WRITE_REPLACE && select_img>_Ndim)
            REPORT_ERROR(ERR_VALUE_INCORRECT,"write: cannot replace image stack is not large enough");
        if(auxI.replaceNsize <1 &&
           (mode==WRITE_REPLACE || mode==WRITE_APPEND))
            REPORT_ERROR(ERR_IO,"write: output file is not an stack");
    }
    else if(!_exists && mode==WRITE_APPEND)
    {
        ;
    }
    else//If new file we are in the WRITE_OVERWRITE mode
    {
        mode=WRITE_OVERWRITE;
    }
    /*
     * SELECT FORMAT
     */

    if(ext_name.contains("spi") || ext_name.contains("xmp") ||
       ext_name.contains("stk") || ext_name.contains("vol"))
        err = writeSPIDER(select_img,isStack,mode);
    else if (ext_name.contains("mrcs"))
        writeMRC(select_img,true,mode);
    else if (ext_name.contains("mrc"))
        writeMRC(select_img,false,mode);
    else if (ext_name.contains("img") || ext_name.contains("hed"))
        writeIMAGIC(select_img,mode);
    else if (ext_name.contains("dm3"))
        writeDM3(select_img,false,mode);
    else if (ext_name.contains("ser"))
        writeTIA(select_img,false,mode);
    else if (ext_name.contains("raw") || ext_name.contains("inf"))
        writeINF(select_img,false,mode);
    else if (ext_name.contains("tif") || ext_name.contains("tiff"))
        writeTIFF(select_img,isStack,mode,imParam);
    else if (ext_name.contains("spe"))
        writeSPE(select_img,isStack,mode);
    else
        err = writeSPIDER(select_img,isStack,mode);

    if ( err < 0 )
    {
        std::cerr << " Filename = " << filename << " Extension= " << ext_name << std::endl;
        REPORT_ERROR(ERR_IO_NOWRITE, "Error writing file");
    }
    //unlock file
}

/** Cast a page of data from type dataType to type Tdest
 *    input pointer  char *
 */
template<typename T>
void Image<T>::castPage2T(char * page, T * ptrDest, DataType datatype, size_t pageSize )
{
    switch (datatype)
    {
    case Unknown_Type:
        REPORT_ERROR(ERR_TYPE_INCORRECT,"ERROR: datatype is Unknown_Type");
    case UChar:
        {
            if (typeid(T) == typeid(unsigned char))
                memcpy(ptrDest, page, pageSize*sizeof(T));
            else
            {
                unsigned char * ptr = (unsigned char *) page;
                for(int i=0; i<pageSize; i++)
                    ptrDest[i]=(T) ptr[i];
            }
            break;
        }
    case SChar:
            {
                if (typeid(T) == typeid(signed char))
            {
                memcpy(ptrDest, page, pageSize*sizeof(T));
                }
                else
                {
                    signed char * ptr = (signed char *) page;
                    for(int i=0; i<pageSize; i++)
                        ptrDest[i]=(T) ptr[i];
                }
            break;
        }
    case UShort:
            {
                if (typeid(T) == typeid(unsigned short))
            {
                memcpy(ptrDest, page, pageSize*sizeof(T));
                }
                else
                {
                    unsigned short * ptr = (unsigned short *) page;
                    for(int i=0; i<pageSize; i++)
                        ptrDest[i]=(T) ptr[i];
                }
            break;
        }
    case Short:
            {
                if (typeid(T) == typeid(short))
            {
                memcpy(ptrDest, page, pageSize*sizeof(T));
                }
                else
                {
                    short * ptr = (short *) page;
                    for(int i=0; i<pageSize; i++)
                        ptrDest[i]=(T) ptr[i];
                }
            break;
        }
    case UInt:
            {
                if (typeid(T) == typeid(unsigned int))
            {
                memcpy(ptrDest, page, pageSize*sizeof(T));
                }
                else
                {
                    unsigned int * ptr = (unsigned int *) page;
                    for(int i=0; i<pageSize; i++)
                        ptrDest[i]=(T) ptr[i];
                }
            break;
        }
    case Int:
            {
                if (typeid(T) == typeid(int))
            {
                memcpy(ptrDest, page, pageSize*sizeof(T));
                }
                else
                {
                    int * ptr = (int *) page;
                    for(int i=0; i<pageSize; i++)
                        ptrDest[i]=(T) ptr[i];
                }
            break;
        }
    case Long:
            {
                if (typeid(T) == typeid(long))
            {
                memcpy(ptrDest, page, pageSize*sizeof(T));
                }
                else
                {
                    long * ptr = (long *) page;
                    for(int i=0; i<pageSize; i++)
                        ptrDest[i]=(T) ptr[i];
                }
            break;
        }
    case Float:
            {
                if (typeid(T) == typeid(float))
            {
                memcpy(ptrDest, page, pageSize*sizeof(T));
                }
                else
                {
                    float * ptr = (float *) page;
                    for(int i=0; i<pageSize; i++)
                        ptrDest[i]=(T) ptr[i];
                }
            break;
        }
    case Double:
            {
                if (typeid(T) == typeid(double))
            {
                memcpy(ptrDest, page, pageSize*sizeof(T));
                }
                else
                {
                    double * ptr = (double *) page;
                    for(int i=0; i<pageSize; i++)
                        ptrDest[i]=(T) ptr[i];
                }
            break;
        }
    default:
            {
                std::cerr<<"Datatype= "<<datatype<<std::endl;
                REPORT_ERROR(ERR_TYPE_INCORRECT," ERROR: cannot cast datatype to T");
                break;
            }
        }

}

/** Cast page from T to datatype
 *  input pointer char *
 */
template<typename T>
void Image<T>::castPage2Datatype(T * srcPtr, char * page, DataType datatype, size_t pageSize )
{
    switch (datatype)
    {
    case Float:
        {
            if (typeid(T) == typeid(float))
            {
                memcpy(page, srcPtr, pageSize*sizeof(T));
            }
            else
            {
                float * ptr = (float *) page;
                for(int i=0; i<pageSize; i++)
                    ptr[i] = (float)srcPtr[i];
            }
            break;
        }
    case Double:
            {
                if (typeid(T) == typeid(double))
            {
                memcpy(page, srcPtr, pageSize*sizeof(T));
                }
                else
                {
                    double * ptr = (double *) page;
                    for(int i=0; i<pageSize; i++)
                        ptr[i] = (double)srcPtr[i];
                }
            break;
        }
    case UShort:
            {
                if (typeid(T) == typeid(unsigned short))
            {
                memcpy(page, srcPtr, pageSize*sizeof(T));
                }
                else
                {
                    unsigned short * ptr = (unsigned short *) page;
                    for(int i=0; i<pageSize; i++)
                        ptr[i] = (unsigned short)srcPtr[i];
                }
            break;
        }
    case UChar:
            {
                if (typeid(T) == typeid(unsigned char))
            {
                memcpy(page, srcPtr, pageSize*sizeof(T));
                }
                else
                {
                    unsigned char * ptr = (unsigned char *) page;
                    for(int i=0; i<pageSize; i++)
                        ptr[i] = (unsigned char)srcPtr[i];
                }
            break;
        }
    default:
            {
                std::cerr<<"outputDatatype= "<<datatype<<std::endl;
                REPORT_ERROR(ERR_TYPE_INCORRECT," ERROR: cannot cast T to outputDatatype");
                break;
            }
        }
}

/** Check file Datatype is same as T type to use mmap.
 */
template<typename T>
bool Image<T>::checkMmapT(DataType datatype)
{

    switch (datatype)
    {
    case Unknown_Type:
        REPORT_ERROR(ERR_TYPE_INCORRECT,"ERROR: datatype is Unknown_Type");
    case UChar:
        {
            if (typeid(T) == typeid(unsigned char))
                return 1;
            else
                return 0;
            break;
        }
    case SChar:
        {
            if (typeid(T) == typeid(signed char))
                return 1;
            else
                return 0;
            break;
        }
    case UShort:
        {
            if (typeid(T) == typeid(unsigned short))
                return 1;
            else
                return 0;
            break;
        }
    case Short:
        {
            if (typeid(T) == typeid(short))
                return 1;
            else
                return 0;
            break;
        }
    case UInt:
        {
            if (typeid(T) == typeid(unsigned int))
                return 1;
            else
                return 0;
            break;
        }
    case Int:
        {
            if (typeid(T) == typeid(int))
                return 1;
            else
                return 0;
            break;
        }
    case Long:
        {
            if (typeid(T) == typeid(long))
                return 1;
            else
                return 0;
            break;
        }
    case Float:
        {
            if (typeid(T) == typeid(float))
                return 1;
            else
                return 0;
            break;
        }
    case Double:
        {
            if (typeid(T) == typeid(double))
                return 1;
            else
                return 0;
            break;
        }
    default:
        {
            std::cerr<<"Datatype= "<<datatype<<std::endl;
            REPORT_ERROR(ERR_TYPE_INCORRECT," ERROR: cannot cast datatype to T");
            break;
        }
    }
    //               int * iTemp = (int*) map;
    //                ptrDest = reinterpret_cast<T*> (iTemp);
}

/** Write an entire page as datatype
 *
 * A page of datasize_n elements T is cast to datatype and written to fimg
 * The memory for the casted page is allocated and freed internally.
 */
template<typename T>
void Image<T>::writePageAsDatatype(FILE * fimg, DataType datatype, size_t datasize_n )
{
    size_t datasize = datasize_n * gettypesize(datatype);
    char * fdata = (char *) askMemory(datasize);
    castPage2Datatype(MULTIDIM_ARRAY(data), fdata, datatype, datasize_n);
    fwrite( fdata, datasize, 1, fimg );
    freeMemory(fdata, datasize);
}

/** Swap an entire page
  * input pointer char *
  */
template<typename T>
void Image<T>::swapPage(char * page, size_t pageNrElements, DataType datatype)
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

/** Read the raw data
  */
template<typename T>
void Image<T>::readData(FILE* fimg, int select_img, DataType datatype, unsigned long pad)
{
    //#define DEBUG
#ifdef DEBUG
    std::cerr<<"entering readdata"<<std::endl;
    std::cerr<<" readData flag= "<<dataflag<<std::endl;
#endif

    if ( dataflag < 1 )
        return;

    // If only half of a transform is stored, it needs to be handled
    if (transform == Hermitian || transform == CentHerm )
        data.setXdim(XSIZE(data)/2 + 1);

    size_t myoffset, readsize, readsize_n, pagemax = 1073741824; //1Gb
    size_t datatypesize=gettypesize(datatype);
    size_t pagesize  =ZYXSIZE(data)*datatypesize;
    size_t haveread_n=0;

    //Multidimarray mmapOn is priority over image mmapOn
    if(data.mmapOn)
        mmapOn = false;

    // Flag to know that data is not going to be mapped although mmapOn is true
    if (mmapOn && !checkMmapT(datatype))
    {
        std::cout << "WARNING: Image Class. File datatype and image declaration not compatible with mmap. Loading into memory." <<std::endl;
        mmapOn = false;
        mFd = -1;
    }

    if (mmapOn)
    {
        if ( NSIZE(data) > 1 )
        {
            REPORT_ERROR(ERR_MMAP,"Image Class::ReadData: mmap with multiple \
                         images file not compatible. Try selecting a unique image.");
        }

        fclose(fimg);

        if ( ( mFd = open(filename.c_str(), O_RDWR, S_IREAD | S_IWRITE) ) == -1 )
            REPORT_ERROR(ERR_IO_NOTOPEN,"Image Class::ReadData: Error opening the image file.");

        char * map;
        mappedSize = pagesize+offset;

        if ( (map = (char*) mmap(0,mappedSize, PROT_READ | PROT_WRITE, MAP_SHARED, mFd, 0)) == (void*) -1 )
            REPORT_ERROR(ERR_MMAP_NOTADDR,"Image Class::ReadData: mmap of image file failed.");
        data.data = reinterpret_cast<T*> (map+offset);
    }
    else
    {
        // Reset select to get the correct offset
        if ( select_img < 0 )
            select_img = 0;

        char* page = NULL;

        // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
        //if memory already allocated use it (no resize allowed)
        data.coreAllocateReuse();
        myoffset = offset + select_img*(pagesize + pad);
        //#define DEBUG

#ifdef DEBUG

        data.printShape();
        printf("DEBUG: Page size: %ld offset= %d \n", pagesize, offset);
        printf("DEBUG: Swap = %d  Pad = %ld  Offset = %ld\n", swap, pad, offset);
        printf("DEBUG: myoffset = %d select_img= %d \n", myoffset, select_img);
#endif

        if (pagesize > pagemax)
            page = (char *) askMemory(pagemax*sizeof(char));
        else
            page = (char *) askMemory(pagesize*sizeof(char));

        fseek( fimg, myoffset, SEEK_SET );
        for ( size_t myn=0; myn<NSIZE(data); myn++ )
        {
            for (size_t myj=0; myj<pagesize; myj+=pagemax )//pagesize size of object
            {
                // Read next page. Divide pages larger than pagemax
                readsize = pagesize - myj;
                if ( readsize > pagemax )
                    readsize = pagemax;
                readsize_n = readsize/datatypesize;

                //Read page from disc
                fread( page, readsize, 1, fimg );
                //swap per page
                if (swap)
                    swapPage(page, readsize, datatype);
                // cast to T per page
                castPage2T(page, MULTIDIM_ARRAY(data) + haveread_n, datatype, readsize_n);
                haveread_n += readsize_n;
            }
            if ( pad > 0 )
                //fread( padpage, pad, 1, fimg);
                fseek( fimg, pad, SEEK_CUR );
        }
        //if ( pad > 0 )
        //    freeMemory(padpage, pad*sizeof(char));
        if ( page > 0 )
            freeMemory(page, pagesize*sizeof(char));

#ifdef DEBUG

        printf("DEBUG img_read_data: Finished reading and converting data\n");
#endif

    }
    return;
}

/** Data access
 *
 * This operator can be used to access the data multidimarray.
 * In this way we could resize an image just by
 * resizing its associated matrix or we could add two images by adding their
 * matrices.
 * @code
 * I().resize(128, 128);
 * I2() = I1() + I2();
 * @endcode
 */
template<typename T>
MultidimArray<T>& Image<T>::operator()()
{
    return data;
}
template<typename T>
const MultidimArray<T>& Image<T>::operator()() const
{
    return data;
}

/** Pixel access
*
* This operator is used to access a pixel within a 2D image. This is a
* logical access, so you could access to negative positions if the image
* has been defined so (see the general explanation for the class).
*
* @code
* std::cout << "Grey level of pixel (-3,-3) of the image = " << I(-3, -3)
* << std::endl;
*
* I(-3, -3) = I(-3, -2);
* @endcode
*/
template<typename T>
T& Image<T>::operator()(int i, int j) const
{
    return A2D_ELEM(data, i, j);
}
/** Set pixel
 * (direct access) needed by swig
 */
template<typename T>
void Image<T>::setPixel(int i, int j, T v)
{
    IMGPIXEL(*this,i,j)=v;
}

/** Get pixel
 * (direct acces) needed by swig
 */
template<typename T>
T Image<T>::getPixel(int i, int j) const
{
    return IMGPIXEL(*this,i,j);
}

/** Voxel access
 *
 * This operator is used to access a voxel within a 3D image. This is a
 * logical access, so you could access to negative positions if the image
 * has been defined so (see the general explanation for the class).
 *
 * @code
 * std::cout << "Grey level of pixel (-3,-3, 1) of the volume = " << I(-3, -3, 1)
 * << std::endl;
 *
 * I(-3, -3, 1) = I(-3, -2, 0);
 * @endcode
 */
template<typename T>
T& Image<T>::operator()(int k, int i, int j) const
{
    return A3D_ELEM(data, k, i, j);
}

/** Get file name
 *
 * @code
 * std::cout << "Image name = " << I.name() << std::endl;
 * @endcode
 */
template<typename T>
const FileName & Image<T>::name() const
{
    return filename;
}

/** Get Image dimensions
 */
template<typename T>
void Image<T>::getDimensions(int &Xdim, int &Ydim, int &Zdim, int &Ndim) const
{
    Xdim = XSIZE(data);
    Ydim = YSIZE(data);
    Zdim = ZSIZE(data);
    Ndim = NSIZE(data);
}

template<typename T>
long unsigned int Image<T>::getSize() const
{
    return NZYXSIZE(data);
}

/** Get Image offset and swap
 */
template<typename T>
void Image<T>::getOffsetAndSwap(unsigned long &_offset, int &_swap) const
{
    _offset = offset;
    _swap = swap;
}

/* Is there label in the individual header */
template<typename T>
bool Image<T>::individualContainsLabel(MDLabel label) const
{
    return MD[0].containsLabel(label);
}

/* Is there label in the main header */
template<typename T>
bool Image<T>::mainContainsLabel(MDLabel label) const
{
    return MDMainHeader.containsLabel(label);
}

/** Get Rot angle
*
* @code
* std::cout << "First Euler angle " << I.rot() << std::endl;
* @endcode
*/
template<typename T>
double Image<T>::rot(const long int n) const
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
template<typename T>
double Image<T>::tilt(const long int n) const
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
template<typename T>
double Image<T>::psi(const long int n) const
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
template<typename T>
double Image<T>::Xoff(const long int n) const
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
template<typename T>
double Image<T>::Yoff(const long int n) const
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
template<typename T>
double Image<T>::Zoff(const long int n) const
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
template<typename T>
double Image<T>::weight(const long int n) const
{
    double dummy = 1;
    MD[n].getValue(MDL_WEIGHT, dummy);
    return dummy;
}

/** Get Flip
*
* @code
* std::cout << "flip= " << flip() << std::endl;
* @endcode
*/
template<typename T>
bool Image<T>::flip(const long int n) const
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
template<typename T>
int Image<T>::dataType() const
{
    int dummy;
    MDMainHeader.getValue(MDL_DATATYPE, dummy);
    return dummy;
}

/** Sampling RateX
*
* @code
* std::cout << "sampling= " << samplingRateX() << std::endl;
* @endcode
*/
template<typename T>
double Image<T>::samplingRateX() const
{
    double dummy = 1.;
    MDMainHeader.getValue(MDL_SAMPLINGRATEX, dummy);
    return dummy;
}

/** Set file name
 */
template<typename T>
void Image<T>::setName(const FileName &_filename)
{
    filename = _filename;
}

/** Set Euler angles in image header
 */
template<typename T>
void Image<T>::setEulerAngles(double rot, double tilt, double psi,
                              long int n)
{
    MD[n].setValue(MDL_ANGLEROT, rot);
    MD[n].setValue(MDL_ANGLETILT, tilt);
    MD[n].setValue(MDL_ANGLEPSI, psi);
}

/** Get Euler angles from image header
 */
template<typename T>
void Image<T>::getEulerAngles(double &rot, double &tilt, double &psi,
                              long int n)
{
    MD[n].getValue(MDL_ANGLEROT, rot);
    MD[n].getValue(MDL_ANGLETILT, tilt);
    MD[n].getValue(MDL_ANGLEPSI, psi);
}

/** Set Rotation angle to image */
template<typename T>
void Image<T>::setRot(double rot, long int n)
{
    MD[n].setValue(MDL_ANGLEROT, rot);
}

/** Set Tilt angle to image */
template<typename T>
void Image<T>::setTilt(double tilt, long int n)
{
    MD[n].setValue(MDL_ANGLETILT, tilt);
}

/** Set Rotation angle to image */
template<typename T>
void Image<T>::setPsi(double psi, long int n)
{
    MD[n].setValue(MDL_ANGLEPSI, psi);
}

/** Set origin offsets in image header
 */
template<typename T>
void Image<T>::setShifts(double xoff, double yoff, double zoff,
                         long int n)
{
    MD[n].setValue(MDL_ORIGINX, xoff);
    MD[n].setValue(MDL_ORIGINY, yoff);
    MD[n].setValue(MDL_ORIGINZ, zoff);
}
/** Get origin offsets from image header
  */
template<typename T>
void Image<T>::getShifts(double &xoff, double &yoff, double &zoff,
                         long int n)
{
    MD[n].getValue(MDL_ORIGINX, xoff);
    MD[n].getValue(MDL_ORIGINY, yoff);
    MD[n].getValue(MDL_ORIGINZ, zoff);
}

/** Set X offset in image header
 */
template<typename T>
void Image<T>::setXoff(double xoff, long int n)
{
    MD[n].setValue(MDL_ORIGINX, xoff);
}

/** Set Y offset in image header
 */
template<typename T>
void Image<T>::setYoff(double yoff, long int n)
{
    MD[n].setValue(MDL_ORIGINY, yoff);
}

/** Set Z offset in image header
 */
template<typename T>
void Image<T>::setZoff(double zoff, long int n)
{
    MD[n].setValue(MDL_ORIGINZ, zoff);
}

/** Set flip in image header
 */
template<typename T>
void Image<T>::setFlip(bool flip, long int n)
{
    MD[n].setValue(MDL_FLIP, flip);
}

/** Set Weight in image header
*/
template<typename T>
void Image<T>::setWeight(double weight, long int n)
{
    MD[n].setValue(MDL_WEIGHT, weight);
}

/** Get geometric transformation matrix from 2D-image header
  */
template<typename T>
Matrix2D< double > Image<T>::getTransformationMatrix(bool only_apply_shifts,
        long int n)
{
    // This has only been implemented for 2D images...
    (*this)().checkDimension(2);

    double phi,psi,theta,xoff,yoff;
    bool flip;
    MD[n].getValue(MDL_ANGLEROT, phi);
    phi = realWRAP(phi, 0., 360.);
    MD[n].getValue(MDL_ANGLETILT, theta);
    theta = realWRAP(theta, 0., 360.);
    MD[n].getValue(MDL_ANGLEPSI, psi);
    psi = realWRAP(psi, 0., 360.);
    MD[n].getValue(MDL_ORIGINX, xoff);
    MD[n].getValue(MDL_ORIGINY, yoff);

    Matrix2D< double > A(3, 3);
    A.initIdentity();

    if (only_apply_shifts)
    {
        Euler_angles2matrix(0., 0., 0., A);
        A(0, 2) = -xoff;
        A(1, 2) = -yoff;
    }
    else
    {
        if (theta == 0.)
        {
            // For untilted images: apply Euler matrix
            Euler_angles2matrix(phi, 0., psi, A);
        }
        else
        {
            // For tilted images: only apply Psi
            // Take another_set into account
            if (theta < 0.)
            {
                theta = -theta;
                psi = realWRAP(psi - 180., -180, 180);
            }
            Euler_angles2matrix(0., 0., psi, A);
        }
        A(0, 2) = -xoff;
        A(1, 2) = -yoff;
    }

    // Also for only_apply_shifts: mirror if necessary!
    MD[n].getValue(MDL_FLIP, flip);

    if (flip)
    {
        A(0, 0) = -A(0, 0);
        A(0, 1) = -A(0, 1);
    }

    return A;
}

/** Sum this object with other file and keep in this object
  */
template<typename T>
void Image<T>::sumWithFile(const FileName &fn)
{
    Image<T> aux;
    aux.read(fn);
    (*this)()+=aux();
}

//test
template<>
void Image< std::complex< double > >::castPage2T(char * page,
        std::complex<double> * ptrDest,
        DataType datatype,
        size_t pageSize)
{

    switch (datatype)
    {
    case ComplexShort:
        {
            std::complex<short> * ptr = (std::complex<short> *) page;
            for(int i=0; i<pageSize;i++)
                ptrDest[i]= std::complex<double> (real(ptr[i]),imag(ptr[i]));
            break;
        }
    case ComplexInt:
            {
                std::complex<int> * ptr = (std::complex<int> *) page;
                for(int i=0; i<pageSize;i++)
                ptrDest[i]= std::complex<double> (real(ptr[i]),imag(ptr[i]));
                break;
            }
        case ComplexFloat:
                {
                    std::complex<float> * ptr = (std::complex<float> *) page;
                    for(int i=0; i<pageSize;i++)
                        ptrDest[i]= std::complex<double> (real(ptr[i]),imag(ptr[i]));
                        break;
                    }
                case ComplexDouble:
                        {
                            memcpy(ptrDest, page, pageSize*sizeof(std::complex<double>));
                            break;
                        }
                    default:
                            {
                                std::cerr<<"Datatype= "<<datatype<<std::endl;
                                REPORT_ERROR(ERR_TYPE_INCORRECT," ERROR: cannot cast datatype to std::complex<double>");
                                break;
                            }
                        }

}

template<>
void Image< std::complex< double > >::castPage2Datatype(std::complex<double> * srcPtr,
        char * page,
        DataType datatype,
        size_t pageSize)
{
    switch (datatype)
    {
    case ComplexShort:
        {
            std::complex<short> * ptr = (std::complex<short> *) page;
            for(int i=0; i<pageSize;i++)
                ptr[i] = (std::complex<short>)srcPtr[i];
            break;
        }
    case ComplexInt:
            {
                std::complex<int> * ptr = (std::complex<int> *) page;
                for(int i=0; i<pageSize;i++)
                ptr[i] = (std::complex<int>)srcPtr[i];
                break;
            }
        case ComplexFloat:
                {
                    std::complex<float> * ptr = (std::complex<float> *) page;
                    for(int i=0; i<pageSize;i++)
                        ptr[i] = (std::complex<float>)srcPtr[i];
                        break;
                    }
                case ComplexDouble:
                        {
                            memcpy(page, srcPtr, pageSize*sizeof(std::complex<double>));
                            break;
                        }
                    default:
                            {
                                std::cerr<<"Datatype= "<<datatype<<std::endl;
                                REPORT_ERROR(ERR_TYPE_INCORRECT," ERROR: cannot cast datatype to complex<double>");
                                break;
                            }
                        }


}

// Get size of datatype
unsigned long  gettypesize(DataType type)
{
    unsigned long   size;

    switch ( type )
    {
    case UChar:
    case SChar:
        size = sizeof(char);
        break;
    case UShort:
    case Short:
        size = sizeof(short);
        break;
    case UInt:
    case Int:
        size = sizeof(int);
        break;
    case Float:
        size = sizeof(float);
        break;
    case Double:
        size = sizeof(double);
        break;
    case ComplexShort:
        size = sizeof(std::complex<short>);
        break;
    case ComplexInt:
        size = sizeof(std::complex<int>);
        break;
    case ComplexFloat:
        size = sizeof(std::complex<float>);
        break;
    case ComplexDouble:
        size = sizeof(std::complex<double>);
        break;
    case Bool:
        size = sizeof(bool);
        break;
    default:
        size = 0;
    }

    return(size);
}
