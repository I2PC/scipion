/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medicine
 * Univ. of California, Los Angeles.
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

#ifndef IMAGE_H
#define IMAGE_H

#include <typeinfo>
#include "funcs.h"
#include "memory.h"
#include "multidim_array.h"
#include "transformations.h"
#include "metadata.h"
#include <fcntl.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

//static std::vector<MDLabel> emptyVector;
//static MetaData emptyMetaData;

typedef enum
{
    NoTransform = 0,        // No transform
    Standard = 1,           // Standard transform: origin = (0,0,0)
    Centered = 2,           // Centered transform: origin = (nx/2,ny/2,nz/2)
    Hermitian = 3,          // Hermitian half: origin = (0,0,0)
    CentHerm = 4            // Centered hermitian: origin = (0,ny/2,nz/2)
} TransformType;

typedef enum
{
    Unknown_Type = 0,       // Undefined data type
    UChar = 1,              // Unsigned character or byte type
    SChar = 2,              // Signed character (for CCP4)
    UShort = 3,             // Unsigned integer (2-byte)
    Short = 4,              // Signed integer (2-byte)
    UInt = 5,       // Unsigned integer (4-byte)
    Int = 6,                // Signed integer (4-byte)
    Long = 7,               // Signed integer (4 or 8 byte, depending on system)
    Float = 8,              // Floating point (4-byte)
    Double = 9,             // Double precision floating point (8-byte)
    ComplexShort = 10,       // Complex two-byte integer (4-byte)
    ComplexInt = 11,        // Complex integer (8-byte)
    ComplexFloat = 12,      // Complex floating point (8-byte)
    ComplexDouble = 13,      // Complex floating point (16-byte)
    Bool = 14    // Boolean (1-byte?)
} DataType;

// Ask memory size of datatype
unsigned long   gettypesize(DataType type);

/** Image Matrix access.
 * @ingroup ImagesSpeedUp
 *
 * This macro does the same as the normal 3D matrix access but in a faster way
 * as no function call is generated.
 *
 * @code
 * VOLMATRIX(V).resize(128, 128, 128);
 *
 * VOLMATRIX(V2) = VOLMATRIX(V1) + VOLMATRIX(V2);
 * @endcode
 */
#define VOLMATRIX(V) ((V).data)

#define SWAPTRIG     65535   // Threshold file z size above which bytes are swapped
// For fast access to pixel values (and for backwards compatibility of the code)
#define IMGPIXEL(I, i, j) A2D_ELEM(((I).data), (i), (j))

/** Voxel access.
 * @ingroup ImagesSpeedUp
 *
 * This macro does the same as the normal voxel access (remember, logical
 * access) but in a faster way as no function call is generated.
 *
 * @code
 * std::cout << "Grey level of voxel (2,-3,-3) of the Volume = " <<
 *     VOLVOXEL(V, 2, -3, -3) << std::endl;
 *
 * VOLVOXEL(I, 2, -3, -3) = VOLVOXEL(I, 2, -3, -2);
 * @endcode
 */
#define VOLVOXEL(V, k, i, j) A3D_ELEM(((V).data), (k), (i), (j))

/** Physical voxel access.
 * @ingroup ImagesSpeedUp
 *
 * The physical voxel access gives you access to a voxel by its physical
 * position and not by its logical one. This access shouldn't be used as a
 * custom, use instead the logical access, but there might be cases in which
 * this access might be interesting. Physical positions start at index 0 in C.
 *
 * @code
 * std::cout << "This is the first voxel stored in the Volume " <<
 *     DIRECT_VOLVOXEL(V, 0, 0, 0) << std::endl;
 * @endcode
 */
#define DIRECT_VOLVOXEL(I, k, i, j) DIRECT_A3D_ELEM(((I).data), (k), (i), (j))

#define DIRECT_IMGPIXEL(I, i, j) DIRECT_A2D_ELEM(((I).data), (i), (j))
#define IMGMATRIX(I) ((I).data)

//dummy vectors for default function inizialization
//static MetaDataContainer emptyMetaDataContainer;
//static std::vector<MDLabel> emptyVector;

/** Template class for images
 * @ingroup Images
 *
 * The image class is the general image handling class.
 *
 */
template<typename T>
class Image
{
public:

    MultidimArray<T>    data;       // The image data array
    // FIXME: why cant this one be private as well?
    MetaData MD;//data for each subimage
    MetaData MDMainHeader;//data for the file
private:
    FileName            filename;    // File name
    int                 dataflag;    // Flag to force reading of the data
    unsigned long       i;           // Current image number (may be > NSIZE)
    unsigned long       offset;      // Data offset
    int                 swap;        // Perform byte swapping upon reading
    TransformType       transform;   // Transform type
    int                 replaceNsize;// Stack size in the replace case
    int                 _exists;     // does target file exists?
    // equal 0 is not exists or not a stack

    bool    mmapOn;  // Mapping when loading from file
    int     mFd;     //Handle the file in reading method and mmap
    size_t       mappedSize;
    /*
    double    min, max; // Limits
    double    avg, std; // Average and standard deviation
    double    smin, smax; // Limits for display
    double    scale;  // Scale of last density conversion operation
    double    shift;  // Shift of last density conversion operation before scaling
    double    resolution; // Resolution limit of data - used for low-pass filtering
    double    ux, uy, uz; // Voxel units (angstrom/pixel edge)
    double    ua, ub, uc; // Unit cell dimensions (angstrom)
    double    alf, bet, gam; // Unit cell angles (radian)
    unsigned int  spacegroup; // Space group

    */

public:
    /** Empty constructor
     *
     * An empty image is created.
     *
     * @code
     * Image<double> I;
     * @endcode
     */
    Image()
    {
        mmapOn = false;
        clear();
    }

    /** Constructor with size
     *
     * A blank image (0.0 filled) is created with the given size. Pay attention
     * to the dimension order: Y and then X.
     *
     * @code
     * Image I(64,64);
     * @endcode
     */
    Image(int Xdim, int Ydim, int Zdim=1, int Ndim=1)
    {
        mmapOn = false;
        clear();
        data.resize(Ndim, Zdim, Ydim, Xdim);
        for (int n = 1; n < Ndim; n++)
            MD.addObject();
    }

    /** Clear.
     * Initialize everything to 0
     */
    void clear()
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
     *
     */
    void clearHeader()
    {
        MDMainHeader.clear();
        MD.clear();
        MD.addObject(); // Each image has at lest one MD object
        MDMainHeader.addObject(); // Each image has at lest one MD object
    }

    /** Check whether image is complex based on T
       */
    bool isComplexT() const
    {
        return ( typeid(T) == typeid(std::complex<double>) ||
                 typeid(T) == typeid(std::complex<float>) );
    }

    /** Check whether image is complex based on transform
         */
    bool isComplex() const
    {
        return !(transform==NoTransform);
    }

    /** Destructor.
     *
     */
    ~Image()
    {
        clear();
    }


    /** Specific read functions for different file formats
    */
#include "rwSPIDER.h"
#include "rwMRC.h"
#include "rwIMAGIC.h"
#include "rwTIA.h"
#include "rwDM3.h"
#include "rwRAW.h"
    //#include "rwTIFF.h"

    /** Is this file an image
     *
     *  Check whether a real-space image can be read
     *
     */
    bool isImage(const FileName &name)
    {
        return !read(name, false);
    }

    /** Is this file a real-valued image
     *
     *  Check whether a real-space image can be read
     *
     */
    bool isRealImage(const FileName &name)
    {
        return (isImage(name) && !isComplex());
    }

    /** Is this file a complex image
     *
     *  Check whether a fourier-space (complex) image can be read
     *
     */
    bool isComplexImage(const FileName &name)
    {
        return (isImage(name) && isComplex());
    }

    /** Rename the image
      */
    void rename (const FileName &name)
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
    int read(const FileName &name, bool readdata=true, int select_img=-1,
             bool apply_geo = false, bool only_apply_shifts = false,
             const MetaData *docFilePtr = NULL,
             std::vector<MDLabel> *activeLabelsPtr = NULL, bool mapData = false)
    {
        //const MetaData &docFile = *docFilePtr;
        //std::vector<MDLabel> &activeLabels = *activeLabelsPtr;

        int err = 0;
        // Check whether to read the data or only the header
        if ( readdata )
            dataflag = 1;
        else
            dataflag = -1;

        // Check whether to map the data or not
        mmapOn = mapData;

        FileName ext_name = name.get_file_format();
        size_t found;
        filename = name;
        found=filename.find_first_of("@");
        if (found!=std::string::npos)
        {
            select_img =  atoi(filename.substr(0, found).c_str());
            filename = filename.substr(found+1) ;
        }

        double imParam = NULL;
        found=filename.find_first_of("%");
        if (found!=std::string::npos)
        {
            imParam =  atof(filename.substr(found+1).c_str());
            filename = filename.substr(0, found) ;
        }

        filename = filename.remove_file_format();

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
        else if (ext_name.contains("raw"))//RAW
            err = readRAW(select_img,false);
        //        else if (ext_name.contains("tif") || ext_name.contains("tiff"))//TIFF
        //            err = readTIFF(select_img,false);
        else
        {
            err = readSPIDER(select_img,true);
        }

        //This implementation does not handle stacks,
        //read in a block

        if (docFilePtr != NULL)
        {
            if (activeLabelsPtr == NULL)
                activeLabelsPtr = docFilePtr->geActiveLabelsAddress();

            std::vector<MDLabel>::iterator strIt;
            double dd;
            std::string ss;
            int ii;
            bool bb;
            std::vector<double> vv;
            //FIXME: This could be done better
            for (strIt = activeLabelsPtr->begin(); strIt != activeLabelsPtr->end(); strIt++)
            {
                switch (MDL::labelType(*strIt))
                {
                case LABEL_DOUBLE:
                    docFilePtr->getValue(*strIt,dd);
                    MD.setValue(*strIt,dd);
                    break;
                case LABEL_STRING:
                    docFilePtr->getValue(*strIt,ss);
                    MD.setValue(*strIt,ss);
                    break;
                case LABEL_INT:
                    docFilePtr->getValue(*strIt,ii);
                    MD.setValue(*strIt,ii);
                    break;
                case LABEL_BOOL:
                    docFilePtr->getValue(*strIt,bb);
                    MD.setValue(*strIt,bb);
                    break;
                case LABEL_VECTOR:
                    docFilePtr->getValue(*strIt,vv);
                    MD.setValue(*strIt,vv);
                    break;
                default:
                    REPORT_ERROR(-1, "Image.read: Unknown label type");
                }
            }//close for activeLabels
        }

        //apply geo has not been defined for volumes
        if(this->data.getDim()>2)
            apply_geo=false;

        //err = readMRC(*this, imgno);
        /*
        if ( filename.contains("#") || ext.empty() )
        err = readRAW(p, select_img);
        else if ( ext.contains("raw") )
        err = readRAW(p, select_img);
        else if ( ext.contains("asc") || ext.contains( "txt") )
        err = readASCII(p);
        else if ( ext.contains("pic") )
        err = readBIORAD(p);
        else if ( ext.contains("brx" ) || ext.contains("brix" ) )
        err = readBRIX(p);
        else if ( ext.contains("dat" ) )
        err = readBrookhavenSTEM(p);
        else if ( ext.contains("ccp") || ext.contains("map") )
        err = readCCP4(p);
        else if ( ext.contains("di") )
        err = readDI(p, select_img);
        else if ( ext.contains("dm") )
        err = readDM(p);
        else if ( ext.contains("omap" ) || ext.contains("dsn6") || ext.contains("dn6") )
        err = readDSN6(p);
        else if ( ext.contains("em") )
        err = readEM(p);
        else if ( ext.contains("pot") )
        err = readGOODFORD(p);
        else if ( ext.contains("grd") )
        err = readGRD(p);
        else if ( ext.contains("hkl") )
        err = readHKL(p);
        else if ( ext.contains("ip") )
        err = readIP(p);
        else if ( ext.contains("jpg") || ext.contains("jpeg") )
        err = readJPEG(p);
        else if ( ext.contains("mif") )
        err = readMIFF(p, select_img);
        else if ( ext.contains("mff") )
        err = readMFF(p);
        else if ( ext.contains("pif") || ext.contains("sf") )
        err = readPIF(p, select_img);
        else if ( ext.contains("bp") || ext.contains("bq") )
        err = readPIC(p);
        else if ( ext.contains("png") )
        err = readPNG(p);
        else if ( ext.contains("spe") )
        err = readSPE(p, select_img);
        else if ( ext.contains("spm") || ext.contains("sup") || ext == "f" )
        err = readSUPRIM(p);
        else if ( ext.contains("tif") )
        err = readTIFF(p, select_img);
        else if ( ext.contains("xpl") || ext.contains("cns") || ext.contains("rfl") )
        err = readXPLOR(p);
        else {
        fprintf(stderr, "Error: File format with extension \"%s\" not supported!\n", ext.c_str());
        err = -1;
            }
        */
        /*
        if ( err < 0 ) {
            REPORT_ERROR(10,"Error reading file");
            }
        */

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
    //so far I assume that we are going to write everything
    //replace and append no valid for stacks
    void write(FileName name="",
               int select_img=-1,
               bool isStack=false,
               int mode=WRITE_OVERWRITE)
    {

        int err = 0;

        if (name == "")
            name = filename;
        /*
                if(
              (NSIZE(data)>1 && mode==WRITE_APPEND ) ||
              (NSIZE(data)>1 && mode==WRITE_REPLACE )
                )
                {
                 REPORT_ERROR(1,"append and replace are not available options for stacks");
                }
        */
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
        filNamePlusExt = filename;
        found=filename.find_first_of(":");
        if ( found!=std::string::npos)
        {
            filename   = filename.substr(0, found);
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
            REPORT_ERROR(1,"write Image ERROR: image is empty!");

        // CHECK FOR INCONSISTENCIES BETWEEN data.xdim and x, etc???
        int Xdim, Ydim, Zdim, Ndim;
        this->getDimensions(Xdim,Ydim, Zdim, Ndim);

        _exists = exists(filename);
        Image<T> auxI;
        replaceNsize=0;//reset replaceNsize in case image is reused
        if(select_img==-1 && mode==WRITE_REPLACE)
            REPORT_ERROR(1,"writeSPIDER: Please specify object to be replaced");
        else if(!_exists && mode==WRITE_REPLACE)
        {
            std:: stringstream replace_number;
            replace_number << select_img;
            REPORT_ERROR(1,(std::string)"Cannot replace object number: "
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
                REPORT_ERROR(1,"write: target and source objects have different size");
            if(mode==WRITE_REPLACE && select_img>_Ndim)
                REPORT_ERROR(1,"write: can not replace image stack is not large enough");
            if(auxI.replaceNsize <1 &&
               (mode==WRITE_REPLACE || mode==WRITE_APPEND))
                REPORT_ERROR(1,"write: output file is not an stack");
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

        if(ext_name.contains("spi") || ext_name.contains("xmp"))
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
        else if (ext_name.contains("raw"))
            writeRAW(select_img,false,mode);
        //        else if (ext_name.contains("tif") || ext_name.contains("tiff"))
        //            writeTIFF(select_img,mode);
        else
            err = writeSPIDER(select_img,isStack,mode);

        if ( err < 0 )
        {
            std::cerr<<" Filename = "<<filename<<" Extension= "<<ext_name<<std::endl;
            REPORT_ERROR(10,"Error writing file");
        }
        //unlock file
    }

    /** Cast a page of data from type dataType to type Tdest
     * @ingroup LittleBigEndian
     *    input pointer  char *
     */
    void castPage2T(char * page, T * ptrDest, DataType datatype, size_t pageSize )
    {
        switch (datatype)
        {
        case Unknown_Type:
            REPORT_ERROR(12,"ERROR: datatype is Unknown_Type");
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
                    REPORT_ERROR(16," ERROR: cannot cast datatype to T");
                    break;
                }
            }

    }

    /** Cast page from T to datatype
     * @ingroup XXX
     *    input pointer  char *
     */
    void castPage2Datatype(T * srcPtr, char * page, DataType datatype, size_t pageSize )

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
        default:
                {
                    std::cerr<<"outputDatatype= "<<datatype<<std::endl;
                    REPORT_ERROR(16," ERROR: cannot cast T to outputDatatype");
                    break;
                }
            }
    }

    /** Check file Datatype is same as T type to use mmap.
     * @ingroup XXX.
     */
    bool checkMmapT(DataType datatype)
    {

        switch (datatype)
        {
        case Unknown_Type:
            REPORT_ERROR(12,"ERROR: datatype is Unknown_Type");
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
                REPORT_ERROR(16," ERROR: cannot cast datatype to T");
                break;
            }
        }
        //               int * iTemp = (int*) map;
        //                ptrDest = reinterpret_cast<T*> (iTemp);
    }

    /** Write an entire page as datatype
     * @ingroup XXX
     *
     * A page of datasize_n elements T is cast to datatype and written to fimg
     * The memory for the casted page is allocated and freed internally.
     *
     */
    void writePageAsDatatype(FILE * fimg, DataType datatype, size_t datasize_n )
    {
        size_t datasize = datasize_n * gettypesize(datatype);
        char * fdata = (char *) askMemory(datasize);
        castPage2Datatype(MULTIDIM_ARRAY(data), fdata, datatype, datasize_n);
        fwrite( fdata, datasize, 1, fimg );
        freeMemory(fdata, datasize);
    }

    /** Swap an entire page
      * @ingroup XXX
      *    input pointer  char *
      */
    void swapPage(char * page, size_t pageNrElements, DataType datatype)
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

    void readData(FILE* fimg, int select_img, DataType datatype, unsigned long pad)
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
        size_t pagemax_n = ROUND(pagemax/datatypesize);
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
                REPORT_ERROR(22,"Image Class::ReadData: mmap with multiple \
                             images file not compatible. Try selecting a unique image.");
            }

            fclose(fimg);

            if ( ( mFd = open(filename.c_str(), O_RDWR, S_IREAD | S_IWRITE) ) == -1 )
                REPORT_ERROR(20,"Image Class::ReadData: Error opening the image file.");

            char * map;
            mappedSize = pagesize+offset;

            if ( (map = (char*) mmap(0,mappedSize, PROT_READ | PROT_WRITE, MAP_SHARED, mFd, 0)) == (void*) -1 )
                REPORT_ERROR(21,"Image Class::ReadData: mmap of image file failed.");
            data.data = reinterpret_cast<T*> (map+offset);
        }
        else
        {
            // Reset select to get the correct offset
            if ( select_img < 0 )
                select_img = 0;

            char* page = NULL;
            char* padpage = NULL;

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

            //if ( pad > 0)
            //    padpage = (char *) askMemory(pad*sizeof(char));
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

     ********* FIXME!!! withx,y,z also being part of image class this resizing would be DANGEROUS!!!

     *
     * @code
     * I().resize(128, 128);
     * I2() = I1() + I2();
     * @endcode
     */
    MultidimArray<T>& operator()()
    {
        return data;
    }
    const MultidimArray<T>& operator()() const
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
    T& operator()(int i, int j) const
    {
        return A2D_ELEM(data, i, j);
    }
    /**
    *  set pixel (direct acces) needed by swig
    */
    void setPixel(int i, int j, T v)
    {
        IMGPIXEL(*this,i,j)=v;
    }

    /**
    *  get pixel (direct acces) needed by swig
    */
    T getPixel(int i, int j)
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
    T& operator()(int k, int i, int j) const
    {
        return A3D_ELEM(data, k, i, j);
    }

    /** Get file name
     *
     * @code
     * std::cout << "Image name = " << I.name() << std::endl;
     * @endcode
     */
    const FileName & name() const
    {
        return filename;
    }

    /** Get Image dimensions
     *
     */
    void getDimensions(int &Xdim, int &Ydim, int &Zdim, int &Ndim) const
    {
        Xdim = XSIZE(data);
        Ydim = YSIZE(data);
        Zdim = ZSIZE(data);
        Ndim = NSIZE(data);
    }

    long unsigned int getSize() const
    {
        return NZYXSIZE(data);
    }

    /** Get Image dimensions
     *
     */
    void getHeaderInfo(unsigned long &_offset, int &_swap) const
    {
        _offset = offset;
        _swap = swap;
    }

    /* Is there label in the individual header */
    bool individualContainsLabel(MDLabel label) const
    {
        return MD.containsLabel(label);
    }

    /* Is there label in the main header */
    bool mainContainsLabel(MDLabel label) const
    {
        return MDMainHeader.containsLabel(label);
    }

    /** Get Rot angle
    *
    * @code
    * std::cout << "First Euler angle " << I.rot() << std::endl;
    * @endcode
    */
    double rot(const long int n = -1) const
    {
        double dummy=0;
        MD.getValue(MDL_ANGLEROT,dummy,n);
        return (dummy);
    }

    /** Get Tilt angle
     *
     * @code
     * std::cout << "Second Euler angle " << I.tilt() << std::endl;
     * @endcode
     */
    double tilt(const long int n = -1) const
    {
        double dummy=0;
        MD.getValue(MDL_ANGLETILT,dummy,n);
        return (dummy);
    }

    /** Get Psi angle
     *
     * @code
     * std::cout << "Third Euler angle " << I.psi() << std::endl;
     * @endcode
     */
    double psi(const long int n = -1) const
    {
        double dummy=0;
        MD.getValue(MDL_ANGLEPSI,dummy,n);
        return (dummy);
    }

    /** Get Xoff
     *
     * @code
     * std::cout << "Origin offset in X " << I.Xoff() << std::endl;
     * @endcode
     */
    double Xoff(const long int n = -1) const
    {
        double dummy=0;
        MD.getValue(MDL_ORIGINX,dummy,n);
        return (dummy);
    }

    /** Get Yoff
     *
     * @code
     * std::cout << "Origin offset in Y " << I.Yoff() << std::endl;
     * @endcode
     */
    double Yoff(const long int n = -1) const
    {
        double dummy=0;
        MD.getValue(MDL_ORIGINY,dummy,n);
        return (dummy);
    }

    /** Get Zoff
     *
     * @code
     * std::cout << "Origin offset in Z " << I.Zoff() << std::endl;
     * @endcode
     */
    double Zoff(const long int n = -1) const
    {
        double dummy=0;
        MD.getValue(MDL_ORIGINZ,dummy,n);
        return (dummy);
    }

    /** Get Weight
    *
    * @code
    * std::cout << "weight= " << I.weight() << std::endl;
    * @endcode
    */
    double weight(const long int n = -1) const
    {
        double dummy=1;
        MD.getValue(MDL_WEIGHT,dummy,n);
        return (dummy);
    }

    /** Get Flip
    *
    * @code
    * std::cout << "flip= " << flip() << std::endl;
    * @endcode
    */
    bool flip(const long int n = -1) const
    {
        bool dummy=false;
        MD.getValue(MDL_FLIP,dummy,n);
        return (dummy);
    }

    /** Sampling RateX
    *
    * @code
    * std::cout << "sampling= " << samplingRateX() << std::endl;
    * @endcode
    */
    bool samplingRateX(const long int n = -1) const
    {
        double dummy=1.;
        MDMainHeader.getValue(MDL_SAMPLINGRATEX,dummy,n);
        return (dummy);
    }

    /** Set file name
     *
     */
    void setName(const FileName &_filename)
    {
        filename = _filename;
    }

    /** Set Euler angles in image header
     *
     */
    void setEulerAngles(double rot, double tilt, double psi,
                        long int n = -1)
    {
        MD.setValue(MDL_ANGLEROT,rot,n);
        MD.setValue(MDL_ANGLETILT,tilt,n);
        MD.setValue(MDL_ANGLEPSI,psi,n);
    }

    /** Get Euler angles from image header
     *
     */
    void getEulerAngles(double &rot, double &tilt, double &psi,
                        long int n = -1)
    {
        MD.getValue(MDL_ANGLEROT,rot,n);
        MD.getValue(MDL_ANGLETILT,tilt,n);
        MD.getValue(MDL_ANGLEPSI,psi,n);
    }

    /** Set Rotation angle to image */
    void setRot(double rot, long int n = -1)
    {
        MD.setValue(MDL_ANGLEROT,rot,n);
    }

    /** Set Tilt angle to image */
    void setTilt(double tilt, long int n = -1)
    {
        MD.setValue(MDL_ANGLETILT,tilt,n);
    }

    /** Set Rotation angle to image */
    void setPsi(double psi, long int n = -1)
    {
        MD.setValue(MDL_ANGLEPSI,psi,n);
    }

    /** Set origin offsets in image header
     *
     */
    void setShifts(double xoff, double yoff, double zoff = 0.,
                   long int n = -1)
    {
        MD.setValue(MDL_ORIGINX,xoff,n);
        MD.setValue(MDL_ORIGINY,yoff,n);
        MD.setValue(MDL_ORIGINZ,zoff,n);
    }
    /** Get origin offsets from image header
         *
         */
    void getShifts(double &xoff, double &yoff, double &zoff = 0.,
                   long int n = -1)
    {
        MD.getValue(MDL_ORIGINX,xoff,n);
        MD.getValue(MDL_ORIGINY,yoff,n);
        MD.getValue(MDL_ORIGINZ,zoff,n);
    }


    /** Set X offset in image header
     */
    void setXoff(double xoff, long int n = -1)
    {
        MD.setValue(MDL_ORIGINX,xoff,n);
    }

    /** Set Y offset in image header
     */
    void setYoff(double yoff, long int n = -1)
    {
        MD.setValue(MDL_ORIGINY,yoff,n);
    }

    /** Set Z offset in image header
     */
    void setZoff(double zoff, long int n = -1)
    {
        MD.setValue(MDL_ORIGINZ,zoff,n);
    }

    /** Set flip in image header
     *
     */
    void setFlip(bool flip, long int n = -1)
    {
        MD.setValue(MDL_FLIP,flip,n);
    }

    /** Set Weight in image header
    */
    void setWeight(double weight, long int n = -1)
    {
        MD.setValue(MDL_WEIGHT,weight, n);
    }

    /** Get geometric transformation matrix from 2D-image header
      */
    Matrix2D< double > getTransformationMatrix(bool only_apply_shifts = false,
            long int n = -1)
    {
        // This has only been implemented for 2D images...
        (*this)().checkDimension(2);

        double phi,psi,theta,xoff,yoff;
        bool flip;
        MD.getValue(MDL_ANGLEROT,phi,n);
        phi = realWRAP(phi, 0., 360.);
        MD.getValue(MDL_ANGLETILT,theta,n);
        theta = realWRAP(theta, 0., 360.);
        MD.getValue(MDL_ANGLEPSI,psi,n);
        psi = realWRAP(psi, 0., 360.);
        MD.getValue(MDL_ORIGINX,xoff,n);
        MD.getValue(MDL_ORIGINY,yoff,n);

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
        MD.getValue(MDL_FLIP,flip,n);
        if (flip)
        {
            A(0, 0) = -A(0, 0);
            A(0, 1) = -A(0, 1);
        }

        return A;
    }

    /** Show image properties
      */
    friend std::ostream& operator<<(std::ostream& o, const Image<T>& I)
    {
        o << "Image type   : ";
        if (I.isComplex())
            o << "Fourier-space image" << std::endl;
        else
            o << "Real-space image" << std::endl;

        o << "Reversed     : ";
        if (I.swap)
            o << "TRUE"  << std::endl;
        else
            o << "FALSE" << std::endl;

        o << "dimensions   : " << NSIZE(I()) << " x " << ZSIZE(I()) << " x " << YSIZE(I()) << " x " << XSIZE(I());
        o << "  (noObjects x slices x rows x columns)" << std::endl;
        if (I.individualContainsLabel(MDL_ANGLEROT))
        {
            o << "Euler angles : " << std::endl;
            o << "  Phi   (rotation around Z axis) = " << I.rot() << std::endl;
            o << "  theta (tilt, second rotation around new Y axis) = " << I.tilt() << std::endl;
            o << "  Psi   (third rotation around new Z axis) = " << I.psi() << std::endl;
        }
        if (I.individualContainsLabel(MDL_ORIGINX))
        {
            o << "Origin Offsets : " << std::endl;
            o << "  Xoff  (origin offset in X-direction) = " << I.Xoff() << std::endl;
            o << "  Yoff  (origin offset in Y-direction) = " << I.Yoff() << std::endl;
            o << "  Zoff  (origin offset in Z-direction) = " << I.Zoff() << std::endl;
        }
        o << "Header size  : " << I.offset << std::endl;
        if (I.individualContainsLabel(MDL_WEIGHT))
            o << "Weight  : " << I.weight() << std::endl;
        if (I.individualContainsLabel(MDL_FLIP))
            o << "Flip    : " << I.flip() << std::endl;
        return o;
    }

    /** Sum this object with other file and keep in this object
      */
    void sumWithFile(const FileName &fn)
    {
    	Image<T> aux;
    	aux.read(fn);
    	(*this)()+=aux();
    }
};


// Special case for complex numbers
template<>
void Image< std::complex< double > >::castPage2T(char * page,
        std::complex<double> * ptrDest,
        DataType datatype,
        size_t pageSize);
template<>
void Image< std::complex< double > >::castPage2Datatype(std::complex< double > * srcPtr,
        char * page,
        DataType datatype,
        size_t pageSize);


#endif
