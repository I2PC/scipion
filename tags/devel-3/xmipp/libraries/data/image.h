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
#include "gcc_version.h"
#include "multidim_array.h"
#include "transformations.h"

typedef enum TransformType {
        NoTransform = 0,        // No transform
        Standard = 1,           // Standard transform: origin = (0,0,0)
        Centered = 2,           // Centered transform: origin = (nx/2,ny/2,nz/2)
        Hermitian = 3,          // Hermitian half: origin = (0,0,0)
        CentHerm = 4            // Centered hermitian: origin = (0,ny/2,nz/2)
} ;

typedef enum DataType {
        Unknown_Type = 0,       // Undefined data type
        UChar = 1,              // Unsigned character or byte type
        SChar = 2,              // Signed character (for CCP4)
        UShort = 3,             // Unsigned integer (2-byte)
        Short = 4,              // Signed integer (2-byte)
        Int = 5,                // Signed integer (4-byte)
        Long = 6,               // Signed integer (4 or 8 byte, depending on system)
        Float = 7,              // Floating point (4-byte)
        Double = 8,             // Double precision floating point (8-byte)
        ComplexShort = 9,       // Complex two-byte integer (4-byte)
        ComplexInt = 10,        // Complex integer (8-byte)
        ComplexFloat = 11,      // Complex floating point (8-byte)
        ComplexDouble = 12      // Complex floating point (16-byte)
} ;

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
#define DIRECT_VOLVOXEL(V, k, i, j) DIRECT_A3D_ELEM(((V).img), (k), (i), (j))

#define DIRECT_IMGPIXEL(I, i, j) DIRECT_A2D_ELEM(((I).data), (i), (j))
#define IMGMATRIX(I) ((I).data)


/** Template class for subimage header information
 * @ingroup Images
 *
 * The header information of subimages is handled in a way similar to Bsoft
 * In this way, a multi-image image can hold assignments for all individual images
 * 
 */
class SubImage
{
public:
    double shiftX, shiftY, shiftZ;        // Origin 
    double angleRot, angleTilt, anglePsi; // Euler angles
    double weight;	                  // Weight
    bool   flip;                          // Apply mirror operation

public:

    /** Empty constructor
     *
     * An empty SubImage is created.
     *
     * @code
     * SubImage I;
     * @endcode
     */
    SubImage()
    {
        emptyInit();
    }

    /** emptyInit.
     * Initialize offsets and angles to 0 and weight to 1.
     */
    void emptyInit()
    {
        shiftX = shiftY = shiftZ = 0.;
        angleRot = angleTilt = anglePsi = 0.;
        weight = 1.;
        flip = false;
    }    

    /** Destructor.
     */
     ~SubImage()
     {
     }
    
};

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

    SubImage*           image;      // Sub-images
    MultidimArray<T>    data;       // The image data array
    // FIXME: why cant this one be private as well?
private:
    FileName            filename;   // File name
    int        			dataflag;	// Flag to force reading of the data
    unsigned long	i;			// Current image number (may be > NSIZE)
    unsigned long	offset; 	// Data offset
    int                swap;       // Perform byte swapping upon reading
    TransformType		transform;  // Transform type
    double				min, max;	// Limits
    double				avg, std;	// Average and standard deviation
    double				smin, smax; // Limits for display
    double				scale;		// Scale of last density conversion operation
    double				shift;		// Shift of last density conversion operation before scaling
    double				resolution; // Resolution limit of data - used for low-pass filtering
    double				ux, uy, uz;	// Voxel units (angstrom/pixel edge)
    double				ua, ub, uc; // Unit cell dimensions (angstrom)
    double				alf, bet, gam;	// Unit cell angles (radian)
    unsigned int		spacegroup;	// Space group


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
        image = NULL;
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
        image = NULL;
        clear();
        data.resize(Ndim, Zdim, Ydim, Xdim);
        image = new SubImage [Ndim];
        if (image == NULL)
            REPORT_ERROR(1001, "Allocate: No space left for subimage structure");
    }

    /** Clear.
     * Initialize everything to 0
     */
    void clear()
    {
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
        min = max = avg = 0.;
        std = -1.;
        scale = shift = 1.;
        resolution = 0.;
        ux = uy = uz = 1.;
        ua = ub = uc = 1.;
        alf, bet, gam = DEG2RAD(90.);
        spacegroup = 1;
        // delete subimage structure if present
        if (image != NULL)
            delete[] image;
        // and put a single empty copy
        image = new SubImage [1];
        if (image == NULL)
            REPORT_ERROR(1001, "Allocate: No space left for subimage structure");
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
     */
     ~Image()
     {
         clear();
     }


     /** Specific read functions for different file formats
     */
     #include "rwSPIDER.h"
     #include "rwMRC.h"

     /** Is this file an image
      *
      *  Check whether a real-space image can be read
      *
      */
     bool isImage(const FileName name)
     {
    	 return !read(name, false);
     }

     /** Is this file a real-valued image
      *
      *  Check whether a real-space image can be read
      *
      */
     bool isRealImage(const FileName name)
     {
         return (isImage(name) && !isComplex());
     }

     /** Is this file a complex image
      *
      *  Check whether a fourier-space (complex) image can be read
      *
      */
     bool isComplexImage(const FileName name)
     {
         return (isImage(name) && isComplex());
     }

     /** Rename the image
       */
     void rename (const FileName name)
     {
    	 filename = name;
     }

    /** General read function
     */
     int read(const FileName name, bool readdata=true, int select_img=-1,
              bool apply_geo = false, bool only_apply_shifts = false)
     {
         int err = 0;

         // Check whether to read the data or only the header
         if ( readdata ) dataflag = 1;
         else dataflag = -1;

         FileName ext_name = name.get_file_format();
         if ( name.contains("#") )
        	 filename = name;
         else
        	 filename = name.before_first_of(":");
//#define DEBUG
#ifdef DEBUG
         std::cerr << "name="<<name <<std::endl;
         std::cerr << "ext= "<<ext_name <<std::endl;
         std::cerr<<" now reading: "<< filename<<" dataflag= "<<dataflag<<std::endl;
#endif

         if (ext_name.contains("mrc"))
             err = readMRC();
         else
             err = readSPIDER(select_img);

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
	else if ( ext.contains("img") || ext.contains("hed") )
		err = readIMAGIC(p, select_img);
	else if ( ext.contains("ip") )
		err = readIP(p);
	else if ( ext.contains("jpg") || ext.contains("jpeg") )
		err = readJPEG(p);
	else if ( ext.contains("mif") )
		err = readMIFF(p, select_img);
	else if ( ext.contains("mff") )
		err = readMFF(p);
	else if ( ext.contains("mrc") )
		err = readMRC(p);
	else if ( ext.contains("pif") || ext.contains("sf") )
		err = readPIF(p, select_img);
	else if ( ext.contains("bp") || ext.contains("bq") )
		err = readPIC(p);
	else if ( ext.contains("png") )
		err = readPNG(p);
	else if ( ext.contains("spe") )
		err = readSPE(p, select_img);
	else if ( ext.contains("spi") )
		err = readSPIDER(p, select_img);
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
                 applyGeometry(3, (*this)(), tmp, A, IS_INV, WRAP);
             }
         }

         // Negative errors are bad.
         return err;

     }

    /** General write function
     */
     void write(FileName name="")
     {
         int err = 0;

         if (name == "")
        	 name = filename;

         FileName ext_name = name.get_file_format();
         filename = name.before_first_of(":");
         filename = filename.before_first_of("#");

#ifdef DEBUG
         std::cerr<<"extension for write= "<<ext_name<<std::endl;
         std::cerr<<"filename= "<<filename<<std::endl;
#endif

         // Check that image is not empty
         if (getSize() < 1)
        	 REPORT_ERROR(1,"write Image ERROR: image is empty!");

          // PERHAPS HERE CHECK FOR INCONSISTENCIES BETWEEN data.xdim and x, etc???
         if (ext_name.contains("mrc"))
         {
        	 err = writeMRC();
         }
         else
             err = writeSPIDER();

         if ( err < 0 ) {
             std::cerr<<" Filename = "<<filename<<" Extension= "<<ext_name<<std::endl;
             REPORT_ERROR(10,"Error writing file");
         }

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
             {
                 memcpy(ptrDest, page, pageSize*sizeof(T));
             }
             else
             {
                 unsigned char * ptr = (unsigned char *) page;
                 for(int i=0; i<pageSize;i++) ptrDest[i]=(T) ptr[i];
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
                 for(int i=0; i<pageSize;i++) ptrDest[i]=(T) ptr[i];
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
                 for(int i=0; i<pageSize;i++) ptrDest[i]=(T) ptr[i];
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
                 for(int i=0; i<pageSize;i++) ptrDest[i]=(T) ptr[i];
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
                 for(int i=0; i<pageSize;i++) ptrDest[i]=(T) ptr[i];
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
                 for(int i=0; i<pageSize;i++) ptrDest[i]=(T) ptr[i];
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
                 for(int i=0; i<pageSize;i++) ptrDest[i]=(T) ptr[i];
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
                 for(int i=0; i<pageSize;i++) ptrDest[i]=(T) ptr[i];
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
                 for(int i=0; i<pageSize;i++)
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
                 for(int i=0; i<pageSize;i++) ptr[i] = (double)srcPtr[i];
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
                  << swap<<" datatypesize= "<<datatypesize<<std::endl;;
#endif

         // Swap bytes if required
         if ( swap == 1 ) {
             if ( datatype >= ComplexShort ) 
                 datatypesize /= 2;
             for ( unsigned long i=0; i<pageNrElements; i+=datatypesize ) 
                 swapbytes(page+i, datatypesize);
         } else if ( swap > 1 ) {
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
        if ( dataflag < 1 ) return;
	
		// If only half of a transform is stored, it needs to be handled
		if (transform == Hermitian || transform == CentHerm )
			data.setXdim(XSIZE(data)/2 + 1);

		// Reset select to get the correct offset
		if ( select_img < 0 ) select_img = 0;

		size_t myoffset, readsize, readsize_n, pagemax = 1073741824; //1Gb
		size_t datatypesize=gettypesize(datatype);
		size_t pagesize  =ZYXSIZE(data)*datatypesize;
		size_t pagemax_n = ROUND(pagemax/datatypesize);
		size_t haveread_n=0;

		char*	page = NULL;
		char*	padpage = NULL;
	
		// Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
		data.coreAllocate();
		myoffset = offset + select_img*(pagesize + pad);

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

		if ( pad > 0) padpage = (char *) askMemory(pad*sizeof(char));
		fseek( fimg, myoffset, SEEK_SET );
		for ( size_t myn=0; myn<NSIZE(data); myn++ )
		{
			for (size_t myj=0; myj<pagesize; myj+=pagemax )
			{
				// Read next page. Divide pages larger than pagemax
				readsize = pagesize - myj;
				if ( readsize > pagemax ) readsize = pagemax;
				readsize_n = readsize/datatypesize;

				//Read page from disc
				fread( page, readsize, 1, fimg );
				//swap per page
				if (swap) swapPage(page, readsize_n, datatype);
				// cast to T per page
				castPage2T(page, MULTIDIM_ARRAY(data) + haveread_n, datatype, readsize_n);
				haveread_n += readsize_n;
			}
			if ( pad > 0 ) fread( padpage, pad, 1, fimg);
		}
		if ( pad > 0 )
			freeMemory(padpage, pad*sizeof(char));
	
#ifdef DEBUG
        printf("DEBUG img_read_data: Finished reading and converting data\n");
#endif

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


    /** Get file name
     *
     * @code
     * std::cout << "Image name = " << I.name() << std::endl;
     * @endcode
     */
     FileName name() const
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
     
     /** Get Rot angle
     *
     * @code
     * std::cout << "First Euler angle " << I.rot() << std::endl;
     * @endcode
     */
    double rot(unsigned long n = 0) const
    {
        return (image[n]).angleRot;
    }

    /** Get Tilt angle
     *
     * @code
     * std::cout << "Second Euler angle " << I.tilt() << std::endl;
     * @endcode
     */
    double tilt(unsigned long n = 0) const
    {
        return (image[n]).angleTilt;
    }

    /** Get Psi angle
     *
     * @code
     * std::cout << "Third Euler angle " << I.psi() << std::endl;
     * @endcode
     */
    double psi(unsigned long n = 0) const
    {
        return (image[n]).anglePsi;
    }

    /** Get Xoff
     *
     * @code
     * std::cout << "Origin offset in X " << I.Xoff() << std::endl;
     * @endcode
     */
    double Xoff(unsigned long n = 0) const
    {
        return (image[n]).shiftX;
    }

    /** Get Yoff
     *
     * @code
     * std::cout << "Origin offset in Y " << I.Yoff() << std::endl;
     * @endcode
     */
    double Yoff(unsigned long n = 0) const
    {
        return (image[n]).shiftY;
    }

    /** Get Zoff
     *
     * @code
     * std::cout << "Origin offset in Z " << I.Zoff() << std::endl;
     * @endcode
     */
    double Zoff(unsigned long n = 0) const
    {
        return (image[n]).shiftZ;
    }

    /** Get Weight
    *
    * @code
    * std::cout << "weight= " << I.weight() << std::endl;
    * @endcode
    */
   double weight(unsigned long n = 0) const
   {
       return (image[n]).weight;
   }

   /** Get Flip
   *
   * @code
   * std::cout << "flip= " << flip() << std::endl;
   * @endcode
   */
   bool flip(unsigned long n = 0) const
   {
	   return (image[n]).flip;
   }

    /** Set file name
     *
     */
     void setName(const FileName _filename)
     {
        filename = _filename;
     }

    /** Set Euler angles in image header
     *
     */
    void setEulerAngles(double rot, double tilt, double psi, 
                        unsigned long n = 0)
    {
        (image[n]).angleRot = rot;
        (image[n]).angleTilt = rot;
        (image[n]).anglePsi = rot;
    }

    /** Set origin offsets in image header
     *
     */
    void setShifts(double xoff, double yoff, double zoff = 0., 
                        unsigned long n = 0)
    {
        (image[n]).shiftX = xoff;
        (image[n]).shiftY = yoff;
        (image[n]).shiftZ = zoff;
    }

    /** Set flip in image header
     *
     */
    void setFlip(bool _flip, unsigned long n = 0)
    {
        (image[n]).flip = _flip;
    }

     /** Set Weight in image header
     *
     */
    void setWeight(bool _weight, unsigned long n = 0)
    {
        (image[n]).weight = _weight;
    }

   /** Get geometric transformation matrix from 2D-image headerq
     */
    Matrix2D< double > getTransformationMatrix(bool only_apply_shifts = false, 
                                               unsigned long n = 0)
    {
        // This has only been implemented for 2D images...
        (*this)().checkDimension(2);
        
        double phi = realWRAP((image[n]).angleRot, 0., 360.);
        double psi = realWRAP((image[n]).anglePsi, -180, 180);
        double theta = realWRAP((image[n]).angleTilt, 0., 360);
        double xoff = (image[n]).shiftX;
        double yoff = (image[n]).shiftY;

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
        if ((image[n]).flip)
        {
            A(0, 0) = -A(0, 0);
            A(0, 1) = -A(0, 1);
        }

        return A;
    }
    
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

    	o << "dimensions   : " << ZSIZE(I()) << " x " << YSIZE(I()) << " x " << XSIZE(I());
    	o << "  (slices x rows x columns)" << std::endl;
    	o << "Euler angles : " << std::endl;
    	o << "  Phi   (rotation around Z axis) = " << I.rot() << std::endl;
    	o << "  theta (tilt, second rotation around new Y axis) = " << I.tilt() << std::endl;
    	o << "  Psi   (third rotation around new Z axis) = " << I.psi() << std::endl;
    	o << "Origin Offsets : " << std::endl;
    	o << "  Xoff  (origin offset in X-direction) = " << I.Xoff() << std::endl;
    	o << "  Yoff  (origin offset in Y-direction) = " << I.Yoff() << std::endl;
    	o << "  Zoff  (origin offset in Z-direction) = " << I.Zoff() << std::endl;
    	o << "Header size  : " << I.offset << std::endl;
    	o << "Weight  : " << I.weight() << std::endl;
    	o << "Flip    : " << I.flip() << std::endl;
    	 return o;
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
