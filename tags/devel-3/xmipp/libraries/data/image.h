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


#define SWAPTRIG     65535   // Threshold file z size above which bytes are swapped


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

    unsigned int        dataflag;	// Flag to force reading of the data
    FileName            filename;       // File name
    unsigned long	x, y, z;        // Image dimensions in X, Y and Z
    unsigned long	n, i;		// Number of images and image number (may be > n)
    unsigned long	px, py, pz; 	// Page dimensions
    unsigned long	offset; 	// Data offset
    DataType		datatype;	// Data type
    TransformType	transform;  	// Transform type
    double		min, max;	// Limits
    double		avg, std;	// Average and standard deviation
    double		smin, smax; 	// Limits for display
    double		scale;		// Scale of last density conversion operation
    double		shift;		// Shift of last density conversion operation before scaling
    double		resolution; 	// Resolution limit of data - used for low-pass filtering
    double		ux, uy, uz;	// Voxel units (angstrom/pixel edge)
    double		ua, ub, uc; 	// Unit cell dimensions (angstrom)
    double		alf, bet, gam;	// Unit cell angles (radian)
    unsigned int	spacegroup;	// Space group
    SubImage*           image;          // Sub-images
    MultidimArray<T>    data;           // The image data array

public:
    /** Empty constructor
     *
     * An empty image is created.
     *
     * @code
     * Image I;
     * @endcode
     */
    Image()
    {
        image = NULL;
        clear();
    }

    /** Constructor with DataType and size
     *
     * A blank image (0.0 filled) is created with the given size. Pay attention
     * to the dimension order: Y and then X.
     *
     * @code
     * Image I(64,64);
     * @endcode
     */
    Image(DataType type, int Xdim, int Ydim, int Zdim=1, int Ndim=1)
    {
        image = NULL;
        clear();
        datatype = type;
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
        datatype = Unknown_Type;
        x = y = z = n = 0;
        px = py = pz = 0;
        filename = "";
        offset = 0;
        transform = NoTransform;
        min = max = avg = 0.;
        std = -1.;
        scale = shift = 1.;
        resolution = 0.;
        ux = uy = uz = 1.;
        ua = ub = uc = 1.;
        alf, bet, gam = DEG2RAD(90.);
        spacegroup = 1;
        if (image != NULL)
            delete[] image;
        image = NULL;
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

    /** General read function
     */
     void read(const FileName name, bool readdata=true, int select_img=-1)
     {
         int err = 0;
         FileName basename, extension;
         
         filename = name;

         if ( readdata ) dataflag = 1;
         /*
           Bimage* 	p = init_img();
           if ( filename.contains("#") )
           p->filename = filename;
           else
           p->filename = clean_filename;
         */

         if (filename.get_extension()=="mrc")
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
	
         if ( err < 0 ) {
             REPORT_ERROR(10,"Error reading file");
         }

     }

    /** General write function
     */
     void write(const FileName name)
     {
         int err = 0;
         FileName basename, extension;
         
         filename = name;

         // PERHAPS HERE CHECK FOR INCONSISTENCIES BETWEEN data.xdim and x, etc???

         if (filename.get_extension()=="mrc")
             //REPORT_ERROR(22,"writeMRC NOT IMPLEMENTED YET");
             err = writeMRC();
         else
             err = writeSPIDER();

         if ( err < 0 ) {
             std::cerr<<" Filename = "<<filename<<" Extension= "<<extension<<std::endl;
             REPORT_ERROR(10,"Error writing file");
         }

     }


/** Cast a page of data from type dataType to type Tdest
 * @ingroup LittleBigEndian
 *    input pointer  char *
 */

     void castPage2T(char * page, T * ptrDest, size_t pageSize )
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

/** Cast page from T to outputDataType
 * @ingroup LittleBigEndian
 *    input pointer  char *
 */
     void castPage2Datatype(T * srcPtr, char * page, DataType outputDataType, size_t pageSize )

     {

         switch (outputDataType)
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
                 for(int i=0; i<pageSize;i++) ptr[i] = (float)srcPtr[i];
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
             std::cerr<<"Datatype= "<<datatype<<std::endl;
             REPORT_ERROR(16," ERROR: cannot cast datatype to T");
             break;
         }
         }
     }

     void swapPage(char * page, size_t pageNrElements, int swap)
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

     void readData(FILE* fimg, int select_img, int swap, unsigned long pad)
     {
#define DEBUG
#ifdef DEBUG
         std::cerr<<"entering readdata"<<std::endl;
         std::cerr<<" readData flag= "<<dataflag<<std::endl;
#endif
        if ( dataflag < 1 ) return;
	
	if ( px < 1 ) px = x;
	if ( py < 1 ) py = y;
	if ( pz < 1 ) pz = z;
	
	// If only half of a transform is stored, it need to be handled
	unsigned long xstore = x;
	unsigned long xpage = px;
	if (transform == Hermitian || transform == CentHerm ) {
            xstore = x/2 + 1;
		if ( px > xstore ) xpage = xstore;
	}
	
	// Reset select to get the correct offset
	if ( select_img < 0 ) select_img = 0;

        size_t myoffset, readsize, readsize_n, pagemax = 1073741824; //1Gb
        size_t datatypesize=gettypesize(datatype);
        size_t pagesize_n=xpage*py*pz;
        size_t pagesize  =pagesize_n*datatypesize;
        size_t pagesize_t=pagesize_n*sizeof(T);
        size_t pagemax_n = ROUND(pagemax/datatypesize);
        size_t haveread_n=0;

 	char*	page = NULL;
	char*	padpage = NULL;

        // Allocate memory for image data
        data.resize(n, z, y, xstore);

#ifdef DEBUG
        data.printShape();
        printf("DEBUG img_read_data: Data size: %ld %ld %ld %ld\n", XSIZE(data), YSIZE(data), ZSIZE(data), NSIZE(data));
        printf("DEBUG img_read_data: Page size: %ld %ld %ld (%ld)\n", px, py, pz, pagesize);
        printf("DEBUG img_read_data: Swap = %d  Pad = %ld  Offset = %ld\n", swap, pad, offset);
#endif 

	if ( XSIZE(data) == px && YSIZE(data) == py && ZSIZE(data) == pz ) { // 3D block
            myoffset = offset + select_img*(pagesize + pad);
#ifdef DEBUG
            std::cerr<<"DEBUG img_read_data: Reading 3D blocks: myoffset = "
                     <<myoffset<<" select_img= "<<select_img<<" pagesize= "<<pagesize<<" offset= "<<offset<<std::endl;
            std::cerr<<"devel3"<<std::endl;
#endif
            if (pagesize > pagemax)
                page = (char *) askMemory(pagemax*sizeof(char));
            else
                page = (char *) askMemory(pagesize*sizeof(char));

            if ( pad > 0) padpage = (char *) askMemory(pad*sizeof(char));
            fseek( fimg, myoffset, SEEK_SET );
            std::cerr<<"NSIZE(data)= "<<NSIZE(data)<<std::endl;
            for ( size_t myn=0; myn<NSIZE(data); myn++ ) 
            {
                for (size_t myj=0; myj<pagesize; myj+=pagemax ) 
                {
                    std::cerr<<"myn= "<<myn<<" myj= "<<myj<<" pagesize= "<<pagesize<<std::endl;

                    // Read next page. Divide pages larger than pagemax 
                    readsize = pagesize - myj;
                    if ( readsize > pagemax ) readsize = pagemax;
                    readsize_n = readsize/datatypesize;
                    std::cerr<<"readsize_n= "<<readsize_n <<" readsize= "<<readsize<<" datatypesize= "<<datatypesize<<std::endl;
               
                    //Read page from disc
                    fread( page, readsize, 1, fimg );
                    std::cerr<<"after fread"<<std::endl;
#ifdef DEBUG
                    float* tt;
                    tt=(float*)page;
                    std::cerr<<"page="<<tt[0]<<std::endl;
 #endif

                    //swap per page
                    if (swap) swapPage(page, readsize_n, swap);
                    std::cerr<<"after swap haveread_n= "<<haveread_n<<std::endl;

                    // cast to T per page
                    castPage2T(page, MULTIDIM_ARRAY(data) + haveread_n, readsize_n);
                    std::cerr<<"after cast"<<std::endl;

#ifdef DEBUG
                    T* ttt;
                    ttt=(T*)(MULTIDIM_ARRAY(data) + haveread_n);
                    std::cerr<<"data="<<ttt[0]<<std::endl;
#endif
                    haveread_n += readsize_n;
               }
               if ( pad > 0 ) fread( padpage, pad, 1, fimg);
            }
            if ( pad > 0 ) 
                freeMemory(padpage, pad*sizeof(char));
	} 
        else 
        {
            REPORT_ERROR(12,"BUG: entering bsoft rwimg code that was meant for unsupported MFF, DSN6 or BRIX format...");
        }
	
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

};

// Special case for complex numbers
template<>
void Image< std::complex< double > >::castPage2T(char * page, 
                                                 std::complex<double> * ptrDest, 
                                                 size_t pageSize);
template<>
void Image< std::complex< double > >::castPage2Datatype(std::complex< double > * srcPtr, 
                                                        char * page, 
                                                        DataType outputDataType, 
                                                        size_t pageSize);


#endif
