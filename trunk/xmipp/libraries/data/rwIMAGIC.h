/*
 * rwIMAGIC.h
 *
 *  Created on: May 17, 2010
 *      Author: roberto
 */
/*
 Base on rwIMAGIC.h
 Header file for reading and writing Image Science's Imagic files
 Format: 2D image file format for the program Imagic (Image Science)
 Author: Bernard Heymann
 Created: 19990424  Modified: 20011030
*/

#ifndef RWIMAGIC_H_
#define RWIMAGIC_H_

#define IMAGICSIZE 1024 // Size of the IMAGIC header for each image

struct IMAGIChead
{             // file header for IMAGIC data
    int imn;          //  0      image location number (1,2,...)
    int ifn;          //  1      # images following
    int ierror;       //  2      error code: error if >0
    int nhfr;         //  3      # header records per image
    int ndate;        //  4      creation day
    int nmonth;       //  5      creation month
    int nyear;        //  6      creation year
    int nhour;        //  7      creation hour
    int nminut;       //  8      creation minute
    int nsec;         //  9      creation second
    int npix2;        // 10      # 4-byte reals in image
    int npixel;       // 11      # image elements
    int ixlp;       // 12      lines per image (Y)
    int iylp;        // 13      pixels per line (X)
    char type[4];      // 14      image type
    int ixold;       // 15      top-left X coordinate
    int iyold;       // 16      top-left Y coordinate
    float avdens;       // 17      average
    float sigma;       // 18      standard deviation
    float varian;       // 19      variance
    float oldavd;      // 20      old average
    float densmax;       // 21      maximum
    float densmin;       // 22      minimum
    //     double sum;       // 23+24  sum of densities
    //     double squares;    // 25+26  sum of squares
    float dummy[4];   // 23-26  dummy place holder
    char lastpr[8];      // 27+28     last program writing file
    char name[80];       // 29-48     image name
    float extra_1[8];   // 49-56     additional parameters
    float eman_alt;   // 57      EMAN: equiv to psi & PFT omega
    float eman_az;    // 58      EMAN: equiv to theta
    float eman_phi;   // 59      EMAN: equiv to phi
    float extra_2[69];   // 60-128     additional parameters
    float euler_alpha;  // 129   Euler angles: psi
    float euler_beta;  // 130       theta
    float euler_gamma;  // 131       phi
    float proj_weight;  // 132   weight of each projection
    float extra_3[66];   // 133-198     additional parameters
    char history[228];      // 199-255   history
} ;

/************************************************************************
@Function: readIMAGIC
@Description:
 Reading an IMAGIC image format.
@Algorithm:
 A 2D file format for the IMAGIC package.
 The header is stored in a separate file with extension ".hed" and
  a fixed size of 1024 bytes per image.
 The image data is stored in a single block in a file with the
  extension ".img".
 Byte order determination: Year and hour values
        must be less than 256*256.
 Data types:     PACK = byte, INTG = short, REAL = float,
        RECO,COMP = complex float.
 Transform type:    Centered (COMP data type)
        RECO is not a transform
 Note that the x and y dimensions are interchanged (actually a display issue).
@Arguments:
 Bimage* p   the image structure.
 int select   image selection in multi-image file (-1 = all images).
@Returns:
 int     error code (<0 means failure).
**************************************************************************/
int  readIMAGIC(int img_select)
{
#define DEBUG
#undef DEBUG
 #ifdef DEBUG
    printf("DEBUG readIMAGIC: Reading Imagic file\n");
#endif
    // get the filename without extension to find the header file
    FileName headername;
    headername = filename.substr(0, filename.find_last_of('.')) + ".hed";
    filename   = filename.substr(0, filename.find_last_of('.')) + ".img";


    // open and read the header file
    FILE  *fhed;
    if ( ( fhed = fopen(headername.c_str(), "r") ) == NULL )
        REPORT_ERROR(1,(std::string)"readIMAGIC: header file " +headername+ " does not exist");
    ;

    IMAGIChead* header = new IMAGIChead;

    if ( fread( header, IMAGICSIZE, 1, fhed ) < 1 )
        REPORT_ERROR(1,(std::string)"readIMAGIC: header file " +headername+ " cannot be read");
    ;

    // Determine byte order and swap bytes if from little-endian machine
    char*   b = (char *) header;
    int    swap = 0;
    unsigned long i, extent = IMAGICSIZE - 916;  // exclude char bytes from swapping
    if ( ( abs(header->nyear) > SWAPTRIG ) || ( header->ixlp > SWAPTRIG ) )
    {
        swap = 1;
        for ( i=0; i<extent; i+=4 )
            if ( i != 56 )          // exclude type string
                swapbytes(b+i, 4);
    }
    int _xDim,_yDim,_zDim;
    unsigned long int _nDim;
    _xDim = (int) header->iylp;
    _yDim = (int) header->ixlp;
    _zDim = (int) 1;
    _nDim = (unsigned long int) header->ifn + 1 ;

    std::stringstream Num;
    std::stringstream Num2;
    if ( img_select > (int)_nDim )
    {
        Num  << img_select;
        Num2 << _nDim;
        REPORT_ERROR(1,(std::string)"readImagic: Image number " + Num.str() +
                     " exceeds stack size " + Num2.str());
    }

    if( img_select > -1)
        _nDim=1;
    data.setDimensions( //setDimensions do not allocate data
        _xDim,
        _yDim,
        _zDim,
        _nDim );
    unsigned long   imgStart=0;
    unsigned long   imgEnd =_nDim;
    if (img_select != -1)
    {
    	imgStart=img_select;
    	imgEnd=img_select+1;
    }

    DataType datatype;

    if ( strstr(header->type,"PACK") )
        datatype = UChar;
    else if ( strstr(header->type,"INTG") )
        datatype = Short;
    else if ( strstr(header->type,"REAL") )
        datatype = Float;
    else if ( strstr(header->type,"RECO") )
    {
        datatype = ComplexFloat; // Complex data
        transform = NoTransform;
    }
    else if ( strstr(header->type,"COMP") )
    {
        datatype = ComplexFloat; // Complex transform data
        transform = Centered;
    }

    // Set min-max values and other statistical values
    if ( header->sigma == 0 && header->varian != 0 )
        header->sigma = sqrt(header->varian);
    if ( header->densmax == 0 && header->densmin == 0 && header->sigma != 0 )
    {
        header->densmin = header->avdens - header->sigma;
        header->densmax = header->avdens + header->sigma;
    }

    MDMainHeader.clear();
    MDMainHeader.addObject();
    MDMainHeader.setValue(MDL_MIN,(double)header->densmin);
    MDMainHeader.setValue(MDL_MAX,(double)header->densmax);
    MDMainHeader.setValue(MDL_AVG,(double)header->avdens);
    MDMainHeader.setValue(MDL_STDDEV,(double)header->sigma);
    MDMainHeader.setValue(MDL_SAMPLINGRATEX,(double)1.);
    MDMainHeader.setValue(MDL_SAMPLINGRATEY,(double)1.);
    MDMainHeader.setValue(MDL_SAMPLINGRATEZ,(double)1.);

    offset = 0;   // separate header file

    unsigned long   Ndim = _nDim, j = 0;

    // View   view;
    char*   hend;
    /*already done in image.h
    if ( img_select > -1 )
{
        if ( img_select >= (int) Ndim )
            REPORT_ERROR(1,(std::string)"readIMAGIC: Image number " + str(img_select)+" exceeds stack size");
        data.setNdim(1);
        Ndim = 1;
        i = img_select;
}
    */


    // Get the header information
    if ( img_select > -1 )
        fseek( fhed, 0, SEEK_SET );
    else
        fseek( fhed, img_select * IMAGICSIZE, SEEK_SET );

    MD.clear();
    for ( i=imgStart; i<imgEnd; i++ )
    //for ( i=0; i<Ndim; i++ )
    {
        if ( fread( header, IMAGICSIZE, 1, fhed ) < 1 )
            return(-2);
        //if ( img_select < 0 || img_select == i )
        {
            hend = (char *) header + extent;
            if ( swap )
                for ( b = (char *) header; b<hend; b+=4 )
                    swapbytes(b, 4);
            MD.addObject();
            MD.setValue(MDL_ORIGINX,  (double)-1. * header->iyold);
            MD.setValue(MDL_ORIGINY,  (double)-1. * header->ixold);
            MD.setValue(MDL_ORIGINZ,  zeroD);
            MD.setValue(MDL_ANGLEROT, (double)-1. * header->euler_alpha);
            MD.setValue(MDL_ANGLETILT,(double)-1. * header->euler_beta);
            MD.setValue(MDL_ANGLEPSI, (double)-1. * header->euler_gamma);
            MD.setValue(MDL_WEIGHT,   (double)oneD);

            j++;
        }
    }

    fclose(fhed);

    delete header;

    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "r") ) == NULL )
        REPORT_ERROR(1,(std::string)"readIMAGIC: image file " +filename+ " does not exist");

    int pad=0;
    readData(fimg, img_select, datatype, pad );

    fclose(fimg);

    return(0);
}

/************************************************************************
@Function: writeIMAGIC
@Description:
 Writing an IMAGIC image format.
@Algorithm:
 A file format for the IMAGIC package.
@Arguments:
 Bimage*    the image structure.
@Returns:
 int     error code (<0 means failure).
**************************************************************************/
int  writeIMAGIC()
{
#ifdef NEVERDEFINED
    //    if ( p->transform != NoTransform )
    //        img_convert_fourier(p, Centered);

    IMAGIChead* header = new IMAGIChead;

    // fill in the file header
    header->nhfr = 1;
    header->npix2 = XSIZE(data)*YSIZE(data);
    header->npixel = header->npix2;
    header->iylp = XSIZE(data);
    header->ixlp = YSIZE(data);

    time_t timer;
    time ( &timer );
    tm* t = localtime(&timer);

    header->ndate = t->tm_mday;
    header->nmonth = t->tm_mon + 1;
    header->nyear = t->tm_year;
    header->nhour = t->tm_hour;
    header->nminut = t->tm_min;
    header->nsec = t->tm_sec;


    // Convert T to datatype
    if ( typeid(T) == typeid(double) ||
         typeid(T) == typeid(float) ||
         typeid(T) == typeid(int) )
        strcpy(header->type,"REAL");
    else if ( typeid(T) == typeid(unsigned char) ||
              typeid(T) == typeid(signed char) )
        strcpy(header->type,"PACK");
    else if ( typeid(T) == typeid(std::complex<float>) ||
              typeid(T) == typeid(std::complex<double>) )
        strcpy(header->type,"COMP");
    else
        REPORT_ERROR(1,"ERROR write IMAGIC image: invalid typeid(T)");

    double aux;

    if (MDMainHeader.firstObject() != MetaData::NO_OBJECTS_STORED)
    {


        if(MDMainHeader.getValue(MDL_MIN,   aux))
            header->densmin = (float)aux;
        if(MDMainHeader.getValue(MDL_MAX,   aux))
            header->densmax = (float)aux;
        if(MDMainHeader.getValue(MDL_AVG,   aux))
            header->avdens   = (float)aux;
        if(MDMainHeader.getValue(MDL_STDDEV,aux))
        {
            header->sigma  = (float)aux;
            header->varian = (float)(aux*aux);
        }
    }

    memcpy(header->lastpr, "Xmipp", 5);
    memcpy(header->name, filename.c_str(), 80);

    // get the filename without extension to find the header file
    FileName headername;
    headername = filename.substr(0, filename.find_last_of('.')) + ".hed";

    FILE        *fhed, *fimg;
    if ( ( fhed = fopen(headername.c_str(), "w") ) == NULL )
        REPORT_ERROR(1,(std::string)"writeIMAGIC: header file " +headername+ " cannot be created");
    if ( ( fhed = fopen(headername.c_str(), "w") ) == NULL )
        REPORT_ERROR(1,(std::string)"writeIMAGIC: data file "   +filename+ " cannot be created");
    MD.firstObject();
    for (unsigned long i=0; i<NSIZE(data); i++ )
    {
        long int next_result = MD.nextObject();
        if (next_result != MetaData::NO_OBJECTS_STORED && next_result != MetaData::NO_MORE_OBJECTS)
        {
            header->imn = i + 1;   // Image number
            header->ifn = NSIZE(data) - i - 1; // Number of images following this one
            if(MD.getValue(MDL_ORIGINX,  aux))
                header->iyold  = (int) (-aux + 0.5);
            if(MD.getValue(MDL_ORIGINY,  aux))
                header->ixold  = (int) (-aux + 0.5);

            if(MD.getValue(MDL_ORIGINZ,  aux))
                header->zoff  =(float)-aux;
            else
                header->zoff  =0.;

            if(MD.getValue(MDL_ANGLEROT, aux))
                header->euler_alpha = header->eman_alt = (float)-aux;
            if(MD.getValue(MDL_ANGLETILT,aux))
                header->euler_beta  = header->eman_az  = (float)-aux;
            if(MD.getValue(MDL_ANGLEPSI, aux))
                header->euler_gamma = header->eman_phi = (float)-aux;
            fwrite( header, IMAGICSIZE, 1, fhed );

            if ( strstr(header->type,"COMP"))
                writePageAsDatatype(fimg, ComplexFloat, datasize_n);
            castPage2Datatype(MULTIDIM_ARRAY(data) + i*datasize_n, fdata, ComplexFloat, datasize_n);
            else if ( strstr(header->type,"REAL"))
                castPage2Datatype(MULTIDIM_ARRAY(data) + i*datasize_n, fdata, Float, datasize_n);
            else if ( typeid(T) == typeid(unsigned char))
                castPage2Datatype(MULTIDIM_ARRAY(data) + i*datasize_n, fdata, UChar, datasize_n);
            else if ( typeid(T) == typeid( char))
                castPage2Datatype(MULTIDIM_ARRAY(data) + i*datasize_n, fdata, SChar, datasize_n);
            else
                REPORT_ERROR(var,"writeIMAGIC : Unknown type");
            fwrite( fdata, datasize, 1, fimg );
        }
    }

    fclose(fhed);

    delete header;


    unsigned long datatypesize = gettypesize(p->datatype);
    unsigned long  datasize = p->x*p->y*p->z*p->n*datatypesize;

    FILE        *fimg;

    if ( p->dataflag )
    {
        if ( ( fimg = fopen(p->filename.c_str(), "w") ) == NULL )
            return(-3);
        fwrite( p->data, datasize, 1, fimg );
        fclose(fimg);
    }
#endif
    return(0);
}

#endif /* RWIMAGIC_H_ */
