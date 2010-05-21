/*
 rwSPIDER.h
 Header file for reading and writing SPIDER files
 Format: 3D image file format for the SPIDER package
 Author: Bernard Heymann
 Created: 19990410  Modified: 20010928
*/

#ifndef RWSPIDER_H
#define RWSPIDER_H

#define SPIDERSIZE 1024 // Minimum size of the SPIDER header (variable)
struct SPIDERhead
{          // file header for SPIDER data
    float nslice;    //  0      slices in volume (image = 1)
    float nrow;         //  1      rows per slice
    float irec;         //  2      # records in file (unused)
    float nhistrec;      //  3      (obsolete)
    float iform;        //  4      file type specifier
    float imami;        //  5      max/min flag (=1 if calculated)
    float fmax;       //  6      maximum
    float fmin;       //  7      minimum
    float av;         //  8      average
    float sig;        //  9      standard deviation (=-1 if not calculated)
    float ihist;        // 10      (obsolete)
    float nsam;         // 11      pixels per row
    float labrec;      // 12      # records in header
    float iangle;       // 13      flag: tilt angles filled
    float phi;        // 14      tilt angles
    float theta;      // 15
    float gamma;      // 16      (=psi)
    float xoff;       // 17      translation
    float yoff;       // 18
    float zoff;       // 19
    float scale;      // 20      scaling
    float labbyt;       // 21      # bytes in header
    float lenbyt;       // 22      record length in bytes (row length)
    float istack;       // 23      indicates stack of images
    float inuse;        // 24      indicates this image in stack is used (not used)
    float maxim;        // 25      max image in stack used
    float imgnum;        // 26      number of current image
    float unused[2]; // 27-28     (unused)
    float kangle;       // 29      flag: additional angles set
    float phi1;       // 30      additional angles
    float theta1;      // 31
    float psi1;       // 32
    float phi2;       // 33
    float theta2;      // 34
    float psi2;       // 35

    double fGeo_matrix[3][3]; // x9 = 72 bytes: Geometric info
    float fAngle1; // angle info

    float fr1;
    float fr2; // lift up cosine mask parameters

    /** Fraga 23/05/97  For Radon transforms **/
    float RTflag; // 1=RT, 2=FFT(RT)
    float Astart;
    float Aend;
    float Ainc;
    float Rsigma; // 4*7 = 28 bytes
    float Tstart;
    float Tend;
    float Tinc; // 4*3 = 12, 12+28 = 40B

    /** Sjors Scheres 17/12/04 **/
    float weight; // For Maximum-Likelihood refinement
    float flip;   // 0=no flipping operation (false), 1=flipping (true)

    char fNada2[576]; // empty 700-76-40=624-40-8= 576 bytes

    char cdat[12];   // 211-213   creation date
    char ctim[8];  // 214-215   creation time
    char ctit[160];  // 216-255   title
} ;


/************************************************************************
@Function: readSPIDER
@Description:
 Reading a SPIDER image file format.
@Algorithm:
 A 3D multi-image format used in electron microscopy.
 Header size:    1024 bytes (not same as data offset!).
 Data offset:    sizeof(float)*x_size*ceil(1024/x_size)
 File format extensions:   .spi
 The identifier is a 4-byte machine stamp:
    1 Big-endian IEEE  17 17 00 00
                2 VAX     34 65 00 00
    3 Cray    -
                4 Little-endian IEEE 68 65 00 00
                5 Convex    85 17 00 00
    6 Fijitsu VP   -
    (Note: not always implemented - so unreliable)
 Byte order determination: File type and third dimension values
        must be less than 256*256.
 Data type:      only float.
 Transform type:    Hermitian
        The x-dimension contains the x-size
        of the full transform
 A multi-image file has a global header followed by a header and data
 for each sub-image.
@Arguments:
 Bimage* p   the image structure.
 int img_select  image selection in multi-image file (-1 = all images).
@Returns:
 int     error code (<0 means failure).
**************************************************************************/
int  readSPIDER(int img_select)
{
//#define DEBUG
#undef DEBUG
#ifdef DEBUG
    printf("DEBUG readSPIDER: Reading Spider file\n");
#endif

    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "r") ) == NULL )
        REPORT_ERROR(1,"rwSPIDER: cannot read image.");

    SPIDERhead* header = (SPIDERhead *) askMemory(sizeof(SPIDERhead));
    if ( fread( header, SPIDERSIZE, 1, fimg ) < 1 )
        REPORT_ERROR(1,"rwSPIDER: cannot allocate memory for header");

    swap = 0;

    // Determine byte order and swap bytes if from different-endian machine
    char*    b = (char *) header;
    int      i, j;
    int      extent = SPIDERSIZE - 180;  // exclude char bytes from swapping
    if ( ( fabs(header->nslice) > SWAPTRIG ) || ( fabs(header->iform) > SWAPTRIG ) ||
         ( fabs(header->nslice) < 1 ) )
    {
#ifdef DEBUG
        fprintf(stderr, "Warning: Swapping header byte order for 4-byte types\n");
#endif

        swap = 1;
        for ( i=0; i<extent; i+=4 )
            swapbytes(b+i, 4);
    }

    // Map the parameters
    data.setDimensions(
        (int) header->nsam,
        (int) header->nrow,
        (int) header->nslice,
        (unsigned long int)1 );
    DataType datatype = Float;
    transform = NoTransform;
    if ( header->iform < 0 )
    {
        transform = Hermitian;
        datatype = ComplexFloat;
        data.setXdim(XSIZE(data) - 2);
    }

    offset = (int) header->labbyt;
    MDMainHeader.clear();
    MDMainHeader.addObject();
    MDMainHeader.setValue(MDL_MIN,(double)header->fmin);
    MDMainHeader.setValue(MDL_MAX,(double)header->fmax);
    MDMainHeader.setValue(MDL_AVG,(double)header->av);
    MDMainHeader.setValue(MDL_STDDEV,(double)header->sig);
    MDMainHeader.setValue(MDL_SAMPLINGRATEX,(double)header->scale);
    MDMainHeader.setValue(MDL_SAMPLINGRATEY,(double)header->scale);
    MDMainHeader.setValue(MDL_SAMPLINGRATEZ,(double)header->scale);

    size_t header_size = offset;
    size_t image_size = header_size + ZYXSIZE(data)*sizeof(float);
    size_t pad = offset;

#ifdef DEBUG

    if (!(ABS(header->istack)<1e-10 || ABS(header->istack - 2)<1e-10))
    {
        std::cerr << "istack="<< header->istack<<std::endl;
        REPORT_ERROR(1,"INVALID ISTACK");
    }
#endif
    if ( header->istack > 0 )
        data.setNdim((int) header->maxim);

    unsigned long   Ndim = NSIZE(data), imgstart = 0;
    unsigned long   imgend = NSIZE(data);
    char*   hend;
    if ( img_select > -1 )
    {
        if ( img_select >= Ndim )
            img_select = Ndim - 1;
        imgstart = img_select;
        imgend = img_select + 1;
        data.setNdim(1);
        Ndim = 1;
        i = img_select;
    }
    MD.clear();
    MD.addObject();
    MD.setValue(MDL_ORIGINX,  (double)header->xoff);
    MD.setValue(MDL_ORIGINY,  (double)header->yoff);
    MD.setValue(MDL_ORIGINZ,  (double)header->zoff);
    MD.setValue(MDL_ANGLEROT, (double)header->phi);
    MD.setValue(MDL_ANGLETILT,(double)header->theta);
    MD.setValue(MDL_ANGLEPSI, (double)header->gamma);
    MD.setValue(MDL_WEIGHT,   (double)header->weight);
    bool baux;
    if(header->flip == 1)
        baux=true;
    else
        baux=false;
    MD.setValue(MDL_FLIP,     baux);


    if ( header->istack > 0 )
    {
        offset += offset;
        for ( i=imgstart; i<imgend; i++ )
        {
            fseek( fimg, header_size + i*image_size, SEEK_SET );
            if ( fread( header, SPIDERSIZE, 1, fimg ) < 1 )
                REPORT_ERROR(3,"rwSPIDER: cannot read multifile header information");
            hend = (char *) header + extent;
            if ( swap )
                for ( b = (char *) header; b<hend; b+=4 )
                    swapbytes(b, 4);
            //j = ( Ndim > 1 )? j = i: 0;
            MD.addObject();//I think first image is readed twice
            //test with spider stack
            MD.setValue(MDL_ORIGINX,header->xoff);
            MD.setValue(MDL_ORIGINY,header->yoff);
            MD.setValue(MDL_ORIGINZ,header->zoff);
            MD.setValue(MDL_ANGLEROT,header->phi);
            MD.setValue(MDL_ANGLETILT,header->theta);
            MD.setValue(MDL_ANGLEPSI,header->gamma);
            MD.setValue(MDL_WEIGHT,header->weight);
            MD.setValue(MDL_FLIP,header->flip);
        }
    }

    freeMemory(header, sizeof(SPIDERhead));

#ifdef DEBUG

    std::cerr<<"DEBUG readSPIDER: header_size = "<<header_size<<" image_size = "<<image_size<<std::endl;
    std::cerr<<"DEBUG readSPIDER: img_select= "<<img_select<<" n= "<<Ndim<<" pad = "<<pad<<std::endl;
#endif

    readData(fimg, img_select, datatype, pad );

    fclose(fimg);

    return(0);

}
/************************************************************************
@Function: writeSPIDER
@Description:
 Writing a SPIDER image file format.
@Algorithm:
 A 3D image format used in electron microscopy.
@Arguments:
 Bimage*    the image structure.
@Returns:
 int     error code (<0 means failure).
**************************************************************************/
int  writeSPIDER()
{
//#define DEBUG
#ifdef DEBUG
    printf("DEBUG writeSPIDER: Writing Spider file\n");
    printf("DEBUG writeSPIDER: File %s\n", filename.c_str());
#endif
#undef DEBUG

    float  lenbyt = sizeof(float)*XSIZE(data);  // Record length (in bytes)
    float  labrec = floor(SPIDERSIZE/lenbyt); // # header records
    if ( fmod(SPIDERSIZE,lenbyt) != 0 )
        labrec++;
    float  labbyt = labrec*lenbyt;   // Size of header in bytes
    offset = (int) labbyt;
    SPIDERhead* header = (SPIDERhead *) askMemory((int)labbyt*sizeof(char));

    // Map the parameters
    header->lenbyt = lenbyt;     // Record length (in bytes)
    header->labrec = labrec;     // # header records
    header->labbyt = labbyt;      // Size of header in bytes

    header->irec = labrec + floor((ZYXSIZE(data)*sizeof(float))/lenbyt + 0.999999); // Total # records
    header->nsam = XSIZE(data);
    header->nrow = YSIZE(data);
    header->nslice = ZSIZE(data);

    // If a transform, then the physical storage in x is only half+1
    size_t xstore = XSIZE(data);
    if ( transform == Hermitian )
    {
        xstore = XSIZE(data)/2 + 1;
        header->nsam = 2*xstore;
    }

#ifdef DEBUG
    printf("DEBUG writeSPIDER: Size: %g %g %g\n", header->nsam, header->nrow, header->nslice);
#endif

    if ( ZSIZE(data) < 2 )
    {
        if ( transform == NoTransform )
            header->iform = 1;     // 2D image
        else
            header->iform = -12 + (int)header->nsam%2;   // 2D Fourier transform
    }
    else
    {
        if ( transform == NoTransform )
            header->iform = 3;     // 3D volume
        else
            header->iform = -22 + (int)header->nsam%2;   // 3D Fourier transform
    }
    double aux;
    bool baux;
    header->imami = 1;
    /**
        std::vector< MDLabel >::iterator strIt;
        for( strIt  = SF.getActiveLabels().begin(); strIt != SF.getActiveLabels().end(); strIt ++ )
        {
            switch ((*strIt)) {
                case MDL_MIN:
                    MDc.getValue(MDL_MIN,aux);    header->fmin = (float)aux;
                    break;
                case MDL_MAX:
                    MDc.getValue(MDL_MAX,aux);    header->fmax = (float)aux;
                    break;
                case MDL_AVG:
                    MDc.getValue(MDL_AVG,aux);    header->av   = (float)aux;
                    break;
                case MDL_STDDEV:
                    MDc.getValue(MDL_STDDEV,aux); header->sig  = (float)aux;
                    break;

        }
    **/

    if (MDMainHeader.firstObject() != MetaData::NO_OBJECTS_STORED)
    {
        if(MDMainHeader.getValue(MDL_MIN,   aux))
            header->fmin = (float)aux;
        if(MDMainHeader.getValue(MDL_MAX,   aux))
            header->fmax = (float)aux;
        if(MDMainHeader.getValue(MDL_AVG,   aux))
            header->av   = (float)aux;
        if(MDMainHeader.getValue(MDL_STDDEV,aux))
            header->sig  = (float)aux;
    }
    // For multi-image files
    if (NSIZE(data) > 1 )
    {
        header->istack = 2;
        header->inuse = -1;
        header->maxim = NSIZE(data);
    }
    else
    {
        header->istack = 0;
        header->inuse = 0;
        header->maxim = 1;
    }

    if (MD.firstObject() != MetaData::NO_OBJECTS_STORED)
    {
        if(MD.getValue(MDL_ORIGINX,  aux))
            header->xoff  =(float)aux;
        if(MD.getValue(MDL_ORIGINY,  aux))
            header->yoff  =(float)aux;
        if(MD.getValue(MDL_ORIGINZ,  aux))
            header->zoff  =(float)aux;
        if(MD.getValue(MDL_ANGLEROT, aux))
            header->phi   =(float)aux;
        if(MD.getValue(MDL_ANGLETILT,aux))
            header->theta =(float)aux;
        if(MD.getValue(MDL_ANGLEPSI, aux))
            header->gamma =(float)aux;
        if(MD.getValue(MDL_WEIGHT,   aux))
            header->weight=(float)aux;
        if(MD.getValue(MDL_FLIP,    baux))
            header->flip  =(float)baux;
    }
    // Set time and date
    time_t timer;
    time ( &timer );
    tm* t = localtime(&timer);
    while ( t->tm_year > 100 )
        t->tm_year -= 100;
    sprintf(header->ctim, "%02d:%02d:%02d", t->tm_hour, t->tm_min, t->tm_sec);
    sprintf(header->cdat, "%02d-%02d-%02d", t->tm_mday, t->tm_mon, t->tm_year);

    size_t datasize, datasize_n;
    datasize_n = xstore*YSIZE(data)*ZSIZE(data);
    if (isComplexT())
        datasize = datasize_n * gettypesize(ComplexFloat);
    else
        datasize = datasize_n * gettypesize(Float);

#ifdef DEBUG

    printf("DEBUG writeSPIDER: Date and time: %s %s\n", header->cdat, header->ctim);
    printf("DEBUG writeSPIDER: Text label: %s\n", header->ctit);
    printf("DEBUG writeSPIDER: Header size: %g\n", header->labbyt);
    printf("DEBUG writeSPIDER: Header records and record length: %g %g\n", header->labrec, header->lenbyt);
    printf("DEBUG writeSPIDER: Data size: %ld\n", datasize);
    printf("DEBUG writeSPIDER: Data offset: %ld\n", offset);
    printf("DEBUG writeSPIDER: File %s\n", filename.c_str());
#endif

    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "w") ) == NULL )
        return(-1);
    fwrite( header, offset, 1, fimg );

    char* fdata = (char *) askMemory(datasize);
    if ( NSIZE(data) == 1 )
    {
        if (isComplexT())
            castPage2Datatype(MULTIDIM_ARRAY(data), fdata, ComplexFloat, datasize_n);
        else
            castPage2Datatype(MULTIDIM_ARRAY(data), fdata, Float, datasize_n);
        fwrite( fdata, datasize, 1, fimg );
    }
    else
    {
        header->istack = 0;
        header->inuse = 0;
        header->maxim = 0;
        for ( size_t i=0; i<NSIZE(data); i++ )
        {
            //header->imgnum = i + 1;
            long int next_result = MD.nextObject();
            if (next_result != MetaData::NO_OBJECTS_STORED && next_result != MetaData::NO_MORE_OBJECTS)
            {
                if(MD.getValue(MDL_ORIGINX,  aux))
                    header->xoff  =(float)aux;
                if(MD.getValue(MDL_ORIGINY,  aux))
                    header->yoff  =(float)aux;
                if(MD.getValue(MDL_ORIGINZ,  aux))
                    header->zoff  =(float)aux;
                if(MD.getValue(MDL_ANGLEROT, aux))
                    header->phi   =(float)aux;
                if(MD.getValue(MDL_ANGLETILT,aux))
                    header->theta =(float)aux;
                if(MD.getValue(MDL_ANGLEPSI, aux))
                    header->gamma =(float)aux;
                if(MD.getValue(MDL_WEIGHT,   aux))
                    header->weight=(float)aux;
                if(MD.getValue(MDL_FLIP,    baux))
                    header->flip  =(float)baux;

            }

            fwrite( header, offset, 1, fimg );
            if (isComplexT())
                castPage2Datatype(MULTIDIM_ARRAY(data) + i*datasize_n, fdata, ComplexFloat, datasize_n);
            else
                castPage2Datatype(MULTIDIM_ARRAY(data) + i*datasize_n, fdata, Float, datasize_n);
            fwrite( fdata, datasize, 1, fimg );
        }
    }

    fclose(fimg);
    freeMemory(fdata, datasize);
    freeMemory(header, (int)labbyt*sizeof(char));

    return(0);
}
#endif

