/*
        Base on rwMRC.h
        Header file for reading and writing MRC files
        Format: 3D crystallographic image file format for the MRC package
        Author: Bernard Heymann
        Created: 19990321       Modified: 20030723
*/

#ifndef RWMRC_H
#define RWMRC_H

#define MRCSIZE    1024 // Minimum size of the MRC header (when nsymbt = 0)

struct MRCheadold
{          // file header for MRC data
    int nx;              //  0   0       image size
    int ny;              //  1   4
    int nz;              //  2   8
    int mode;            //  3           0=uchar,1=short,2=float
    int nxStart;         //  4           unit cell offset
    int nyStart;         //  5
    int nzStart;         //  6
    int mx;              //  7           unit cell size in voxels
    int my;              //  8
    int mz;              //  9
    float a;             // 10   40      cell dimensions in A
    float b;             // 11
    float c;             // 12
    float alpha;         // 13           cell angles in degrees
    float beta;          // 14
    float gamma;         // 15
    int mapc;            // 16           column axis
    int mapr;            // 17           row axis
    int maps;            // 18           section axis
    float amin;          // 19           minimum density value
    float amax;          // 20   80      maximum density value
    float amean;         // 21           average density value
    int ispg;            // 22           space group number
    int nsymbt;          // 23           bytes used for sym. ops. table
    float extra[29];     // 24           user-defined info
    float xOrigin;       // 53           phase origin in pixels
    float yOrigin;       // 54
    int nlabl;           // 55           number of labels used
    char labels[10][80]; // 56-255       10 80-character labels
} ;

struct MRChead
{             // file header for MRC data
    int nx;              //  0   0       image size
    int ny;              //  1   4
    int nz;              //  2   8
    int mode;            //  3           0=char,1=short,2=float
    int nxStart;         //  4           unit cell offset
    int nyStart;         //  5
    int nzStart;         //  6
    int mx;              //  7           unit cell size in voxels
    int my;              //  8
    int mz;              //  9
    float a;             // 10   40      cell dimensions in A
    float b;             // 11
    float c;             // 12
    float alpha;         // 13           cell angles in degrees
    float beta;          // 14
    float gamma;         // 15
    int mapc;            // 16           column axis
    int mapr;            // 17           row axis
    int maps;            // 18           section axis
    float amin;          // 19           minimum density value
    float amax;          // 20   80      maximum density value
    float amean;         // 21           average density value
    int ispg;            // 22           space group number
    int nsymbt;          // 23           bytes used for sym. ops. table
    float extra[25];     // 24           user-defined info
    float xOrigin;       // 49           phase origin in pixels
    float yOrigin;       // 50
    float zOrigin;       // 51
    char map[4];         // 52       identifier for map file ("MAP ")
    char machst[4];      // 53           machine stamp
    float arms;          // 54       RMS deviation
    int nlabl;           // 55           number of labels used
    char labels[800];    // 56-255       10 80-character labels
} ;


// I/O prototypes
int readMRC(int img_select,bool isStack=false)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readMRC: Reading MRC file\n");
#endif

    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "r") ) == NULL )
        return(-1);

    MRChead*        header = (MRChead *) askMemory(sizeof(MRChead));
    if ( fread( header, MRCSIZE, 1, fimg ) < 1 )
        return(-2);

    // Determine byte order and swap bytes if from little-endian machine
    swap = 0;
    char*       b = (char *) header;
    int         i;
    if ( ( abs( header->mode ) > SWAPTRIG ) || ( abs(header->nz) > SWAPTRIG ) )
    {
#ifdef DEBUG
        fprintf(stderr, "Warning: Swapping header byte order for 4-byte types\n");
#endif

        swap = 1;
        int     extent = MRCSIZE - 800; // exclude labels from swapping
        for ( i=0; i<extent; i+=4 )
            swapbytes(b+i, 4);
    }

    // Convert VAX floating point types if necessary
    if ( header->amin > header->amax )
        REPORT_ERROR(1,"readMRC: amin > max: VAX floating point conversion unsupported");
    int _xDim,_yDim,_zDim;
    unsigned long int _nDim;
    _xDim = (int) header->nx;
    _yDim = (int) header->ny;
    _zDim = (int) header->nz;
    _nDim = 1;
    if(isStack)
    {
        _nDim = (unsigned long int) _zDim;
        _zDim = 1;
        replaceNsize=_nDim;
        std::stringstream Num;
        std::stringstream Num2;
        if ( img_select > (int)_nDim )
        {
            Num  << img_select;
            Num2 << _nDim;
            REPORT_ERROR(1,(std::string)"readImagic: Image number " + Num.str() +
                         " exceeds stack size " + Num2.str());
        }
    }
    else
        replaceNsize=0;

    // Map the parameters
    if (isStack && img_select==-1)
        data.setDimensions(_xDim, _yDim, 1, _nDim);
    else if(isStack && img_select!=-1)
        data.setDimensions(_xDim, _yDim, 1, 1);
    else
        data.setDimensions(_xDim, _yDim, _zDim,1);

    unsigned long   imgStart=0;
    unsigned long   imgEnd =_nDim;
    if (img_select != -1)
    {
        imgStart=img_select;
        imgEnd=img_select+1;
    }

    DataType datatype;
    switch ( header->mode%5 )
    {
    case 0:
        datatype = SChar;
        break;
    case 1:
        datatype = Short;
        break;
    case 2:
        datatype = Float;
        break;
    case 3:
        datatype = ComplexShort;
        break;
    case 4:
        datatype = ComplexFloat;
        break;
    default:
        datatype = UChar;
        break;
    }
    offset = MRCSIZE + header->nsymbt;
    if ( header->mode%5 > 2 && header->mode%5 < 5 )
    {
        transform = CentHerm;
        fseek(fimg, 0, SEEK_END);
        if ( ftell(fimg) > offset + 0.8*ZYXSIZE(data)*gettypesize(datatype) )
            data.setXdim(2 * (XSIZE(data) - 1));
        if ( header->mx%2 == 1 )
            data.setXdim(XSIZE(data) + 1);     // Quick fix for odd x-size maps
    }

    MDMainHeader.removeObjects();
    MDMainHeader.addObject();
    MDMainHeader.setValue(MDL_MIN,(double)header->amin);
    MDMainHeader.setValue(MDL_MAX,(double)header->amax);
    MDMainHeader.setValue(MDL_AVG,(double)header->amean);
    MDMainHeader.setValue(MDL_STDDEV,(double)header->arms);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    if ( header->mx && header->a!=0)//ux
        MDMainHeader.setValue(MDL_SAMPLINGRATEX,(double)header->a/header->mx);
    if ( header->my && header->b!=0)//yx
        MDMainHeader.setValue(MDL_SAMPLINGRATEY,(double)header->b/header->my);
    if ( header->mz && header->c!=0)//zx
        MDMainHeader.setValue(MDL_SAMPLINGRATEZ,(double)header->c/header->mz);

    int Ndim = NSIZE(data);
    if (isStack && dataflag<0)   // Don't read the individual header
    {                            // and the data if not necessary
        fclose(fimg);
        return 0;
    }


    MD.removeObjects();
    for ( i=imgStart; i<imgEnd; i++ )
        //for(int i=0;i< Ndim;i++)
    {
        MD.addObject();
        double aux;
        if(MDMainHeader.getValue(MDL_SAMPLINGRATEX,aux))
        {
            aux = -header->xOrigin/aux;
            MD.setValue(MDL_ORIGINX, aux);
        }
        if(MDMainHeader.getValue(MDL_SAMPLINGRATEY,aux))
        {
            aux = -header->yOrigin/aux;
            MD.setValue(MDL_ORIGINY, aux);
        }

        if(MDMainHeader.getValue(MDL_SAMPLINGRATEZ,aux))
        {
            aux = -header->zOrigin/aux;
            MD.setValue(MDL_ORIGINZ, aux);
        }
        MD.setValue(MDL_ANGLEROT, zeroD);
        MD.setValue(MDL_ANGLETILT,zeroD);
        MD.setValue(MDL_ANGLEPSI, zeroD);
        MD.setValue(MDL_WEIGHT,   oneD);
        MD.setValue(MDL_FLIP,     falseb);
    }

    //#define DEBUG
#ifdef DEBUG
    MDMainHeader.write("/dev/stderr");
    MD.write("/dev/stderr");
#endif

    freeMemory(header, sizeof(MRChead));
    readData(fimg, img_select, datatype, 0);

    if ( !mmapOn )
    	fclose(fimg);

    return(0);
}

int writeMRC(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE)
{
    /*
        if ( transform != NoTransform )
            img_convert_fourier(p, CentHerm);
    */
    MRChead*        header = (MRChead *) askMemory(sizeof(MRChead));

    // Map the parameters
    strncpy(header->map, "MAP ", 4);
    // FIXME TO BE DONE WITH rwCCP4!!
    //set_CCP4_machine_stamp(header->machst);
    int Xdim = XSIZE(data);
    int Ydim = YSIZE(data);
    int Zdim = ZSIZE(data);
    int Ndim = NSIZE(data);
    unsigned long   imgStart=0;
    unsigned long   imgEnd =Ndim;
    if (img_select != -1)
    {
        imgStart=img_select;
        imgEnd=img_select+1;
    }
    if (mode = WRITE_APPEND)
    {
        imgStart=0;
        imgEnd=1;
    }
    header->nx = Xdim;
    header->ny = Ydim;
    header->nz = Zdim;
    if ( transform == CentHerm )
        header->nx = Xdim/2 + 1;        // If a transform, physical storage is nx/2 + 1

    // Convert T to datatype
    if ( typeid(T) == typeid(double) ||
         typeid(T) == typeid(float) ||
         typeid(T) == typeid(int) )
        header->mode = 2;
    else if ( typeid(T) == typeid(unsigned char) ||
              typeid(T) == typeid(signed char) )
        header->mode = 0;
    else if ( typeid(T) == typeid(std::complex<float>) ||
              typeid(T) == typeid(std::complex<double>) )
        header->mode = 4;
    else
        REPORT_ERROR(1,"ERROR write MRC image: invalid typeid(T)");

    //Set this to zero till we decide if we want to update it
    header->mx = 0;//(int) (ua/ux + 0.5);
    header->my = 0;//(int) (ub/uy + 0.5);
    header->mz = 0;//(int) (uc/uz + 0.5);
    header->mapc = 1;
    header->mapr = 2;
    header->maps = 3;
    double aux,aux2;

    header->a = 0.;// ua;
    header->b = 0.;// ub;
    header->c = 0.;// uc;

    if (MDMainHeader.firstObject() != MetaData::NO_OBJECTS_STORED)
    {
        if(MDMainHeader.getValue(MDL_MIN, aux))
            header->amin  = (float)aux;
        if(MDMainHeader.getValue(MDL_MAX, aux))
            header->amax  = (float)aux;
        if(MDMainHeader.getValue(MDL_AVG, aux))
            header->amean = (float)aux;
        if(MDMainHeader.getValue(MDL_STDDEV, aux))
            header->arms  = (float)aux;
        if(MDMainHeader.getValue(MDL_ORIGINX, aux))
            header->nxStart = (int)(aux-0.5);
        if (MDMainHeader.getValue(MDL_SAMPLINGRATEX,aux2))//header is init to zero
            header->xOrigin = (float)(aux*aux2);
        if (MDMainHeader.getValue(MDL_ORIGINY, aux))
            header->nyStart = (int)(aux-0.5);
        if (MDMainHeader.getValue(MDL_SAMPLINGRATEY,aux2))//header is init to zero
            header->yOrigin = (float)(aux*aux2);
        if (MDMainHeader.getValue(MDL_ORIGINZ, aux))
            header->nzStart = (int)(aux-0.5);
        if (MDMainHeader.getValue(MDL_SAMPLINGRATEZ,aux2))//header is init to zero
            header->zOrigin = (float)(aux*aux2);
    }

    header->nsymbt = 0;
    header->nlabl = 10; // or zero?
    //strncpy(header->labels, p->label.c_str(), 799);

    offset = MRCSIZE + header->nsymbt;
    size_t datasize, datasize_n;
    datasize_n = Xdim*Ydim*Zdim;
    if (isComplexT())
        datasize = datasize_n * gettypesize(ComplexFloat);
    else
        datasize = datasize_n * gettypesize(Float);

//#define DEBUG
#ifdef DEBUG

    printf("DEBUG rwMRC: Offset = %ld,  Datasize_n = %ld\n", offset, datasize_n);
#endif

    // For multi-image files
    if (mode == WRITE_APPEND && isStack)
    {
        header->nz = replaceNsize +1;
    }
    //else header-> is correct

    //locking
    struct flock fl;
    int fd;

    fl.l_type   = F_WRLCK;  /* F_RDLCK, F_WRLCK, F_UNLCK    */
    fl.l_whence = SEEK_SET; /* SEEK_SET, SEEK_CUR, SEEK_END */
    fl.l_start  = 0;        /* Offset from l_whence         */
    fl.l_len    = 0;        /* length, 0 = to EOF           */
    fl.l_pid    = getpid(); /* our PID                      */

    FILE        *fimg;
    if (mode==WRITE_OVERWRITE || (!_exists && mode==WRITE_APPEND))//open in overwrite mode
    {
        if ( ( fimg = fopen(filename.c_str(), "w") ) == NULL )
            REPORT_ERROR(1,(std::string)"Cannot create file " + filename);
    }
    else //open in append mode
    {
        if ( ( fimg = fopen(filename.c_str(), "r+") ) == NULL )
            REPORT_ERROR(1,(std::string)"Cannot create file " + filename);
    }

    //BLOCK
    fl.l_type   = F_WRLCK;
    fcntl(fileno(fimg), F_SETLKW, &fl); /* locked */

    // Write header
    if(mode==WRITE_OVERWRITE || mode==WRITE_APPEND)
        fwrite( header, MRCSIZE, 1, fimg );
    freeMemory(header, sizeof(MRChead) );


    /*
    // Write 3D map
    if ( typeid(T) == typeid(double) ||
            typeid(T) == typeid(float) ||
            typeid(T) == typeid(int) )
{
        writePageAsDatatype(fimg, Float, datasize_n);
}
    else if ( typeid(T) == typeid(unsigned char) ||
              typeid(T) == typeid(signed char) )
        writePageAsDatatype(fimg, SChar, datasize_n);
    else if ( typeid(T) == typeid(std::complex<float>) ||
              typeid(T) == typeid(std::complex<double>) )
        writePageAsDatatype(fimg, ComplexFloat, datasize_n);
    else
        REPORT_ERROR(1,"ERROR write MRC image: invalid typeid(T)");
    */
    //write only once, ignore select_img
    char* fdata = (char *) askMemory(datasize);
    //think about writing in several chucks

    if ( NSIZE(data) == 1 && mode==WRITE_OVERWRITE)
    {
        if (isComplexT())
            castPage2Datatype(MULTIDIM_ARRAY(data), fdata, ComplexFloat, datasize_n);
        else
            castPage2Datatype(MULTIDIM_ARRAY(data), fdata, Float, datasize_n);
        fwrite( fdata, datasize, 1, fimg );
    }
    else
    {
        if(mode==WRITE_APPEND)
            fseek( fimg, 0, SEEK_END);
        else if(mode==WRITE_REPLACE)
            fseek( fimg,offset + (datasize)*img_select, SEEK_SET);
        for ( size_t i =imgStart; i<imgEnd; i++ )
        {
            if (isComplexT())
                castPage2Datatype(MULTIDIM_ARRAY(data) + i*datasize_n, fdata, ComplexFloat, datasize_n);
            else
                castPage2Datatype(MULTIDIM_ARRAY(data) + i*datasize_n, fdata, Float, datasize_n);
            fwrite( fdata, datasize, 1, fimg );
        }
    }
    //I guess I do not need to unlock since we are going to close the file
    fclose(fimg);

    freeMemory(fdata, datasize);

    return(0);
}
#endif
