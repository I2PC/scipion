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

#define IMAGICSIZE 1024 // Size of the IMAGIC header for each image

///@defgroup Imagic Imagic File format
///@ingroup ImageFormats

/** Imagic Header
  * @ingroup Imagic
*/
struct IMAGIChead
{             // file header for IMAGIC data
    int imn;          //  0      image location number (1,2,...)
    int ifn;          //  1      # images following, only of importance in the first location
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
/** Imagic reader
  * @ingroup Imagic
*/
int  ImageBase::readIMAGIC(size_t select_img)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG readIMAGIC: Reading Imagic file\n");
#endif

    IMAGIChead* header = new IMAGIChead;

    if ( fread( header, IMAGICSIZE, 1, fhed ) < 1 )
        REPORT_ERROR(ERR_IO_NOREAD,(String)"readIMAGIC: header file of " + filename + " cannot be read");

    // Determine byte order and swap bytes if from little-endian machine
    if ( (swap = (( abs(header->nyear) > SWAPTRIG ) || ( header->ixlp > SWAPTRIG ))) )
        swapPage((char *) header, IMAGICSIZE - 916, DT_Float); // IMAGICSIZE - 916 is to exclude labels from swapping

    DataType datatype;

    if ( strstr(header->type,"PACK") )
        datatype = DT_UChar;
    else if ( strstr(header->type,"INTG") )
        datatype = DT_Short;
    else if ( strstr(header->type,"REAL") )
        datatype = DT_Float;
    else if ( strstr(header->type,"RECO") )
    {
        datatype = DT_CFloat; // Complex data
        transform = NoTransform;
    }
    else if ( strstr(header->type,"COMP") )
    {
        datatype = DT_CFloat; // Complex transform data
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

    offset = 0;   // separate header file

    MDMainHeader.setValue(MDL_MIN,(double)header->densmin);
    MDMainHeader.setValue(MDL_MAX,(double)header->densmax);
    MDMainHeader.setValue(MDL_AVG,(double)header->avdens);
    MDMainHeader.setValue(MDL_STDDEV,(double)header->sigma);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_X,(double)1.);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_Y,(double)1.);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_Z,(double)1.);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    int _xDim,_yDim,_zDim;
    size_t _nDim;
    _xDim = (int) header->iylp;
    _yDim = (int) header->ixlp;
    _zDim = (int) 1;
    _nDim = (size_t) header->ifn + 1 ;

    if ( select_img > _nDim )
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, formatString("readImagic: Image number %lu exceeds stack size %lu", select_img, _nDim));

    if( select_img != ALL_IMAGES )
        _nDim = 1;

    replaceNsize = _nDim;
    setDimensions(_xDim, _yDim, _zDim, _nDim );

    if (dataMode == HEADER || (dataMode == _HEADER_ALL && _nDim > 1)) // Stop reading if not necessary
    {
        delete header;
        return 0;
    }

    // Get the header information
    fseek( fhed, IMG_INDEX(select_img) * IMAGICSIZE, SEEK_SET );

    MD.clear();
    MD.resize(_nDim,MDL::emptyHeader);
    double daux=1.;
    for ( size_t i = 0; i < _nDim; ++i )
    {
        if ( fread( header, IMAGICSIZE, 1, fhed ) < 1 )
            return(-2);
        {
            if ( swap )
                swapPage((char *) header, IMAGICSIZE - 916, DT_Float);

            if (dataMode == _HEADER_ALL || dataMode == _DATA_ALL)
            {
                MD[i].setValue(MDL_SHIFT_X,  (double)-1. * header->ixold);
                MD[i].setValue(MDL_SHIFT_Y,  (double)-1. * header->iyold);
                MD[i].setValue(MDL_SHIFT_Z,  0.);
                MD[i].setValue(MDL_ANGLE_ROT, (double)-1. * header->euler_alpha);
                MD[i].setValue(MDL_ANGLE_TILT,(double)-1. * header->euler_beta);
                MD[i].setValue(MDL_ANGLE_PSI, (double)-1. * header->euler_gamma);
                MD[i].setValue(MDL_WEIGHT,   1.);
                MD[i].setValue(MDL_SCALE, daux);
            }
        }
    }
    delete header;

    if (dataMode < DATA)   // Don't read the individual header and the data if not necessary
        return 0;

    size_t pad = 0;
    readData(fimg, select_img, datatype, pad );

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
/** Imagic Writer
  * @ingroup Imagic
*/
int  ImageBase::writeIMAGIC(size_t select_img, int mode, const String &bitDepth, bool adjust)
{
#undef DEBUG
    //#define DEBUG
#ifdef DEBUG
    printf("DEBUG writeIMAGIC: Reading Imagic file\n");
#endif

    IMAGIChead header;

    // Cast T to datatype without convert data
    DataType wDType, myTypeID = myT();
    CastWriteMode castMode = CW_CAST;

    if (bitDepth == "")
    {
        switch(myTypeID)
        {
        case DT_Double:
        case DT_Float:
        case DT_Int:
        case DT_UInt:
            wDType = DT_Float;
            strncpy(header.type,"REAL", sizeof(header.type));
            break;
        case DT_UShort:
            castMode = CW_CONVERT;
        case DT_Short:
            wDType = DT_Short;
            strncpy(header.type,"INTG", sizeof(header.type));
            break;
        case DT_SChar:
            castMode = CW_CONVERT;
        case DT_UChar:
            wDType = DT_UChar;
            strncpy(header.type,"PACK", sizeof(header.type));
            break;
        case DT_CFloat:
        case DT_CDouble:
            wDType = DT_CFloat;
            strncpy(header.type,"COMP", sizeof(header.type));
            break;
        default:
            wDType = DT_Unknown;
            REPORT_ERROR(ERR_TYPE_INCORRECT, "ERROR: Unsupported data type by IMAGIC format.");
        }
    }
    else //Convert to other data type
    {
        // Default Value
        wDType = (bitDepth == "default") ? DT_Float : datatypeRAW(bitDepth);

        switch (wDType)
        {
        case DT_UChar:
            strncpy(header.type,"PACK", sizeof(header.type));
            castMode = (adjust)? CW_ADJUST : CW_CONVERT;
            break;
        case DT_Short:
            strncpy(header.type,"INTG", sizeof(header.type));
            castMode = (adjust)? CW_ADJUST : CW_CONVERT;
            break;
        case DT_Float:
            (strncpy)(header.type,"REAL", sizeof(header.type));
            break;
        case DT_CFloat:
            strncpy(header.type,"COMP", sizeof(header.type));
            break;
        default:
            REPORT_ERROR(ERR_TYPE_INCORRECT,"ERROR: incorrect IMAGIC bits depth value.");
        }
    }

    if (mmapOnWrite)
    {
        MDMainHeader.setValue(MDL_DATATYPE,(int) wDType);
        if (!checkMmapT(wDType))
        {
            if (dataMode < DATA && castMode == CW_CAST) // This means ImageGeneric wants to know which DataType must use in mapFile2Write
                return 0;
            else //Mapping is an extra. When not available, go on and do not report an error.
            {
                /* In this case we cannot map the file because required and feasible datatypes are
                 * not compatible. Then we denote to MapFile2Write the same incoming datatype to
                 * keep using this Image object as usual, without mapping on write.
                 */
                mmapOnWrite = false;
                dataMode = DATA;
                MDMainHeader.setValue(MDL_DATATYPE,(int) myTypeID);

                // In case Image size great then, at least, map the multidimarray
                if (mdaBase->nzyxdim*gettypesize(myTypeID) > tiff_map_min_size)
                    mdaBase->setMmap(true);

                // Allocate memory for image data (Assume xdim, ydim, zdim and ndim are already set
                //if memory already allocated use it (no resize allowed)
                mdaBase->coreAllocateReuse();

                return 0;
            }
        }
        else
            dataMode = DATA;
    }

    size_t Xdim, Ydim, Zdim, Ndim;
    getDimensions(Xdim, Ydim, Zdim, Ndim);

    if (Zdim > 1)
        REPORT_ERROR(ERR_MULTIDIM_DIM, "writeIMAGIC: Imagic format does not support volumes.");

    size_t datasize, datasize_n;
    datasize_n = (size_t)Xdim*Ydim*Zdim;
    datasize = datasize_n * gettypesize(wDType);

    // fill in the file header
    header.nhfr = 1;
    header.npix2 = Xdim*Ydim;
    header.npixel = header.npix2;
    header.iylp = Xdim;
    header.ixlp = Ydim;

    time_t timer;
    time ( &timer );
    tm* t = localtime(&timer);

    header.ndate = t->tm_mday;
    header.nmonth = t->tm_mon + 1;
    header.nyear = t->tm_year;
    header.nhour = t->tm_hour;
    header.nminut = t->tm_min;
    header.nsec = t->tm_sec;

    double aux;

    if (!MDMainHeader.empty())
    {
#define SET_MAIN_HEADER_VALUE(field, label)  MDMainHeader.getValueOrDefault(label, aux, 0.); header.field = (float)aux
        SET_MAIN_HEADER_VALUE(densmin, MDL_MIN);
        SET_MAIN_HEADER_VALUE(densmax, MDL_MAX);
        SET_MAIN_HEADER_VALUE(avdens, MDL_AVG);
        SET_MAIN_HEADER_VALUE(sigma, MDL_STDDEV);
        header.varian = header.sigma*header.sigma;
    }

    memcpy(header.lastpr, "Xmipp", 5);
    memcpy(header.name, filename.c_str(), 80);

    size_t  imgStart = IMG_INDEX(select_img);

    header.ifn = replaceNsize - 1 ;
    header.imn = 1;

    if ( mode == WRITE_APPEND )
    {
        imgStart = replaceNsize;
        header.ifn = replaceNsize + Ndim - 1 ;
    }
    else if( mode == WRITE_REPLACE && imgStart + Ndim > replaceNsize)
        header.ifn = imgStart + Ndim - 1;
    else if (Ndim > replaceNsize)
        header.ifn = Ndim - 1;

    /*
     * BLOCK HEADER IF NEEDED
     */
    FileLock flockHead, flockImg;
    flockHead.lock(fhed);
    flockImg.lock(fimg);

    if (replaceNsize == 0) // Header written first time
    {
        if ( swapWrite )
        {
            IMAGIChead headTemp = header;
            swapPage((char *) &headTemp, IMAGICSIZE - 916, DT_Float);
            fwrite( &headTemp, IMAGICSIZE, 1, fhed );
        }
        fwrite( &header, IMAGICSIZE, 1, fhed );
    }
    else if( header.ifn + 1 > (int)replaceNsize && imgStart > 0 ) // Update number of images when needed
    {
        fseek( fhed, sizeof(int), SEEK_SET);
        if ( swapWrite )
        {
            int ifnswp = header.ifn;
            swapPage((char *) &ifnswp, SIZEOF_INT, DT_Int);
            fwrite(&(ifnswp),SIZEOF_INT,1,fhed);
        }
        else
            fwrite(&(header.ifn),SIZEOF_INT,1,fhed);
    }

    // Jump to the selected imgStart position
    fseek(fimg, datasize   * imgStart, SEEK_SET);
    fseek(fhed, IMAGICSIZE * imgStart, SEEK_SET);

    std::vector<MDRow>::iterator it = MD.begin();

    for (size_t i = 0; i < Ndim; ++i, ++it)
    {
        header.iyold=header.ixold=0;
        header.euler_alpha=header.euler_beta=header.euler_gamma=0.;

        // Write the individual image header
        if (it != MD.end() && (dataMode == _HEADER_ALL || dataMode == _DATA_ALL))
        {
#define SET_HEADER_VALUEInt(field, label)  it->getValueOrDefault((label), (aux), 0); header.field = -(int)(aux)
#define SET_HEADER_VALUEDouble(field, label)  it->getValueOrDefault((label), (aux), 0.); header.field = -(float)(aux)

            SET_HEADER_VALUEInt(ixold, MDL_SHIFT_X);
            SET_HEADER_VALUEInt(iyold, MDL_SHIFT_Y);
            SET_HEADER_VALUEDouble(euler_alpha, MDL_ANGLE_ROT);
            SET_HEADER_VALUEDouble(euler_beta, MDL_ANGLE_TILT);
            SET_HEADER_VALUEDouble(euler_gamma, MDL_ANGLE_PSI);
        }
        // Update index number of image
        header.imn = imgStart + i + 1;

        if ( swapWrite )
            swapPage((char *) &header, IMAGICSIZE - 916, DT_Float);
        fwrite( &header, IMAGICSIZE, 1, fhed );

        if (dataMode >= DATA)
        {
            if (mmapOnWrite && Ndim == 1) // Can map one image at a time only
            {
                mappedOffset = ftell(fimg);
                mappedSize = mappedOffset + datasize;
                fseek(fimg, datasize-1, SEEK_CUR);
                fputc(0, fimg);
            }
            else
                writeData(fimg, i*datasize_n, wDType, datasize_n, castMode);
        }
        else
            fseek(fimg, datasize, SEEK_CUR);
    }

    //Unlock
    flockHead.unlock();
    flockImg.unlock();

    if (mmapOnWrite)
        mmapFile();

    return(0);
}
