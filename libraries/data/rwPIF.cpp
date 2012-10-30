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

#define PIFHEADERSIZE 512 // size of EM file header

int  ImageBase::readPIF(size_t select_img)
{
    PIFMainHeader  mainHeader;
    if ( fread( &mainHeader, PIFHEADERSIZE, 1, fimg ) != 1 )
        REPORT_ERROR(ERR_IO_NOREAD, formatString("rwPIF: cannot read Spider main header from file %s"
                     ". Error message: %s", filename.c_str() ,strerror(errno)));

    // Check Machine endianess
    bool isLE = IsLittleEndian();

    // Determine byte order and swap bytes if from different-endian machine
    swap = (isLE == mainHeader.endianNess);

    if (swap)
        swapPage((char *) &mainHeader, PIFHEADERSIZE, DT_Int);


    // Setting image dimensions
    ArrayDim aDim;

    aDim.ndim = mainHeader.numImages;
    aDim.zdim = mainHeader.nz;
    aDim.ydim = mainHeader.ny;
    aDim.xdim = mainHeader.nx;

    replaceNsize = aDim.ndim;


    offset = (size_t) header->labbyt;
    DataType datatype  = DT_Float;

    MDMainHeader.setValue(MDL_MIN,(double)header->fmin);
    MDMainHeader.setValue(MDL_MAX,(double)header->fmax);
    MDMainHeader.setValue(MDL_AVG,(double)header->av);
    MDMainHeader.setValue(MDL_STDDEV,(double)header->sig);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_X,(double)header->scale);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_Y,(double)header->scale);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_Z,(double)header->scale);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    bool isStack = ( header->istack > 0 );
    int _xDim,_yDim,_zDim;
    size_t _nDim, _nDimSet;
    _xDim = (int) header->nsam;
    _yDim = (int) header->nrow;
    _zDim = (int) header->nslice;
    _nDim = (isStack)? header->maxim : 1;

    if (_xDim < 1 || _yDim < 1 || _zDim < 1 || _nDim < 1)
        REPORT_ERROR(ERR_IO_NOTFILE,formatString("Invalid Spider file:  %s", filename.c_str()));

    replaceNsize = _nDim;

    /************
     * BELLOW HERE DO NOT USE HEADER BUT LOCAL VARIABLES
     */

    // Map the parameters, REad the whole object (-1) or a slide
    // Only handle stacks of images not of volumes
    if(!isStack)
        _nDimSet = 1;
    else
        _nDimSet = (select_img == ALL_IMAGES) ? _nDim : 1;

    setDimensions(_xDim, _yDim, _zDim, _nDimSet);

    //image is in stack? and set right initial and final image
    size_t header_size = offset;

    if ( isStack)
    {
        if ( select_img > _nDim )
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, formatString("readSpider (%s): Image number %lu exceeds stack size %lu" ,filename.c_str(),select_img, _nDim));
        offset += offset;
    }

    if (dataMode == HEADER || (dataMode == _HEADER_ALL && _nDimSet > 1)) // Stop reading if not necessary
    {
        delete header;
        return 0;
    }

    size_t datasize_n  = _xDim*_yDim*_zDim;
    size_t image_size  = header_size + datasize_n*sizeof(float);
    size_t pad         = (size_t) header->labbyt;
    size_t   imgStart = IMG_INDEX(select_img);
    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : _nDim;
    size_t   img_seek = header_size + imgStart * image_size;
    char*   hend;

    MD.clear();
    MD.resize(imgEnd - imgStart,MDL::emptyHeader);
    double daux;

    //std::cerr << formatString("DEBUG_JM: header_size: %10lu, datasize_n: %10lu, image_size: %10lu, imgStart: %10lu, img_seek: %10lu",
    //    header_size, datasize_n, image_size, imgStart, img_seek) <<std::endl;

    for (size_t n = 0, i = imgStart; i < imgEnd; ++i, ++n, img_seek += image_size )
    {
        if (fseek( fimg, img_seek, SEEK_SET ) != 0)//fseek return 0 on success
            REPORT_ERROR(ERR_IO, formatString("rwSPIDER: error seeking %lu for read image %lu", img_seek, i));

        // std::cerr << formatString("DEBUG_JM: rwSPIDER: seeking %lu for read image %lu", img_seek, i) <<std::endl;

        if(isStack)
        {
            if ( fread( header, SPIDERSIZE, 1, fimg ) != 1 )
                REPORT_ERROR(ERR_IO_NOREAD, formatString("rwSPIDER: cannot read Spider image %lu header", i));
            if ( swap )
                swapPage((char *) header, SPIDERSIZE - 180, DT_Float);
        }
        if (dataMode == _HEADER_ALL || dataMode == _DATA_ALL)
        {
            daux = (double)header->xoff;
            MD[n].setValue(MDL_SHIFT_X, daux);
            daux = (double)header->yoff;
            MD[n].setValue(MDL_SHIFT_Y, daux);
            daux = (double)header->zoff;
            MD[n].setValue(MDL_SHIFT_Z, daux);
            daux = (double)header->phi;
            MD[n].setValue(MDL_ANGLE_ROT, daux);
            daux = (double)header->theta;
            MD[n].setValue(MDL_ANGLE_TILT, daux);
            daux = (double)header->gamma;
            MD[n].setValue(MDL_ANGLE_PSI, daux);
            daux = (double)header->weight;
            MD[n].setValue(MDL_WEIGHT, daux);
            bool baux = (header->flip == 1);
            MD[n].setValue(MDL_FLIP, baux);
            daux = (double) header->scale;
            if (daux==0.)
                daux=1.0;
            MD[n].setValue(MDL_SCALE, daux);
        }
    }
    delete header;

    if (dataMode < DATA)   // Don't read  data if not necessary but read the header
        return 0;

#ifdef DEBUG

    std::cerr<<"DEBUG readSPIDER: header_size = "<<header_size<<" image_size = "<<image_size<<std::endl;
    std::cerr<<"DEBUG readSPIDER: select_img= "<<select_img<<" n= "<<Ndim<<" pad = "<<pad<<std::endl;
#endif
    //offset should point to the begin of the data
    readData(fimg, select_img, datatype, pad );

    return(0);
}
