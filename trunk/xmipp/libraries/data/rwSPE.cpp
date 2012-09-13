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

///@defgroup SPE Princeton Instruments File Format
///@ingroup ImageFormats

// I/O prototypes
/** SPE Reader
  * @ingroup SPE
*/
int ImageBase::readSPE(size_t select_img,bool isStack)
{

    int _xDim,_yDim,_zDim;
    size_t _nDim;

    short int aux;
    fseek(fimg,42,SEEK_SET);
    xmippFREAD(&aux, sizeof(short int), 1, fimg, swap );
    _xDim = aux;
    fseek(fimg,656,SEEK_SET);
    xmippFREAD(&aux, sizeof(short int), 1, fimg, swap );
    _yDim = aux;

    _zDim = 1;
    _nDim = 1;

    // Map the parameters
    setDimensions(_xDim, _yDim, _zDim, _nDim);

    offset = 4100;
    DataType datatype = DT_UShort;

    MDMainHeader.setValue(MDL_SAMPLINGRATE_X,(double) -1);
    MDMainHeader.setValue(MDL_SAMPLINGRATE_Y,(double) -1);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    if (dataMode==HEADER) // Stop reading if not necessary
        return 0;

    size_t   imgStart = IMG_INDEX(select_img);
    size_t   imgEnd = (select_img != ALL_IMAGES) ? imgStart + 1 : _nDim;

    MD.clear();
    MD.resize(imgEnd - imgStart,MDL::emptyHeader);
    if( dataMode < DATA )
        return 0;

    size_t pad = 0;

    readData(fimg, select_img, datatype, pad);

    return(0);
}
