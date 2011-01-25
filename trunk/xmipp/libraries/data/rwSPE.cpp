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

#include "image_base.h"

///@defgroup SPE Princeton Instruments File Format
///@ingroup ImageFormats

// I/O prototypes
/** SPE Reader
  * @ingroup SPE
*/
int ImageBase::readSPE(int img_select,bool isStack)
{

    int _xDim,_yDim,_zDim;
    unsigned long int _nDim;

    short int aux;
    fseek(fimg,42,SEEK_SET);
    xmippFREAD(&aux, sizeof(short int), 1, fimg, swap );
    _xDim = aux;
    fseek(fimg,656,SEEK_SET);
    xmippFREAD(&aux, sizeof(short int), 1, fimg, swap );
    _yDim = aux;

    _zDim = (int) 1;
    _nDim = (int) 1;

    // Map the parameters
    setDimensions(_xDim, _yDim, _zDim, _nDim);

    unsigned long   imgStart=0;
    unsigned long   imgEnd =_nDim;
    if (img_select != -1)
    {
        imgStart=img_select;
        imgEnd=img_select+1;
    }

    DataType datatype = UShort;

    MDMainHeader.setValue(MDL_SAMPLINGRATEX,(double) -1);
    MDMainHeader.setValue(MDL_SAMPLINGRATEY,(double) -1);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    if( dataflag == -2 )
        return 0;

    MD.clear();
    MD.resize(imgEnd - imgStart);
    for ( i = imgStart; i < imgEnd; ++i )
    {
        MD[i-imgStart].setValue(MDL_ORIGINX, zeroD);
        MD[i-imgStart].setValue(MDL_ORIGINY, zeroD);
        MD[i-imgStart].setValue(MDL_ORIGINZ,  zeroD);
        MD[i-imgStart].setValue(MDL_ANGLEROT, zeroD);
        MD[i-imgStart].setValue(MDL_ANGLETILT,zeroD);
        MD[i-imgStart].setValue(MDL_ANGLEPSI, zeroD);
        MD[i-imgStart].setValue(MDL_WEIGHT,   oneD);
        MD[i-imgStart].setValue(MDL_FLIP,     falseb);
        MD[i-imgStart].setValue(MDL_SCALE,    oneD);
    }

    offset = 4100;
    size_t pad = 0;

    readData(fimg, img_select, datatype, pad);

    return(0);
}
