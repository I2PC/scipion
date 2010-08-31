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

#ifndef RWSPE_H_
#define RWSPE_H_

///@defgroup SPE Princeton Instruments CCD camera
///@ingroup ImageFormats

// I/O prototypes
/** SPE Reader
  * @ingroup SPE
*/
int readSPE(int img_select,bool isStack=false)
{

    int xDim,yDim,zDim, depth;
    unsigned long int nDim;

    FILE        *fimg;
    if ( ( fimg = fopen(filename.c_str(), "r") ) == NULL )
        REPORT_ERROR(ERR_IO_NOTOPEN,"readSPE: error opening image file.");

    short int aux;
    fseek(fimg,42,SEEK_SET);
    xmippFREAD(&aux, sizeof(short int), 1, fimg, swap );
    xDim = aux;
    fseek(fimg,656,SEEK_SET);
    xmippFREAD(&aux, sizeof(short int), 1, fimg, swap );
    yDim = aux;

    zDim = (int) 1;
    nDim = (int) 1;

    // Map the parameters
    data.setDimensions(xDim, yDim, zDim, nDim);

    unsigned long   imgStart=0;
    unsigned long   imgEnd =nDim;
    if (img_select != -1)
    {
        imgStart=img_select;
        imgEnd=img_select+1;
    }

    DataType datatype = UShort;

    MDMainHeader.removeObjects();
    MDMainHeader.setColumnFormat(false);
    MDMainHeader.addObject();
    MDMainHeader.setValue(MDL_SAMPLINGRATEX,(double) -1);
    MDMainHeader.setValue(MDL_SAMPLINGRATEY,(double) -1);
    MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

    if( dataflag == -2 )
    {
        fclose(fimg);
        return 0;
    }

    MD.removeObjects();
    for ( i=imgStart; i<imgEnd; i++ )
    {
        MD.addObject();
        MD.setValue(MDL_ORIGINX, zeroD);
        MD.setValue(MDL_ORIGINY, zeroD);
        MD.setValue(MDL_ORIGINZ,  zeroD);
        MD.setValue(MDL_ANGLEROT, zeroD);
        MD.setValue(MDL_ANGLETILT,zeroD);
        MD.setValue(MDL_ANGLEPSI, zeroD);
        MD.setValue(MDL_WEIGHT,   oneD);
        MD.setValue(MDL_FLIP,     falseb);
    }

    offset = 4100;
    size_t pad = 0;

    readData(fimg, img_select, datatype, pad);

    if ( !mmapOn )
        fclose(fimg);

    return(0);
}

/** SPE Writer
  * @ingroup SPE
*/
int writeSPE(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE)
{
    REPORT_ERROR(ERR_IMG_NOWRITE, "writeSPE is not implemented.");
    return(-1);
}
#endif /* RWSPE_H_ */
