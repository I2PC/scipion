/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es  (2007)
 *             Joaquin Oton            joton@cnb.csic.es (2010)
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
#ifndef IMAGE_CONVERT_H_
#define IMAGE_CONVERT_H_

#include "xmipp_program.h"

typedef enum
{
    MD2MD,
    MD2VOL,
    VOL2MD
} ImageConv;

class ProgConvImg: public XmippMetadataProgram
{
private:
    std::string  type;       // Type of output conversion
    std::string  depth;
    ImageGeneric imIn, *imOut;
    MDRow        row;
    ImageConv    convMode;
    CastWriteMode  castMode;
    size_t        k;
    int         writeMode;
    bool        appendToStack;
    bool        swap;

protected:
    void defineParams();
    void readParams();
    void preProcess();
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);
    void finishProcessing();
    void show();
};//class ProgConvImg

#endif /* IMAGE_CONVERT_H_ */
