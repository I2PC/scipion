/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include <data/args.h>
#include <data/geometry.h>
#include <data/xmipp_program.h>

class FlipParameters: public XmippMetadataProgram
{
public:
    bool flipX;
    bool flipY;
    bool flipZ;
    ImageGeneric img, imgOut;

    void defineParams()
    {
        each_image_produces_an_output = true;
        XmippMetadataProgram::defineParams();
        addUsageLine("Flip/Mirror around one of the main axis an image/volume");
        addExampleLine("mirror around x axis:", false);
        addExampleLine("xmipp_transform_mirror   -i 128.mrc -o 128x.spi --flipX");
        addParamsLine("--flipX     : flip around X axis");
        addParamsLine("or --flipY     : flip around Y axis");
        addParamsLine("or --flipZ     : flip around Z axis");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        flipX = checkParam("--flipX");
        flipY = checkParam("--flipY");
        flipZ = checkParam("--flipZ");
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
    {
        img.read(fnImg);
        if (flipX)
            img().selfReverseX();
        if (flipY)
            img().selfReverseY();
        if (flipZ)
            img().selfReverseZ();
        img.write(fnImgOut);
    }

};


int main(int argc, char **argv)
{
    FlipParameters program;
    program.read(argc, argv);
    return program.tryRun();
}

