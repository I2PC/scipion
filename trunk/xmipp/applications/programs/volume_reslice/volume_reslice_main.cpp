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

#include <data/program.h>
#include <data/transformations.h>

class ProgResliceVol : public XmippProgram
{
private:

    FileName fnImgIn, fnImgOut;
    ImageGeneric imgIn;   // Input Volume
    ImageGeneric imgOut;  // Output Volume
    String      face;
    Matrix2D<double> R;


    void defineParams()
    {
        addUsageLine("Reslice a volume aligning any X or Y direction with Z axis.");
        addUsageLine("The output volume also has the dimensions rotated, so all the information in the volume is kept.");
        addSeeAlsoLine("transform_geometry");

        addParamsLine("  -i <input_file>    : Input volume file");
        addParamsLine("  [-o <output_file>] : Output volume file");
        addParamsLine("  --face <new_front_face>  : Face to be brought to front.");
        addParamsLine("     where <new_front_face>");
        addParamsLine("                 top    : Align -Y axis to Z axis, rotating 90 degrees around X axis");
        addParamsLine("                 left   : Align -X axis to Z axis, rotating -90 degrees around Y axis");
        addParamsLine("                 bottom : Align Y axis to Z axis, rotating -90 degrees around X axis");
        addParamsLine("                 right  : Align X axis to Z axis, rotating 90 degrees around Y axis");
    }

    void readParams()
    {
        fnImgIn = getParam("-i");
        fnImgOut = getParam("-o");
        face = getParam("--face");
    }


    void run()
    {
        imgIn.readMapped(fnImgIn);
        imgIn().setXmippOrigin();

        int xDim, yDim, zDim, XdimOut, yDimOut, zDimOut;
        imgIn.getDimensions(xDim, yDim, zDim);

        Matrix1D<double> xyz(3);
        xyz.initZeros();

        if (face == "top" || face == "bottom")
        {
            XdimOut = xDim;
            yDimOut = zDim;
            zDimOut = yDim;

            YY(xyz) = (face == "bottom") ? 1 : -1;
        }
        else if (face == "left" || face == "right")
        {
            XdimOut = zDim;
            yDimOut = yDim;
            zDimOut = xDim;

            XX(xyz) = (face == "right") ? 1 : -1;
        }

        alignWithZ(xyz, R);


        // Create output file
        imgOut.setDatatype(imgIn.getDatatype());
        imgOut.mapFile2Write(XdimOut, yDimOut, zDimOut, fnImgOut, fnImgIn == fnImgOut);
        imgOut().setXmippOrigin();

        applyGeometry(NEAREST, imgOut(), imgIn(), R, IS_NOT_INV, false, 0.);

        imgOut.write(fnImgOut);

    }


}
;

int main(int argc, char **argv)
{
    ProgResliceVol program;
    program.read(argc, argv);
    return program.tryRun();
}
