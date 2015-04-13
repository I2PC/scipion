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

#include <data/xmipp_program.h>
#include <data/transformations.h>
#include <data/matrix2d.h>


class ProgAngularRotate: public XmippProgram
{
protected:
    FileName   fnIn, fnOut;
    MetaData   mdIn, mdOut;
    Matrix2D<double> R, T, B;
    bool       writeMatrix;

    void defineParams()
    {
        addUsageLine("Calculate the Euler angles after applying a rotation to a geometry set from a metadata.");
        addKeywords("angular, geometry, rotate");
        addSeeAlsoLine("transform_geometry");

        addParamsLine("  -i <metadata>   : Input metadata file");
        addParamsLine(" [-o <metadata>]  : Output metadata file, if not provided result will overwrite input file");
        addParamsLine("    --ang <angle>                  : Rotate in plane");
        addParamsLine(" or --euler <rot> <tilt> <psi>     : Rotate with these Euler angles");
        addParamsLine(" or --alignZ <x> <y> <z>           : Align (x,y,z) with Z axis");
        addParamsLine(" or --axis <ang> <x=0> <y=0> <z=1> : Rotate <ang> degrees around (x,y,z)");
        addParamsLine("[--write_matrix]                    : Print transformation matrix to screen");

        addExampleLine("Modify the geometry of a selfile applying a rotation of 0.1 degrees around (x,y,z)=(1,1,1) direction:", false);
        addExampleLine("xmipp_angular_rotate -i selfile.sel -o newselfile.sel --axis 0.1 1 1 1");

    }

    void readParams()
    {
        fnIn = getParam("-i");
        fnOut = (checkParam("-o"))? getParam("-o") : fnIn;

        mdIn.read(fnIn);
        int dim;

        if (mdIn.containsLabel(MDL_ANGLE_PSI) || mdIn.containsLabel(MDL_ANGLE_TILT))
            dim = 3;
        else
            dim = 2;

        R.initIdentity(dim + 1);
        B.initIdentity(dim + 1);
        Matrix1D<double> xyz(3);

        if (checkParam("--ang"))
            Euler_angles2matrix(0, 0,getDoubleParam("--ang"), R, true);
        else if (checkParam("--axis"))
        {
            double ang = getDoubleParam("--axis", 0);
            XX(xyz) = getDoubleParam("--axis", 1); //axis x
            YY(xyz) = getDoubleParam("--axis", 2); //y
            ZZ(xyz) = getDoubleParam("--axis", 3);//z
            rotation3DMatrix(ang, xyz, R, true);
        }
        else if (checkParam("--euler"))
        {
            XX(xyz) = getDoubleParam("--euler", 0); //rot
            YY(xyz) = getDoubleParam("--euler", 1); //tilt
            ZZ(xyz) = getDoubleParam("--euler", 2);//psi
            Euler_angles2matrix(XX(xyz), YY(xyz), ZZ(xyz), R, true);
        }
        else if (checkParam("--alignZ"))
        {
            XX(xyz) = getDoubleParam("--alignZ", 0); //rot
            YY(xyz) = getDoubleParam("--alignZ", 1); //tilt
            ZZ(xyz) = getDoubleParam("--alignZ", 2);//psi
            alignWithZ(xyz, R);
        }

        writeMatrix = checkParam("--write_matrix");
    }

    void run()
    {
        MDRow input;
        bool containsScale=mdIn.containsLabel(MDL_SCALE);
        FOR_ALL_OBJECTS_IN_METADATA(mdIn)
        {
            mdIn.getRow(input, __iter.objId);//Get geometric transformation for image
            geo2TransformationMatrix(input, B);

            T = R * B;

            transformationMatrix2Geo(T, input);
            if (!containsScale)
            	input.setValue(MDL_SCALE,1.0);
            mdOut.addRow(input);

            if (writeMatrix)
                std::cout << T << std::endl;
        }

        mdOut.write(fnOut);
    }
};

int main(int argc, char **argv)
{
    ProgAngularRotate program;
    program.read(argc, argv);
    return program.tryRun();
} //main
