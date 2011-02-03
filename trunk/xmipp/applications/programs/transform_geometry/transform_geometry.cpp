/***************************************************************************
 *
 * Authors:     Carlos Oscar (           coss@cnb.csic.es (2000)
 *              Roberto Marabini (added fourier option)
 *              Jose Miguel de la Rosa Trevin (fusion of shift, rotate and scale)
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

#include <data/image.h>
#include <data/args.h>
#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <data/fftw.h>
#include <data/program.h>
#include <data/transformations.h>
#include <data/matrix2d.h>


class ProgTransformGeometry: public XmippMetadataProgram
{
protected:
    int             zdim, ydim, xdim, splineDegree;
    bool            applyGeo, inverse, wrap, doResize;
    Matrix2D<double> R, T, S;
    Matrix1D<double>          shiftV, rotV;
    double scaleX, scaleY, scaleZ;

    void defineParams()
    {
        each_image_produces_an_output = true;
        XmippMetadataProgram::defineParams();
        //usage
        addUsageLine("Apply geometric transformations to images.");
        addUsageLine("This program allows you to shift (providing 2 or 3 translations),");
        addUsageLine("rotate (2 or 3 Euler angles) and scale a group of images/volumes.");
        addUsageLine("By default, if rotation angles and translations are stored in metada");
        addUsageLine("it will be read, also the original images will be preserved if possible.");
        //examples

        //params
        addParamsLine("== Scale ==");
        addParamsLine("[--scale <scale_type>]                :Perform scaling");
        addParamsLine("         where <scale_type>");
        addParamsLine("             factor <n=1>           : Scaling factor, 0.5 halves and 2 doubles");
        addParamsLine("                                    : This option doesn't cause a change on image dimensions.");
        addParamsLine("             dim <x> <y=x> <z=x>    : New x,y and z dimensions");
        addParamsLine("                                    : This cause a change on image dimensions.");
        addParamsLine(" alias -s;");

        addParamsLine("== Rotate ==");
        addParamsLine("[--rotate <rotation_type>]   :Perform rotation");
        addParamsLine("         where <rotation_type>");
        addParamsLine("             euler <rot> <tilt> <psi>     : Rotate with these Euler angles");
        addParamsLine("             alignZ <x> <y> <z>           : Align (x,y,z) with Z axis");
        addParamsLine("             axis <ang> <x=0> <y=0> <z=1> : Rotate <ang> degrees around (x,y,z)");

        addParamsLine("== Shift ==");
        addParamsLine("[--shift <x> <y=x> <z=x>]    : Shift by x, y and z");

        addParamsLine("== Other options ==");
        addParamsLine(" [--interp <interpolation_type=spline>] : Interpolation type to be used. ");
        addParamsLine("      where <interpolation_type>");
        addParamsLine("        spline          : Use spline interpolation");
        addParamsLine("        linear          : Use bilinear/trilinear interpolation");
        addParamsLine("        fourier <thr=1> : Use padding/windowing in Fourier Space (only for scale in 2D)");
        addParamsLine("[--inverse]                         : Apply inverse transformations");
        addParamsLine("[--apply_geo]                       : By default, the original images are preserved");
        addParamsLine("                                    : and the alignment information is stored in metadata");
        addParamsLine("[--dont_wrap]                       : By default, the image/volume is wrapped");
        addParamsLine("[--write_matrix]                    : Print transformation matrix to screen");

    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        applyGeo = checkParam("--apply_geo");
        inverse = checkParam("--inverse");
        wrap = ! checkParam("--dont_wrap");
        String degree = getParam("--interp");
        if (degree == "spline")
            splineDegree = BSPLINE3;
        else if (degree == "linear")
            splineDegree = LINEAR;
        else
            splineDegree = -1;//Fourier interpolation, only for scale
    }

    void preProcess()
    {
        unsigned long n;
        ImgSize(mdIn, xdim, ydim , zdim, n);
        R.initZeros(4, 4);
        S.initIdentity(4);
        T.initZeros(4, 4);
        shiftV.resizeNoCopy(3);
        rotV.resizeNoCopy(3);

        if (checkParam("--scale"))
        {
            if (String(getParam("--scale"))=="dim")
            {
                //Calculate scale factor from images sizes and given dimensions
                //this approach assumes that all images have equal size
                applyGeo = doResize = true;
                int x = getIntParam("--scale", 1);
                int y = String(getParam("--scale", 2)) == "x" ? x : getIntParam("--scale", 2);
                scaleX = (double)x / xdim;
                scaleY = (double)y / ydim;
                xdim = x;
                ydim = y;
                if (zdim > 1)
                {
                    int z = String(getParam("--scale", 3)) == "x" ? x : getIntParam("--scale", 3);
                    scaleZ = (double)z / zdim;
                    zdim = z;
                }
            }
            else //scale factor case
            {
                scaleX = scaleY = scaleZ = getDoubleParam("--factor");
            }
            //calculate scale matrix
            scaleMatrix(S, scaleX, scaleY, scaleZ);
        }

        if (checkParam("--rotate"))
        {
            //In any case read following 3 values
            Matrix1D<double> xyz;
            XX(xyz) = getDoubleParam("--rotate", 1); //rot
            YY(xyz) = getDoubleParam("--rotate", 2); //tilt
            ZZ(xyz) = getDoubleParam("--rotate", 3);//psi

            if (String(getParam("--rotate")) == "euler")
            {
                Euler_angles2matrix(XX(xyz), YY(xyz), ZZ(xyz), R, true);
            }
            else if (String(getParam("--rotate")) == "alignZ")
            {
                alignWithZ(xyz, R);
            }
            else
            {
                double ang = getDoubleParam("--rotate", 1);
                XX(xyz) = getDoubleParam("--rotate", 2); //rot
                YY(xyz) = getDoubleParam("--rotate", 3); //tilt
                ZZ(xyz) = getDoubleParam("--rotate", 4);//psi
                rotationMatrix(ang, xyz, R);
            }
        }
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        Image<double> img, imgOut;
        img.read(fnImg); //Doesn't apply geo to work with original images
        img().setXmippOrigin();
        //MultidimArray<double> out;

        double angRot = 0., angTilt = 0., angPsi = 0.;

        //try to values from metadata
        Matrix1D<double> xyz(3);
        Matrix2D<double> A(4, 4);

        //------Scale-----------------------
        if (checkParam("--scale"))
        {
            if (doResize)
                imgOut().resizeNoCopy(zdim, ydim, xdim);

            A = S;
            std::cerr << "scaling...matrix: " << S << std::endl;
        }

        //-------Shift-------------------
        mdIn.getValue(MDL_SHIFTX, XX(xyz), objId);
        mdIn.getValue(MDL_SHIFTY, YY(xyz), objId);
        mdIn.getValue(MDL_SHIFTZ, ZZ(xyz), objId);

        if (checkParam("--shift"))
        {
            std::cerr << "shifting..." << std::endl;
            double x = getDoubleParam("--shift", 0);
            XX(xyz) += x;
            YY(xyz) += (String(getParam("--shift", 1)) == "x") ? x : getDoubleParam("--shift", 1);
            ZZ(xyz) += (String(getParam("--shift", 2)) == "x") ? x : getDoubleParam("--shift", 2);
        }
        translationMatrix(T, XX(xyz), YY(xyz), ZZ(xyz));

        //------Rotate----------------------
        mdIn.getValue(MDL_ANGLEROT, XX(xyz), objId);
        mdIn.getValue(MDL_ANGLETILT, YY(xyz), objId);
        mdIn.getValue(MDL_ANGLEPSI, ZZ(xyz), objId);

        Matrix2D<double> R2(R);
        Euler_angles2matrix(XX(xyz), YY(xyz), ZZ(xyz), R2, true);

        if (checkParam("--rotate"))
        {
            std::cerr << "rotating..." << std::endl;
            A = R2 * R;
        }

        if (applyGeo)
        {
            //A = S * R2 * R * T;
            applyGeometry(splineDegree, imgOut(), img(), A, inverse, wrap);
        }
        else
        {
            mdIn.setValue(MDL_IMAGE, fnImgOut, objId);
            mdIn.setValue(MDL_SHIFTX, XX(xyz), objId);
            mdIn.setValue(MDL_SHIFTY, YY(xyz), objId);
            mdIn.setValue(MDL_SHIFTZ, ZZ(xyz), objId);
        }
        imgOut.write(fnImgOut);
    }

    void finishProcessing()
    {
        // To keep the mdIn info we overwrite the mdOut done by XmippMetadataProgram
        if (!applyGeo)
            mdOut = mdIn;
        XmippMetadataProgram::finishProcessing();
    }

public:
    /** Constructor */
    ProgTransformGeometry()
{}
}
; //end of class ProgScale


int main(int argc, char **argv)
{
    ProgTransformGeometry program;
    program.read(argc, argv);
    program.tryRun();
} //main

