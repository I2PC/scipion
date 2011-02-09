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

#define INTERP_FOURIER -1

class ProgTransformGeometry: public XmippMetadataProgram
{
protected:
    int             zdim, ydim, xdim, splineDegree, dim;
    bool            applyTransform, inverse, wrap, doResize, isVol;
    Matrix2D<double> R, T, S, A, B;
    Matrix1D<double>          shiftV, rotV, scaleV;
    MDRow            input, transformation;
    Image<double> img, imgOut;

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
        addParamsLine("== Transformations ==");
        addParamsLine("[--rotate <rotation_type>]   :Perform rotation");
        addParamsLine("         where <rotation_type>");
        addParamsLine("             ang <angle>     : Rotate an image (positive angle values is a clockwise rotation)");
        addParamsLine("             euler <rot> <tilt> <psi>     : Rotate with these Euler angles");
        addParamsLine("             alignZ <x> <y> <z>           : Align (x,y,z) with Z axis");
        addParamsLine("             axis <ang> <x=0> <y=0> <z=1> : Rotate <ang> degrees around (x,y,z)");
        addParamsLine("[--scale <scale_type>]                :Perform scaling");
        addParamsLine("         where <scale_type>");
        addParamsLine("             factor <n=1>           : Scaling factor, 0.5 halves and 2 doubles");
        addParamsLine("                                    : This option doesn't cause a change on image dimensions.");
        addParamsLine("             dim <x> <y=x> <z=x>    : New x,y and z dimensions");
        addParamsLine("                                    : This cause a change on image dimensions.");
        addParamsLine(" alias -s;");
        addParamsLine("[--shift <x> <y=0> <z=0>]    : Shift by x, y and z");
        addParamsLine("[--flip]                                : Flip images, only valid for 2D");
        addParamsLine("== Other options ==");
        addParamsLine(" [--interp <interpolation_type=spline>] : Interpolation type to be used. ");
        addParamsLine("      where <interpolation_type>");
        addParamsLine("        spline          : Use spline interpolation");
        addParamsLine("        linear          : Use bilinear/trilinear interpolation");
        addParamsLine("        fourier <thr=1> : Use padding/windowing in Fourier Space (only for scale in 2D)");
        addParamsLine("            requires --scale;");
        addParamsLine("[--inverse]                         : Apply inverse transformations");
        addParamsLine("[--apply_transform]                 : By default, the original images are preserved");
        addParamsLine("                                    : and the alignment information is stored in metadata");
        addParamsLine("[--dont_wrap]                       : By default, the image/volume is wrapped");
        addParamsLine("[--write_matrix]                    : Print transformation matrix to screen");

    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        applyTransform = checkParam("--apply_transform");
        inverse = checkParam("--inverse");
        wrap = ! checkParam("--dont_wrap");
        String degree = getParam("--interp");
        if (degree == "spline")
            splineDegree = BSPLINE3;
        else if (degree == "linear")
            splineDegree = LINEAR;
        else
        {
            splineDegree = INTERP_FOURIER;

        }
    }


    void preProcess()
    {
        unsigned long n;
        ImgSize(mdIn, xdim, ydim , zdim, n);
        //If zdim greater than 1, is a volume and should apply transform
        dim = (isVol = (zdim > 1)) ? 3 : 2;
        //Check that fourier interpolation is only for scale in 2d
        if (splineDegree == INTERP_FOURIER &&
            (checkParam("--shift") || checkParam("--rotate") || isVol))
            REPORT_ERROR(ERR_ARG_INCORRECT, "Fourier interpolation is only allowed for scale in 2D");

        applyTransform = isVol || mdIn.size() == 1;
        S.initIdentity(dim + 1);
        R.initIdentity(dim + 1);
        T.initIdentity(dim + 1);
        A.initIdentity(dim + 1);

        imgOut().resizeNoCopy(zdim, ydim, xdim);
        scaleV.resizeNoCopy(dim);
        shiftV.resizeNoCopy(dim);
        rotV.resizeNoCopy(dim);

        if (checkParam("--shift"))
        {
            XX(shiftV) = getDoubleParam("--shift", 0);
            YY(shiftV) = getDoubleParam("--shift", 1);
            if (isVol)
            {
                ZZ(shiftV) = getDoubleParam("--shift", 2);
                translation3DMatrix(shiftV, T);
            }
            else
                translation2DMatrix(shiftV, T);
        }

        if (checkParam("--scale"))
        {
            if (STR_EQUAL(getParam("--scale" ), "dim"))
            {
                //Calculate scale factor from images sizes and given dimensions
                //this approach assumes that all images have equal size
                applyTransform = doResize = true;
                double oxdim = xdim, oydim = ydim, ozdim = zdim;
                xdim = getIntParam("--scale", 1);
                ydim = STR_EQUAL(getParam("--scale", 2), "x") ? xdim : getIntParam("--scale", 2);

                XX(scaleV) = (double)xdim / oxdim;
                YY(scaleV) = (double)ydim / oydim;

                if (isVol)
                {
                    zdim = STR_EQUAL(getParam("--scale", 3), "x") ? xdim : getIntParam("--scale", 3);
                    ZZ(scaleV) = (double)zdim / ozdim;
                }
                imgOut().resizeNoCopy(zdim, ydim, xdim);
            }
            else //scale factor case
            {
                double factor = getDoubleParam("--scale", 1);
                //Some extra validations for factor
                if (factor <= 0)
                    REPORT_ERROR(ERR_VALUE_INCORRECT,"Factor must be a positive number");
                scaleV.initConstant(factor);
            }

            if (isVol)
                scale3DMatrix(scaleV, S);
            else
            {
                dMij(S, 0, 0) = XX(scaleV);
                dMij(S, 1, 1) = YY(scaleV);
            }
        }

        if (checkParam("--rotate"))
        {
            if (STR_EQUAL(getParam("--rotate"), "ang"))
            {
                if (isVol)
                    REPORT_ERROR(ERR_PARAM_INCORRECT, "The selected rotation option is only valid for images");
                rotation2DMatrix(getDoubleParam("--rotate", 1), R, true);
            }
            else
            {
                if (!isVol)
                    REPORT_ERROR(ERR_PARAM_INCORRECT, "The selected rotation option is only valid for volumes");

                //In any case read following 3 values, leave euler angles in xyz
                Matrix1D<double> xyz(3);
                XX(xyz) = getDoubleParam("--rotate", 1); //rot
                YY(xyz) = getDoubleParam("--rotate", 2); //tilt
                ZZ(xyz) = getDoubleParam("--rotate", 3);//psi

                if (STR_EQUAL(getParam("--rotate"), "euler"))
                {
                    Euler_angles2matrix(XX(xyz), YY(xyz), ZZ(xyz), R, true);
                }
                else if (STR_EQUAL(getParam("--rotate"), "alignZ"))
                {
                    alignWithZ(xyz, R);
                }
                else
                {
                    double ang = getDoubleParam("--rotate", 1);
                    XX(xyz) = getDoubleParam("--rotate", 2); //rot
                    YY(xyz) = getDoubleParam("--rotate", 3); //tilt
                    ZZ(xyz) = getDoubleParam("--rotate", 4);//psi
                    rotation3DMatrix(ang, xyz, R, true);
                }
            }
        }

        A = S * T * R;
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        img.read(fnImg); //Doesn't apply geo to work with original images
        img().setXmippOrigin();
        B.initIdentity(dim + 1);

        if (!isVol)
        {
            mdIn.getRow(input, objId);//Get geometric transformation for image
            geo2TransformationMatrix(input, B);
        }
        //MultidimArray<double> out;
        std::cerr << "A: " << std::endl << A;
        std::cerr << "B: " << std::endl << B;
        T = A * B;
        if (checkParam("--write_matrix"))
            std::cout << T << std::endl;

        if (applyTransform)
        {
            applyGeometry(splineDegree, imgOut(), img(), T, !inverse, wrap);
        }
        else
        {
            if (inverse)
                T = T.inv();
            transformationMatrix2Geo(T, input);
            mdOut.setRow(input, newId);
            imgOut() = img();
        }
        imgOut.write(fnImgOut);
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

    return 0;
} //main

