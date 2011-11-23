/***************************************************************************
 * Authors:     AUTHOR_NAME (aerey@cnb.csic.es)
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
#include "transform_geometry.h"

ProgTransformGeometry::ProgTransformGeometry()
{}
ProgTransformGeometry::~ProgTransformGeometry()
{}

void ProgTransformGeometry::defineParams()
{
    each_image_produces_an_output = true;
    save_metadata_stack = true;
    keep_input_columns = true;
    temporaryOutput = false;
    XmippMetadataProgram::defineParams();
    //usage
    addUsageLine("Apply geometric transformations to images. You can shift, rotate and scale");
    addUsageLine("+a group of images/volumes. By default the geometric transformations will be");
    addUsageLine("+read from a metadata, if provided. Also, the original images will be preserved");
    addUsageLine("+if possible, trying to write out transformations to the output metadata.");
    addUsageLine("+If output is not specified, the original images will be overwritten.");
    //keywords
    addKeywords("transform, geometry, shift, rotate, scale, flip, pyramid, fourier");
    //examples
    addExampleLine("To simply apply the transformations in a metadata to the images:", false);
    addExampleLine("xmipp_transform_geometry -i mD1.doc --apply_transform");
    addExampleLine("Shift a volume by 10, 5 and -10 in x,y and z and dont wrapping", false);
    addExampleLine("xmipp_transform_geometry -i a.vol --shift 10 5 -10 -o b.vol --dont_wrap");
    addExampleLine("Scale a group of images to half size, not modifying image dimensions neither original image files", false);
    addExampleLine("xmipp_transform_geometry -i images.doc --scale factor 0.5 --oroot halvedOriginal");
    addExampleLine("Scaling (from 128 to 64) using Fourier (this will apply tranformations)", false);
    addExampleLine("xmipp_transform_geometry -i images.doc --scale dim 128 --oroot halvedFourierDim");
    addExampleLine("For pyramid scaling(reducing image size 4 times):", false);
    addExampleLine("xmipp_transform_geometry -i images.doc --scale pyramid -2 --oroot halvedPyramidDim");
    //params
    addParamsLine("== Transformations ==");
    addParamsLine("[--rotate <rotation_type>]   :Perform rotation");
    addParamsLine("         where <rotation_type>");
    addParamsLine("             ang <angle>     : Rotate an image (positive angle values is a clockwise rotation)");
    addParamsLine("             euler <rot> <tilt> <psi>     : Rotate with these Euler angles");
    addParamsLine("             alignZ <x> <y> <z>           : Align (x,y,z) with Z axis");
    addParamsLine("             axis <ang> <x=0> <y=0> <z=1> : Rotate <ang> degrees around (x,y,z)");
    addParamsLine("[--scale <scale_type>]                :Perform scaling");
    addParamsLine("                                    : All types of scaling produces a resize of images, except =factor_no_resize=");
    addParamsLine("         where <scale_type>");
    addParamsLine("             factor <n=1>           : Scaling factor, 0.5 halves and 2 doubles. (resize)");
    addParamsLine("             factor_no_resize <n=1> : Scaling factor, 0.5 halves and 2 doubles. (noresize)");
    addParamsLine("             dim <x> <y=x> <z=x>    : New x,y and z dimensions");
    addParamsLine("             fourier <x> <y=x> <thr=1> : Use padding/windowing in Fourier Space");
    addParamsLine("             requires --disable_metadata;");
    addParamsLine("             pyramid <levels=1>    : Use positive value to expand and negative to reduce");
    addParamsLine("             requires --disable_metadata;");
    addParamsLine(" alias -s;");
    addParamsLine("[--shift <x> <y=0> <z=0>]    : Shift by x, y and z");
    addParamsLine("[--flip]                                : Flip images, only valid for 2D");
    addParamsLine("== Other options ==");
    addParamsLine(" [--disable_metadata]               : Do not apply information stored in metadata file");
    addParamsLine(" alias -d;");
    addParamsLine(" [--interp <interpolation_type=spline>] : Interpolation type to be used. ");
    addParamsLine("      where <interpolation_type>");
    addParamsLine("        spline          : Use spline interpolation");
    addParamsLine("        linear          : Use bilinear/trilinear interpolation");
    addParamsLine("[--inverse]                         : Apply inverse transformations");
    addParamsLine("[--apply_transform]                 : By default, the original images are preserved");
    addParamsLine("                                    : and the alignment information is stored in metadata");
    addParamsLine("[--dont_wrap]                       : By default, the image/volume is wrapped");
    addParamsLine("[--write_matrix]                    : Print transformation matrix to screen");

}

void ProgTransformGeometry::readParams()
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
    flip = checkParam("--flip");
    scale_type = SCALE_NONE;
    disableMetadata = checkParam("--disable_metadata");
}


void ProgTransformGeometry::preProcess()
{
    //If zdimOut greater than 1, is a volume and should apply transform
    dim = (isVol = (zdimOut > 1)) ? 3 : 2;
    //Check that fourier interpolation is only for scale in 2d
    if (splineDegree == INTERP_FOURIER &&
        (checkParam("--shift") || checkParam("--rotate") || isVol))
        REPORT_ERROR(ERR_ARG_INCORRECT, "Fourier interpolation is only allowed for scale in 2D");

    applyTransform = applyTransform || isVol || mdIn.size() == 1;
    S.initIdentity(dim + 1);
    R.initIdentity(dim + 1);
    T.initIdentity(dim + 1);
    A.initIdentity(dim + 1);

    scaleV.resizeNoCopy(dim);
    scaleV.initConstant(1.);
    shiftV.resizeNoCopy(dim);
    rotV.resizeNoCopy(dim);

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
                XX(xyz) = getDoubleParam("--rotate", 2); //axis x
                YY(xyz) = getDoubleParam("--rotate", 3); //y
                ZZ(xyz) = getDoubleParam("--rotate", 4);//z
                rotation3DMatrix(ang, xyz, R, true);
            }
        }
    }

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
    	if (fn_out.empty() && oroot.empty())
    	{
    		fn_out.initRandom(10);
    		fn_out=fn_out.addExtension(fn_in.getExtension());
    		temporaryOutput=true;
    	}

        if (STR_EQUAL(getParam("--scale" ), "dim") || STR_EQUAL(getParam("--scale" ), "factor"))
        {
            //Calculate scale factor from images sizes and given dimensions
            //this approach assumes that all images have equal size

            double oxdim = xdimOut, oydim = ydimOut, ozdim = zdimOut;
            double factor;
            if (STR_EQUAL(getParam("--scale" ), "dim") )
            {
                xdimOut = getIntParam("--scale", 1);
                ydimOut = STR_EQUAL(getParam("--scale", 2), "x") ? xdimOut : getIntParam("--scale", 2);
                XX(scaleV) = (double)xdimOut / oxdim;
                YY(scaleV) = (double)ydimOut / oydim;
            }
            else
            {
                factor = getDoubleParam("--scale", 1);
                xdimOut = (int) (xdimOut * factor);
                ydimOut = (int) (ydimOut * factor);
                XX(scaleV) = factor;
                YY(scaleV) = factor;
            }


            //if scale factor is large splines s not the way to go, print a warning
            if( fabs(1.0-XX(scaleV)) > 0.1 )
            {
                std::cerr << "WARNING: Do not apply large scale factor using B-splines "
                << " choose fourier option ()"
                << std::endl;
            }
            if (isVol)
            {
                if (STR_EQUAL(getParam("--scale" ), "dim") )
                {
                    zdimOut = STR_EQUAL(getParam("--scale", 3), "x")
                              ? xdimOut : getIntParam("--scale", 3);
                    ZZ(scaleV) = (double)zdimOut / ozdim;
                }
                else
                {
                    zdimOut = (int) (zdimOut * factor);
                    ZZ(scaleV) = factor;
                }
            }
        }
        else if (STR_EQUAL(getParam("--scale"), "fourier"))
        {
            if (isVol)
                REPORT_ERROR(ERR_PARAM_INCORRECT, "The 'fourier' scaling type is only valid for images");
            int oxdim = xdimOut, oydim = ydimOut;
            scale_type = SCALE_FOURIER;

            xdimOut = getIntParam("--scale", 1);
            ydimOut = STR_EQUAL(getParam("--scale", 2), "x") ? xdimOut : getIntParam("--scale", 2);
            fourier_threads = getIntParam("--scale", 3);
            //Do not think this is true
            //            if (oxdim < xdimOut || oydim < ydimOut)
            //                REPORT_ERROR(ERR_PARAM_INCORRECT, "The 'fourier' scaling type can only be used for reducing size");
        }
        else if (STR_EQUAL(getParam("--scale"), "pyramid"))
        {
            scale_type = SCALE_PYRAMID_EXPAND;
            pyramid_level = getIntParam("--scale", 1);
            double scale_factor = (double)(pow(2.0, pyramid_level));
            xdimOut *= scale_factor;
            ydimOut *= scale_factor;
            if (isVol)
                zdimOut *= scale_factor;
            if (pyramid_level < 0)
            {
                pyramid_level *= -1; //change sign, negative means reduce operation
                scale_type = SCALE_PYRAMID_REDUCE;
            }
        }
        else //scale factor_no_resize case
        {
            double factor = getDoubleParam("--scale", 1);
            //Some extra validations for factor
            if (factor <= 0)
                REPORT_ERROR(ERR_VALUE_INCORRECT,"Factor must be a positive number");
            scaleV.initConstant(factor);
        }

        if (scale_type != SCALE_FACTOR)
            applyTransform = true;

        if (isVol)
            scale3DMatrix(scaleV, S);
        else
        {
            dMij(S, 0, 0) = XX(scaleV);
            dMij(S, 1, 1) = YY(scaleV);
        }
    }

    A = S * T * R;
    if (flip)
    {
        MAT_ELEM(A, 0, 0) *= -1.;
        MAT_ELEM(A, 0, 1) *= -1.;
    }
    if (inverse)
        A = A.inv();
}

void ProgTransformGeometry::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{

    B.initIdentity(dim + 1);

    if (!disableMetadata)
        geo2TransformationMatrix(rowOut, B);

    T = A * B;

    if (checkParam("--write_matrix"))
        std::cout << T << std::endl;

    if (applyTransform || fnImg != fnImgOut)
    {
        img.read(fnImg);
        img().setXmippOrigin();
    }

    if (applyTransform)
    {
        imgOut.setDatatype(img.getDatatype());

        /* if special scale do not resize output just rotate image */
        if (scale_type != SCALE_PYRAMID_EXPAND &&
            scale_type != SCALE_PYRAMID_REDUCE &&
            scale_type != SCALE_FOURIER )
            imgOut().resize(1, zdimOut, ydimOut, xdimOut, false);
        else
        	imgOut().resize(img(),false);
        imgOut().setXmippOrigin();
        //        if (only_scale)
        //         T.initIdentity();
        applyGeometry(splineDegree, imgOut(), img(), T, IS_NOT_INV, wrap, 0.);

        //imgOut.write(fnImgOut + ".before");

        switch (scale_type)
        {
        case SCALE_PYRAMID_EXPAND:
            selfPyramidExpand(splineDegree, imgOut(), pyramid_level);
            break;
        case SCALE_PYRAMID_REDUCE:
            selfPyramidReduce(splineDegree, imgOut(), pyramid_level);
            break;
        case SCALE_FOURIER:
            selfScaleToSizeFourier(ydimOut, xdimOut,imgOut(), fourier_threads);
            break;
        }
        imgOut.write(fnImgOut);
    }
    else
    {
        transformationMatrix2Geo(T, rowOut);
        if (fnImg != fnImgOut )
            img.write(fnImgOut);
    }
}

void ProgTransformGeometry::postProcess()
{
	if (temporaryOutput)
		std::rename(fn_out.c_str(),fn_in.c_str());
}
