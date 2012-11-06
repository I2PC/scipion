/***************************************************************************
 * Authors:     Joaquin Oton    (joton@cnb.csic.es)
 *              J.M. De la Rosa (jmdelarosa@cnb.csic.es)
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
#include "image_resize.h"

ProgImageResize::ProgImageResize()
{}

ProgImageResize::~ProgImageResize()
{}

void ProgImageResize::defineParams()
{
    each_image_produces_an_output = true;
    save_metadata_stack = true;
    keep_input_columns = true;
    temporaryOutput = false;
    XmippMetadataProgram::defineParams();
    //usage
    addUsageLine("Resize image or volume dimensions.");
    addUsageLine("+You can select different ways to resize your image/volumes.");
    //keywords
    addKeywords("transform, dimension, pyramid, fourier, scale");
    //examples

    //params
    addParamsLine("--factor <n=0.5>                 : Resize a factor of dimensions, 0.5 halves and 2 doubles.");
    addParamsLine(" alias -n;");
    addParamsLine("or --dim <x> <y=x> <z=x>         : New x,y and z dimensions");
    addParamsLine(" alias -d;");
    addParamsLine("or --fourier <x> <y=x> <thr=1>   : Use padding/windowing in Fourier Space to resize");
    addParamsLine(" alias -f;");
    addParamsLine("or --pyramid <levels=1>          : Use positive value to expand and negative to reduce");
    addParamsLine(" alias -p;");
    addParamsLine(" [--interp <interpolation_type=spline>] : Interpolation type to be used. ");
    addParamsLine("      where <interpolation_type>");
    addParamsLine("        spline          : Use spline interpolation");
    addParamsLine("        linear          : Use bilinear/trilinear interpolation");
}

void ProgImageResize::readParams()
{
    XmippMetadataProgram::readParams();
    String degree = getParam("--interp");
    if (degree == "spline")
        splineDegree = BSPLINE3;
    else if (degree == "linear")
        splineDegree = LINEAR;
    scale_type = RESIZE_NONE;
    apply_geo = !checkParam("--dont_apply_geo");
}


void ProgImageResize::preProcess()
{
    double factor;

    //If zdimOut greater than 1, is a volume and should apply transform
    dim = (isVol = (zdimOut > 1)) ? 3 : 2;
    //Check that fourier interpolation is only for scale in 2d
    if (splineDegree == INTERP_FOURIER &&
        (checkParam("--shift") || checkParam("--rotate") || isVol))
        REPORT_ERROR(ERR_ARG_INCORRECT, "Fourier interpolation is only allowed for scale in 2D");

    scaleV.resizeNoCopy(dim);
    scaleV.initConstant(1.);

    if (fn_out.empty() && oroot.empty())
    {
        fn_out=fn_in+"_tmp";
        fn_out=fn_out.addExtension(fn_in.getExtension());
        temporaryOutput=true;
    }

    if (checkParam("--factor") || checkParam("--dim"))
    {
        //Calculate scale factor from images sizes and given dimensions
        //this approach assumes that all images have equal size

        double oxdim = xdimOut, oydim = ydimOut, ozdim = zdimOut;

        if (checkParam("--dim"))
        {
            xdimOut = getIntParam("--dim", 0);
            ydimOut = STR_EQUAL(getParam("--dim", 1), "x") ? xdimOut : getIntParam("--dim", 1);
            XX(scaleV) = (double)xdimOut / oxdim;
            YY(scaleV) = (double)ydimOut / oydim;
        }
        else
        {
            factor = getDoubleParam("--factor");
            if (factor <= 0)
                REPORT_ERROR(ERR_VALUE_INCORRECT,"Resize factor must be a positive number");
            xdimOut = (int) (xdimOut * factor);
            ydimOut = (int) (ydimOut * factor);
            XX(scaleV) = factor;
            YY(scaleV) = factor;
        }

        //if scale factor is large splines s not the way to go, print a warning
        if( fabs(1.0-XX(scaleV)) > 0.1 )
        {
            reportWarning("Do not apply large scale factor using B-splines "
                          "choose fourier option.");
        }

        if (isVol)
        {
            if (checkParam("--dim"))
            {
                zdimOut = STR_EQUAL(getParam("--dim", 2), "x")
                          ? xdimOut : getIntParam("--dim", 2);
                ZZ(scaleV) = (double)zdimOut / ozdim;
            }
            else
            {
                zdimOut = (int) (zdimOut * factor);
                ZZ(scaleV) = factor;
            }
        }
        scale_type = RESIZE_FACTOR;
    }
    else if (checkParam("--fourier"))
    {
        if (isVol)
            REPORT_ERROR(ERR_PARAM_INCORRECT, "The 'fourier' scaling type is only valid for images");
        int oxdim = xdimOut, oydim = ydimOut;
        scale_type = RESIZE_FOURIER;

        xdimOut = getIntParam("--fourier", 0);
        ydimOut = STR_EQUAL(getParam("--fourier", 1), "x") ? xdimOut : getIntParam("--scale", 1);
        fourier_threads = getIntParam("--fourier", 2);
        //Do not think this is true
        //            if (oxdim < xdimOut || oydim < ydimOut)
        //                REPORT_ERROR(ERR_PARAM_INCORRECT, "The 'fourier' scaling type can only be used for reducing size");
    }
    else if (checkParam("--pyramid"))
    {
        scale_type = RESIZE_PYRAMID_EXPAND;
        pyramid_level = getIntParam("--pyramid", 0);
        factor = (double)(pow(2.0, pyramid_level));
        xdimOut *= factor;
        ydimOut *= factor;
        if (isVol)
            zdimOut *= factor;
        if (pyramid_level < 0)
        {
            pyramid_level *= -1; //change sign, negative means reduce operation
            scale_type = RESIZE_PYRAMID_REDUCE;
        }
    }
}

void ProgImageResize::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{

//    if (!apply_geo)
//        geo2TransformationMatrix(rowOut, B);

    img.read(fnImg);
    img().setXmippOrigin();

    imgOut.setDatatype(img.getDatatype());
    imgOut().resize(1, zdimOut, ydimOut, xdimOut, false);
    imgOut().setXmippOrigin();

    switch (scale_type)
    {
    case RESIZE_PYRAMID_EXPAND:
        selfPyramidExpand(splineDegree, imgOut(), pyramid_level);
        break;
    case RESIZE_PYRAMID_REDUCE:
        selfPyramidReduce(splineDegree, imgOut(), pyramid_level);
        break;
    case RESIZE_FOURIER:
        selfScaleToSizeFourier(ydimOut, xdimOut,imgOut(), fourier_threads);
        break;
    }
    imgOut.write(fnImgOut);
}

void ProgImageResize::postProcess()
{
    if (temporaryOutput)
        std::rename(fn_out.c_str(),fn_in.c_str());
}
