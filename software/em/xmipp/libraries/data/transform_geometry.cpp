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
    allow_apply_geo = true;
    mdVol = false;
    XmippMetadataProgram::defineParams();
    //usage
    addUsageLine("Apply geometric transformations to images. You can shift, rotate and scale a group of images/volumes.");
    addUsageLine("+ By default the geometric transformations will be read from a metadata, if provided. Also, ");
    addUsageLine("+ the original images will be preserved if possible, trying to write out transformations to ");
    addUsageLine("+ the output metadata. If output is not specified, the original files will be overwritten.");
    addUsageLine("+ When output param -o is a stack file, associated metadata file is also generated using same rootname.");
    addUsageLine("+ If geometrical transformation are applied, involved labels are reset in output metadata.");
    //keywords
    addKeywords("transform, geometry, shift, rotate, scale, flip");
    //params
    addParamsLine("== Transformations ==");
    addParamsLine("[--rotate <ang=0>]   : Inplane rotation in 2D images.");
    addParamsLine("               : Positive angle is a clockwise rotation");
    addParamsLine("[--rotate_volume <rotation_type>]   : Rotation of volumes.");
    addParamsLine("         where <rotation_type>");
    addParamsLine("             euler <rot> <tilt> <psi>     : Rotate with these Euler angles (ZYZ convention)");
    addParamsLine("             matrix <r11> <r12> <r13> <r21> <r22> <r23> <r31> <r32> <r33>    : 3x3 rotation matrix, first index is row");
    addParamsLine("             alignZ <x> <y> <z>           : Align (x,y,z) with Z axis");
    addParamsLine("             axis <ang> <x=0> <y=0> <z=1> : Rotate <ang> degrees around (x,y,z)");
    addParamsLine("             icosahedral <from> <to>      : Rotate an icosahedral volume from i1 to i2, for instance");
    addParamsLine("                                          : valid symmetries are i1, i2, i3 and i4");
    addParamsLine("[--scale <factor=1>]                      : Perfom Scaling. Factor 0.5 halves and 2 doubles");
    addParamsLine(" alias -s;");
    addParamsLine("[--shift <x=0> <y=0> <z=0>]               : Shift by x, y and z");
    addParamsLine("[--flip]                                  : Flip images, only valid for 2D");
    addParamsLine("[--matrix <...>]                 : Apply directly the matrix transformation");
    addParamsLine("== Other options ==");
    addParamsLine(" [--interp <interpolation_type=spline>] : Interpolation type to be used. ");
    addParamsLine("      where <interpolation_type>");
    addParamsLine("        spline          : Use spline interpolation");
    addParamsLine("        linear          : Use bilinear/trilinear interpolation");
    addParamsLine("[--inverse]                         : Apply inverse transformations");
    addParamsLine("[--apply_transform]                 : By default, the original images are preserved");
    addParamsLine("                                    : and the alignment information is stored in metadata");
    addParamsLine("[--dont_wrap]                       : By default, the image/volume is wrapped");
    addParamsLine("[--write_matrix]                    : Print transformation matrix to screen");
    //examples
    addExampleLine("Write a metadata with geometrical transformations keeping the reference to original images:", false);
    addExampleLine("xmipp_transform_geometry -i mD1.xmd --shift 2 3 4 --scale 1.2 --rotate 23 -o newGeo.xmd");
    addExampleLine("Write a metadata with geometrical transformations copying original images to new stack file:", false);
    addExampleLine("xmipp_transform_geometry -i mD1.xmd --shift 2 3 4 --scale 1.2 --rotate 23 -o newGeo.stk");
    addExampleLine("Apply geometrical transformations to images ,reset involved labels and save in new metadata and stack files:", false);
    addExampleLine("xmipp_transform_geometry -i mD1.xmd --shift 2 3 4 --scale 1.2 --rotate 23 -o newGeo.xmd --apply_transform");
    addExampleLine("To simply apply the transformations in a metadata to the images:", false);
    addExampleLine("xmipp_transform_geometry -i mD1.xmd --apply_transform");
    addExampleLine("Shift a volume by 10, 5 and -10 in x,y and z and do not wrap", false);
    addExampleLine("xmipp_transform_geometry -i a.vol --shift 10 5 -10 -o b.vol --dont_wrap");
    addExampleLine("Scale a group of images to half size, not modifying image dimensions neither original image files", false);
    addExampleLine("xmipp_transform_geometry -i images.xmd --scale 0.5 -o halvedOriginal.xmd");

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

    /** In most cases output "-o" is a metadata with the new geometry keeping 
     *  the names of input images so we set the flags to keep the same image 
     *  names in the output metadata
     */
    if ( !checkParam("--oroot") && fn_out.hasMetadataExtension())
    {
        if ( input_is_metadata && !applyTransform )
            each_image_produces_an_output = !(produces_a_metadata = true);
        else /** If "-o" is a metadata but we are writing output images, 
              *  -o can only be a stack if --oroot is no passed, and 
              *  then MD is generated automatically **/
            fn_out = fn_out.replaceExtension("stk");
    }
    else if ( !checkParam("--oroot") && !checkParam("-o") )
        produces_a_metadata = true;
}


void ProgTransformGeometry::calculateRotationMatrix()
{
    const char * rotateType = getParam("--rotate_volume");
    Matrix1D<double> xyz(3);
    if (!STR_EQUAL(rotateType, "icosahedral"))
    {   // params are char for icosahedral option
        XX(xyz) = getDoubleParam("--rotate_volume", 1); //rot
        YY(xyz) = getDoubleParam("--rotate_volume", 2); //tilt
        ZZ(xyz) = getDoubleParam("--rotate_volume", 3); //psi
    }

    if (STR_EQUAL(rotateType, "euler"))
        Euler_angles2matrix(XX(xyz), YY(xyz), ZZ(xyz), R, true);
    else if (STR_EQUAL(rotateType, "matrix"))
    {
        R.initZeros(4,4);
        MAT_ELEM(R, 3, 3) = 1;
        for (int i = 0; i < 9; ++i)
        {
            int r = i / 3;
            int c = i % 3;
            MAT_ELEM(R, r, c) =  getDoubleParam("--rotate_volume", i+1);
        }
    }
    else if (STR_EQUAL(rotateType, "alignZ"))
        alignWithZ(xyz, R);
    else if (STR_EQUAL(rotateType, "icosahedral"))
    {
        const char * icoFrom = getParam("--rotate_volume", 1);
        const char * icoTo = getParam("--rotate_volume", 2);
        rotation3DMatrixFromIcoOrientations(icoFrom, icoTo, R);
    }
    else
    {
        double ang = getDoubleParam("--rotate_volume", 1);
        XX(xyz) = getDoubleParam("--rotate_volume", 2); //axis x
        YY(xyz) = getDoubleParam("--rotate_volume", 3); //y
        ZZ(xyz) = getDoubleParam("--rotate_volume", 4);//z
        rotation3DMatrix(ang, xyz, R, true);
    }
}

void ProgTransformGeometry::preProcess()
{
    //If zdimOut greater than 1, is a volume and should apply transform
    dim = (isVol = (zdimOut > 1)) ? 3 : 2;

    applyTransform = applyTransform || isVol || mdInSize == 1;
    R.initIdentity(dim + 1);
    A.initIdentity(dim + 1);

    MDRow rowGeo;
    rowGeo.setValue(MDL_SHIFT_X, getDoubleParam("--shift", 0));
    rowGeo.setValue(MDL_SHIFT_Y, getDoubleParam("--shift", 1));
    rowGeo.setValue(MDL_SCALE, getDoubleParam("--scale"));

    if (isVol)
    {
        if (checkParam("--rotate_volume"))
            calculateRotationMatrix();
        else
        	mdVol = true;
        rowGeo.setValue(MDL_SHIFT_Z, getDoubleParam("--shift", 2));

    }
    else
        rowGeo.setValue(MDL_ANGLE_PSI, getDoubleParam("--rotate"));

    geo2TransformationMatrix(rowGeo, A);

    A = A * R;

    if (flip)
    {
        MAT_ELEM(A, 0, 0) *= -1.;
        MAT_ELEM(A, 0, 1) *= -1.;
        if (dim == 3)
            MAT_ELEM(A, 0, 2) *= -1.;
    }
    if (inverse)
        A = A.inv();
}

void ProgTransformGeometry::processImage(const FileName &fnImg,
                                         const FileName &fnImgOut,
                                         const MDRow &rowIn,
                                         MDRow &rowOut)
{

    if (checkParam("--matrix"))
    {
      // In this case we are directly reading the transformation matrix
      // from the arguments passed
      matrixStr = getParam("--matrix");
      string2TransformationMatrix(matrixStr, T);
    }
    else
    {
      B.initIdentity(dim + 1);


    if (apply_geo || mdVol)
        geo2TransformationMatrix(rowOut, B);

      T = A * B;
    }

    if (checkParam("--write_matrix"))
        std::cerr << T << std::endl;

    if (applyTransform || fnImg != fnImgOut)
        img.read(fnImg);

    if (applyTransform)
    {
        img().setXmippOrigin();
        imgOut.setDatatype(img.getDatatype());
        imgOut().resize(1, zdimOut, ydimOut, xdimOut, false);
        imgOut().setXmippOrigin();
        applyGeometry(splineDegree, imgOut(), img(), T, IS_NOT_INV, wrap, 0.);
        imgOut.write(fnImgOut);
        rowOut.resetGeo(false);
    }
    else
    {
        transformationMatrix2Geo(T, rowOut);
        if (fnImg != fnImgOut )
            img.write(fnImgOut);
    }
}
