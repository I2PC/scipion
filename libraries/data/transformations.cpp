/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Sjors H.W. Scheres
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

#include "transformations.h"
#include "filters.h"

void geo2TransformationMatrix(const MDRow &imageGeo, Matrix2D<double> &A,
                              bool only_apply_shifts)
{
    // This has only been implemented for 2D images...
    double psi = 0, shiftX = 0., shiftY = 0., scale = 1.;
    bool flip = false;

    imageGeo.getValue(MDL_ANGLE_PSI, psi);
    imageGeo.getValue(MDL_SHIFT_X, shiftX);
    imageGeo.getValue(MDL_SHIFT_Y, shiftY);
    imageGeo.getValue(MDL_SCALE, scale);
    imageGeo.getValue(MDL_FLIP, flip);

    psi = realWRAP(psi, 0., 360.);

    int dim = A.Xdim() - 1;
    //This check the case when matrix A is not initialized with correct size
    if (dim < 2 || dim > 3)
    {
        dim = 3;
        A.resizeNoCopy(dim + 1, dim + 1);
    }

    if (only_apply_shifts)
        A.initIdentity();
    else if (dim == 2) //2D geometry
        rotation2DMatrix(psi, A, true);
    else if (dim == 3)//3D geometry
    {
        double rot = 0., tilt = 0., shiftZ = 0.;
        imageGeo.getValue(MDL_ANGLE_ROT, rot);
        imageGeo.getValue(MDL_ANGLE_TILT, tilt);
        imageGeo.getValue(MDL_SHIFT_Z, shiftZ);
        Euler_angles2matrix(rot, tilt, psi, A, true);
        dMij(A, 2, dim) = shiftZ;
    }
    dMij(A, 0, dim) = shiftX;
    dMij(A, 1, dim) = shiftY;

    if (scale != 1.)
    {
        if (dim == 2)
        {
            M3x3_BY_CT(A, A, scale);
        }
        else if (dim == 3)
        {
            M4x4_BY_CT(A, A, scale);
        }
        dMij(A, dim, dim) = 1.;
    }

    if (flip)
    {
        dMij(A, 0, 0) *= -1.;
        dMij(A, 0, 1) *= -1.;
        if (dim == 3)
            dMij(A, 0, 2) *= -1.;
    }
}

void transformationMatrix2Parameters2D(const Matrix2D<double> &A, bool &flip,
                                       double &scale, double &shiftX, double &shiftY, double &psi)
{
    //Calculate determinant for getting flip
    flip = ((dMij(A, 0, 0) * dMij(A, 1, 1) - dMij(A, 0, 1) * dMij(A, 1, 0) ) < 0);
    int sgn = flip ? -1 : 1;
    double cosine = sgn * dMij(A, 0, 0), sine = sgn * dMij(A, 0, 1);
    double scale2 = cosine * cosine +  sine * sine;
    scale = sqrt(scale2);
    double invScale = 1 / scale;
    shiftX = dMij(A, 0, 2) * invScale;
    shiftY = dMij(A, 1, 2) * invScale;
    psi = RAD2DEG(atan2(sine, cosine));
}

void transformationMatrix2Parameters3D(const Matrix2D<double> &A, bool &flip, double &scale,
                                       double &shiftX, double &shiftY, double &shiftZ,
                                       double &rot, double &tilt, double &psi)
{
    scale =  sqrt(dMij(A,2,0)*dMij(A,2,0) \
                  + dMij(A,2,1)*dMij(A,2,1)\
                  + dMij(A,2,2)*dMij(A,2,2) );
    double invScale = 1./ scale;

    Matrix2D<double> tmpMatrix(4,4);
    M4x4_BY_CT(tmpMatrix, A, invScale);

    Matrix2D<double> eulerMatrix(3,3);

    FOR_ALL_ELEMENTS_IN_MATRIX2D(eulerMatrix)
    dMij(eulerMatrix,i,j) = dMij(tmpMatrix, i, j);
    //check determinant if -1 then flip = true
    flip = tmpMatrix.det3x3() < 0;
    if (flip)
    {
        dMij(eulerMatrix, 0, 0) *= -1.;
        dMij(eulerMatrix, 0, 1) *= -1.;
        dMij(eulerMatrix, 0, 2) *= -1.;
    }
    Euler_matrix2angles(eulerMatrix, rot, tilt, psi);

    shiftX = dMij(tmpMatrix,0,3);
    shiftY = dMij(tmpMatrix,1,3);
    shiftZ = dMij(tmpMatrix,2,3);
}

#define ADD_IF_EXIST_NONZERO(label, value) if (imageGeo.containsLabel(label) || !XMIPP_EQUAL_ZERO(value))\
                                                  imageGeo.setValue(label, value);
void transformationMatrix2Geo(const Matrix2D<double> &A, MDRow & imageGeo)
{
    bool flip = false;
    double scale = 1 , shiftX = 0, shiftY = 0, psi = 0, shiftZ = 0, rot = 0, tilt = 0;

    int dim = A.Xdim() - 1;

    if (dim == 2)
        transformationMatrix2Parameters2D(A, flip, scale, shiftX, shiftY, psi);
    else if (dim == 3)
    {
        transformationMatrix2Parameters3D(A, flip, scale, shiftX, shiftY, shiftZ, rot,tilt, psi);
        ADD_IF_EXIST_NONZERO(MDL_ANGLE_ROT, rot);
        ADD_IF_EXIST_NONZERO(MDL_ANGLE_TILT, tilt);
        ADD_IF_EXIST_NONZERO(MDL_SHIFT_Z, shiftZ);
    }

    ADD_IF_EXIST_NONZERO(MDL_ANGLE_PSI, psi);
    ADD_IF_EXIST_NONZERO(MDL_SHIFT_X, shiftX);
    ADD_IF_EXIST_NONZERO(MDL_SHIFT_Y, shiftY);

    if (imageGeo.containsLabel(MDL_SCALE) || !XMIPP_EQUAL_REAL(scale, 1.))
        imageGeo.setValue(MDL_SCALE, scale);

    if (imageGeo.containsLabel(MDL_FLIP) || flip)
        imageGeo.setValue(MDL_FLIP, flip);
}

/* Rotation 2D ------------------------------------------------------------- */
void rotation2DMatrix(double ang, Matrix2D< double > &result, bool homogeneous)
{
    double cosine, sine;

    ang = DEG2RAD(ang);
    cosine = cos(ang);
    sine = sin(ang);

    if (homogeneous)
    {
        if (MAT_XSIZE(result)!=3 || MAT_YSIZE(result)!=3)
            result.resizeNoCopy(3,3);
        dMij(result,0, 2) = 0;
        dMij(result,1, 2) = 0;
        dMij(result,2, 0) = 0;
        dMij(result,2, 1) = 0;
        dMij(result,2, 2) = 1;
    }
    else
        if (MAT_XSIZE(result)!=2 || MAT_YSIZE(result)!=2)
            result.resizeNoCopy(2,2);
    dMij(result,0, 0) = cosine;
    dMij(result,0, 1) = sine;
    dMij(result,1, 0) = -sine;
    dMij(result,1, 1) = cosine;
}

/* Translation 2D ---------------------------------------------------------- */
void translation2DMatrix(const Matrix1D<double> &v,
                         Matrix2D< double > &result,
                         bool inverse)
{
    if (VEC_XSIZE(v) != 2)
        REPORT_ERROR(ERR_MATRIX_SIZE, "Translation2D_matrix: vector is not in R2");

    result.initIdentity(3);
    if (inverse)
    {
        dMij(result,0, 2) = -XX(v);
        dMij(result,1, 2) = -YY(v);
    }
    else
    {
        dMij(result,0, 2) = XX(v);
        dMij(result,1, 2) = YY(v);
    }

}

/* Rotation 3D around the system axes -------------------------------------- */
void rotation3DMatrix(double ang, char axis, Matrix2D< double > &result,
                      bool homogeneous)
{
    if (homogeneous)
    {
        result.initZeros(4,4);
        dMij(result,3, 3) = 1;
    }
    else
        result.initZeros(3,3);

    double cosine, sine;
    ang = DEG2RAD(ang);
    cosine = cos(ang);
    sine = sin(ang);

    switch (axis)
    {
    case 'Z':
        dMij(result,0, 0) = cosine;
        dMij(result,0, 1) = sine;
        dMij(result,1, 0) = -sine;
        dMij(result,1, 1) = cosine;
        dMij(result,2, 2) = 1;
        break;
    case 'Y':
        dMij(result,0, 0) = cosine;
        dMij(result,0, 2) = sine;
        dMij(result,2, 0) = -sine;
        dMij(result,2, 2) = cosine;
        dMij(result,1, 1) = 1;
        break;
    case 'X':
        dMij(result,1, 1) = cosine;
        dMij(result,1, 2) = sine;
        dMij(result,2, 1) = -sine;
        dMij(result,2, 2) = cosine;
        dMij(result,0, 0) = 1;
        break;
    default:
        REPORT_ERROR(ERR_VALUE_INCORRECT, "rotation3DMatrix: Unknown axis");
    }
}

/* Align a vector with Z axis */
void alignWithZ(const Matrix1D<double> &axis, Matrix2D<double>& result,
                bool homogeneous)
{
    if (axis.size() != 3)
        REPORT_ERROR(ERR_MATRIX_SIZE, "alignWithZ: Axis is not in R3");
    if (homogeneous)
    {
        result.initZeros(4,4);
        dMij(result,3, 3) = 1;
    }
    else
        result.initZeros(3,3);
    Matrix1D<double>  Axis(axis);
    Axis.selfNormalize();

    // Compute length of the projection on YZ plane
    double proj_mod = sqrt(YY(Axis) * YY(Axis) + ZZ(Axis) * ZZ(Axis));
    if (proj_mod > XMIPP_EQUAL_ACCURACY)
    {   // proj_mod!=0
        // Build Matrix result, which makes the turning axis coincident with Z
        dMij(result,0, 0) = proj_mod;
        dMij(result,0, 1) = -XX(Axis) * YY(Axis) / proj_mod;
        dMij(result,0, 2) = -XX(Axis) * ZZ(Axis) / proj_mod;
        dMij(result,1, 0) = 0;
        dMij(result,1, 1) = ZZ(Axis) / proj_mod;
        dMij(result,1, 2) = -YY(Axis) / proj_mod;
        dMij(result,2, 0) = XX(Axis);
        dMij(result,2, 1) = YY(Axis);
        dMij(result,2, 2) = ZZ(Axis);
    }
    else
    {
        // I know that the Axis is the X axis, EITHER POSITIVE OR NEGATIVE!!
        dMij(result,0, 0) = 0;
        dMij(result,0, 1) = 0;
        dMij(result,0, 2) = (XX(Axis) > 0)? -1 : 1;
        dMij(result,1, 0) = 0;
        dMij(result,1, 1) = 1;
        dMij(result,1, 2) = 0;
        dMij(result,2, 0) = (XX(Axis) > 0)? 1 : -1;
        dMij(result,2, 1) = 0;
        dMij(result,2, 2) = 0;
    }
}

/* Rotation 3D around any axis -------------------------------------------- */
void rotation3DMatrix(double ang, const Matrix1D<double> &axis,
                      Matrix2D<double> &result, bool homogeneous)
{
#ifdef NEVERDEFINED
    // Compute a matrix which makes the turning axis coincident with Z
    // And turn around this axis
    Matrix2D<double> A, R;
    alignWithZ(axis, A, homogeneous);
    rotation3DMatrix(ang, 'Z', R, homogeneous);
    result = A.transpose() * R * A;
#else
    // http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    if (homogeneous)
        result.initIdentity(4);
    else
        result.initIdentity(3);
    double s,c;
    sincos(-DEG2RAD(ang),&s,&c);
    double c1=1-c;
    double x=XX(axis);
    double y=YY(axis);
    double z=ZZ(axis);
    double xy=x*y;
    double xz=x*z;
    double yz=y*z;
    double x2=x*x;
    double y2=y*y;
    double z2=z*z;
    dMij(result,0,0)=c+x2*c1;
    dMij(result,0,1)=xy*c1-z*s;
    dMij(result,0,2)=xz*c1+y*s;
    dMij(result,1,0)=xy*c1+z*s;
    dMij(result,1,1)=c+y2*c1;
    dMij(result,1,2)=yz*c1-x*s;
    dMij(result,2,0)=xz*c1-y*s;
    dMij(result,2,1)=yz*c1+x*s;
    dMij(result,2,2)=c+z2*c1;
#endif
}

/* Translation 3D ---------------------------------------------------------- */
void translation3DMatrix(const Matrix1D<double> &v, Matrix2D<double> &result, bool inverse)
{
    if (VEC_XSIZE(v) != 3)
        REPORT_ERROR(ERR_MATRIX_SIZE, "Translation3D_matrix: vector is not in R3");

    result.initIdentity(4);
    if (inverse)
    {
        dMij(result,0, 3) = -XX(v);
        dMij(result,1, 3) = -YY(v);
        dMij(result,2, 3) = -ZZ(v);
    }
    else
    {
        dMij(result,0, 3) = XX(v);
        dMij(result,1, 3) = YY(v);
        dMij(result,2, 3) = ZZ(v);
    }
}

/* Scale 3D ---------------------------------------------------------------- */
void scale3DMatrix(const Matrix1D<double> &sc, Matrix2D<double>& result,
                   bool homogeneous)
{
    if (VEC_XSIZE(sc) != 3)
        REPORT_ERROR(ERR_MATRIX_SIZE, "Scale3D_matrix: vector is not in R3");

    if (homogeneous)
    {
        result.initZeros(4,4);
        dMij(result,3, 3) = 1;
    }
    else
        result.initZeros(3,3);
    dMij(result,0, 0) = XX(sc);
    dMij(result,1, 1) = YY(sc);
    dMij(result,2, 2) = ZZ(sc);
}

// Special case for complex numbers
template<>
void applyGeometry(int SplineDegree,
                   MultidimArray< std::complex<double> >& V2,
                   const MultidimArray< std::complex<double> >& V1,
                   const Matrix2D< double > &A, bool inv,
                   bool wrap, std::complex<double> outside)
{

    if (SplineDegree > 1)
    {
        MultidimArray<double> re, im, rotre, rotim;
        MultidimArray<std::complex<double> > oneImg;
        double outre, outim;
        re.resize(ZSIZE(V1), YSIZE(V1), XSIZE(V1));
        im.resize(ZSIZE(V1), YSIZE(V1), XSIZE(V1));
        outre = outside.real();
        outim = outside.imag();
        oneImg=V1;
        Complex2RealImag(MULTIDIM_ARRAY(oneImg),
                         MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                         MULTIDIM_SIZE(oneImg));
        applyGeometry(SplineDegree, rotre, re, A, inv, wrap, outre);
        applyGeometry(SplineDegree, rotim, im, A, inv, wrap, outim);
        V2.resize(oneImg);
        RealImag2Complex(MULTIDIM_ARRAY(rotre), MULTIDIM_ARRAY(rotim),
                         MULTIDIM_ARRAY(V2), MULTIDIM_SIZE(re));
    }
    else
    { //FIXME I do not think you want to recall your self
        REPORT_ERROR(ERR_NOT_IMPLEMENTED,"I do not think you want to recall your self");
        applyGeometry(SplineDegree, V2, V1, A, inv, wrap, outside);
    }
}

// Special case for complex numbers
template<>
void selfApplyGeometry(int Splinedegree,
                       MultidimArray< std::complex<double> > &V1,
                       const Matrix2D<double> &A, bool inv,
                       bool wrap, std::complex<double> outside)
{
    MultidimArray<std::complex<double> > aux = V1;
    applyGeometry(Splinedegree, V1, aux, A, inv, wrap, outside);
}

void applyGeometry(int SplineDegree,
                   MultidimArrayGeneric &V2,
                   const MultidimArrayGeneric &V1,
                   const Matrix2D< double > &A, bool inv,
                   bool wrap, double outside)
{
#define APPLYGEO(type)  applyGeometry(SplineDegree,(*(MultidimArray<type>*)(V2.im)), (*(MultidimArray<type>*)(V1.im)), A, inv, wrap, (type) outside);
    SWITCHDATATYPE(V1.datatype, APPLYGEO)
#undef APPLYGEO

}

// Special case for complex arrays
void produceSplineCoefficients(int SplineDegree,
                               MultidimArray< double > &coeffs,
                               const MultidimArray< std::complex<double> > &V1)
{
    // TODO Implement
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"Spline coefficients of a complex matrix is not implemented.");
}


void selfScaleToSize(int SplineDegree,
                     MultidimArrayGeneric &V1,
                     int Xdim, int Ydim, int Zdim)
{
#define SELFSCALETOSIZE(type) selfScaleToSize(SplineDegree,MULTIDIM_ARRAY_TYPE(V1,type),Xdim,Ydim,Zdim);
    SWITCHDATATYPE(V1.datatype,SELFSCALETOSIZE)
#undef SELFSCALETOSIZE
}

void scaleToSize(int SplineDegree,
                        MultidimArrayGeneric &V2,
                        const MultidimArrayGeneric &V1,int Xdim, int Ydim, int Zdim)
{
  if (V1.datatype != V2.datatype)
    REPORT_ERROR(ERR_PARAM_INCORRECT, "scaleToSize: MultidimArrayGeneric requires same datatype");
#define SCALETOSIZE(type) scaleToSize(SplineDegree,MULTIDIM_ARRAY_TYPE(V2,type),MULTIDIM_ARRAY_TYPE(V1,type),Xdim,Ydim,Zdim);
    SWITCHDATATYPE(V1.datatype, SCALETOSIZE)
#undef SCALETOSIZE
}

// Special case for complex arrays
void scaleToSize(int SplineDegree,
                 MultidimArray< std::complex<double> > &V2,
                 const MultidimArray< std::complex<double> > &V1,
                 int Xdim, int Ydim, int Zdim)
{
    if (SplineDegree > 1)
    {
        MultidimArray< double > re, im, aux;
        MultidimArray<std::complex<double> > oneImg;

        re.resize(ZSIZE(V1), YSIZE(V1), XSIZE(V1));
        im.resize(ZSIZE(V1), YSIZE(V1), XSIZE(V1));

        oneImg=V1;
        Complex2RealImag(MULTIDIM_ARRAY(oneImg),
                         MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                         MULTIDIM_SIZE(oneImg));
        aux = re;
        scaleToSize(SplineDegree, re, aux, Ydim, Xdim, Zdim);
        aux = im;
        scaleToSize(SplineDegree, im, aux, Ydim, Xdim, Zdim);
        RealImag2Complex(MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                         MULTIDIM_ARRAY(V2), MULTIDIM_SIZE(re));
    }
    else
        scaleToSize(SplineDegree, V2, V1, Xdim, Ydim, Zdim);

}

// Special case for complex arrays
void selfScaleToSize(int SplineDegree,
                     MultidimArray< std::complex<double> > &V1,
                     int Xdim, int Ydim, int Zdim)
{
    MultidimArray<std::complex<double> > aux;
    scaleToSize(SplineDegree, V1, aux, Xdim, Ydim, Zdim);
}

/** Same as template version but for MultidimArrayGeneric */
void selfPyramidReduce(int SplineDegree,
                       MultidimArrayGeneric &V1,
                       int levels)
{
#define SELFPYRAMIDREDUCE(type) selfPyramidReduce(SplineDegree, *((MultidimArray<type>*)(V1.im)), levels);
    SWITCHDATATYPE(V1.datatype,SELFPYRAMIDREDUCE);
#undef SELFPYRAMIDREDUCE
}

/** Same as previous but for MultidimArrayGeneric */
void selfPyramidExpand(int SplineDegree,
                       MultidimArrayGeneric &V1,
                       int levels)
{
#define SELFPYRAMIDEXPAND(type) selfPyramidExpand(SplineDegree, *((MultidimArray<type>*)(V1.im)), levels);
    SWITCHDATATYPE(V1.datatype,SELFPYRAMIDEXPAND);
#undef SELFPYRAMIDEXPAND
}

void pyramidExpand(int SplineDegree,
                   MultidimArrayGeneric &V2,
                   const MultidimArrayGeneric &V1,
                   int levels)
{
  if (V1.datatype != V2.datatype)
    REPORT_ERROR(ERR_PARAM_INCORRECT, "pyramidExpand: MultidimArrayGeneric requires same datatype");
#define PYRAMIDEXPAND(type) pyramidExpand(SplineDegree,MULTIDIM_ARRAY_TYPE(V2,type),MULTIDIM_ARRAY_TYPE(V1,type), levels);
    SWITCHDATATYPE(V1.datatype, PYRAMIDEXPAND)
#undef PYRAMIDEXPAND
}

void pyramidReduce(int SplineDegree,
                   MultidimArrayGeneric &V2,
                   const MultidimArrayGeneric &V1,
                   int levels)
{
  if (V1.datatype != V2.datatype)
    REPORT_ERROR(ERR_PARAM_INCORRECT, "pyramidReduce: MultidimArrayGeneric requires same datatype");
#define PYRAMIDREDUCE(type) pyramidReduce(SplineDegree,MULTIDIM_ARRAY_TYPE(V2,type),MULTIDIM_ARRAY_TYPE(V1,type), levels);
    SWITCHDATATYPE(V1.datatype, PYRAMIDREDUCE)
#undef PYRAMIDREDUCE
}

/** Interpolates the value of the 3D matrix M at the point (x,y,z) knowing
* that this image is a set of B-spline coefficients. And making the diff
* of x, such->  V=sum(Coef diff(Bx) By Bz)
* Only for BSplines of degree 3!!
* @ingroup VolumesMemory
*
* (x,y,z) are in logical coordinates.
*/
double interpolatedElementBSplineDiffX(MultidimArray<double> &vol, double x, double y, double z,
                                       int SplineDegree)
{
    int SplineDegree_1 = SplineDegree - 1;
    double aux;

    // Logical to physical
    z -= STARTINGZ(vol);
    y -= STARTINGY(vol);
    x -= STARTINGX(vol);

    int l1 = CEIL(x - SplineDegree_1);
    int l2 = l1 + SplineDegree;

    int m1 = CEIL(y - SplineDegree_1);
    int m2 = m1 + SplineDegree;

    int n1 = CEIL(z - SplineDegree_1);
    int n2 = n1 + SplineDegree;

    double zyxsum = 0.0;
    int Zdim=(int)ZSIZE(vol);
    int Ydim=(int)YSIZE(vol);
    int Xdim=(int)XSIZE(vol);
    for (int n = n1; n <= n2; n++)
    {
        int equivalent_n=n;
        if      (n<0)
            equivalent_n=-n-1;
        else if (n>=Zdim)
            equivalent_n=2*Zdim-n-1;
        double yxsum = 0.0;
        for (int m = m1; m <= m2; m++)
        {
            int equivalent_m=m;
            if      (m<0)
                equivalent_m=-m-1;
            else if (m>=Ydim)
                equivalent_m=2*Ydim-m-1;
            double xsum = 0.0;
            for (int l = l1; l <= l2; l++)
            {
                double xminusl = x - (double) l;
                int equivalent_l=l;
                if      (l<0)
                    equivalent_l=-l-1;
                else if (l>=Xdim)
                    equivalent_l=2*Xdim-l-1;
                double Coeff = (double) DIRECT_A3D_ELEM(vol,
                                                        equivalent_n,equivalent_m,equivalent_l);
                switch (SplineDegree)
                {
                case 2:
                    xsum += Coeff * Bspline02(xminusl);
                    break;
                case 3:
                    BSPLINE03DIFF1(aux,xminusl);
                    xsum += Coeff * aux;
                    break;
                case 4:
                    xsum += Coeff * Bspline04(xminusl);
                    break;
                case 5:
                    xsum += Coeff * Bspline05(xminusl);
                    break;
                case 6:
                    xsum += Coeff * Bspline06(xminusl);
                    break;
                case 7:
                    xsum += Coeff * Bspline07(xminusl);
                    break;
                case 8:
                    xsum += Coeff * Bspline08(xminusl);
                    break;
                case 9:
                    xsum += Coeff * Bspline09(xminusl);
                    break;
                }
            }

            double yminusm = y - (double) m;
            switch (SplineDegree)
            {
            case 2:
                yxsum += xsum * Bspline02(yminusm);
                break;
            case 3:
                BSPLINE03(aux,yminusm);
                yxsum += xsum * aux;
                break;
            case 4:
                yxsum += xsum * Bspline04(yminusm);
                break;
            case 5:
                yxsum += xsum * Bspline05(yminusm);
                break;
            case 6:
                yxsum += xsum * Bspline06(yminusm);
                break;
            case 7:
                yxsum += xsum * Bspline07(yminusm);
                break;
            case 8:
                yxsum += xsum * Bspline08(yminusm);
                break;
            case 9:
                yxsum += xsum * Bspline09(yminusm);
                break;
            }
        }

        double zminusn = z - (double) n;
        switch (SplineDegree)
        {
        case 2:
            zyxsum += yxsum * Bspline02(zminusn);
            break;
        case 3:
            BSPLINE03(aux,zminusn);
            zyxsum += yxsum * aux;
            break;
        case 4:
            zyxsum += yxsum * Bspline04(zminusn);
            break;
        case 5:
            zyxsum += yxsum * Bspline05(zminusn);
            break;
        case 6:
            zyxsum += yxsum * Bspline06(zminusn);
            break;
        case 7:
            zyxsum += yxsum * Bspline07(zminusn);
            break;
        case 8:
            zyxsum += yxsum * Bspline08(zminusn);
            break;
        case 9:
            zyxsum += yxsum * Bspline09(zminusn);
            break;
        }
    }

    return zyxsum;
}

/** Interpolates the value of the 3D matrix M at the point (x,y,z) knowing
 * that this image is a set of B-spline coefficients. And making the diff
 * of y, such->  V=sum(Coef Bx diff(By) Bz)
 * Only for BSplines of degree 3!!
 * @ingroup VolumesMemory
 *
 * (x,y,z) are in logical coordinates.
 */
double interpolatedElementBSplineDiffY(MultidimArray<double> &vol, double x, double y, double z,
                                       int SplineDegree)
{
    int SplineDegree_1 = SplineDegree - 1;
    double aux;

    // Logical to physical
    z -= STARTINGZ(vol);
    y -= STARTINGY(vol);
    x -= STARTINGX(vol);

    int l1 = CEIL(x - SplineDegree_1);
    int l2 = l1 + SplineDegree;

    int m1 = CEIL(y - SplineDegree_1);
    int m2 = m1 + SplineDegree;

    int n1 = CEIL(z - SplineDegree_1);
    int n2 = n1 + SplineDegree;

    double zyxsum = 0.0;
    int Zdim=(int)ZSIZE(vol);
    int Ydim=(int)YSIZE(vol);
    int Xdim=(int)XSIZE(vol);
    for (int n = n1; n <= n2; n++)
    {
        int equivalent_n=n;
        if      (n<0)
            equivalent_n=-n-1;
        else if (n>=Zdim)
            equivalent_n=2*Zdim-n-1;
        double yxsum = 0.0;
        for (int m = m1; m <= m2; m++)
        {
            int equivalent_m=m;
            if      (m<0)
                equivalent_m=-m-1;
            else if (m>=Ydim)
                equivalent_m=2*Ydim-m-1;
            double xsum = 0.0;
            for (int l = l1; l <= l2; l++)
            {
                double xminusl = x - (double) l;
                int equivalent_l=l;
                if      (l<0)
                    equivalent_l=-l-1;
                else if (l>=Xdim)
                    equivalent_l=2*Xdim-l-1;
                double Coeff = (double) DIRECT_A3D_ELEM(vol,
                                                        equivalent_n,equivalent_m,equivalent_l);
                double aux;
                switch (SplineDegree)
                {
                case 2:
                    xsum += Coeff * Bspline02(xminusl);
                    break;
                case 3:
                    BSPLINE03(aux,xminusl);
                    xsum += Coeff * aux;
                    break;
                case 4:
                    xsum += Coeff * Bspline04(xminusl);
                    break;
                case 5:
                    xsum += Coeff * Bspline05(xminusl);
                    break;
                case 6:
                    xsum += Coeff * Bspline06(xminusl);
                    break;
                case 7:
                    xsum += Coeff * Bspline07(xminusl);
                    break;
                case 8:
                    xsum += Coeff * Bspline08(xminusl);
                    break;
                case 9:
                    xsum += Coeff * Bspline09(xminusl);
                    break;
                }
            }

            double yminusm = y - (double) m;
            switch (SplineDegree)
            {
            case 2:
                yxsum += xsum * Bspline02(yminusm);
                break;
            case 3:
                BSPLINE03DIFF1(aux,yminusm);
                yxsum += xsum * aux;
                break;
            case 4:
                yxsum += xsum * Bspline04(yminusm);
                break;
            case 5:
                yxsum += xsum * Bspline05(yminusm);
                break;
            case 6:
                yxsum += xsum * Bspline06(yminusm);
                break;
            case 7:
                yxsum += xsum * Bspline07(yminusm);
                break;
            case 8:
                yxsum += xsum * Bspline08(yminusm);
                break;
            case 9:
                yxsum += xsum * Bspline09(yminusm);
                break;
            }
        }

        double zminusn = z - (double) n;
        switch (SplineDegree)
        {
        case 2:
            zyxsum += yxsum * Bspline02(zminusn);
            break;
        case 3:
            BSPLINE03(aux,zminusn);
            zyxsum += yxsum * aux;
            break;
        case 4:
            zyxsum += yxsum * Bspline04(zminusn);
            break;
        case 5:
            zyxsum += yxsum * Bspline05(zminusn);
            break;
        case 6:
            zyxsum += yxsum * Bspline06(zminusn);
            break;
        case 7:
            zyxsum += yxsum * Bspline07(zminusn);
            break;
        case 8:
            zyxsum += yxsum * Bspline08(zminusn);
            break;
        case 9:
            zyxsum += yxsum * Bspline09(zminusn);
            break;
        }
    }

    return zyxsum;
}

/** Interpolates the value of the 3D matrix M at the point (x,y,z) knowing
 * that this image is a set of B-spline coefficients. And making the diff
 * of z, such->  V=sum(Coef Bx By diff(Bz))
 * Only for BSplines of degree 3!!
 * @ingroup VolumesMemory
 *
 * (x,y,z) are in logical coordinates.
 */
double interpolatedElementBSplineDiffZ(MultidimArray<double> &vol, double x, double y, double z,
                                       int SplineDegree)
{
    int SplineDegree_1 = SplineDegree - 1;
    double aux;

    // Logical to physical
    z -= STARTINGZ(vol);
    y -= STARTINGY(vol);
    x -= STARTINGX(vol);

    int l1 = CEIL(x - SplineDegree_1);
    int l2 = l1 + SplineDegree;

    int m1 = CEIL(y - SplineDegree_1);
    int m2 = m1 + SplineDegree;

    int n1 = CEIL(z - SplineDegree_1);
    int n2 = n1 + SplineDegree;

    double zyxsum = 0.0;
    int Zdim=(int)ZSIZE(vol);
    int Ydim=(int)YSIZE(vol);
    int Xdim=(int)XSIZE(vol);
    for (int n = n1; n <= n2; n++)
    {
        int equivalent_n=n;
        if      (n<0)
            equivalent_n=-n-1;
        else if (n>=Zdim)
            equivalent_n=2*Zdim-n-1;
        double yxsum = 0.0;
        for (int m = m1; m <= m2; m++)
        {
            int equivalent_m=m;
            if      (m<0)
                equivalent_m=-m-1;
            else if (m>=Ydim)
                equivalent_m=2*Ydim-m-1;
            double xsum = 0.0;
            for (int l = l1; l <= l2; l++)
            {
                double xminusl = x - (double) l;
                int equivalent_l=l;
                if      (l<0)
                    equivalent_l=-l-1;
                else if (l>=Xdim)
                    equivalent_l=2*Xdim-l-1;
                double Coeff = (double) DIRECT_A3D_ELEM(vol,
                                                        equivalent_n,equivalent_m,equivalent_l);
                double aux;
                switch (SplineDegree)
                {
                case 2:
                    xsum += Coeff * Bspline02(xminusl);
                    break;
                case 3:
                    BSPLINE03(aux,xminusl);
                    xsum += Coeff * aux;
                    break;
                case 4:
                    xsum += Coeff * Bspline04(xminusl);
                    break;
                case 5:
                    xsum += Coeff * Bspline05(xminusl);
                    break;
                case 6:
                    xsum += Coeff * Bspline06(xminusl);
                    break;
                case 7:
                    xsum += Coeff * Bspline07(xminusl);
                    break;
                case 8:
                    xsum += Coeff * Bspline08(xminusl);
                    break;
                case 9:
                    xsum += Coeff * Bspline09(xminusl);
                    break;
                }
            }

            double yminusm = y - (double) m;
            switch (SplineDegree)
            {
            case 2:
                yxsum += xsum * Bspline02(yminusm);
                break;
            case 3:
                BSPLINE03(aux,yminusm);
                yxsum += xsum * aux;
                break;
            case 4:
                yxsum += xsum * Bspline04(yminusm);
                break;
            case 5:
                yxsum += xsum * Bspline05(yminusm);
                break;
            case 6:
                yxsum += xsum * Bspline06(yminusm);
                break;
            case 7:
                yxsum += xsum * Bspline07(yminusm);
                break;
            case 8:
                yxsum += xsum * Bspline08(yminusm);
                break;
            case 9:
                yxsum += xsum * Bspline09(yminusm);
                break;
            }
        }

        double zminusn = z - (double) n;
        switch (SplineDegree)
        {
        case 2:
            zyxsum += yxsum * Bspline02(zminusn);
            break;
        case 3:
            BSPLINE03DIFF1(aux,zminusn);
            zyxsum += yxsum * aux;
            break;
        case 4:
            zyxsum += yxsum * Bspline04(zminusn);
            break;
        case 5:
            zyxsum += yxsum * Bspline05(zminusn);
            break;
        case 6:
            zyxsum += yxsum * Bspline06(zminusn);
            break;
        case 7:
            zyxsum += yxsum * Bspline07(zminusn);
            break;
        case 8:
            zyxsum += yxsum * Bspline08(zminusn);
            break;
        case 9:
            zyxsum += yxsum * Bspline09(zminusn);
            break;
        }
    }

    return zyxsum;
}

void radiallySymmetrize(const MultidimArray<double>& img, MultidimArray<double> &radialImg)
{
	Matrix1D<int> center(2);
	center.initZeros();
	MultidimArray<int> distance, radial_count;
	MultidimArray<double> radial_mean;
	int dim;
	radialAveragePrecomputeDistance(img, center, distance, dim);
	fastRadialAverage(img, distance, dim, radial_mean, radial_count);

	radialImg.initZeros(img);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(radialImg)
    {
        int d=DIRECT_MULTIDIM_ELEM(distance,n);
        DIRECT_MULTIDIM_ELEM(radialImg,n)=A1D_ELEM(radial_mean,d);
    }
}
