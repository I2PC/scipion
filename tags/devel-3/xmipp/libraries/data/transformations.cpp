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


/* Rotation 2D ------------------------------------------------------------- */
Matrix2D<double> rotation2DMatrix(double ang)
{
    Matrix2D<double> result(3, 3);
    rotation2DMatrix(ang, result);
    return result;
}

void rotation2DMatrix(double ang, Matrix2D< double > &result)
{
    double cosine, sine;

    ang = DEG2RAD(ang);
    cosine = cos(ang);
    sine = sin(ang);

    result(0, 0) = cosine;
    result(0, 1) = -sine;
    result(0, 2) = 0;

    result(1, 0) = sine;
    result(1, 1) = cosine;
    result(1, 2) = 0;

    result(2, 0) = 0;
    result(2, 1) = 0;
    result(2, 2) = 1;
}

/* Translation 2D ---------------------------------------------------------- */
Matrix2D<double> translation2DMatrix(const Matrix1D<double> &v)
{
    if (v.size() != 2)
        REPORT_ERROR(1002, "Translation2D_matrix: vector is not in R2");

    Matrix2D<double> result(3, 3);

    result.initIdentity();
    result(0, 2) = XX(v);
    result(1, 2) = YY(v);

    return result;
}

/* Rotation 3D around the system axes -------------------------------------- */
Matrix2D<double> rotation3DMatrix(double ang, char axis)
{
    Matrix2D<double> result(4, 4);
    double cosine, sine;

    ang = DEG2RAD(ang);
    cosine = cos(ang);
    sine = sin(ang);

    result.initZeros();
    result(3, 3) = 1;
    switch (axis)
    {
    case 'Z':
        result(0, 0) = cosine;
        result(0, 1) = -sine;
        result(1, 0) = sine;
        result(1, 1) = cosine;
        result(2, 2) = 1;
        break;
    case 'Y':
        result(0, 0) = cosine;
        result(0, 2) = -sine;
        result(2, 0) = sine;
        result(2, 2) = cosine;
        result(1, 1) = 1;
        break;
    case 'X':
        result(1, 1) = cosine;
        result(1, 2) = -sine;
        result(2, 1) = sine;
        result(2, 2) = cosine;
        result(0, 0) = 1;
        break;
    default:
        REPORT_ERROR(1105, "rotation3DMatrix: Unknown axis");
    }
    return result;
}

/* Align a vector with Z axis */
Matrix2D<double> alignWithZ(const Matrix1D<double> &axis)
{
    Matrix1D<double>  Axis;
    Matrix2D<double>  A(4, 4);

    if (axis.size() != 3)
        REPORT_ERROR(1002, "alignWithZ: Axis is not in R3");

    // Copy axis and compute length of the projection on YZ plane
    Axis = axis;
    Axis.selfNormalize();
    double proj_mod = sqrt(YY(Axis) * YY(Axis) + ZZ(Axis) * ZZ(Axis));

    A(3, 3) = 1;
    if (proj_mod > XMIPP_EQUAL_ACCURACY)
    { // proj_mod!=0
        // Build Matrix A, which makes the turning axis coincident with Z
        A(0, 0) = proj_mod;
        A(0, 1) = -XX(Axis) * YY(Axis) / proj_mod;
        A(0, 2) = -XX(Axis) * ZZ(Axis) / proj_mod;
        A(1, 0) = 0;
        A(1, 1) = ZZ(Axis) / proj_mod;
        A(1, 2) = -YY(Axis) / proj_mod;
        A(2, 0) = XX(Axis);
        A(2, 1) = YY(Axis);
        A(2, 2) = ZZ(Axis);
    }
    else
    {
        // I know that the Axis is the X axis
        A(0, 0) = 0;
        A(0, 1) = 0;
        A(0, 2) = -1;
        A(1, 0) = 0;
        A(1, 1) = 1;
        A(1, 2) = 0;
        A(2, 0) = 1;
        A(2, 1) = 0;
        A(2, 2) = 0;
    }
    return A;
}

/* Rotation 3D around any axis -------------------------------------------- */
Matrix2D<double> rotation3DMatrix(double ang, const Matrix1D<double> &axis)
{
    // Compute a matrix which makes the turning axis coincident with Z
    // And turn around this axis
    Matrix2D<double> A = alignWithZ(axis);
    return A.transpose() * rotation3DMatrix(ang, 'Z') * A;
}

/* Translation 3D ---------------------------------------------------------- */
Matrix2D<double> translation3DMatrix(const Matrix1D<double> &v)
{
    if (v.size() != 3)
        REPORT_ERROR(1002, "Translation3D_matrix: vector is not in R3");

    Matrix2D<double> result(4, 4);

    result.initIdentity();
    result(0, 3) = XX(v);
    result(1, 3) = YY(v);
    result(2, 3) = ZZ(v);

    return result;
}

/* Scale 3D ---------------------------------------------------------------- */
Matrix2D<double> scale3DMatrix(const Matrix1D<double> &sc)
{
    if (sc.size() != 3)
        REPORT_ERROR(1002, "Scale3D_matrix: vector is not in R3");

    Matrix2D<double> result(4, 4);

    result.initIdentity();
    result(0, 0) = XX(sc);
    result(1, 1) = YY(sc);
    result(2, 2) = ZZ(sc);

    return result;
}


// Special case for complex numbers
void applyGeometry(int SplineDegree, 
                   MultidimArray< std::complex<double> > &V2,
                   const MultidimArray< std::complex<double> > &V1,
                   const Matrix2D<double> &A, bool inv, 
                   bool wrap, std::complex<double> outside, unsigned long n)
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
        V1.getImage(n, oneImg);
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
        applyGeometry(SplineDegree, V2, V1, A, inv, wrap, outside, n);
        

}

// Special case for complex numbers
void selfApplyGeometry(int Splinedegree, 
                       MultidimArray< std::complex<double> > &V1,
                       const Matrix2D<double> &A, bool inv, 
                       bool wrap, std::complex<double> outside, unsigned long n)
{
    MultidimArray<std::complex<double> > aux = V1;
    applyGeometry(Splinedegree, V1, aux, A, inv, wrap, outside, n);
}

// Special case for complex arrays
void produceSplineCoefficients(int SplineDegree, 
                               MultidimArray< double > &coeffs,
                               const MultidimArray< std::complex<double> > &V1,  
                               unsigned long n)
{
    // TODO Implement
    REPORT_ERROR(222,"Spline coefficients of a complex matrix is not implemented.");
}

void produceImageFromSplineCoefficients(int SplineDegree, 
                                        MultidimArray< double >& img, 
                                        const MultidimArray< double > &coeffs)
{

    img.initZeros(ZSIZE(coeffs), YSIZE(coeffs), XSIZE(coeffs));
    STARTINGX(img) = STARTINGX(coeffs);
    STARTINGY(img) = STARTINGY(coeffs);
    STARTINGZ(img) = STARTINGZ(coeffs);

    int Status;
    MultidimArray< double > aux;
    typeCast(coeffs, aux); // This will create a single volume!
        
    ChangeBasisVolume(MULTIDIM_ARRAY(aux), MULTIDIM_ARRAY(img),
                      XSIZE(coeffs), YSIZE(coeffs), ZSIZE(coeffs),
                      BasicSpline, CardinalSpline, SplineDegree,
                      MirrorOnBounds, DBL_EPSILON, &Status);
    if (Status)
        REPORT_ERROR(1201, "Error in ImageFromSplineCoefficients...");
    
}


// Special case for complex arrays
void scaleToSize(int SplineDegree, 
                 MultidimArray< std::complex<double> > &V2,
                 const MultidimArray< std::complex<double> > &V1,
                 int Xdim, int Ydim, int Zdim,
                 unsigned long n)
{
    if (SplineDegree > 1)
    {
        MultidimArray< double > re, im, aux;
        MultidimArray<std::complex<double> > oneImg;

        re.resize(ZSIZE(V1), YSIZE(V1), XSIZE(V1));
        im.resize(ZSIZE(V1), YSIZE(V1), XSIZE(V1));

        V1.getImage(n, oneImg);
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
        scaleToSize(SplineDegree, V2, V1, Xdim, Ydim, Zdim, n);
    
}
