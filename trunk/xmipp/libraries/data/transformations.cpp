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
void rotationMatrix(double ang, Matrix2D< double > &result, bool homogeneous)
{
    double cosine, sine;

    ang = DEG2RAD(ang);
    cosine = cos(ang);
    sine = sin(ang);

    if (homogeneous)
    {
    	result.initZeros(4,4);
    	MAT_ELEM(result,3,3)=1;
    }
    else
    	result.initZeros(3,3);
    MAT_ELEM(result,0, 0) = cosine;
    MAT_ELEM(result,0, 1) = -sine;

    MAT_ELEM(result,1, 0) = sine;
    MAT_ELEM(result,1, 1) = cosine;

    MAT_ELEM(result,2, 2) = 1;
}

/* Translation 2D ---------------------------------------------------------- */
void translationMatrix(Matrix2D< double > &result, double xshift, double yshift,
                       double zshift)
{
    result.initIdentity(4);
    MAT_ELEM(result,0, 3) = xshift;
    MAT_ELEM(result,1, 3) = yshift;
    MAT_ELEM(result,2, 3) = zshift;
}

/* Rotation 3D around the system axes -------------------------------------- */
void rotationMatrix(double ang, char axis, Matrix2D< double > &result, bool homogeneous)
{
	if (homogeneous)
	{
		result.initZeros(4,4);
		MAT_ELEM(result,3,3)=1;
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
        MAT_ELEM(result,0, 0) = cosine;
        MAT_ELEM(result,0, 1) = -sine;
        MAT_ELEM(result,1, 0) = sine;
        MAT_ELEM(result,1, 1) = cosine;
        MAT_ELEM(result,2, 2) = 1;
        break;
    case 'Y':
        MAT_ELEM(result,0, 0) = cosine;
        MAT_ELEM(result,0, 2) = -sine;
        MAT_ELEM(result,1, 1) = 1;
        MAT_ELEM(result,2, 0) = sine;
        MAT_ELEM(result,2, 2) = cosine;
        break;
    case 'X':
        MAT_ELEM(result,0, 0) = 1;
        MAT_ELEM(result,1, 1) = cosine;
        MAT_ELEM(result,1, 2) = -sine;
        MAT_ELEM(result,2, 1) = sine;
        MAT_ELEM(result,2, 2) = cosine;
        break;
    default:
        REPORT_ERROR(ERR_VALUE_INCORRECT, "rotation3DMatrix: Unknown axis");
    }
}

/* Align a vector with Z axis */
void alignWithZ(const Matrix1D<double> &axis, Matrix2D<double>& result,
		bool homogeneous)
{
#ifndef RELEASE_MODE
    if (axis.size() != 3)
        REPORT_ERROR(ERR_MATRIX_SIZE, "alignWithZ: Axis is not in R3");
#endif

    if (homogeneous)
    {
    	result.initZeros(4,4);
    	MAT_ELEM(result,3, 3) = 1;
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
        MAT_ELEM(result,0, 0) = proj_mod;
        MAT_ELEM(result,0, 1) = -XX(Axis) * YY(Axis) / proj_mod;
        MAT_ELEM(result,0, 2) = -XX(Axis) * ZZ(Axis) / proj_mod;
        MAT_ELEM(result,1, 0) = 0;
        MAT_ELEM(result,1, 1) = ZZ(Axis) / proj_mod;
        MAT_ELEM(result,1, 2) = -YY(Axis) / proj_mod;
        MAT_ELEM(result,2, 0) = XX(Axis);
        MAT_ELEM(result,2, 1) = YY(Axis);
        MAT_ELEM(result,2, 2) = ZZ(Axis);
    }
    else
    {
        // I know that the Axis is the X axis
        MAT_ELEM(result,0, 0) = 0;
        MAT_ELEM(result,0, 1) = 0;
        MAT_ELEM(result,0, 2) = -1;
        MAT_ELEM(result,1, 0) = 0;
        MAT_ELEM(result,1, 1) = 1;
        MAT_ELEM(result,1, 2) = 0;
        MAT_ELEM(result,2, 0) = 1;
        MAT_ELEM(result,2, 1) = 0;
        MAT_ELEM(result,2, 2) = 0;
    }
}

/* Rotation 3D around any axis -------------------------------------------- */
void rotationMatrix(double ang, const Matrix1D<double> &axis, Matrix2D<double> &result,
		bool homogeneous)
{
    // Compute a matrix which makes the turning axis coincident with Z
    // And turn around this axis
    Matrix2D<double> A,R;
    alignWithZ(axis,A,homogeneous);
    rotationMatrix(ang, 'Z', R,homogeneous);
    result=A.transpose() * R * A;
}

/* Scale 3D ---------------------------------------------------------------- */
void scaleMatrix(Matrix2D< double > &result, double scaleX, double scaleY, double scaleZ)
{
    result.initZeros(4,4);
    MAT_ELEM(result,3, 3) = 1;
    MAT_ELEM(result,0, 0) = scaleX;
    MAT_ELEM(result,1, 1) = scaleY;
    MAT_ELEM(result,2, 2) = scaleZ;
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
        applyGeometry(SplineDegree, V2, V1, A, inv, wrap, outside);
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

// Special case for complex arrays
void produceSplineCoefficients(int SplineDegree,
                               MultidimArray< double > &coeffs,
                               const MultidimArray< std::complex<double> > &V1)
{
    // TODO Implement
    REPORT_ERROR(ERR_NOT_IMPLEMENTED,"Spline coefficients of a complex matrix is not implemented.");
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
        REPORT_ERROR(ERR_UNCLASSIFIED, "Error in ImageFromSplineCoefficients...");

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

/** Interpolates the value of the 3D matrix M at the point (x,y,z) knowing
* that this image is a set of B-spline coefficients. And making the diff
* of x, such->  V=sum(Coef diff(Bx) By Bz)
* ��Only for BSplines of degree 3!!
* @ingroup VolumesMemory
*
* (x,y,z) are in logical coordinates.
*/
double interpolatedElementBSplineDiffX(MultidimArray<double> &vol, double x, double y, double z,
                                       int SplineDegree)
{
    int SplineDegree_1 = SplineDegree - 1;

    // Logical to physical
    z -= STARTINGZ(vol);
    y -= STARTINGY(vol);
    x -= STARTINGX(vol);

    int lmax = XSIZE(vol);
    int mmax = YSIZE(vol);
    int nmax = ZSIZE(vol);

    int l1 = CEIL(x - SplineDegree_1);
    int l2 = l1 + SplineDegree;

    int m1 = CEIL(y - SplineDegree_1);
    int m2 = m1 + SplineDegree;

    int n1 = CEIL(z - SplineDegree_1);
    int n2 = n1 + SplineDegree;

    double zyxsum = 0.0;
    for (int n = n1; n <= n2; n++)
    {
        int equivalent_n=n;
        if      (n<0)
            equivalent_n=-n-1;
        else if (n>=ZSIZE(vol))
            equivalent_n=2*ZSIZE(vol)-n-1;
        double yxsum = 0.0;
        for (int m = m1; m <= m2; m++)
        {
            int equivalent_m=m;
            if      (m<0)
                equivalent_m=-m-1;
            else if (m>=YSIZE(vol))
                equivalent_m=2*YSIZE(vol)-m-1;
            double xsum = 0.0;
            for (int l = l1; l <= l2; l++)
            {
                double xminusl = x - (double) l;
                int equivalent_l=l;
                if      (l<0)
                    equivalent_l=-l-1;
                else if (l>=XSIZE(vol))
                    equivalent_l=2*XSIZE(vol)-l-1;
                double Coeff = (double) DIRECT_A3D_ELEM(vol,
                                                        equivalent_n,equivalent_m,equivalent_l);
                switch (SplineDegree)
                {
                case 2:
                    xsum += Coeff * Bspline02(xminusl);
                    break;
                case 3:
                    xsum += Coeff * Bspline03Diff1(xminusl);
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
                yxsum += xsum * Bspline03(yminusm);
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
            zyxsum += yxsum * Bspline03(zminusn);
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
 * ��Only for BSplines of degree 3!!
 * @ingroup VolumesMemory
 *
 * (x,y,z) are in logical coordinates.
 */
double interpolatedElementBSplineDiffY(MultidimArray<double> &vol, double x, double y, double z,
                                       int SplineDegree)
{
    int SplineDegree_1 = SplineDegree - 1;

    // Logical to physical
    z -= STARTINGZ(vol);
    y -= STARTINGY(vol);
    x -= STARTINGX(vol);

    int lmax = XSIZE(vol);
    int mmax = YSIZE(vol);
    int nmax = ZSIZE(vol);

    int l1 = CEIL(x - SplineDegree_1);
    int l2 = l1 + SplineDegree;

    int m1 = CEIL(y - SplineDegree_1);
    int m2 = m1 + SplineDegree;

    int n1 = CEIL(z - SplineDegree_1);
    int n2 = n1 + SplineDegree;

    double zyxsum = 0.0;
    for (int n = n1; n <= n2; n++)
    {
        int equivalent_n=n;
        if      (n<0)
            equivalent_n=-n-1;
        else if (n>=ZSIZE(vol))
            equivalent_n=2*ZSIZE(vol)-n-1;
        double yxsum = 0.0;
        for (int m = m1; m <= m2; m++)
        {
            int equivalent_m=m;
            if      (m<0)
                equivalent_m=-m-1;
            else if (m>=YSIZE(vol))
                equivalent_m=2*YSIZE(vol)-m-1;
            double xsum = 0.0;
            for (int l = l1; l <= l2; l++)
            {
                double xminusl = x - (double) l;
                int equivalent_l=l;
                if      (l<0)
                    equivalent_l=-l-1;
                else if (l>=XSIZE(vol))
                    equivalent_l=2*XSIZE(vol)-l-1;
                double Coeff = (double) DIRECT_A3D_ELEM(vol,
                                                        equivalent_n,equivalent_m,equivalent_l);
                switch (SplineDegree)
                {
                case 2:
                    xsum += Coeff * Bspline02(xminusl);
                    break;
                case 3:
                    xsum += Coeff * Bspline03(xminusl);
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
                yxsum += xsum * Bspline03Diff1(yminusm);
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
            zyxsum += yxsum * Bspline03(zminusn);
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
 * ��Only for BSplines of degree 3!!
 * @ingroup VolumesMemory
 *
 * (x,y,z) are in logical coordinates.
 */
double interpolatedElementBSplineDiffZ(MultidimArray<double> &vol, double x, double y, double z,
                                       int SplineDegree)
{
    int SplineDegree_1 = SplineDegree - 1;

    // Logical to physical
    z -= STARTINGZ(vol);
    y -= STARTINGY(vol);
    x -= STARTINGX(vol);

    int lmax = XSIZE(vol);
    int mmax = YSIZE(vol);
    int nmax = ZSIZE(vol);

    int l1 = CEIL(x - SplineDegree_1);
    int l2 = l1 + SplineDegree;

    int m1 = CEIL(y - SplineDegree_1);
    int m2 = m1 + SplineDegree;

    int n1 = CEIL(z - SplineDegree_1);
    int n2 = n1 + SplineDegree;

    double zyxsum = 0.0;
    for (int n = n1; n <= n2; n++)
    {
        int equivalent_n=n;
        if      (n<0)
            equivalent_n=-n-1;
        else if (n>=ZSIZE(vol))
            equivalent_n=2*ZSIZE(vol)-n-1;
        double yxsum = 0.0;
        for (int m = m1; m <= m2; m++)
        {
            int equivalent_m=m;
            if      (m<0)
                equivalent_m=-m-1;
            else if (m>=YSIZE(vol))
                equivalent_m=2*YSIZE(vol)-m-1;
            double xsum = 0.0;
            for (int l = l1; l <= l2; l++)
            {
                double xminusl = x - (double) l;
                int equivalent_l=l;
                if      (l<0)
                    equivalent_l=-l-1;
                else if (l>=XSIZE(vol))
                    equivalent_l=2*XSIZE(vol)-l-1;
                double Coeff = (double) DIRECT_A3D_ELEM(vol,
                                                        equivalent_n,equivalent_m,equivalent_l);
                switch (SplineDegree)
                {
                case 2:
                    xsum += Coeff * Bspline02(xminusl);
                    break;
                case 3:
                    xsum += Coeff * Bspline03(xminusl);
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
                yxsum += xsum * Bspline03(yminusm);
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
            zyxsum += yxsum * Bspline03Diff1(zminusn);
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
