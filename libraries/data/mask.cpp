/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "mask.h"
#include "xmipp_program.h"

/*---------------------------------------------------------------------------*/
/* Multidim Masks                                                                  */
/*---------------------------------------------------------------------------*/
void RaisedCosineMask(MultidimArray<double> &mask,
                      double r1, double r2, int mode, double x0, double y0, double z0)
{
    double K = PI / (r2 - r1);
    double r1_2=r1*r1;
    double r2_2=r2*r2;
    for (int k=STARTINGZ(mask); k<=FINISHINGZ(mask); ++k)
    {
    	double k2=(k - z0) * (k - z0);
        for (int i=STARTINGY(mask); i<=FINISHINGY(mask); ++i)
        {
        	double i2_k2=k2+(i - y0) * (i - y0);
            for (int j=STARTINGX(mask); j<=FINISHINGX(mask); ++j)
            {
            	double r2=i2_k2+(j - x0) * (j - x0);
                if (r2 <= r1_2)
                    A3D_ELEM(mask, k, i, j) = 1;
                else if (r2 < r2_2)
                {
                    double r=sqrt(r2);
                    A3D_ELEM(mask, k, i, j) = 0.5*(1 + cos(K * (r - r1)));
                }
                else
                    A3D_ELEM(mask, k, i, j) = 0;
                if (mode == OUTSIDE_MASK)
                    A3D_ELEM(mask, k, i, j) = 1 - A3D_ELEM(mask, k, i, j);
            }
        }
    }
}

void RaisedCrownMask(MultidimArray<double> &mask,
                     double r1, double r2, double pix_width, int mode,
                     double x0, double y0, double z0)
{
    RaisedCosineMask(mask, r1 - pix_width, r1 + pix_width, OUTSIDE_MASK, x0, y0, z0);
    MultidimArray<double> aux;
    aux.resize(mask);
    RaisedCosineMask(aux, r2 - pix_width, r2 + pix_width, INNER_MASK, x0, y0, z0);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        A3D_ELEM(mask, k, i, j) *= A3D_ELEM(aux, k, i, j);
        if (mode == OUTSIDE_MASK)
            A3D_ELEM(mask, k, i, j) = 1 - A3D_ELEM(mask, k, i, j);
    }
}

void KaiserMask(MultidimArray<double> &mask, double delta, double Deltaw)
{
    // Convert Deltaw from a frequency normalized to 1, to a freq. normalized to PI
    Deltaw *= PI;

    // Design Kaiser selfWindow
    double A = -20 * log10(delta);
    int    M = CEIL((A - 8) / (2.285 * Deltaw));
    double beta;
    if (A > 50)
        beta = 0.1102 * (A - 8.7);
    else if (A >= 21)
        beta = 0.5842 * pow(A - 21, 0.4) + 0.07886 * (A - 21);
    else
        beta = 0;

    // "Draw" Kaiser selfWindow
    if (YSIZE(mask)==1 && ZSIZE(mask)==1)
    {
        // 1D
        mask.resize(2*M + 1);
    }
    else if (ZSIZE(mask)==1)
    {
        // 2D
        mask.resize(2*M + 1, 2*M + 1);
    }
    else
    {
        // 3D
        mask.resize(2*M + 1, 2*M + 1, 2*M + 1);
    }

    mask.setXmippOrigin();
    double iI0Beta = 1.0 / bessi0(beta);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        double r = sqrt((double)(i * i + j * j + k * k));
        if (r <= M)
            A3D_ELEM(mask, k, i, j) = bessi0(beta * sqrt(1 - (r / M) * (r / M))) * iI0Beta;
    }

}

void SincMask(MultidimArray<double> &mask,
              double omega, int mode, double x0, double y0, double z0)
{
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        double r = sqrt( (k - z0) * (k - z0) + (i - y0) * (i - y0) + (j - x0) * (j - x0) );
        A3D_ELEM(mask, k, i, j) = omega/PI * SINC(omega/PI * r);
        if (mode == OUTSIDE_MASK)
            A3D_ELEM(mask, k, i, j) = 1 - A3D_ELEM(mask, k, i, j);
    }
}

void SincKaiserMask(MultidimArray<double> &mask,
                    double omega, double delta, double Deltaw)
{
    MultidimArray<double> kaiser;
    kaiser.resizeNoCopy(mask);
    KaiserMask(kaiser, delta, Deltaw);
    mask.resize(kaiser);
    mask.setXmippOrigin();
    SincMask(mask, omega*PI, INNER_MASK);
    mask *= kaiser;
}

void BlackmanMask(MultidimArray<double> &mask, int mode,
                  double x0, double y0, double z0)
{
    double Xdim2 = XMIPP_MAX(1, (XSIZE(mask) - 1) * (XSIZE(mask) - 1));
    double Ydim2 = XMIPP_MAX(1, (YSIZE(mask) - 1) * (YSIZE(mask) - 1));
    double Zdim2 = XMIPP_MAX(1, (ZSIZE(mask) - 1) * (ZSIZE(mask) - 1));
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        double r = sqrt((k - z0) * (k - z0) / Zdim2 + (i - y0) * (i - y0) / Xdim2 + (j - x0) * (j - x0) / Ydim2);
        A3D_ELEM(mask, k, i, j) = 0.42 + 0.5 * cos(2 * PI * r) + 0.08 * cos(4 * PI * r);
        if (mode == OUTSIDE_MASK)
            A3D_ELEM(mask, k, i, j) = 1 - A3D_ELEM(mask, k, i, j);
    }
}

void SincBlackmanMask(MultidimArray<double> &mask,
                      double omega, double power_percentage, int mode,
                      double x0, double y0, double z0)
{
    MultidimArray<double> blackman;

    int N = CEIL(1 / omega * CEIL(-1 / 2 + 1 / (PI * (1 - power_percentage / 100))));

    if (ZSIZE(mask)==1)
    {
        // 2D
        mask.resize(N, N);
        blackman.resize(N, N);
    }
    else
    {
        // 3D
        mask.resize(N, N, N);
        blackman.resize(N, N, N);
    }
    mask.setXmippOrigin();
    SincMask(mask,omega,INNER_MASK,x0,y0,z0);
    blackman.setXmippOrigin();
    BlackmanMask(blackman);
    mask *= blackman;
}

void BinaryCircularMask(MultidimArray<int> &mask,
                        double radius, int mode, double x0, double y0, double z0)
{
    mask.initZeros();
    double radius2 = radius * radius;
    for (int k=STARTINGZ(mask); k<=FINISHINGZ(mask); ++k)
    {
        double diff=k - z0;
        double z2=diff*diff;
        for (int i=STARTINGY(mask); i<=FINISHINGY(mask); ++i)
        {
            diff=i - y0;
            double z2y2=z2+diff*diff;
            for (int j=STARTINGX(mask); j<=FINISHINGX(mask); ++j)
            {
                diff=j - x0;
                double r2 = z2y2+diff*diff;
                if (r2 <= radius2 && mode == INNER_MASK)
                    A3D_ELEM(mask, k, i, j) = 1;
                else if (r2 >= radius2 && mode == OUTSIDE_MASK)
                    A3D_ELEM(mask, k, i, j) = 1;
            }
        }
    }
}

void BlobCircularMask(MultidimArray<double> &mask,
                      double r1, blobtype blob, int mode,
                      double x0, double y0, double z0)
{
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        double r = sqrt((k - z0) * (k - z0) + (i - y0) * (i - y0) + (j - x0) * (j - x0));
        if (mode == INNER_MASK)
        {
            if (r <= r1)
                A3D_ELEM(mask, k, i, j) = 1;
            else
                A3D_ELEM(mask, k, i, j) = blob_val(r-r1, blob);
        }
        else
        {
            if (r >= r1)
                A3D_ELEM(mask, k, i, j) = 1;
            else
                A3D_ELEM(mask, k, i, j) = blob_val(r1-r, blob);
        }
    }

}

void BinaryCrownMask(MultidimArray<int> &mask,
                     double R1, double R2, int mode, double x0, double y0, double z0)
{
    mask.initZeros();
    double R12 = R1 * R1;
    double R22 = R2 * R2;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        double r2 = (k - z0) * (k - z0) + (i - y0) * (i - y0) + (j - x0) * (j - x0);
        bool in_crown = (r2 >= R12 && r2 <= R22);
        if (in_crown  && mode == INNER_MASK)
            A3D_ELEM(mask, k, i, j) = 1;
        else if (!in_crown && mode == OUTSIDE_MASK)
            A3D_ELEM(mask, k, i, j) = 1;
    }
}

void BinaryTubeMask(MultidimArray<int> &mask,
                     double R1, double R2, double H, int mode, double x0, double y0, double z0)
{
    mask.initZeros();
    double R12 = R1 * R1;
    double R22 = R2 * R2;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        double r2 = (i - y0) * (i - y0) + (j - x0) * (j - x0);
        bool in_tube = (r2 >= R12 && r2 <= R22 && ABS(k)<H);
        if (in_tube  && mode == INNER_MASK)
            A3D_ELEM(mask, k, i, j) = 1;
        else if (!in_tube && mode == OUTSIDE_MASK)
            A3D_ELEM(mask, k, i, j) = 1;
    }
}

void BlobCrownMask(MultidimArray<double> &mask,
                   double r1, double r2, blobtype blob, int mode,
                   double x0, double y0, double z0)
{
    MultidimArray<double> aux;
    aux.resize(mask);
    if (mode == INNER_MASK)
    {
        BlobCircularMask(mask, r1, blob,
                         OUTSIDE_MASK, x0, y0, z0);
        BlobCircularMask(aux, r2, blob,
                         INNER_MASK, x0, y0, z0);
        FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
        {
            A3D_ELEM(mask, k, i, j) *= A3D_ELEM(aux, k, i, j);
        }
    }
    else
    {
        BlobCircularMask(mask, r1, blob,
                         INNER_MASK, x0, y0, z0);
        BlobCircularMask(aux, r2, blob,
                         OUTSIDE_MASK, x0, y0, z0);
        FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
        {
            A3D_ELEM(mask, k, i, j) += A3D_ELEM(aux, k, i, j);
        }
    }

}

void BinaryFrameMask(MultidimArray<int> &mask,
                     int Xrect, int Yrect, int Zrect, int mode, double x0, double y0, double z0)
{
    mask.initZeros();
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        bool in_frame =
            (j >= x0 + FIRST_XMIPP_INDEX(Xrect)) && (j <= x0 + LAST_XMIPP_INDEX(Xrect)) &&
            (i >= y0 + FIRST_XMIPP_INDEX(Yrect)) && (i <= y0 + LAST_XMIPP_INDEX(Yrect)) &&
            (k >= z0 + FIRST_XMIPP_INDEX(Zrect)) && (k <= z0 + LAST_XMIPP_INDEX(Zrect));

        if ((in_frame  && mode == INNER_MASK) || (!in_frame && mode == OUTSIDE_MASK))
            A3D_ELEM(mask, k, i, j) = 1;
    }
}

void GaussianMask(MultidimArray<double> &mask,
                  double sigma, int mode, double x0, double y0, double z0)
{
    double sigma2 = sigma * sigma;

    double K = 1. / pow(sqrt(2.*PI)*sigma,(double)mask.getDim());

    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        double r2 = (k - z0) * (k - z0) + (i - y0) * (i - y0) + (j - x0) * (j - x0);
        A3D_ELEM(mask, k, i, j) = K * exp(-0.5 * r2 / sigma2);

        if (mode == OUTSIDE_MASK)
            A3D_ELEM(mask, k, i, j) = 1 - A3D_ELEM(mask, k, i, j);
    }
}

/*---------------------------------------------------------------------------*/
/* 2D Masks                                                                  */
/*---------------------------------------------------------------------------*/

#define DWTCIRCULAR2D_BLOCK(s,quadrant) \
    SelectDWTBlock(s, mask, quadrant, \
                   XX(corner1),XX(corner2),YY(corner1),YY(corner2)); \
    V2_PLUS_V2(center,corner1,corner2); \
    V2_BY_CT(center,center,0.5); \
    FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(corner1,corner2) { \
        double r2=(XX(r)-XX(center))*(XX(r)-XX(center))+ \
                  (YY(r)-YY(center))*(YY(r)-YY(center)); \
        A2D_ELEM(mask,YY(r),XX(r))=(r2<=radius2); \
    }
void BinaryDWTCircularMask2D(MultidimArray<int> &mask, double radius,
                             int smin, int smax, const std::string &quadrant)
{
    double radius2 = radius * radius / (4 * (smin + 1));
    mask.initZeros();
    for (int s = smin; s <= smax; s++)
    {
        Matrix1D<int> corner1(2), corner2(2), r(2);
        Matrix1D<double> center(2);
        if (quadrant == "xx")
        {
            DWTCIRCULAR2D_BLOCK(s, "01");
            DWTCIRCULAR2D_BLOCK(s, "10");
            DWTCIRCULAR2D_BLOCK(s, "11");
        }
        else
            DWTCIRCULAR2D_BLOCK(s, quadrant);
        radius2 /= 4;
    }
}

void SeparableSincKaiserMask2D(MultidimArray<double> &mask,
                               double omega, double delta, double Deltaw)
{
    // Convert Deltaw from a frequency normalized to 1, to a freq. normalized to PI
    Deltaw *= PI;
    omega *= PI;

    // Design Kaiser selfWindow
    double A = -20 * log10(delta);
    double M = CEIL((A - 8) / (2.285 * Deltaw));
    double beta;
    if (A > 50)
        beta = 0.1102 * (A - 8.7);
    else if (A >= 21)
        beta = 0.5842 * pow(A - 21, 0.4) + 0.07886 * (A - 21);
    else
        beta = 0;

    // "Draw" Separable Kaiser Sinc selfWindow
    mask.resize((size_t)(2*M + 1),(size_t)(2*M + 1));
    mask.setXmippOrigin();
    double iI0Beta = 1.0 / bessi0(beta);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
    {
        mask(i, j) = omega/PI * SINC(omega/PI * i) * omega/PI * SINC(omega/PI * j) *
                     bessi0(beta * sqrt((double)(1 - (i / M) * (i / M)))) * iI0Beta *
                     bessi0(beta * sqrt((double)(1 - (j / M) * (j / M)))) * iI0Beta;
    }
}

void mask2D_4neig(MultidimArray<int> &mask, int value, int center)
{
    mask.resize(3, 3);
    mask.initZeros();
    mask(0, 1) = mask(1, 0) = mask(1, 2) = mask(2, 1) = value;
    mask(1, 1) = center;

}
void mask2D_8neig(MultidimArray<int> &mask, int value1, int value2, int center)
{
    mask.resize(3, 3);
    mask.initZeros();
    mask(0, 1) = mask(1, 0) = mask(1, 2) = mask(2, 1) = value1;
    mask(0, 0) = mask(0, 2) = mask(2, 0) = mask(2, 2) = value2;
    mask(1, 1) = center;

}

/*---------------------------------------------------------------------------*/
/* 3D Masks                                                                  */
/*---------------------------------------------------------------------------*/

#define DWTSPHERICALMASK_BLOCK(s,quadrant) \
    SelectDWTBlock(s, mask, quadrant, \
                   XX(corner1),XX(corner2),YY(corner1),YY(corner2), \
                   ZZ(corner1),ZZ(corner2)); \
    V3_PLUS_V3(center,corner1,corner2); \
    V3_BY_CT(center,center,0.5); \
    FOR_ALL_ELEMENTS_IN_ARRAY3D_BETWEEN(corner1,corner2) { \
        double r2=(XX(r)-XX(center))*(XX(r)-XX(center))+ \
                  (YY(r)-YY(center))*(YY(r)-YY(center))+ \
                  (ZZ(r)-ZZ(center))*(ZZ(r)-ZZ(center)); \
        A3D_ELEM(mask,ZZ(r),YY(r),XX(r))=(r2<=radius2); \
    }
void BinaryDWTSphericalMask3D(MultidimArray<int> &mask, double radius,
                              int smin, int smax, const std::string &quadrant)
{
    mask.initZeros();
    double radius2 = radius * radius / (4 * (smin + 1));
    for (int s = smin; s <= smax; s++)
    {
        Matrix1D<int> corner1(3), corner2(3), r(3);
        Matrix1D<double> center(3);
        if (quadrant == "xxx")
        {
            DWTSPHERICALMASK_BLOCK(s, "001");
            DWTSPHERICALMASK_BLOCK(s, "010");
            DWTSPHERICALMASK_BLOCK(s, "011");
            DWTSPHERICALMASK_BLOCK(s, "100");
            DWTSPHERICALMASK_BLOCK(s, "101");
            DWTSPHERICALMASK_BLOCK(s, "110");
            DWTSPHERICALMASK_BLOCK(s, "111");
        }
        else
            DWTSPHERICALMASK_BLOCK(s, quadrant);
        radius2 /= 4;
    }
}


void BinaryCylinderMask(MultidimArray<int> &mask,
                        double R, double H, int mode, double x0, double y0, double z0)
{
    mask.initZeros();
    double R2 = R * R;
    double H_2 = H / 2;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        double r2 = (i - y0) * (i - y0) + (j - x0) * (j - x0);
        int in_cyilinder = (r2 <= R2 && ABS(k - z0) <= H_2);
        if (in_cyilinder  && mode == INNER_MASK)
            A3D_ELEM(mask, k, i, j) = 1;
        else if (!in_cyilinder && mode == OUTSIDE_MASK)
            A3D_ELEM(mask, k, i, j) = 1;
    }
}

void BinaryConeMask(MultidimArray<int> &mask, double theta, int mode,bool centerOrigin)
{

    int halfX, halfY,halfZ;
    int minX,minY,minZ;
    int maxX,maxY,maxZ;
    int ipp,jpp,kpp;

    halfX = mask.xdim/2;
    halfY = mask.ydim/2;
    halfZ = mask.zdim/2;

    minX = -halfX;
    minY = -halfY;
    minZ = -halfZ;

    maxX = (int)((mask.xdim-0.5)/2);
    maxY = (int)((mask.ydim-0.5)/2);
    maxZ = (int)((mask.zdim-0.5)/2);

    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        if (centerOrigin)
        {
            kpp = intWRAP (k+halfZ,minZ,maxZ);
            ipp = intWRAP (i+halfY,minY,maxY);
            jpp = intWRAP (j+halfX,minX,maxX);
        }
        else
        {
            kpp=k;
            ipp=i;
            jpp=j;
        }

        double rad = tan(PI * theta / 180.) * (double)k;
        if ((double)(i*i + j*j) < (rad*rad))
            A3D_ELEM(mask, kpp, ipp, jpp) = 0;
        else
            A3D_ELEM(mask, kpp, ipp, jpp) = 1;
        if (mode == OUTSIDE_MASK)
            A3D_ELEM(mask, kpp, ipp, jpp) = 1 - A3D_ELEM(mask, kpp, ipp, jpp);
    }
}

void BinaryWedgeMask(MultidimArray<int> &mask, double theta0, double thetaF,
                     const Matrix2D<double> &A, bool centerOrigin)
{
    int halfX, halfY,halfZ;
    int minX,minY,minZ;
    int maxX,maxY,maxZ;
    int ipp,jpp,kpp;

    halfX = mask.xdim/2;
    halfY = mask.ydim/2;
    halfZ = mask.zdim/2;

    minX = -halfX;
    minY = -halfY;
    minZ = -halfZ;

    maxX = (int)((mask.xdim-0.5)/2);
    maxY = (int)((mask.ydim-0.5)/2);
    maxZ = (int)((mask.zdim-0.5)/2);


    double xp, zp;
    double tg0, tgF, limx0, limxF;

    tg0 = -tan(PI * (-90. - thetaF) / 180.);
    tgF = -tan(PI * (90. - theta0) / 180.);
    if (ABS(tg0) < XMIPP_EQUAL_ACCURACY)
        tg0=0.;
    if (ABS(tgF) < XMIPP_EQUAL_ACCURACY)
        tgF=0.;
    // ROB: A=A.inv(); A no const
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
    {
        double di=(double)i;
        double dj=(double)j;
        double dk=(double)k;
        xp = MAT_ELEM(A, 0, 0) * dj + MAT_ELEM(A, 0, 1) * di + MAT_ELEM(A, 0, 2) * dk;
        zp = MAT_ELEM(A, 2, 0) * dj + MAT_ELEM(A, 2, 1) * di + MAT_ELEM(A, 2, 2) * dk;

        if (centerOrigin)
        {
            kpp = intWRAP (k+halfZ,minZ,maxZ);
            ipp = intWRAP (i+halfY,minY,maxY);
            jpp = intWRAP (j+halfX,minX,maxX);
        }
        else
        {
            kpp=k;
            ipp=i;
            jpp=j;
        }

        limx0 = tg0 * zp;// + 0.5;
        limxF = tgF * zp;// + 0.5;
        if (zp >= 0)
        {
            if (xp <= limx0 || xp >= limxF)
            {
                A3D_ELEM(mask, kpp, ipp, jpp) = 1;
            }
            else
                A3D_ELEM(mask, kpp, ipp, jpp) = 0;
        }
        else
        {
            if (xp <= limxF || xp >= limx0)
            {
                A3D_ELEM(mask, kpp, ipp, jpp) = 1;
            }
            else
                A3D_ELEM(mask, kpp, ipp, jpp) = 0;
        }
    }
}


void mask3D_6neig(MultidimArray<int> &mask, int value, int center)
{
    mask.resize(3, 3, 3);
    mask.initZeros();
    mask(1, 1, 1) = center;
    mask(1, 1, 0) = mask(1, 1, 2) = mask(0, 1, 1) = mask(2, 1, 1) = mask(1, 0, 1) = mask(1, 2, 1) = value;

}

void mask3D_18neig(MultidimArray<int> &mask, int value1, int value2, int center)
{
    mask.resize(3, 3, 3);
    mask.initZeros();
    mask(1, 1, 1) = center;
    //Face neighbors
    mask(1, 1, 0) = mask(1, 1, 2) = mask(0, 1, 1) = mask(2, 1, 1) = mask(1, 0, 1) = mask(1, 2, 1) = value1;
    //Edge neighbors
    mask(0, 1, 0) = mask(0, 0, 1) = mask(0, 1, 2) = mask(0, 2, 1) = value2;
    mask(1, 0, 0) = mask(1, 2, 0) = mask(1, 0, 2) = mask(1, 2, 2) = value2;
    mask(2, 1, 0) = mask(2, 0, 1) = mask(2, 1, 2) = mask(2, 2, 1) = value2;


}
void mask3D_26neig(MultidimArray<int> &mask, int value1, int value2, int value3,
                   int center)
{
    mask.resize(3, 3, 3);
    mask.initZeros();
    mask(1, 1, 1) = center;
    //Face neighbors
    mask(1, 1, 0) = mask(1, 1, 2) = mask(0, 1, 1) = mask(2, 1, 1) = mask(1, 0, 1) = mask(1, 2, 1) = value1;
    //Edge neighbors
    mask(0, 1, 0) = mask(0, 0, 1) = mask(0, 1, 2) = mask(0, 2, 1) = value2;
    mask(1, 0, 0) = mask(1, 2, 0) = mask(1, 0, 2) = mask(1, 2, 2) = value2;
    mask(2, 1, 0) = mask(2, 0, 1) = mask(2, 1, 2) = mask(2, 2, 1) = value2;
    //Vertex neighbors
    mask(0, 0, 0) = mask(0, 0, 2) = mask(0, 2, 0) = mask(0, 2, 2) = value3;
    mask(2, 0, 0) = mask(2, 0, 2) = mask(2, 2, 0) = mask(2, 2, 2) = value3;

}

/*---------------------------------------------------------------------------*/
/* Mask Type                                                                 */
/*---------------------------------------------------------------------------*/
// Constructor -------------------------------------------------------------
Mask::Mask(int _allowed_data_types)
{
    clear();
    allowed_data_types = _allowed_data_types;
}

// Default values ----------------------------------------------------------
void Mask::clear()
{
    type = NO_MASK;
    mode = INNER_MASK;
    H = R1 = R2 = sigma = 0;
    imask.clear();
    dmask.clear();
    allowed_data_types = 0;
    fn_mask = "";
    x0 = y0 = z0 = 0;
}

// Resize ------------------------------------------------------------------
void Mask::resize(size_t Xdim)
{
    switch (datatype())
    {
    case INT_MASK:
        imask.resize(Xdim);
        imask.setXmippOrigin();
        break;
    case DOUBLE_MASK:
        dmask.resize(Xdim);
        dmask.setXmippOrigin();
        break;
    }
}

void Mask::resize(size_t Ydim, size_t Xdim)
{
    switch (datatype())
    {
    case INT_MASK:
        imask.resize(Ydim, Xdim);
        imask.setXmippOrigin();
        break;
    case DOUBLE_MASK:
        dmask.resize(Ydim, Xdim);
        dmask.setXmippOrigin();
        break;
    }
}

void Mask::resize(size_t Zdim, size_t Ydim, size_t Xdim)
{
    switch (datatype())
    {
    case INT_MASK:
        imask.resize(Zdim, Ydim, Xdim);
        imask.setXmippOrigin();
        break;
    case DOUBLE_MASK:
        dmask.resize(Zdim, Ydim, Xdim);
        dmask.setXmippOrigin();
        break;
    }
}

//#ifdef NEVER
// Read from command lines -------------------------------------------------
void Mask::read(int argc, const char **argv)
{
    int i = paremeterPosition(argc, argv, "--center");
    if (i != -1)
    {
        if (i + 3 >= argc)
            REPORT_ERROR(ERR_UNCLASSIFIED, "Mask: Not enough parameters after -center");
        x0 = textToFloat(argv[i+1]);
        y0 = textToFloat(argv[i+2]);
        z0 = textToFloat(argv[i+3]);
    }
    else
    {
        x0 = y0 = z0 = 0;
    }

    i = paremeterPosition(argc, argv, "--mask");
    if (i == -1)
    {
        clear();
        return;
    }
    if (i + 1 >= argc)
        REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: --mask with no mask_type");
    // Circular mask ........................................................
    if (strcmp(argv[i+1], "circular") == 0)
    {
        if (i + 2 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: circular mask with no radius");
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        R1 = textToFloat(argv[i+2]);
        if (R1 < 0)
        {
            mode = INNER_MASK;
            R1 = ABS(R1);
        }
        else if (R1 > 0)
            mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: circular mask with radius 0");
        type = BINARY_CIRCULAR_MASK;
        // Circular DWT mask ....................................................
    }
    else if (strcmp(argv[i+1], "DWT_circular") == 0)
    {
        if (i + 5 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: DWT circular mask with not enough parameters");
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        R1 = ABS(textToFloat(argv[i+2]));
        smin = textToInteger(argv[i+3]);
        smax = textToInteger(argv[i+4]);
        quadrant = argv[i+5];
        type = BINARY_DWT_CIRCULAR_MASK;
        // Rectangular mask .....................................................
    }
    else if (strcmp(argv[i+1], "rectangular") == 0)
    {
        if (i + 3 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: rectangular mask needs at least two dimensions");
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        Xrect = textToInteger(argv[i+2]);
        Yrect = textToInteger(argv[i+3]);
        if (i + 4 < argc)
        {
            Zrect = textToInteger(argv[i+4]);
            if (argv[i+4][0] != '-')
                Zrect = ABS(Zrect);
        }
        else
            Zrect = 0;
        if (Xrect < 0 && Yrect < 0 && Zrect <= 0)
        {
            mode = INNER_MASK;
            Xrect = ABS(Xrect);
            Yrect = ABS(Yrect);
            Zrect = ABS(Zrect);
        }
        else if (Xrect > 0 && Yrect > 0 && Zrect >= 0)
            mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for rectangle");
        type = BINARY_FRAME_MASK;
        // Cone mask ............................................................
    }
    else if (strcmp(argv[i+1], "cone") == 0)
    {
        if (i + 2 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: cone mask needs one angle");
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        R1 = textToFloat(argv[i+2]);
        if (R1 < 0)
        {
            mode = INNER_MASK;
            R1 = ABS(R1);
        }
        else
            mode = OUTSIDE_MASK;
        type = BINARY_CONE_MASK;
        // Wedge mask ............................................................
    }
    else if (strcmp(argv[i+1], "wedge") == 0)
    {
        if (i + 3 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: wedge mask needs two angles");
        if (!(allowed_data_types & DOUBLE_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        R1 = textToFloat(argv[i+2]);
        R2 = textToFloat(argv[i+3]);
        type = BINARY_WEDGE_MASK;
        // Crown mask ...........................................................
    }
    else if (strcmp(argv[i+1], "crown") == 0)
    {
        if (i + 3 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: crown mask needs two radii");
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        R1 = textToFloat(argv[i+2]);
        R2 = textToFloat(argv[i+3]);
        if (R1 < 0 && R2 < 0)
        {
            mode = INNER_MASK;
            R1 = ABS(R1);
            R2 = ABS(R2);
        }
        else if (R1 > 0 && R2 > 0)
            mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for crown");
        type = BINARY_CROWN_MASK;
        // Cylinder mask ........................................................
    }
    else if (strcmp(argv[i+1], "cylinder") == 0)
    {
        if (i + 3 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: cylinder mask needs a radius and a height");
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        R1 = textToFloat(argv[i+2]);
        H = textToFloat(argv[i+3]);
        if (R1 < 0 && H < 0)
        {
            mode = INNER_MASK;
            R1 = ABS(R1);
            H = ABS(H);
        }
        else if (R1 > 0 && H > 0)
            mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for cylinder");
        type = BINARY_CYLINDER_MASK;
        // Gaussian mask ........................................................
    }
    else if (strcmp(argv[i+1], "tube") == 0)
    {
        if (i + 4 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: tube mask needs two radii and a height");
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        R1 = textToFloat(argv[i+2]);
        R2 = textToFloat(argv[i+3]);
        H = textToFloat(argv[i+4]);
        if (R1 < 0 && R2 < 0 && H<0)
        {
            mode = INNER_MASK;
            R1 = ABS(R1);
            R2 = ABS(R2);
            H=ABS(H);
        }
        else if (R1 > 0 && R2 > 0 && H>0)
            mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for crown");
        type = BINARY_TUBE;
        // Cylinder mask ........................................................
    }
    else if (strcmp(argv[i+1], "gaussian") == 0)
    {
        if (i + 2 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: gaussian mask needs a sigma");
        if (!(allowed_data_types & DOUBLE_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: continuous masks are not allowed");
        sigma = textToFloat(argv[i+2]);
        if (sigma < 0)
        {
            mode = INNER_MASK;
            sigma = ABS(sigma);
        }
        else
            mode = OUTSIDE_MASK;
        type = GAUSSIAN_MASK;
        // Raised cosine mask ...................................................
    }
    else if (strcmp(argv[i+1], "raised_cosine") == 0)
    {
        if (i + 3 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: raised_cosine mask needs two radii");
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: continuous masks are not allowed");
        R1 = textToFloat(argv[i+2]);
        R2 = textToFloat(argv[i+3]);
        if (R1 < 0 && R2 < 0)
        {
            mode = INNER_MASK;
            R1 = ABS(R1);
            R2 = ABS(R2);
        }
        else if (R1 > 0 && R2 > 0)
            mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for raised_cosine");
        type = RAISED_COSINE_MASK;
        // Raised crown mask ....................................................
    }
    else if (strcmp(argv[i+1], "raised_crown") == 0)
    {
        if (i + 4 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: raised_crown mask needs two radii & a width");
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: continuous masks are not allowed");
        R1 = textToFloat(argv[i+2]);
        R2 = textToFloat(argv[i+3]);
        pix_width = textToFloat(argv[i+4]);
        if (R1 < 0 && R2 < 0)
        {
            mode = INNER_MASK;
            R1 = ABS(R1);
            R2 = ABS(R2);
        }
        else if (R1 > 0 && R2 > 0)
            mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for raised_cosine");
        type = RAISED_CROWN_MASK;
        // Blob circular mask ....................................................
    }
    else if (strcmp(argv[i+1], "blob_circular") == 0)
    {
        if (i + 3 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: blob_circular mask needs one radius and a width");
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: continuous masks are not allowed");
        R1 = textToFloat(argv[i+2]);
        double aux = textToFloat(argv[i+3]);
        blob_radius= ABS(aux);
        if (aux < 0)
            mode = INNER_MASK;
        else if (aux > 0)
            mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for blob_circular");
        type = BLOB_CIRCULAR_MASK;
        blob_order= textToInteger(getParameter(argc, argv, "-m", "2"));
        blob_alpha= textToFloat(getParameter(argc, argv, "-a", "10.4"));

        // Raised crown mask ....................................................
    }
    else if (strcmp(argv[i+1], "blob_crown") == 0)
    {
        if (i + 4 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: blob_crown mask needs two radii and a with");
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: continuous masks are not allowed");
        R1 = textToFloat(argv[i+2]);
        R2 = textToFloat(argv[i+3]);
        double aux = textToFloat(argv[i+4]);
        blob_radius= ABS(aux);
        if (aux < 0)
            mode = INNER_MASK;
        else if (aux > 0)
            mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for blob_crown");
        type = BLOB_CROWN_MASK;
        blob_order= textToInteger(getParameter(argc, argv, "-m", "2"));
        blob_alpha= textToFloat(getParameter(argc, argv, "-a", "10.4"));

        // Blackman mask ........................................................
    }
    else if (strcmp(argv[i+1], "blackman") == 0)
    {
        mode = INNER_MASK;
        type = BLACKMAN_MASK;
        // Sinc mask ............................................................
    }
    else if (strcmp(argv[i+1], "sinc") == 0)
    {
        if (i + 2 >= argc)
            REPORT_ERROR(ERR_ARG_MISSING, "MaskProgram: sinc mask needs a frequency");
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        omega = textToFloat(argv[i+2]);
        if (omega < 0)
        {
            mode = INNER_MASK;
            omega = ABS(omega);
        }
        else
            mode = OUTSIDE_MASK;
        type = SINC_MASK;
    }
    else if (strcmp(argv[i+1], "binary_file") == 0)
    {
        fn_mask = argv[i+2];
        type = READ_BINARY_MASK;
    }
    else if (strcmp(argv[i+1], "real_file") == 0)
    {
        fn_mask = argv[i+2];
        type = READ_REAL_MASK;
    }
    else
        REPORT_ERROR(ERR_ARG_INCORRECT,"Incorrect mask type");
}
//#endif
// Show --------------------------------------------------------------------
void Mask::show() const
{
#define SHOW_MODE \
    if (mode==INNER_MASK) std::cout << "   mode=INNER MASK\n"; \
    else                  std::cout << "   mode=OUTER MASK\n";
#define SHOW_CENTER \
    std::cout << "   (x0,y0,z0)=(" << x0 << "," << y0 << "," << z0 << ")\n";
    switch (type)
    {
    case NO_MASK:
        std::cout << "Mask type: No mask\n";
        break;
    case BINARY_CIRCULAR_MASK:
        std::cout << "Mask type: Binary circular\n"
        << "   R=" << R1 << std::endl;
        SHOW_MODE;
        SHOW_CENTER;
        break;
    case BINARY_DWT_CIRCULAR_MASK:
        std::cout << "Mask type: Binary DWT circular\n"
        << "   R=" << R1 << std::endl
        << "   smin=" << smin << std::endl
        << "   smax=" << smax << std::endl
        << "   quadrant=" << quadrant << std::endl;
        break;
    case BINARY_CROWN_MASK:
        std::cout << "Mask type: Binary crown\n"
        << "   R1=" << R1 << std::endl
        << "   R2=" << R2 << std::endl;
        SHOW_MODE;
        SHOW_CENTER;
        break;
    case BINARY_CYLINDER_MASK:
        std::cout << "Mask type: Cylinder\n"
        << "   R1=" << R1 << std::endl
        << "   H=" << H << std::endl;
        SHOW_MODE;
        SHOW_CENTER;
        break;
    case BINARY_TUBE:
    	std::cout << "Mask type: Tube\n"
        << "   R1=" << R1 << std::endl
        << "   R2=" << R2 << std::endl
        << "   H=" << H << std::endl;
        SHOW_MODE;
        SHOW_CENTER;
        break;
    case BINARY_FRAME_MASK:
        std::cout << "Mask type: Frame\n"
        << "   Xrect=" << Xrect << std::endl
        << "   Yrect=" << Yrect << std::endl;
        SHOW_MODE;
        SHOW_CENTER;
        break;
    case GAUSSIAN_MASK:
        std::cout << "Mask type: Gaussian\n"
        << "   sigma=" << sigma << std::endl;
        SHOW_MODE;
        SHOW_CENTER;
        break;
    case RAISED_COSINE_MASK:
        std::cout << "Mask type: Raised cosine\n"
        << "   R1=" << R1 << std::endl
        << "   R2=" << R2 << std::endl;
        SHOW_MODE;
        SHOW_CENTER;
        break;
    case RAISED_CROWN_MASK:
        std::cout << "Mask type: Raised crown\n"
        << "   R1=" << R1 << std::endl
        << "   R2=" << R2 << std::endl
        << "   pixwidth=" << pix_width << std::endl;
        SHOW_MODE;
        SHOW_CENTER;
        break;
    case BLOB_CIRCULAR_MASK:
        std::cout << "Mask type: Blob circular\n"
        << "   R1=" << R1 << std::endl
        << "   blob radius=" << blob_radius << std::endl
        << "   blob order="  << blob_order  << std::endl
        << "   blob alpha="  << blob_alpha  << std::endl;
        SHOW_MODE;
        SHOW_CENTER;
        break;
    case BLOB_CROWN_MASK:
        std::cout << "Mask type: Blob crown\n"
        << "   R1=" << R1 << std::endl
        << "   R2=" << R2 << std::endl
        << "   blob radius=" << blob_radius << std::endl
        << "   blob order="  << blob_order  << std::endl
        << "   blob alpha="  << blob_alpha  << std::endl;
        SHOW_MODE;
        SHOW_CENTER;
        break;
    case BLACKMAN_MASK:
        std::cout << "Mask type: Blackman\n";
        SHOW_MODE;
        SHOW_CENTER;
        break;
    case SINC_MASK:
        std::cout << "Mask type: Sinc\n"
        << "   w=" << omega << std::endl;
        SHOW_MODE;
        SHOW_CENTER;
        break;
    default:
        std::cout << "Mask type: Read from disk\n"
        << "   File=" << fn_mask << std::endl;
        break;
    }
}

// Usage -------------------------------------------------------------------
void Mask::usage() const
{
    std::cerr << "Mask usage:\n";
    std::cerr << "   [-center <x0=0> <y0=0> <z0=0>]: Center of the mask\n";
    if (allowed_data_types & INT_MASK)
        std::cerr << "   [-mask circular <R>       : circle/sphere mask\n"
        << "                               if R>0 => outside R\n"
        << "                               if R<0 => inside  R\n"
        << "   [-mask DWT_circular <R> <smin> <smax>: circle/sphere mask\n"
        << "                               smin and smax define the scales\n"
        << "                               to be kept\n"
        << "   |-mask rectangular <Xrect> <Yrect> [<Zrect>]: 2D or 3D rectangle\n"
        << "                               if X,Y,Z > 0 => outside rectangle\n"
        << "                               if X,Y,Z < 0 => inside rectangle\n"
        << "   |-mask crown <R1> <R2>    : 2D or 3D crown\n"
        << "                               if R1,R2 > 0 => outside crown\n"
        << "                               if R1,R2 < 0 => inside crown\n"
        << "   |-mask cylinder <R> <H>   : 2D circle or 3D cylinder\n"
        << "                               if R,H > 0 => outside cylinder\n"
        << "                               if R,H < 0 => inside cylinder\n"
        << "   |-mask tube <R1> <R2> <H> : 2D or 3D tube\n"
        << "                               if R1,R2,H > 0 => outside tube\n"
        << "                               if R1,R2,H < 0 => inside tube\n"
        << "   |-mask cone <theta>       : 3D cone (parallel to Z) \n"
        << "                               if theta > 0 => outside cone\n"
        << "                               if theta < 0 => inside cone\n"
        << "   |-mask wedge <th0> <thF>  : 3D missing-wedge mask for data \n"
        << "                               collected between tilting angles \n"
        << "                               th0 and thF (around the Y-axis) \n"
        << "   |-mask <binary file>      : Read from file\n"
        ;
    if (allowed_data_types & DOUBLE_MASK)
        std::cerr << "   |-mask gaussian <sigma>   : 2D or 3D gaussian\n"
        << "                               if sigma > 0 => outside gaussian\n"
        << "                               if sigma < 0 => inside gaussian\n"
        << "   |-mask raised_cosine <R1> <R2>: 2D or 3D raised_cosine\n"
        << "                               if R1,R2 > 0 => outside sphere\n"
        << "                               if R1,R2 < 0 => inside sphere\n"
        << "   |-mask raised_crown <R1> <R2> <pixwidth>: 2D or 3D raised_crown\n"
        << "                               if R1,R2 > 0 => outside crown\n"
        << "                               if R1,R2 < 0 => inside crown\n"
        << "   |-mask blob_circular <R1> <blob_radius>: 2D or 3D blob circular\n"
        << "                               if blob_radius > 0 => outside blob\n"
        << "                               if blob_radius < 0 => inside blob\n"
        << "   |-mask blob_crown <R1> <R2> <blob_radius>: 2D or 3D blob_crown\n"
        << "                               if blob_radius > 0 => outside crown\n"
        << "                               if blob_radius < 0 => inside crown\n"
        << "   [ -m <blob_order=2>       : Order of blob\n"
        << "   [ -a <blob_alpha=10.4>    : Alpha of blob\n"
        << "   |-mask blackman           : 2D or 3D Blackman mask\n"
        << "                               always inside blackman\n"
        << "   |-mask sinc <w>]          : 2D or 3D sincs\n"
        << "                               if w > 0 => outside sinc\n"
        << "                               if w < 0 => inside sinc\n"
        ;
}

// Write -------------------------------------------------------------------
void Mask::write_mask(const FileName &fn)
{
    Image<double> img;
    if (datatype() == INT_MASK)
        typeCast(imask, img());
    else if (datatype() == DOUBLE_MASK)
        img()=dmask;
    img.write(fn);
}


void Mask::defineParams(XmippProgram * program, int allowed_data_types, const char* prefix, const char* comment, bool moreOptions)
{
    char tempLine[256], tempLine2[256];

    char advanced=' ';
    if (moreOptions)
    	advanced='+';
    if(prefix == NULL)
        sprintf(tempLine, "  [--mask%c <mask_type=circular>] ",advanced);
    else
        sprintf(tempLine,"%s --mask%c <mask_type=circular> ", prefix,advanced);
    if (comment != NULL)
        sprintf(tempLine2, "%s : %s", tempLine, comment);
    else
    	strcpy(tempLine2,tempLine);

    program->addParamsLine(tempLine2);
    program->addParamsLine("        where <mask_type> ");
    // program->addParamsLine("== INT MASK ==");
    if (allowed_data_types & INT_MASK)
    {
        program->addParamsLine("         circular <R>        : circle/sphere mask");
        program->addParamsLine("                :if R>0 => outside R");
        program->addParamsLine("                :if R<0 => inside  R");
        program->addParamsLine("         DWT_circular <R> <smin> <smax>: circle/sphere mask");
        program->addParamsLine("                : smin and smax define the scales to be kept");
        program->addParamsLine("         rectangular <Xrect> <Yrect> <Zrect=-1>: 2D or 3D rectangle");
        program->addParamsLine("                                         :if X,Y,Z > 0 => outside rectangle");
        program->addParamsLine("                                          :if X,Y,Z < 0 => inside rectangle");
        program->addParamsLine("         crown <R1> <R2>    : 2D or 3D crown");
        program->addParamsLine("                                          :if R1,R2 > 0 => outside crown");
        program->addParamsLine("                                          :if R1,R2 < 0 => inside crown");
        program->addParamsLine("         cylinder <R> <H>   : 2D circle or 3D cylinder");
        program->addParamsLine("                                         :if R,H > 0 => outside cylinder");
        program->addParamsLine("                                          :if R,H < 0 => inside cylinder");
        program->addParamsLine("         tube <R1> <R2> <H> : 3D tube");
        program->addParamsLine("                                          :if R1,R2,H > 0 => outside tube");
        program->addParamsLine("                                          :if R1,R2,H < 0 => inside tube");
        program->addParamsLine("         cone <theta>       : 3D cone (parallel to Z) ");
        program->addParamsLine("                                          :if theta > 0 => outside cone");
        program->addParamsLine("                                         :if theta < 0 => inside cone");
        program->addParamsLine("         wedge <th0> <thF>  : 3D missing-wedge mask for data ");
        program->addParamsLine("                                          :collected between tilting angles ");
        program->addParamsLine("                                          :th0 and thF (around the Y-axis) ");
        program->addParamsLine("         binary_file <binary_file>      : Read from file and cast to binary");
    }
    //program->addParamsLine("== DOUBLE MASK ==");
    if (allowed_data_types & DOUBLE_MASK)
    {
        program->addParamsLine("         real_file <float_file>       : Read from file and do not cast");
        program->addParamsLine("         gaussian <sigma>   : 2D or 3D gaussian");
        program->addParamsLine("                              :if sigma > 0 => outside gaussian");
        program->addParamsLine("                              : if sigma < 0 => inside gaussian");
        program->addParamsLine("         raised_cosine <R1> <R2>: 2D or 3D raised_cosine");
        program->addParamsLine("                              : if R1,R2 > 0 => outside sphere");
        program->addParamsLine("                              : if R1,R2 < 0 => inside sphere");
        program->addParamsLine("         raised_crown <R1> <R2> <pixwidth>: 2D or 3D raised_crown");
        program->addParamsLine("                              : if R1,R2 > 0 => outside sphere");
        program->addParamsLine("                              : if R1,R2 < 0 => inside sphere");
        program->addParamsLine("         blob_circular <R1> <blob_radius>: 2D or 3D blob circular");
        program->addParamsLine("                              : if blob_radius > 0 => outside sphere");
        program->addParamsLine("                              : if blob_radius < 0 => inside sphere");
        program->addParamsLine("         blob_crown <R1> <R2> <blob_radius>: 2D or 3D blob_crown");
        program->addParamsLine("                              : if blob_radius > 0 => outside sphere");
        program->addParamsLine("                              : if blob_radius < 0 => inside sphere");
        program->addParamsLine("         blackman           : 2D or 3D Blackman mask");
        program->addParamsLine("                             :  always inside blackman");
        program->addParamsLine("         sinc <w>          : 2D or 3D sincs");
        program->addParamsLine("                             :  if w > 0 => outside sinc");
        program->addParamsLine("                             :  if w < 0 => inside sinc");
        sprintf(tempLine, "   [ -m%c <blob_order=2>]       : Order of blob",advanced);
        program->addParamsLine(tempLine);
        sprintf(tempLine, "   [ -a%c <blob_alpha=10.4>]    : Alpha of blob",advanced);
        program->addParamsLine(tempLine);
    }
    sprintf(tempLine, "   [--center%c <x0=0> <y0=0> <z0=0>]: mask center",advanced);
    program->addParamsLine(tempLine);
}

void Mask::readParams(XmippProgram * program)
{
    x0 = y0 = z0 = 0;
    x0 = program->getDoubleParam("--center",0);
    y0 = program->getDoubleParam("--center",1);
    z0 = program->getDoubleParam("--center",2);
    mask_type = program->getParam("--mask");

    /* Circular mask ........................................................*/
    if (mask_type == "circular")
    {
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        R1 = program->getDoubleParam("--mask","circular");
        if (R1 < 0)
        {
            Mask::mode = INNER_MASK;
            R1 = ABS(R1);
        }
        else if (R1 > 0)
            Mask::mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: circular mask with radius 0");
        type = BINARY_CIRCULAR_MASK;
    }
    /*// Circular DWT mask ....................................................*/
    else if (mask_type == "DWT_circular")
    {
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        R1 = program->getDoubleParam("--mask","DWT_circular",0);
        smin = program->getIntParam("--mask","DWT_circular",1);
        smax = program->getIntParam("--mask","DWT_circular",2);
        quadrant = program->getParam("--mask","DWT_circular",3);
        type = BINARY_DWT_CIRCULAR_MASK;
    }
    /*// Rectangular mask .....................................................*/
    else if (mask_type == "rectangular")
    {
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        Xrect = program->getIntParam("--mask","rectangular",0);
        Yrect = program->getIntParam("--mask","rectangular",1);
        Zrect = program->getIntParam("--mask","rectangular",2);
        if (Xrect < 0 && Yrect < 0 && Zrect <= 1)
        {
            Mask::mode = INNER_MASK;
            Xrect = ABS(Xrect);
            Yrect = ABS(Yrect);
            Zrect = ABS(Zrect);
        }
        else if (Xrect > 0 && Yrect > 0 && Zrect >= -1)
            Mask::mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for rectangle");
        type = BINARY_FRAME_MASK;
    }
    /*// Cone mask ............................................................*/
    else if (mask_type == "cone")
    {
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        R1 = program->getDoubleParam("--mask","cone");
        if (R1 < 0)
        {
            Mask::mode = INNER_MASK;
            R1 = ABS(R1);
        }
        else
            Mask::mode = OUTSIDE_MASK;
        type = BINARY_CONE_MASK;
    }
    /*// Wedge mask ............................................................*/
    else if (mask_type == "wedge")
    {
        if (!(allowed_data_types & DOUBLE_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");
        R1 = program->getDoubleParam("--mask","wedge",0);
        R2 = program->getDoubleParam("--mask","wedge",1);
        type = BINARY_WEDGE_MASK;
    }
    /*// Crown mask ...........................................................*/
    else if (mask_type == "crown")
    {
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");

        R1 = program->getDoubleParam("--mask","crown",0);
        R2 = program->getDoubleParam("--mask","crown",1);

        if (R1 < 0 && R2 < 0)
        {
            Mask::mode = INNER_MASK;
            R1 = ABS(R1);
            R2 = ABS(R2);
        }
        else if (R1 > 0 && R2 > 0)
            Mask::mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for crown");

        type = BINARY_CROWN_MASK;
    }
    /*// Cylinder mask ........................................................*/
    else if (mask_type == "cylinder")
    {
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");

        R1 = program->getDoubleParam("--mask","cylinder",0);
        H = program->getDoubleParam("--mask","cylinder",1);

        if (R1 < 0 && H < 0)
        {
            Mask::mode = INNER_MASK;
            R1 = ABS(R1);
            H = ABS(H);
        }
        else if (R1 > 0 && H > 0)
            Mask::mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for cylinder");

        type = BINARY_CYLINDER_MASK;
    }
    /*// Crown mask ...........................................................*/
    else if (mask_type == "tube")
    {
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");

        R1 = program->getDoubleParam("--mask","tube",0);
        R2 = program->getDoubleParam("--mask","tube",1);
        H = program->getDoubleParam("--mask","tube",2);

        if (R1 < 0 && R2 < 0 && H<0)
        {
            Mask::mode = INNER_MASK;
            R1 = ABS(R1);
            R2 = ABS(R2);
            H=ABS(H);
        }
        else if (R1 > 0 && R2 > 0 && H>0)
            Mask::mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for tube");

        type = BINARY_TUBE;
    }
    /*// Gaussian mask ........................................................*/
    else if (mask_type == "gaussian")
    {
        if (!(allowed_data_types & DOUBLE_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: continuous masks are not allowed");

        sigma = program->getDoubleParam("--mask","gaussian");

        if (sigma < 0)
        {
            Mask::mode = INNER_MASK;
            sigma = ABS(sigma);
        }
        else
            Mask::mode = OUTSIDE_MASK;

        type = GAUSSIAN_MASK;
    }
    /*// Raised cosine mask ...................................................*/
    else if (mask_type == "raised_cosine")
    {
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: continuous masks are not allowed");

        R1 = program->getDoubleParam("--mask","raised_cosine",0);
        R2 = program->getDoubleParam("--mask","raised_cosine",1);

        if (R1 < 0 && R2 < 0)
        {
            Mask::mode = INNER_MASK;
            R1 = ABS(R1);
            R2 = ABS(R2);
        }
        else if (R1 > 0 && R2 > 0)
            Mask::mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for raised_cosine");

        type = RAISED_COSINE_MASK;
    }
    /*// Raised crown mask ....................................................*/
    else if (mask_type == "raised_crown")
    {
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: continuous masks are not allowed");

        R1 = program->getDoubleParam("--mask","raised_crown",0);
        R2 = program->getDoubleParam("--mask","raised_crown",1);
        pix_width = program->getDoubleParam("--mask","raised_crown",2);

        if (R1 < 0 && R2 < 0)
        {
            Mask::mode = INNER_MASK;
            R1 = ABS(R1);
            R2 = ABS(R2);
        }
        else if (R1 > 0 && R2 > 0)
            Mask::mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for raised_cosine");

        type = RAISED_CROWN_MASK;
    }
    /*// Blob circular mask ....................................................*/
    else if (mask_type == "blob_circular")
    {
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: continuous masks are not allowed");

        R1 = program->getDoubleParam("--mask","blob_circular",0);
        double aux = program->getDoubleParam("--mask","blob_circular",1);
        blob_radius= ABS(aux);

        if (aux < 0)
            Mask::mode = INNER_MASK;
        else if (aux > 0)
            Mask::mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for blob_circular");

        type = BLOB_CIRCULAR_MASK;
        blob_order= program->getIntParam("-m");
        blob_alpha= program->getIntParam("-a");
    }
    /*// Raised crown mask ....................................................*/
    else if (mask_type == "blob_crown")
    {
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: continuous masks are not allowed");

        R1 = program->getDoubleParam("--mask","blob_crown",0);
        R2 = program->getDoubleParam("--mask","blob_crown",1);
        double aux = program->getDoubleParam("--mask","blob_crown",2);
        blob_radius= ABS(aux);

        if (aux < 0)
            Mask::mode = INNER_MASK;
        else if (aux > 0)
            Mask::mode = OUTSIDE_MASK;
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: cannot determine mode for blob_crown");

        type = BLOB_CROWN_MASK;
        blob_order= program->getIntParam("-m");
        blob_alpha= program->getDoubleParam("-a");
    }
    /*// Blackman mask ........................................................*/
    else if (mask_type == "blackman")
    {
        Mask::mode = INNER_MASK;
        type = BLACKMAN_MASK;
    }
    /*// Sinc mask ............................................................*/
    else if (mask_type == "sinc")
    {
        if (!(allowed_data_types & INT_MASK))
            REPORT_ERROR(ERR_ARG_INCORRECT, "MaskProgram: binary masks are not allowed");

        omega = program->getDoubleParam("--mask","sinc");

        if (omega < 0)
        {
            Mask::mode = INNER_MASK;
            omega = ABS(omega);
        }
        else
            Mask::mode = OUTSIDE_MASK;

        type = SINC_MASK;
    }
    else if (mask_type == "binary_file")
    {
        fn_mask = program->getParam("--mask","binary_file");
        type = READ_BINARY_MASK;
    }
    else if (mask_type == "real_file")
    {
        fn_mask = program->getParam("--mask","real_file");
        type = READ_REAL_MASK;
    }
    else
    	REPORT_ERROR(ERR_DEBUG_IMPOSIBLE,"You should never see  this meessage, Mask::readParams ");
}
// Generate mask --------------------------------------------------------
void Mask::generate_mask(bool apply_geo)
{
    ImageGeneric img;
    Matrix2D<double> AA(4, 4);
    AA.initIdentity();
    blobtype blob;
    if (type==BLOB_CIRCULAR_MASK || type==BLOB_CROWN_MASK)
    {
        blob.radius = blob_radius;
        blob.order = blob_order;
        blob.alpha = blob_alpha;
    }

    switch (type)
    {
    case NO_MASK:
        imask.initConstant(1);
        break;
    case BINARY_CIRCULAR_MASK:
        BinaryCircularMask(imask, R1, mode, x0, y0, z0);
        break;
    case BINARY_DWT_CIRCULAR_MASK:
        BinaryDWTCircularMask2D(imask, R1, smin, smax, quadrant);
        break;
    case BINARY_DWT_SPHERICAL_MASK:
        BinaryDWTSphericalMask3D(imask, R1, smin, smax, quadrant);
        break;
    case BINARY_CROWN_MASK:
        BinaryCrownMask(imask, R1, R2, mode, x0, y0, z0);
        break;
    case BINARY_CYLINDER_MASK:
        BinaryCylinderMask(imask, R1, H, mode, x0, y0, z0);
        break;
    case BINARY_TUBE:
        BinaryTubeMask(imask, R1, R2, H,  mode, x0, y0, z0);
        break;
    case BINARY_FRAME_MASK:
        BinaryFrameMask(imask, Xrect, Yrect, Zrect, mode, x0, y0, z0);
        break;
    case BINARY_CONE_MASK:
        BinaryConeMask(imask, R1, mode);
        break;
    case BINARY_WEDGE_MASK:
        BinaryWedgeMask(imask, R1, R2, AA);
        break;
    case GAUSSIAN_MASK:
        GaussianMask(dmask, sigma, mode, x0, y0, z0);
        break;
    case RAISED_COSINE_MASK:
        RaisedCosineMask(dmask, R1, R2, mode, x0, y0, z0);
        break;
    case RAISED_CROWN_MASK:
        RaisedCrownMask(dmask, R1, R2, pix_width, mode, x0, y0, z0);
        break;
    case BLOB_CIRCULAR_MASK:
        BlobCircularMask(dmask, R1, blob, mode, x0, y0, z0);
        break;
    case BLOB_CROWN_MASK:
        BlobCrownMask(dmask, R1, R2, blob, mode, x0, y0, z0);
        break;
    case BLACKMAN_MASK:
        BlackmanMask(dmask, mode, x0, y0, z0);
        break;
    case SINC_MASK:
        SincMask(dmask, omega, mode, x0, y0, z0);
        break;
    case READ_BINARY_MASK:
        img.readMapped(fn_mask);
        img().getImage(imask);
        imask.setXmippOrigin();
        break;
    case READ_REAL_MASK://ROB
        img.readMapped(fn_mask);
        img().getImage(dmask);
        dmask.setXmippOrigin();
        break;
    default:
        REPORT_ERROR(ERR_VALUE_INCORRECT, "MaskProgram::generate_mask: Unknown mask type :"
                     + integerToString(type));
    }

    if (apply_geo)
    {
        switch (datatype())
        {
        case INT_MASK:
            if (ZSIZE(imask) > 1)
                REPORT_ERROR(ERR_NOT_IMPLEMENTED,"Error: apply_geo only implemented for 2D masks");
            apply_geo_binary_2D_mask(imask, mask_geo);
            break;
        case DOUBLE_MASK:
            if (ZSIZE(dmask) > 1)
                REPORT_ERROR(ERR_NOT_IMPLEMENTED,"Error: apply_geo only implemented for 2D masks");
            apply_geo_cont_2D_mask(dmask, mask_geo);
            break;
        }
    }
}


/*---------------------------------------------------------------------------*/
/* Mask tools                                                                */
/*---------------------------------------------------------------------------*/

// Apply geometric transformation to a binary mask ========================
void apply_geo_binary_2D_mask(MultidimArray<int> &mask,
                              const Matrix2D<double> &A)
{
    MultidimArray<double> tmp;
    tmp.resize(mask);
    typeCast(mask, tmp);
    double outside = DIRECT_A2D_ELEM(tmp, 0, 0);
    MultidimArray<double> tmp2;
    tmp2 = tmp;
    // Instead of IS_INV for images use IS_NOT_INV for masks!
    applyGeometry(1, tmp, tmp2, A, IS_NOT_INV, DONT_WRAP, outside);
    // The type cast gives strange results here, using round instead
    //typeCast(tmp, mask);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(mask)
    {
        dAij(mask,i,j)=ROUND(dAij(tmp,i,j));
        //        std::cout << "i, j = " << i << "," << j <<std::endl;
    }
}

// Apply geometric transformation to a continuous mask =====================
void apply_geo_cont_2D_mask(MultidimArray<double> &mask,
                            const Matrix2D<double> &A)
{
    double outside = DIRECT_A2D_ELEM(mask, 0, 0);
    MultidimArray<double> tmp = mask;
    // Instead of IS_INV for images use IS_NOT_INV for masks!
    applyGeometry(1, tmp, mask, A, IS_NOT_INV, DONT_WRAP, outside);
}

int count_with_mask(const MultidimArray<int> &mask,
                    const MultidimArray< std::complex<double> > &m, int mode, double th1, double th2)
{
    SPEED_UP_tempsInt;
    int N = 0;
    FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(mask, m)
    if (A2D_ELEM(mask, i, j))
        switch (mode)
        {
        case (COUNT_ABOVE):
                        if (abs(A3D_ELEM(m, k, i, j)) >= th1)
                            N++;
            break;
        case (COUNT_BELOW):
                        if (abs(A3D_ELEM(m, k, i, j)) <= th1)
                            N++;
            break;
        case (COUNT_BETWEEN):
                        if (abs(A3D_ELEM(m, k, i, j)) >= th1 && abs(A3D_ELEM(m, k, i, j)) <= th2)
                            N++;
            break;
        }
    return N;
}

void rangeAdjust_within_mask(const MultidimArray<double> *mask,
                             const MultidimArray<double> &m1, MultidimArray<double> &m2)
{
    Matrix2D<double> A(2, 2);
    A.initZeros();
    Matrix1D<double> b(2);
    b.initZeros();
    SPEED_UP_tempsInt;
    // Compute Least squares solution
    if (mask == NULL)
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(m1, m2)
        {
            A(0, 0) += m2(k, i, j) * m2(k, i, j);
            A(0, 1) += m2(k, i, j);
            A(1, 1) += 1;
            b(0)   += m1(k, i, j) * m2(k, i, j);
            b(1)   += m1(k, i, j);
        }
        A(1, 0) = A(0, 1);
    }
    else
    {
        FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(*mask, m2)
        {
            if ((*mask)(k, i, j))
            {
                A(0, 0) += m2(k, i, j) * m2(k, i, j);
                A(0, 1) += m2(k, i, j);
                A(1, 1) += 1;
                b(0)   += m1(k, i, j) * m2(k, i, j);
                b(1)   += m1(k, i, j);
            }
        }
        A(1, 0) = A(0, 1);
    }
    b = A.inv() * b;

    // Apply to m2
    FOR_ALL_ELEMENTS_IN_ARRAY3D(m2) m2(k, i, j) = b(0) * m2(k, i, j) + b(1);
}

/************** Program Mask implementation ***********************/
/* Define Parameters ----------------------------------------------------------------- */
void ProgMask::defineParams()
{
    each_image_produces_an_output = true;
    save_metadata_stack = true;
    keep_input_columns = true;
    XmippMetadataProgram::defineParams();
    Mask::defineParams(this);

    addUsageLine("Create or Apply a mask. Count pixels/voxels within a mask");
    addUsageLine("+ ");
    addUsageLine("+You do not need to give the dimensions of the mask but you simply provide ");
    addUsageLine("+an example of image/volume you are going to apply the mask to, then the dimensions ");
    addUsageLine("+are taken from this file and the mask is created. In the creation of the mask, ");
    addUsageLine("+a file with the mask is written to disk but it is not applied to the input file.");
    addUsageLine("+ ");
    addUsageLine("+You can generate blank images/volumes with the size of the sample one if you do not ");
    addUsageLine("+supply any mask type.You may also apply masks without having to generate the corresponding");
    addUsageLine("+files (but you also can save them)");
    addUsageLine("+ ");
    addUsageLine("+This utility also allows you to count the number of pixels/voxels in an image/volume");
    addUsageLine("+which are inside a given mask and whose value is below|above or both some threshold.");

    addExampleLine("Sample at circular mask inside radius 72:", false);
    addExampleLine("xmipp_transform_mask  -i reference.vol -o output_volume.vol --mask circular -72");
    addExampleLine("As above but save mask:", false);
    addExampleLine("xmipp_transform_mask  -i reference.vol --create_mask  output_mask.vol --mask circular -25");
    addExampleLine("Mask and overwrite a selection file:", false);
    addExampleLine("xmipp_transform_mask  -i t7_10.sel --mask circular -72");
    addExampleLine("Mask using rectangular mask:", false);
    addExampleLine("xmipp_transform_mask -i singleImage.spi -o salida20.spi --mask rectangular -10 -10");

    addParamsLine("   [--create_mask <output_mask_file>]  : Don't apply and save mask");
    addParamsLine("   [--count_above <th>]                : Voxels within mask >= th");
    addParamsLine("   [--count_below <th>]                : Voxels within mask <= th");
    addParamsLine("   [--substitute <val=\"0\">]  : Value outside the mask: userProvidedValue|min|max|avg");

}

/* Read Parameters ----------------------------------------------------------------- */
void ProgMask::readParams()
{
    XmippMetadataProgram::readParams();
    mask.readParams(this);

    count_above  = checkParam("--count_above");
    if (count_above)
        th_above  = getDoubleParam("-count_above");
    count_below  = checkParam("--count_below");
    if (count_below)
        th_below  = getDoubleParam("--count_below");
    create_mask  = checkParam("--create_mask");
    if (create_mask)
        fn_mask  = getParam("--create_mask");
    //mask.read(argc, argv);
    str_subs_val = getParam("--substitute");
    count = count_below || count_above;
}


/* Preprocess ------------------------------------------------------------- */
void ProgMask::preProcess()
{
//    int max_length = 0;

    if (create_mask && input_is_stack)
        REPORT_ERROR(ERR_MD_NOOBJ, "Mask: Cannot create a mask for a selection file\n");

    // Initialize progress bar
    time_config();
    if (!count && input_is_stack)
        init_progress_bar(mdInSize);
//    else
//        max_length = maxFileNameLength(mdIn);
}

/* Postprocess ------------------------------------------------------------- */
void ProgMask::postProcess()
{
    if (!count)
        progress_bar(mdInSize);
}

/* Process image ------------------------------------------------------------- */
void ProgMask::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    static size_t imageCount = 0;
    ++imageCount;
    Image<double> image;
    image.readApplyGeo(fnImg, rowIn);
    image().setXmippOrigin();

    // Generate mask
    if (ZSIZE(image()) > 1)
        apply_geo=false;
    if (apply_geo)
    {
        if (mask.x0 + mask.y0 != 0.)
            REPORT_ERROR(ERR_ARG_INCORRECT, "Mask: -center option cannot be combined with apply_geo; use -dont_apply_geo");
        else
            // Read geometric transformation from the image and store for mask
            image.getTransformationMatrix(mask.mask_geo);
    }
    mask.generate_mask(image());

    // Apply mask
    if (!create_mask)
    {
        if      (str_subs_val=="min")
            subs_val=image().computeMin();
        else if (str_subs_val=="max")
            subs_val=image().computeMax();
        else if (str_subs_val=="avg")
            subs_val=image().computeAvg();
        else
            subs_val=textToFloat(str_subs_val);

        mask.apply_mask(image(), image(), subs_val, apply_geo);
        if (!count)
            image.write(fnImgOut);
    }
    else
        mask.write_mask(fn_mask);

    // Count
    if (count)
    {
        if (mask.datatype() == INT_MASK)
        {
            int count;
            std::string elem_type="pixels";
            if (ZSIZE(image())>1)
                elem_type="voxels";
            if      (count_above && !count_below)
            {
                std::cout << stringToString(fn_in,max_length)
                << " number of " << elem_type << " above " << th_above;
                count=count_with_mask_above(mask.get_binary_mask(),
                                            image(),th_above);
            }
            else if (count_below && !count_above)
            {
                std::cout << stringToString(fn_in,max_length)
                << " number of " << elem_type << " below " << th_below;
                count=count_with_mask_below(mask.get_binary_mask(),
                                            image(),th_below);
            }
            else if (count_below && count_above)
            {
                std::cout << stringToString(fn_in,max_length)
                << " number of " << elem_type << " above " << th_above
                << " and below " << th_below << " = ";
                count=count_with_mask_between(mask.get_binary_mask(),
                                              image(),th_above,th_below);
            }
            std::cout << " = " << count << std::endl;
        }
        else
            std::cerr << "Cannot count pixels with a continuous mask\n";
    }

    if (imageCount % 25 == 0 && !count)
        progress_bar(imageCount);
}

