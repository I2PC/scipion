/***************************************************************************
 *
 * Authors:     Slavica JONIC (slavica.jonic@impmc.jussieu.fr, slavica.jonic@a3.epfl.ch)
 *              Jean-Noel PIOCHE (jnp95@hotmail.com)
 *
 * Biomedical Imaging Group, EPFL (Lausanne, Suisse).
 * Structures des Assemblages Macromoléculaires, IMPMC UMR 7590 (Paris, France).
 * IUT de Reims-Chlons-Charleville (Reims, France).
 *
 * Last modifications by JNP the 27/05/2009 15:52:56
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

#include "project_real_shears.h"
#include <external/bilib/headers/getputd.h>
#include <external/bilib/headers/changebasis.h>
#include <external/bilib/configs.h>

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503130808472633361816000000000e-16
#endif

RealShearsInfo::RealShearsInfo(const MultidimArray<double> &V)
{
    volume=&V;
    if (XSIZE(V)!=XSIZE(V) || XSIZE(V)!=ZSIZE(V))
        REPORT_ERROR(ERR_MULTIDIM_DIM, "The volume must be cubic");
    Xdim=XSIZE(*volume);

    int    Status = !ERROR;
    MultidimArray<double> planeCoef, inputPlane, inputRow;
    Coef_x.resizeNoCopy(Xdim,Xdim,Xdim);
    Coef_y.resizeNoCopy(Xdim,Xdim,Xdim);
    Coef_z.resizeNoCopy(Xdim,Xdim,Xdim);
    planeCoef.resizeNoCopy(Xdim,Xdim);
    inputPlane.resizeNoCopy(Xdim,Xdim);
    inputRow.resizeNoCopy(Xdim);

    for (long l = 0; l<Xdim; l++)
    {
        for (long m = 0; m < Xdim; m++)
        {
            CopyDoubleToDouble(MULTIDIM_ARRAY(V),             Xdim, Xdim, Xdim,  l, 0L,  m,
                               MULTIDIM_ARRAY(inputRow),        1L, Xdim,   1L, 0L, 0L, 0L,
                               1L, Xdim, 1L);
            memcpy(&DIRECT_A2D_ELEM(inputPlane,m,0),
                   MULTIDIM_ARRAY(inputRow),Xdim*sizeof(double));
        }

        ChangeBasisVolume(MULTIDIM_ARRAY(inputPlane), MULTIDIM_ARRAY(planeCoef),
                          Xdim, Xdim, 1L, CardinalSpline,
                          BasicSpline, 3L, FiniteCoefficientSupport, DBL_EPSILON, &Status);

        CopyDoubleToDouble(MULTIDIM_ARRAY(planeCoef), Xdim, Xdim,   1L, 0L, 0L, 0L,
                           MULTIDIM_ARRAY(Coef_x),    Xdim, Xdim, Xdim, 0L, 0L,  l,
                           Xdim, Xdim, 1L);
    }

    for (long l = 0; l < Xdim; l++)
    {
        for (long m = 0; m < Xdim; m++)
        {
            CopyDoubleToDouble(MULTIDIM_ARRAY(V),             Xdim, Xdim, Xdim, 0L,  l,  m,
                               MULTIDIM_ARRAY(inputRow),      Xdim,   1L,   1L, 0L, 0L, 0L,
                               Xdim, 1L, 1L);
            CopyDoubleToDouble(MULTIDIM_ARRAY(inputRow),   Xdim,   1L, 1L, 0L, 0L, 0L,
                               MULTIDIM_ARRAY(inputPlane), Xdim, Xdim, 1L, 0L,  m, 0L,
                               Xdim, 1L, 1L);
        }

        ChangeBasisVolume(MULTIDIM_ARRAY(inputPlane), MULTIDIM_ARRAY(planeCoef),
                          Xdim, Xdim, 1L, CardinalSpline,
                          BasicSpline, 3L, FiniteCoefficientSupport, DBL_EPSILON, &Status);
        CopyDoubleToDouble(MULTIDIM_ARRAY(planeCoef), Xdim, Xdim,   1L, 0L, 0L, 0L,
                           MULTIDIM_ARRAY(Coef_y),    Xdim, Xdim, Xdim, 0L, 0L,  l,
                           Xdim, Xdim, 1L);
    }

    for (long l = 0L; l < Xdim; l++)
    {
        for (long m = 0L; m < Xdim; m++)
        {
            CopyDoubleToDouble(MULTIDIM_ARRAY(V),             Xdim, Xdim, Xdim, 0L, m, l,
                               MULTIDIM_ARRAY(inputRow),      Xdim,   1L,   1L, 0L, 0L, 0L,
                               Xdim, 1L, 1L);
            CopyDoubleToDouble(MULTIDIM_ARRAY(inputRow),   Xdim,   1L, 1L, 0L, 0L, 0L,
                               MULTIDIM_ARRAY(inputPlane), Xdim, Xdim, 1L, 0L,  m, 0L,
                               Xdim, 1L, 1L);
        }

        ChangeBasisVolume(MULTIDIM_ARRAY(inputPlane), MULTIDIM_ARRAY(planeCoef),
                          Xdim, Xdim, 1L, CardinalSpline,
                          BasicSpline, 3L, FiniteCoefficientSupport, DBL_EPSILON, &Status);
        CopyDoubleToDouble(MULTIDIM_ARRAY(planeCoef), Xdim, Xdim,   1L, 0L, 0L, 0L,
                           MULTIDIM_ARRAY(Coef_z),    Xdim, Xdim, Xdim, 0L, 0L,  l,
                           Xdim, Xdim, 1L);
    }

    Ac.initIdentity(4);
    Acinv.initIdentity(4);
    double halfSize=Xdim/2;
    MAT_ELEM(Ac,0,3)=MAT_ELEM(Ac,1,3)=MAT_ELEM(Ac,2,3)=halfSize;
    MAT_ELEM(Acinv,0,3)=MAT_ELEM(Acinv,1,3)=MAT_ELEM(Acinv,2,3)=-halfSize;
}

/// Transforms angles from (Ry, Rz, Ry) to (Rx, Ry, Rz) system. Returns possible error.
void convertAngles(double &phi, double &theta, double &psi)
{
    Matrix2D<double> E;
    Euler_angles2matrix(phi,theta,psi,E,false);
    double A00=MAT_ELEM(E,0,0);
    double A02=MAT_ELEM(E,0,2);
    double A10=MAT_ELEM(E,1,0);
    double A12=MAT_ELEM(E,1,2);
    double A20=MAT_ELEM(E,2,0);
    double A21=MAT_ELEM(E,2,1);
    double A22=MAT_ELEM(E,2,2);

    double abs_cosay = sqrt(A22*A22+A21*A21);

    //Results angles
    double ax, ay, az;

    if(abs_cosay > 0.0)
    {
        double sign_cosay;

        ax = atan2(-A21, A22);
        az = atan2(-A10, A00);

        if(fabs(cos(ax)) == 0.0)
            sign_cosay = SGN(-A21/sin(ax));
        else
        {
            if (cos(ax) > 0.0)
                sign_cosay = SGN(A22);
            else
                sign_cosay = -SGN(A22);
        }

        ay  = atan2(A20, sign_cosay * abs_cosay);
    }
    else
    {
        //Let's consider the matrix as a rotation around Z
        if (SGN(A20) > 0.0)
        {
            ax = 0.0;
            ay  = PI/2.0;
            az = atan2(A12, -A02);
        }
        else
        {
            ax = 0.0;
            ay  = -PI/2.0;
            az = atan2(-A12, A02);
        }
    }

    phi   = ax;
    theta = ay;
    psi   = az;
}

//-------------------------------------------------------------------------
// Second level of projection
void projectionRealShears2(MultidimArray<double> &CoefVolume,
                           double absscale,
                           Matrix2D<double> &Binv,
                           Matrix1D<double> &BinvCscaled,
                           int *arr,
                           MultidimArray<double> &projection)
{
    double  K[4], Arg[4], Ki[4];
    int Xdim=XSIZE(projection);
    double *ptrProjection=MULTIDIM_ARRAY(projection);
    size_t CC1 = Xdim * Xdim;

    Ki[0]=dMn(Binv,3);
    Ki[1]=dMn(Binv,7);
    Ki[2]=dMn(Binv,11);
    Ki[3]=dMn(Binv,15);
    for (int i = 0; i < Xdim; i++)
    {
        K[0]=Ki[0];
        K[1]=Ki[1];
        K[2]=Ki[2];
        K[3]=Ki[3];
        for (int n = 0; n < Xdim; n++)
        {
            double Proj = 0.0;
            for (int ksi = 0; ksi < Xdim; ksi++)
            {
                size_t CC2 = CC1 * ksi;
                double sc = (double) ksi - K[arr[0]];
                Arg[0]=VEC_ELEM(BinvCscaled,0)*sc+K[0];
                Arg[1]=VEC_ELEM(BinvCscaled,1)*sc+K[1];
                Arg[2]=VEC_ELEM(BinvCscaled,2)*sc+K[2];
                Arg[3]=VEC_ELEM(BinvCscaled,3)*sc+K[3];
                double g = Arg[arr[1]];
                double h = Arg[arr[2]];

                int l1 = (int)ceil(g-2.0);
                int l2 = l1 + 3;
                int m1 = (int)ceil(h-2.0);
                int m2 = m1 + 3;
                double columns = 0.0;
                double aux;
                for (int m = m1; m <= m2; m++)
                {
                    if (m < Xdim && m > -1)
                    {
                        size_t CC3 = CC2 + Xdim * m;
                        double rows = 0.0;
                        for (int l = l1; l <= l2; l++)
                        {
                            if ((l < Xdim && l > -1))
                            {
                                double gminusl = g - (double) l;
                                BSPLINE03(aux,gminusl);
                                double Coeff = DIRECT_A1D_ELEM(CoefVolume,CC3 + l);
                                rows += Coeff * aux;
                            }
                        }
                        double hminusm = h - (double) m;
                        BSPLINE03(aux,hminusm);
                        columns +=  rows * aux;
                    }
                }
                Proj += columns;
            }
            *ptrProjection++ = Proj*absscale;
            K[0]+=dMn(Binv,0);
            K[1]+=dMn(Binv,4);
            K[2]+=dMn(Binv,8);
            K[3]+=dMn(Binv,12);
        }
        Ki[0]+=dMn(Binv,1);
        Ki[1]+=dMn(Binv,5);
        Ki[2]+=dMn(Binv,9);
        Ki[3]+=dMn(Binv,13);
    }
}

//-------------------------------------------------------------------------
// First level of projection
void projectionRealShears1(RealShearsInfo &Data,
                           double phi, double theta, double psi,
                           double shiftX, double shiftY,
                           MultidimArray<double> &projection)
{
    Matrix2D<double> R;
    R.initIdentity(4);
    MAT_ELEM(R,0,3)=shiftX;
    MAT_ELEM(R,1,3)=shiftY;
    Matrix2D<double> B=Data.Ac*R;
    rotation3DMatrix(-RAD2DEG(phi),'X',R,true);
    B=B*R;
    rotation3DMatrix(RAD2DEG(theta),'Y',R,true);
    B=B*R;
    rotation3DMatrix(-RAD2DEG(psi),'Z',R,true);
    B=B*R;
    B=B*Data.Acinv;
    Matrix2D<double> Binv;
    B.inv(Binv);
    Matrix1D<double> BinvC;
    Binv.getCol(2,BinvC);

    double  scale_x, scale_y, scale_z;
    if (XX(BinvC) != 0.0)
        scale_x = 1.0 / XX(BinvC);
    else
        scale_x = 1e38;

    if (YY(BinvC) != 0.0)
        scale_y = 1.0 / YY(BinvC);
    else
        scale_y = 1e38;

    if (ZZ(BinvC) != 0.0)
        scale_z = 1.0 / ZZ(BinvC);
    else
        scale_z = 1e38;
    double m_x = fabs(scale_x);
    double m_y = fabs(scale_y);
    double m_z = fabs(scale_z);

    int     arr[3];
    double minm = m_x;
    double scale = scale_x;
    arr[0] = 0;
    arr[1] = 1;
    arr[2] = 2;
    MultidimArray<double> *Coef_xyz = &Data.Coef_x;
    if (m_y < minm)
    {
        minm = m_y;
        scale = scale_y;
        arr[0] = 1;
        arr[1] = 0;
        arr[2] = 2;
        Coef_xyz = &Data.Coef_y;
    }
    if (m_z < minm)
    {
        minm = m_z;
        scale = scale_z;
        arr[0] = 2;
        arr[1] = 0;
        arr[2] = 1;
        Coef_xyz = &Data.Coef_z;
    }
    Matrix1D<double> BinvCscaled=BinvC;
    BinvCscaled*=scale;

    projectionRealShears2(*Coef_xyz, minm, Binv, BinvCscaled, arr, projection);
}

void projectVolume(RealShearsInfo &Data, Projection &P, int Ydim, int Xdim,
                    double rot, double tilt, double psi, double shiftX, double shiftY)
{
    P.reset(Data.Xdim,Data.Xdim);
    double Phi=-psi;
    double Theta=-tilt;
    double Psi=-rot;
    convertAngles(Phi, Theta, Psi);
    projectionRealShears1(Data, Phi, Theta, Psi, shiftX, shiftY, P());
    if (Ydim!=Data.Xdim || Xdim!=Data.Xdim)
        P().selfWindow(FIRST_XMIPP_INDEX(Ydim),FIRST_XMIPP_INDEX(Xdim),
                       LAST_XMIPP_INDEX(Ydim),LAST_XMIPP_INDEX(Xdim));
}
