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

#include "projection_real_shears.h"

/// Transforms angles from (Ry, Rz, Ry) to (Rx, Ry, Rz) system. Returns possible error.
void angles_transcription(Matrix1D<double> &angles)
{
    double phi = VEC_ELEM(angles,0);
    double theta = VEC_ELEM(angles,1);
    double psi = VEC_ELEM(angles,2);

    Matrix2D<double> E;
    Euler_angles2matrix(-RAD2DEG(phi),RAD2DEG(theta),-RAD2DEG(psi),E,false);
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

        if(abs(cos(ax)) == 0.0)
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

    VEC_ELEM(angles,0) = ax;
    VEC_ELEM(angles,1) = ay;
    VEC_ELEM(angles,2) = az;
}

//-----------------------------------------------------------------------------------------------
///Computes one iteration of the projection. The resulting projection is into the pointer parameter called "Projection".\n
///Returns possible error.
void do_compute_projection   (int Xdim,
                              double *CoefVolume,
                              double absscale,
                              double *Binv,
                              double *BinvCscaled,
                              int *arr,
                              double *Projection)
{
    long    i, n, l, l1, l2, m, m1, m2, ksi, CC1, CC2, CC3, row, column, index;
    double  X[4], K[4], Arg[4];
    double  Proj, sc, g, h, rows, columns, Coeff, Difference;
    double  gminusl, hminusm;

    CC1 = Xdim * Xdim;

    X[2]=0.0;
    X[3]=1.0;
    for (i = 0L; i < Xdim; i++)
    {
        X[1]=i;
        for (n = 0L; n < Xdim; n++)
        {
            X[0]=n;
            MatrixTimesVector(Binv, X, K, 4L, 4L);
            Proj = 0.0;
            for (ksi = 0L; ksi < Xdim; ksi++)
            {
                CC2 = CC1 * ksi;
                sc = (double) ksi - K[arr[0]];
                for (int ii=0; ii<4; ++ii)
                    Arg[ii]=BinvCscaled[ii]*sc+K[ii];
                g = Arg[arr[1]];
                h = Arg[arr[2]];

                double aux=g-2.0;
                l1 = CEIL(aux);
                l2 = l1 + 3L;
                aux=h-2.0;
                m1 = CEIL(aux);
                m2 = m1 + 3L;
                columns = 0.0;
                for (m = m1; m <= m2; m++)
                {
                    if (m < Xdim && m > -1L)
                    {
                        CC3 = CC2 + Xdim * m;
                        rows = 0.0;
                        for (l = l1; l <= l2; l++)
                        {
                            if ((l < Xdim && l > -1L))
                            {
                                gminusl = g - (double) l;
                                double aux;
                                BSPLINE03(aux,gminusl);
                                Coeff = (double) CoefVolume[CC3 + l];
                                rows += Coeff * aux;
                            }
                        }
                        hminusm = h - (double) m;
                        double aux;
                        BSPLINE03(aux,hminusm);
                        columns +=  rows * aux;
                    }
                }
                Proj += columns;
            }
            Proj *= absscale;

            *Projection++ = Proj;
        }
    }
}/* End of do_compute_projection */

//-----------------------------------------------------------------------------------------------
///Computes projection. The resulting projection is into the pointer parameter called "Projection".\n
///Returns possible error.
void Compute_projection(const VolumeStruct &Data,
						const Matrix1D<double> &angles,
						const Matrix1D<double> &shifts,
                       double *Coef_x,
                       double *Coef_y,
                       double *Coef_z,
                       const Matrix2D<double> &RightOperHlp,
                       const Matrix2D<double> &LeftOperHlp,
                       double *projection)
{
    int     Status=!ERROR, arr[3];
    double  scale, scale_x, scale_y, scale_z, m_x, m_y, m_z, minm;

    double  psi, theta, phi, Sinphi, Cosphi, Sinpsi, Cospsi, Sintheta, Costheta;
    psi   = VEC_ELEM(angles,0);
    theta = VEC_ELEM(angles,1);
    phi   = VEC_ELEM(angles,2);
    sincos(phi,&Sinphi,&Cosphi);
    sincos(theta,&Sintheta,&Costheta);
    sincos(psi,&Sinpsi,&Cospsi);

    Matrix2D<double> At;
    translation3DMatrix(shifts,At);

    Matrix2D<double> Help1=At*LeftOperHlp;
    Matrix2D<double> R;
    rotation3DMatrix(RAD2DEG(phi),'X',R,true);

    double *Help2 = (double *)malloc((size_t) 16L * sizeof(double));
    if (Help2 == (double *)NULL)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::Compute_projection: "
                     "ERROR - Not enough memory for Help2");

    if (MatrixMultiply(MATRIX2D_ARRAY(Help1), MATRIX2D_ARRAY(R), Help2, 4L, 4L, 4L) == ERROR)
        REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::Compute_projection: "
                     "Error returned by MatrixMultiply");

    if (GetIdentitySquareMatrix(MATRIX2D_ARRAY(R), 4L) == ERROR)
        REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::Compute_projection: "
                     "Error returned by GetIdentitySquareMatrix");

    double *hlp = MATRIX2D_ARRAY(R);
    *hlp = Costheta;
    hlp += (ptrdiff_t)2L;
    *hlp = Sintheta;
    hlp += (ptrdiff_t)6L;
    *hlp = - Sintheta;
    hlp += (ptrdiff_t)2L;
    *hlp = Costheta;

    double *Help3 = (double *)malloc((size_t) 16L * sizeof(double));
    if (Help3 == (double *)NULL)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::Compute_projection: "
                     "ERROR - Not enough memory for Help3");

    if (MatrixMultiply(Help2, MATRIX2D_ARRAY(R), Help3, 4L, 4L, 4L) == ERROR)
        REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::Compute_projection: "
                     "Error returned by MatrixMultiply");

    if (GetIdentitySquareMatrix(MATRIX2D_ARRAY(R), 4L) == ERROR)
        REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::Compute_projection: "
                     "Error returned by GetIdentitySquareMatrix");

    hlp = MATRIX2D_ARRAY(R);
    *hlp++ = Cospsi;
    *hlp = - Sinpsi;
    hlp += (ptrdiff_t)3L;
    *hlp++ = Sinpsi;
    *hlp = Cospsi;

    double *Help4 = (double *)malloc((size_t) 16L * sizeof(double));
    if (Help4 == (double *)NULL)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::Compute_projection: "
                     "ERROR - Not enough memory for Help4");

    if (MatrixMultiply(Help3, MATRIX2D_ARRAY(R), Help4, 4L, 4L, 4L) == ERROR)
        REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::Compute_projection: "
                     "Error returned by MatrixMultiply");

    Matrix2D<double> B;
    B.resizeNoCopy(4,4);
    if (MatrixMultiply(Help4, MATRIX2D_ARRAY(RightOperHlp), MATRIX2D_ARRAY(B), 4L, 4L, 4L) == ERROR)
        REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::Compute_projection: "
                     "Error returned by MatrixMultiply");
    free(Help4);

    double *Binv = (double *)malloc((size_t) 16L * sizeof(double));
    if (Binv == (double *)NULL)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::Compute_projection: "
                     "ERROR - Not enough memory for Binv");

    if (SquareMatrixInvertGauss(MATRIX2D_ARRAY(B), Binv, 4L, DBL_EPSILON, &Status) == ERROR)
        REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::Compute_projection: "
                     "Error returned by SquareMatrixInvertGauss");

    double *C1 = (double *)malloc((size_t) 4L * sizeof(double));
    if (C1 == (double *)NULL)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::Compute_projection: "
                     "ERROR - Not enough memory for C1");

    hlp = C1;
    *hlp++ = 0.0;
    *hlp++ = 0.0;
    *hlp++ = 1.0;
    *hlp = 0.0;

    double *BinvC = (double *)malloc((size_t) 4L * sizeof(double));
    if (BinvC == (double *)NULL)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::Compute_projection: "
                     "ERROR - Not enough memory for BinvC");

    if (MatrixTimesVector(Binv, C1, BinvC, 4L, 4L) == ERROR)
        REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::Compute_projection: "
                     "Error returned by MatrixTimesVector");

    if (BinvC[0] != 0.0)
        scale_x = 1.0 / BinvC[0];
    else
        scale_x = DBL_MAX;

    if (BinvC[1] != 0.0)
        scale_y = 1.0 / BinvC[1];
    else
        scale_y = DBL_MAX;

    if (BinvC[2] != 0.0)
        scale_z = 1.0 / BinvC[2];
    else
        scale_z = DBL_MAX;

    m_x = fabs(scale_x);
    m_y = fabs(scale_y);
    m_z = fabs(scale_z);

    double *BinvCscaled = (double *)malloc((size_t) 4L * sizeof(double));
    if (BinvCscaled == (double *)NULL)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::Compute_projection: "
                     "ERROR - Not enough memory for BinvCscaled");

    minm = m_x;
    scale = scale_x;
    if (VectorScale(BinvC, BinvCscaled, scale_x, 4L) == ERROR)
        REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::Compute_projection: "
                     "Error returned by VectorScale");
    int Xdim=XSIZE(*Data.volume);
    arr[0] = 0;
    arr[1] = 1;
    arr[2] = 2;
    double *Coef_xyz = Coef_x;
    if (m_y < minm)
    {
        minm = m_y;
        scale = scale_y;
        if (VectorScale(BinvC, BinvCscaled, scale_y, 4L) == ERROR)
            REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::Compute_projection: "
                         "Error returned by VectorScale");
        arr[0] = 1;
        arr[1] = 0;
        arr[2] = 2;
        Coef_xyz = Coef_y;
    }
    if (m_z < minm)
    {
        minm = m_z;
        scale = scale_z;
        VectorScale(BinvC, BinvCscaled, scale_z, 4L);
        arr[0] = 2;
        arr[1] = 0;
        arr[2] = 1;
        Coef_xyz = Coef_z;
    }

    free(BinvC);
    free(C1);

    do_compute_projection(Xdim, Coef_xyz, minm, Binv, BinvCscaled, arr, projection);

    free(BinvCscaled);
    free(Binv);
    free(Help2);
    free(Help3);
}/* End of Compute_projection */

///Main compute function. Returns possible error.
void do_one_projection(VolumeStruct &Data2)
{
    int    Status = !ERROR;
    long    i, m, n, l;
    double    lambda;
    double    *VolumeCoef, *InputVolume, *InputVolumePlane;
    double    *InputVolumeRow, *Coef_x, *Coef_y, *Coef_z;

    angles_transcription(Data2.angles);

    int Xdim = XSIZE(*Data2.volume);

    InputVolume = MULTIDIM_ARRAY(*Data2.volume);

    AllocateVolumeDouble( &Coef_x, Xdim, Xdim, Xdim, &Status);
    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for Coef_x");

    AllocateVolumeDouble( &Coef_y, Xdim, Xdim, Xdim, &Status);
    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for Coef_y");

    AllocateVolumeDouble( &Coef_z, Xdim, Xdim, Xdim, &Status);
    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for Coef_z");


    AllocateVolumeDouble( &VolumeCoef, Xdim, Xdim, 1L, &Status);
    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for VolumeCoef");

    AllocateVolumeDouble( &InputVolumePlane, Xdim, Xdim, 1L, &Status);
    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for InputVolumePlane");

    AllocateVolumeDouble( &InputVolumeRow, 1L, Xdim, 1L, &Status);
    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for InputVolumePlane");

    for (l = 0L; l < Xdim; l++)
    {

        for (m = 0L; m < Xdim; m++)
        {
            if (CopyDoubleToDouble(InputVolume, Xdim, Xdim, Xdim, l,    0L,    m,
                                   InputVolumeRow, 1L, Xdim, 1L, 0L, 0L, 0L, 1L, Xdim, 1L) == ERROR)
                REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::ROUT_project_execute: "
                             "Error returned by CopyDoubleToDouble");
            for (i = 0L; i < Xdim; i++)
            {
                InputVolumePlane[Xdim * m + i] = InputVolumeRow[i];
            }
        }

        ChangeBasisVolume(InputVolumePlane, VolumeCoef, Xdim, Xdim, 1L, CardinalSpline,
                          BasicSpline, 3L, FiniteCoefficientSupport, DBL_EPSILON, &Status);
        if (Status == ERROR)
            REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::ROUT_project_execute: "
                         "ERROR");

        if (CopyDoubleToDouble(VolumeCoef, Xdim, Xdim, 1L, 0L,    0L,    0L,
                               Coef_x, Xdim, Xdim, Xdim, 0L, 0L, l, Xdim, Xdim, 1L) == ERROR)
            REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::ROUT_project_execute: "
                         "Error returned by CopyDoubleToDouble");
    }

    FreeVolumeDouble(&VolumeCoef);
    FreeVolumeDouble(&InputVolumePlane);
    FreeVolumeDouble(&InputVolumeRow);

    AllocateVolumeDouble( &VolumeCoef, Xdim, Xdim, 1L, &Status);
    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for VolumeCoef");

    AllocateVolumeDouble( &InputVolumePlane, Xdim, Xdim, 1L, &Status);
    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for InputVolumePlane");

    AllocateVolumeDouble( &InputVolumeRow, Xdim, 1L, 1L, &Status);
    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for InputVolumePlane");

    for (l = 0L; l < Xdim; l++)
    {
        for (m = 0L; m < Xdim; m++)
        {
            if (CopyDoubleToDouble(InputVolume, Xdim, Xdim, Xdim, 0L, l,    m,
                                   InputVolumeRow, Xdim, 1L, 1L, 0L, 0L, 0L, Xdim, 1L, 1L) == ERROR)
                REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::ROUT_project_execute: "
                             "Error returned by CopyDoubleToDouble");

            if (CopyDoubleToDouble(InputVolumeRow, Xdim, 1L, 1L, 0L,    0L,    0L,
                                   InputVolumePlane, Xdim, Xdim, 1L, 0L, m, 0L, Xdim, 1L, 1L) == ERROR)
                REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::ROUT_project_execute: "
                             "Error returned by CopyDoubleToDouble");

        }

        ChangeBasisVolume((double*)InputVolumePlane, (double*)VolumeCoef, Xdim, Xdim, 1L, CardinalSpline,
                          BasicSpline, 3L, FiniteCoefficientSupport, DBL_EPSILON, &Status);
        if (Status == ERROR)
            REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::ROUT_project_execute: "
                         "ERROR");

        if (CopyDoubleToDouble(VolumeCoef, Xdim, Xdim, 1L, 0L,    0L,    0L,
                               Coef_y, Xdim, Xdim, Xdim, 0L, 0L, l, Xdim, Xdim, 1L) == ERROR)
            REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::ROUT_project_execute: "
                         "Error returned by CopyDoubleToDouble");
    }

    FreeVolumeDouble(&VolumeCoef);
    FreeVolumeDouble(&InputVolumePlane);
    FreeVolumeDouble(&InputVolumeRow);

    AllocateVolumeDouble( &VolumeCoef, Xdim, Xdim, 1L, &Status);
    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for VolumeCoef");

    AllocateVolumeDouble( &InputVolumePlane, Xdim, Xdim, 1L, &Status);
    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for InputVolumePlane");

    AllocateVolumeDouble( &InputVolumeRow, Xdim, 1L, 1L, &Status);
    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for InputVolumePlane");

    for (l = 0L; l < Xdim; l++)
    {
        for (m = 0L; m < Xdim; m++)
        {
            if (CopyDoubleToDouble(InputVolume, Xdim, Xdim, Xdim, 0L, m,    l,
                                   InputVolumeRow, Xdim, 1L, 1L, 0L, 0L, 0L, Xdim, 1L, 1L) == ERROR)
                REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::ROUT_project_execute: "
                             "Error returned by CopyDoubleToDouble");

            if (CopyDoubleToDouble(InputVolumeRow, Xdim, 1L, 1L, 0L,    0L,    0L,
                                   InputVolumePlane, Xdim, Xdim, 1L, 0L, m, 0L, Xdim, 1L, 1L) == ERROR)
                REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::ROUT_project_execute: "
                             "Error returned by CopyDoubleToDouble");
        }

        ChangeBasisVolume(InputVolumePlane, VolumeCoef, Xdim, Xdim, 1L, CardinalSpline,
                          BasicSpline, 3L, FiniteCoefficientSupport, DBL_EPSILON, &Status);
        if (Status == ERROR)
            REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::ROUT_project_execute: "
                         "ERROR");

        if (CopyDoubleToDouble(VolumeCoef, Xdim, Xdim, 1L, 0L,    0L,    0L,
                               Coef_z, Xdim, Xdim, Xdim, 0L, 0L, l, Xdim, Xdim, 1L) == ERROR)
            REPORT_ERROR(ERR_NUMERICAL, "Projection_real_shears::ROUT_project_execute: "
                         "Error returned by CopyDoubleToDouble");
    }

    FreeVolumeDouble(&VolumeCoef);
    FreeVolumeDouble(&InputVolumePlane);
    FreeVolumeDouble(&InputVolumeRow);

    if (Status == ERROR)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Projection_real_shears::ROUT_project_execute: "
                     "ERROR - Not enough memory for Projection");

    Matrix2D<double> Ac, Acinv;
    Ac.initIdentity(4);
    Acinv.initIdentity(4);
    double halfSize=XSIZE(*Data2.volume)/2;
    MAT_ELEM(Ac,0,3)=MAT_ELEM(Ac,1,3)=MAT_ELEM(Ac,2,3)=halfSize;
    MAT_ELEM(Acinv,0,3)=MAT_ELEM(Acinv,1,3)=MAT_ELEM(Acinv,2,3)=-halfSize;

    Compute_projection(Data2, Data2.angles, Data2.shifts,
    			       Coef_x, Coef_y, Coef_z, Acinv, Ac,
                       MULTIDIM_ARRAY(Data2.projection()));

    FreeVolumeDouble(&Coef_x);
    FreeVolumeDouble(&Coef_y);
    FreeVolumeDouble(&Coef_z);
}

///Parameters reading. Note that all parameters are required.
void Projection_real_shears::read(int argc, char **argv)
{
    fn_proj_param = getParameter(argc, argv, "-i");
    fn_sel_file   = getParameter(argc, argv, "-o", "");
    display = !checkParameter(argc, argv, "-quiet");
}

///Description of the projection_real_shears function.
void Projection_real_shears::usage()
{
    printf("\nUsage:\n\n");
    printf("projection_real_shears -i <Parameters File> \n"
           "                      [-o <sel_file>]\n"
           "                      [-quiet]\n");
    printf("\tWhere:\n"
           "\t<Parameters File>:  File containing projection parameters\n"
           "\t                    Note that only the top-four parameters lines are read\n"
           "\t                    Check the manual for a description of the parameters\n"
           "\t<sel_file>:         This is a selection file with all the generated\n"
           "\t                    projections.\n");
}

//-----------------------------------------------------------------------------------------------
///Reads the projection parameters of the input file and inserts it into Projection_Parameters fields.\n
///This is an "overloaded" function in order to use translation parameters.
void Projection_real_shears::read(const FileName &fn_proj_param)
{
    FILE    *fh_param;
    char    line[201];
    int     lineNo = 0;

    if ((fh_param = fopen(fn_proj_param.c_str(), "r")) == NULL)
        REPORT_ERROR(ERR_IO_NOTOPEN,
                     (std::string)"Projection_real_shears::read: There is a problem "
                     "opening the file " + fn_proj_param);

    while (fgets(line, 200, fh_param) != NULL)
    {
        if (line[0] == 0)
            continue;
        if (line[0] == '#')
            continue;
        if (line[0] == '\n')
            continue;
        switch (lineNo)
        {
        case 0: //***** Line 1 *****
            //Volume file
            fnPhantom = firstWord(line);

            if (!exists(fnPhantom))
                REPORT_ERROR(ERR_IO_NOTEXIST, (std::string)"Projection_real_shears::read: "
                             "file " + fnPhantom + " doesn't exist");

            lineNo = 1;
            break;
        case 1: //***** Line 2 *****
            fnProjectionSeed = firstWord(line);

            char    *auxstr;

            auxstr = nextToken();
            if (auxstr != NULL)
                starting =
                    textToInteger(auxstr);

            fn_projection_extension = nextWord();

            lineNo = 2;
            break;
        case 2: //***** Line 3 *****
            proj_Xdim = textToInteger(firstToken(line));
            proj_Ydim = proj_Xdim ;

            lineNo = 3;
            break;
        case 3: //***** Line 4 *****
            // Angle file
            fn_angle = firstWord(line);

            if (!exists(fn_angle))
                REPORT_ERROR(ERR_IO_NOTEXIST, (std::string)"Projection_real_shears::read: "
                             "file " + fn_angle + " doesn't exist");

            lineNo = 4;
            break;
        default:
            break;
        } // switch end
    } // while end

    if (lineNo != 4) //If all parameters was not read
        REPORT_ERROR(ERR_PARAM_MISSING, (std::string)"Projection_real_shears::read: I "
                     "couldn't read all parameters from file " + fn_proj_param);
    fclose(fh_param);
}

///Writes the projection file obtained. Returns possibles errors.
void Projection_real_shears::write_projection_file(int numFile)
{
    FileName fn_proj;
    fn_proj.compose(fnProjectionSeed, numFile, fn_projection_extension);
    SF.setValue(MDL_IMAGE,fn_proj,SF.addObject());
    Data->projection.write(fn_proj);
    if(display)
        progress_bar(numFile);
}

///Reads a DocLine and fills Data fields. Returns possibles errors.
void Projection_real_shears::read_a_DocLine(size_t objId)
{
    double rot     ;
    double tilt    ;
    double psi     ;
    double shiftX=0;
    double shiftY=0;

    DF.getValue(MDL_SHIFTX,shiftX,objId);
    DF.getValue(MDL_SHIFTY,shiftY,objId);
    DF.getValue(MDL_ANGLEROT,rot,objId);
    DF.getValue(MDL_ANGLETILT,tilt,objId);
    DF.getValue(MDL_ANGLEPSI,psi,objId);

    VECTOR_R3(Data->angles,-DEG2RAD(psi),-DEG2RAD(tilt),-DEG2RAD(rot));
    VECTOR_R2(Data->shifts,shiftX,shiftY);
}

///////////////////////// MAIN INSTRUCTION FOR MPI ////////////////////////////////
void project_Volume(VolumeStruct &Data, Projection &P, int Ydim, int Xdim,
                    double rot, double tilt, double psi)
{
    // Prepare Data Structure
    VECTOR_R3(Data.angles,-DEG2RAD(psi),-DEG2RAD(tilt),-DEG2RAD(rot));
    Data.shifts.initZeros();
    do_one_projection(Data);
    P=Data.projection;
    P().selfWindow(FIRST_XMIPP_INDEX(Ydim),FIRST_XMIPP_INDEX(Xdim),
                   LAST_XMIPP_INDEX(Ydim),LAST_XMIPP_INDEX(Xdim));
}

VolumeStruct::VolumeStruct(const MultidimArray<double> &V)
{
    volume=&V;
    if (XSIZE(V)!=YSIZE(V) || XSIZE(V)!=ZSIZE(V))
        REPORT_ERROR(ERR_MULTIDIM_DIM, "The volume must be cubic");
    projection.reset(YSIZE(V),XSIZE(V));
    projection().setXmippOrigin();
    shifts.initZeros(3);
    angles.initZeros(3);
}

///Does start instructions. Returns possibles errors.
void Projection_real_shears::start_to_process()
{
    read(fn_proj_param);
    DF.read(fn_angle);
    V.read(fnPhantom);
    Data=new VolumeStruct(V());

    //Configure the time for the progress bar
    if(display)
    {
        time_config();
        init_progress_bar(DF.size());
    }
}

///Does finish instructions. Returns possibles errors.
void Projection_real_shears::finish_to_process()
{
    if(display)
        progress_bar(DF.size());
    if (fn_sel_file == "") //If the name of the output file is not specified
        fn_sel_file = "sel"+fnProjectionSeed+".sel";
    SF.write(fn_sel_file);
    delete Data;
}

//-------------------------------------- Main function ----------------------------------------
void Projection_real_shears::ROUT_project_real_shears()
{
    start_to_process();
    int num_file = starting;
    FOR_ALL_OBJECTS_IN_METADATA(DF)
    {
        read_a_DocLine(__iter.objId);
        do_one_projection(*Data);
        write_projection_file(num_file++);
    }
    finish_to_process();
}
