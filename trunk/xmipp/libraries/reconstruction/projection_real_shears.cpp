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
void projectionRealShears3   (MultidimArray<double> &CoefVolume,
                              double absscale,
                              double *Binv,
                              double *BinvCscaled,
                              int *arr,
                              MultidimArray<double> &projection)
{
    long    i, n, l, l1, l2, m, m1, m2, ksi, row, column, index;
    double  X[4], K[4], Arg[4];
    double  Proj, sc, g, h, rows, columns, Coeff;
    double  gminusl, hminusm;

    int Xdim=XSIZE(projection);
    double *ptrProjection=MULTIDIM_ARRAY(projection);
    size_t CC1 = Xdim * Xdim;

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
                size_t CC2 = CC1 * ksi;
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
                        size_t CC3 = CC2 + Xdim * m;
                        rows = 0.0;
                        for (l = l1; l <= l2; l++)
                        {
                            if ((l < Xdim && l > -1L))
                            {
                                gminusl = g - (double) l;
                                double aux;
                                BSPLINE03(aux,gminusl);
                                Coeff = DIRECT_A1D_ELEM(CoefVolume,CC3 + l);
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

            *ptrProjection++ = Proj;
        }
    }
}

//-----------------------------------------------------------------------------------------------
///Computes projection. The resulting projection is into the pointer parameter called "Projection".\n
///Returns possible error.
void projectionRealShears2(const Matrix1D<double> &angles,
                           const Matrix1D<double> &shifts,
                           MultidimArray<double> &Coef_x,
                           MultidimArray<double> &Coef_y,
                           MultidimArray<double> &Coef_z,
                           const Matrix2D<double> &LeftOperHlp,
                           const Matrix2D<double> &RightOperHlp,
                           MultidimArray<double> &projection)
{
    double phi   = VEC_ELEM(angles,0);
    double theta = VEC_ELEM(angles,1);
    double psi   = VEC_ELEM(angles,2);

    Matrix2D<double> R;
    translation3DMatrix(shifts,R);
    Matrix2D<double> B=LeftOperHlp*R;
    rotation3DMatrix(-RAD2DEG(phi),'X',R,true);
    B=B*R;
    rotation3DMatrix(RAD2DEG(theta),'Y',R,true);
    B=B*R;
    rotation3DMatrix(-RAD2DEG(psi),'Z',R,true);
    B=B*R;
    B=B*RightOperHlp;
    Matrix2D<double> Binv;
    B.inv(Binv);
    Matrix1D<double> BinvC;
    Binv.getCol(2,BinvC);

    double  scale_x, scale_y, scale_z;
    if (XX(BinvC) != 0.0)
        scale_x = 1.0 / XX(BinvC);
    else
        scale_x = DBL_MAX;

    if (YY(BinvC) != 0.0)
        scale_y = 1.0 / YY(BinvC);
    else
        scale_y = DBL_MAX;

    if (ZZ(BinvC) != 0.0)
        scale_z = 1.0 / ZZ(BinvC);
    else
        scale_z = DBL_MAX;
    double m_x = fabs(scale_x);
    double m_y = fabs(scale_y);
    double m_z = fabs(scale_z);

    int     arr[3];
    double minm = m_x;
    double scale = scale_x;
    arr[0] = 0;
    arr[1] = 1;
    arr[2] = 2;
    MultidimArray<double> &Coef_xyz = Coef_x;
    if (m_y < minm)
    {
        minm = m_y;
        scale = scale_y;
        arr[0] = 1;
        arr[1] = 0;
        arr[2] = 2;
        Coef_xyz = Coef_y;
    }
    if (m_z < minm)
    {
        minm = m_z;
        scale = scale_z;
        arr[0] = 2;
        arr[1] = 0;
        arr[2] = 1;
        Coef_xyz = Coef_z;
    }
    Matrix1D<double> BinvCscaled=BinvC;
    BinvCscaled*=scale;

    projectionRealShears3(Coef_xyz, minm, MATRIX2D_ARRAY(Binv),
                          MATRIX1D_ARRAY(BinvCscaled), arr, projection);
}

///Main compute function. Returns possible error.
void projectionRealShears1(VolumeStruct &Data2, Projection &P)
{
    int    Status = !ERROR;

    angles_transcription(Data2.angles);
    int Xdim = Data2.Xdim;

    MultidimArray<double> Coef_x, Coef_y, Coef_z, planeCoef, inputPlane, inputRow;
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
            CopyDoubleToDouble(MULTIDIM_ARRAY(*Data2.volume), Xdim, Xdim, Xdim,  l, 0L,  m,
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
            CopyDoubleToDouble(MULTIDIM_ARRAY(*Data2.volume), Xdim, Xdim, Xdim, 0L,  l,  m,
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
            CopyDoubleToDouble(MULTIDIM_ARRAY(*Data2.volume), Xdim, Xdim, Xdim, 0L, m, l,
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

    Matrix2D<double> Ac, Acinv;
    Ac.initIdentity(4);
    Acinv.initIdentity(4);
    double halfSize=XSIZE(*Data2.volume)/2;
    MAT_ELEM(Ac,0,3)=MAT_ELEM(Ac,1,3)=MAT_ELEM(Ac,2,3)=halfSize;
    MAT_ELEM(Acinv,0,3)=MAT_ELEM(Acinv,1,3)=MAT_ELEM(Acinv,2,3)=-halfSize;

    P.reset(Data2.Xdim,Data2.Xdim);
    projectionRealShears2(Data2.angles, Data2.shifts, Coef_x, Coef_y, Coef_z,
                          Ac, Acinv, P());
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

void project_Volume(VolumeStruct &Data, Projection &P, int Ydim, int Xdim,
                    double rot, double tilt, double psi)
{
    VECTOR_R3(Data.angles,-psi,-tilt,-rot);
    Data.shifts.initZeros();
    P.reset(Data.Xdim,Data.Xdim);
    projectionRealShears1(Data,P);
    if (Ydim!=Data.Xdim || Xdim!=Data.Xdim)
        P().selfWindow(FIRST_XMIPP_INDEX(Ydim),FIRST_XMIPP_INDEX(Xdim),
                       LAST_XMIPP_INDEX(Ydim),LAST_XMIPP_INDEX(Xdim));
}

VolumeStruct::VolumeStruct(const MultidimArray<double> &V)
{
    volume=&V;
    if (XSIZE(V)!=XSIZE(V) || XSIZE(V)!=ZSIZE(V))
        REPORT_ERROR(ERR_MULTIDIM_DIM, "The volume must be cubic");
    Xdim=XSIZE(*volume);
    shifts.initZeros(3);
    angles.initZeros(3);
}

///Does start instructions. Returns possibles errors.
void Projection_real_shears::produceSideInfo()
{
    read(fn_proj_param);
    DF.read(fn_angle);
    V.read(fnPhantom);
    Data=new VolumeStruct(V());
}

//-------------------------------------- Main function ----------------------------------------
void Projection_real_shears::ROUT_project_real_shears()
{
    produceSideInfo();
    if(display)
        init_progress_bar(DF.size());
    int num_file = starting;
    Projection P;
    MetaData SFout;
    FileName fn_proj;
    FOR_ALL_OBJECTS_IN_METADATA(DF)
    {
        // Get parameters
        double rot, tilt, psi, shiftX=0, shiftY=0;
        DF.getValue(MDL_SHIFTX,shiftX,__iter.objId);
        DF.getValue(MDL_SHIFTY,shiftY,__iter.objId);
        DF.getValue(MDL_ANGLEROT,rot,__iter.objId);
        DF.getValue(MDL_ANGLETILT,tilt,__iter.objId);
        DF.getValue(MDL_ANGLEPSI,psi,__iter.objId);
        VECTOR_R3(Data->angles,-psi,-tilt,-rot);
        VECTOR_R2(Data->shifts,shiftX,shiftY);

        // Project
        projectionRealShears1(*Data,P);

        // Write projection
        fn_proj.compose(fnProjectionSeed, num_file, fn_projection_extension);
        SFout.setValue(MDL_IMAGE,fn_proj,SFout.addObject());
        P.write(fn_proj);
        if(display)
            progress_bar(num_file);
        num_file++;
    }

    if(display)
        progress_bar(DF.size());
    if (fn_sel_file == "") //If the name of the output file is not specified
        fn_sel_file = "sel"+fnProjectionSeed+".sel";
    SFout.write(fn_sel_file);
    delete Data;
}
