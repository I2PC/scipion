/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 1999 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 *
 *****************************************************************************/

#include <strstream>

#include <data/args.h>
#include <data/geometry.h>
#include <interface/aph.h>

#include "spots2realspace3d.h"
#include "art_crystal.h"

#define twoPI 2*M_PI

void Spot2RealSpace3D_Parameters::read_from_file(const FileName &fnprm)
{
    FILE *fh_param;
    if ((fh_param = fopen(fnprm.c_str(), "r")) == NULL)
        REPORT_ERROR(3005,
                     (string)"Spot2RealSpace3D_Parameters::read: There is a problem "
                     "opening the file " + fnprm);
    try
    {
        fnaph_in = getParameter(fh_param, "APH file", 0, NULL,
                             3007, "Spot2RealSpace3D_Parameters::read: APH File not found");
        fn_out = getParameter(fh_param, "Output file", 0, NULL,
                           3007, "Spot2RealSpace3D_Parameters::read: Output File not found");
        NoCells.resize(3);
        string aux_str = getParameter(fh_param, "Output cells no", 0, "1 1 1");
        istrstream *is = NULL;
        is = new istrstream(aux_str.c_str());
        try
        {
            *is >> XX(NoCells) >> YY(NoCells) >> ZZ(NoCells);
        }
        catch (...)
        {
            REPORT_ERROR(3007,
                         "Spot2RealSpace3D_Parameters::read: Cannot read Output cells no");
        }

        Phase_Shift.resize(3);
        aux_str = getParameter(fh_param, "Phase Shift", 0, "0.0 0.0 0.0");
        delete is;
        is = new istrstream(aux_str.c_str());
        try
        {
            *is >> XX(Phase_Shift) >> YY(Phase_Shift) >> ZZ(Phase_Shift);
        }
        catch (...)
        {
            REPORT_ERROR(3007,
                         "Spot2RealSpace3D_Parameters::read: Cannot read Phase Shift");
        }

        KeepContrast = textToInteger(getParameter(fh_param, "Keep Contrast", 0, "1"));
        Celldim.resize(3);
        aux_str = getParameter(fh_param, "Cell Size");
        delete is;
        is = new istrstream(aux_str.c_str());
        try
        {
            *is >> XX(Celldim) >> YY(Celldim) >> ZZ(Celldim);
        }
        catch (...)
        {
            REPORT_ERROR(3007,
                         "Spot2RealSpace3D_Parameters::read: Cannot read Cell Size");
        }
        delete is;
//      Scale_Factor=textToFloat(getParameter(fh_param,"Scale Factor",0,"-1"));
    }
    catch (Xmipp_error XE)
    {
        cout << XE << endl;
        REPORT_ERROR(3007, (string)"There is an error reading " + fnprm);
    }

    fclose(fh_param);
}


/* Show parameters --------------------------------------------------------- */
ostream& operator << (ostream &o, const Spot2RealSpace3D_Parameters &prm)
{
    o << "APH Input            : " << prm.fnaph_in << endl
    << "Output Spider File   : " << prm.fn_out << endl
    << "Ouput Layout (x,y,z) : " << prm.NoCells.transpose() << endl
    << "Phase_Shift (x,y,z)  : " << prm.Phase_Shift.transpose()  << endl
    << "Cell Size (x,y,z)    : " << prm.Celldim .transpose() << endl
    << "Keep Contrast        : " << prm.KeepContrast << endl;
    return o;
}

/* DFT^-1 ------------------------------------------------------------------ */
void IDFT_3D(const Matrix3D< complex<double> > &FT, Matrix3D<double> &V1)
{
    double x, y, z;
    int s_iz, s_iy, s_ix, e_iz, e_iy, e_ix;
    int s_kz, s_ky, s_kx, e_kz, e_ky, e_kx;
    double x_dim, y_dim, z_dim;
    int X_dim, Y_dim, Z_dim;
    long ix, iy, iz, kx, ky, kz;
    double argx, argy, argz;
    double cosarg, sinarg;
    int dir = 1;//if dir = -1 forward FT, if dir=1 backward FT

    s_iz = V1.startingZ();
    s_iy = V1.startingY();
    s_ix = V1.startingX();
    e_iz = V1.finishingZ();
    e_iy = V1.finishingY();
    e_ix = V1.finishingX();

    s_kz = FT.startingZ();
    s_ky = FT.startingY();
    s_kx = FT.startingX();
    e_kz = FT.finishingZ();
    e_ky = FT.finishingY();
    e_kx = FT.finishingX();

    z_dim = (double)V1.sliceNumber();
    y_dim = (double)V1.rowNumber();
    x_dim = (double)V1.colNumber();

    time_config();
    init_progress_bar(-s_iz + e_iz + 1);


    for (iz = s_iz;iz <= e_iz;iz++)
    {
        progress_bar(iz - s_iz);
        argz = /*- dir * */ 2.0 * PI * iz / z_dim;
        for (iy = s_iy;iy <= e_iy;iy++)
        {
            argy = /* - dir * */2.0 * PI * iy / y_dim;
            for (ix = s_ix;ix <= e_ix;ix++)
            {
                argx = /*- dir * */ 2.0 * PI * ix / x_dim;
                for (kz = s_kz;kz <= e_kz;kz++)
                {
                    for (ky = s_ky;ky <= e_ky;ky++)
                    {
                        for (kx = s_kx;kx <= e_kx;kx++)
                        {
                            if (VOL_ELEM(FT, kz, ky, kx) == (complex<double>)0.0)
                                continue;
                            cosarg = cos(kz * argz
                                         + ky * argy
                                         + kx * argx);
                            sinarg = sin(kz * argz
                                         + ky * argy
                                         + kx * argx);
                            VOL_ELEM(V1, iz, iy, ix) += real(VOL_ELEM(FT, kz, ky, kx)) * cosarg +
                                                        /* if dir use "-"*/ imag(VOL_ELEM(FT, kz, ky, kx)) * sinarg;
//img part                  I[iz][iy][ix] += (x1[kx][ky][kx] * sinarg +
//                          y1[kz][ky][kx] * cosarg);
                        }
                    }
                }//kx,ky,kz
            }//ix
        }//iy
    }//iz
    progress_bar(-s_iz + e_iz + 1);
#ifdef NEVEEVER
    /* Copy the data back */
    double 3m;
    3m = z_dim * y_dim * x_dim;
    if (dir == 1)
    {
        FOR_ALL_ELEMENTS_IN_MATRIX3D(V1())
        {
            x1[iz][iy][ix] = x2[iz][iy][ix] / 3m;
            y1[iz][iy][ix] = y2[iz][iy][ix] / 3m;
        }
    }
    else
    {
        for (i = 0;i < m;i++)
        {
            x1[iz][iy][ix] = x2[iz][iy][ix];
            y1[iz][iy][ix] = y2[iz][iy][ix];
        }
    }
#endif

}
//void Spot2RealSpace3D_Parameters::produce_SideInfo() {
//    FOR_ALL_ELEMENTS_IN_MATRIX3D(aph_file.spots_abs) {//l,k,h
//       if(aph_file.spots_abs(k,i,j)!=0.)
//         cout << j <<" " << i << " " << k << " " <<aph_file.spots_abs(k,i,j) << endl;
//    }
//}
/* Main routine for transforming ------------------------------------------- */
void ROUT_Spots2RealSpace_3D(Spot2RealSpace3D_Parameters &prm,
                             VolumeXmipp &V1)
{
    prm.aph_file.read_from_prepmklcf(prm.fnaph_in);
    cout << prm;
    // Apply phase shift
    double phase_shift = (prm.KeepContrast) ? 0 : 180;
    FOR_ALL_ELEMENTS_IN_MATRIX3D(prm.aph_file.spots_arg)
    {
        VOL_ELEM(prm.aph_file.spots_arg, k, i, j)  = (double) /*DEG2RAD */
                (VOL_ELEM(prm.aph_file.spots_arg, k, i, j) + phase_shift +
                 k * ZZ(prm.Phase_Shift) +
                 i * YY(prm.Phase_Shift) +
                 j * XX(prm.Phase_Shift));
    }

    // Resize and symmetrize Fourier Transform
    int lmin = STARTINGZ(prm.aph_file.spots_abs);
    int lmax = FINISHINGZ(prm.aph_file.spots_abs);
    int kmin = STARTINGY(prm.aph_file.spots_abs);
    int kmax = FINISHINGY(prm.aph_file.spots_abs);
    int hmin = STARTINGX(prm.aph_file.spots_abs);
    int hmax = FINISHINGX(prm.aph_file.spots_abs);

    int lsize = MAX(ABS(lmin), ABS(lmax));
    int ksize = MAX(ABS(kmin), ABS(kmax));
    int hsize = MAX(ABS(hmin), ABS(hmax));

    Matrix3D< complex<double> > FT;
    FT.initZeros(2*lsize + 1, 2*ksize + 1, 2*hsize + 1);
    STARTINGZ(FT) = -lsize;
    STARTINGY(FT) = -ksize;
    STARTINGX(FT) = -hsize;

    switch (prm.aph_file.Space_Group)
    {
    case(1):
                    cout << "Detected P1 symmetry" << endl;
        symmetrize_P1(FT, prm);
        break;
    case(90)://P4212
                    cout << "Detected P4212 symmetry" << endl;
        symmetrize_P4212(FT, prm);
        break;
    default:
        cerr << "\nHORROR: Symmetry not implemented!!!" << endl;
        exit(false);
        break;

    }//switch end



    // Create volumen
    if (XX(prm.Celldim) == -1 ||
        YY(prm.Celldim) == -1 ||
        ZZ(prm.Celldim) == -1
       )
{
        ZZ(prm.Celldim) = ZSIZE(FT);
        YY(prm.Celldim) = XSIZE(FT);
        XX(prm.Celldim) = YSIZE(FT);
    }
    V1().initZeros(ZZ(prm.Celldim),
                    YY(prm.Celldim),
                    XX(prm.Celldim));
    V1().setXmippOrigin();


    IDFT_3D(FT, V1());
//   V1.write(prm.fn_out);
//   V1() *= prm.Scale_Factor;
//ADD HERE ANY OTHER FACTOR, MAY BE AS OPTION

    // Any repetition?
    int s_k = STARTINGZ(V1()); int e_k = FINISHINGZ(V1());
    int s_i = STARTINGY(V1()); int e_i = FINISHINGY(V1());
    int s_j = STARTINGX(V1()); int e_j = FINISHINGX(V1());
    V1().resize(ZZ(prm.Celldim)*ZZ(prm.NoCells),
                YY(prm.Celldim)*YY(prm.NoCells),
                XX(prm.Celldim)*XX(prm.NoCells));
    int s_K = STARTINGZ(V1()); int e_K = FINISHINGZ(V1());
    int s_I = STARTINGY(V1()); int e_I = FINISHINGY(V1());
    int s_J = STARTINGX(V1()); int e_J = FINISHINGX(V1());


    for (int K = s_K, k = s_k; K <= e_K; K++, k++)
    {
        if (k > e_k) k = s_k;
        for (int I = s_I, i = s_i; I <= e_I; I++, i++)
        {
            if (i > e_i) i = s_i;
            for (int J = s_J, j = s_j; J <= e_J; J++, j++)
            {
                if (j > e_j) j = s_j;
                VOLVOXEL(V1, K, I, J) = VOLVOXEL(V1, k, i, j);
            }//j
        }//i
    }//k
    // Save Image
    V1.write(prm.fn_out);

}

void symmetrize_P1(Matrix3D< complex<double> > &FT,
                   Spot2RealSpace3D_Parameters &prm)
{
    FOR_ALL_ELEMENTS_IN_MATRIX3D(prm.aph_file.spots_abs)
    {
        VOL_ELEM(FT, k, i, j) = polar(VOL_ELEM(prm.aph_file.spots_abs, k, i, j),
                                      DEG2RAD(VOL_ELEM(prm.aph_file.spots_arg, k, i, j)));
        VOL_ELEM(FT, -k, -i, -j) = conj(VOL_ELEM(FT, k, i, j));
    }
}
void symmetrize_P4212(Matrix3D< complex<double> > &FT,
                      Spot2RealSpace3D_Parameters &prm)
{
    int asymmh, asymmk, asymml;   /* Reflection equivalent in the asymm. unit. */
    int ip1, ip2;              /* To bring the phase from the asymmetric unit. */
    int spec;                       /* Indicates if the reflection is special. */
    int  iptest;/* Indicates the phase of reflection in case of being special. */
    double phase, amplitude;
    int horder = STARTINGX(FT);
    int korder  = STARTINGY(FT);
    int lorder  = STARTINGZ(FT);
    for (int l = lorder; l <= -lorder; l++)
    {
        for (int k = korder; k <= -korder; k++)
        {
            for (int h = horder; h <= -horder; h++)
            {
                asymmh = h;
                asymmk = k;
                asymml = l;
                /* Computing the reflection equivalent in the asymmetric unit. */
                AsymmUnitP4212(&asymmh, &asymmk, &asymml, &ip1, &ip2, &spec, &iptest);
//cout << "\n (" << l <<"," << k << "," << h <<")->";
//cout << " (" << asymml <<"," << asymmk << "," << asymmh <<")= ";
                if (VOL_ELEM(prm.aph_file.spots_abs, asymml, asymmk, asymmh) == 0)
                    continue;
                /* The amplitude is the same. */
                amplitude = VOL_ELEM(prm.aph_file.spots_abs, asymml, asymmk, asymmh);
                /* Bringing the phase from the asymmetric unit. */
                phase = VOL_ELEM(prm.aph_file.spots_arg, asymml, asymmk, asymmh) * ip1 - ip2;
                phase = fmod(phase, 360.0);
                if (phase > 180.0)  phase -= 360.0;
                if (phase < -180.0) phase += 360.0;
                /* Imposing phase from the symmetry. */
                // This line is for imposing the symmetry. I am going to keep it
                // commented since the simmetry should has been imposed before.
                // if(impose && spec) refptr->ap.phs = iptest;

                VOL_ELEM(FT, l, k, h) = polar(amplitude, DEG2RAD(phase));
//cout << "FT" <<VOL_ELEM(FT, l,k,h);
            }//for h
        }//for k
    }//for l

}
/*----------------------------------------------------------------------------*/
/*
    Brings the input reflection  H K L into the asymmetric unit according to
    the symmetry P4212. The routine returns the new values of the indexes
    H' K' L' within the asymmetric unit, as well as a flag indicating whether
    the reflection is special. And, in such a case, returns its phase value
    (0 or 90 degrees) corresponding to its character real or imaginary.

    The asymmetric unit in P4212 involves H,K,Z >=0 and H <= K. When H < 0,
    or H==0 and K < 0, or H==K==0 and Z < 0 the Conjugate Symmetry Property
    of the DFT is applied in order to bring the reflection into the H,K >= 0
    and make easier the subsequent processes.

                                                                              */
void AsymmUnitP4212(int *ih, int *ik, int *il, int *ip1, int *ip2,
                    int *spec, int *iptest)

{
    static int matrix[4][8] =
        {
            {
                -1, 0, 0, -1, -1,   0,   0, -1
            },
            { 1, 0, 0, -1,  1, 180, 180, -1},
            { 0, 1, 1,  0, -1,   0,   0,  1},
            { 0, 1, 1,  0,  1,   0,   0, -1}
        };

    static int gomatrix[7] =
        {
            1, 2, 1, 3, 1, 2, 1
        };

    int pass, index;

    /* Initialization of ip1 and ip2. */
    *ip1 = 1;
    *ip2 = 0;

    /* The asymmetric unit in P4212 involves H,K,Z >=0 and H <= K.
       In the cases in which H<0 or (H=0 & K<0) or (H=K=0 & Z<0) the program
       generates the corresponding reflection by applying the Conjugate Symmetry.
       This is got by multiplying by matrix[0]. */
    if (*ih < 0               ||              /* H < 0 */
        *ih == 0 && *ik < 0   ||              /* H == 0 && K < 0 */
        *ih == 0 && *ik == 0 && *il < 0)      /* H == 0 && K == 0 && Z < 0 */
        MatrixMult(matrix[0], ih, ik, il, ip1, ip2);

    /* Bringing the current reflection into the asymmetric unit. */
    for (pass = 0; pass < 2;)
    {
        index = 0;
        if (*ik >= 0) index ++;
        if (*il >= 0) index += 2;
        if (*ih < abs(*ik)) index += 4;

        /* <index> classifies the reflection by its indices.
           <gomatrix> indicates which matrix will bring the reflection into the
           unique asymmetric unit for a given index.

                       index      K>=0      Z>=0    |K|>=|H|
                       -----------------------------------
                         0         NO        NO        NO
                         1         YES       NO        NO
                         2         NO        YES       NO
                         3         YES       YES       NO
                         4         NO        NO        YES
                         5         YES       NO        YES
                         6         NO        YES       YES
                         7         YES       YES       YES

           P622 in the highest symmetry and its asymmetric unit is only. index = 7.

                                                                                  */

        if (index < 7)
        {
            MatrixMult(matrix[gomatrix[index]], ih, ik, il, ip1, ip2);

            if (gomatrix[index] < 3) continue;
        }

        if (*ih == 0 && *ik < 0)
            MatrixMult(matrix[0], ih, ik, il, ip1, ip2);

        pass++;
    }

    /* After reflections have been placed into the asymmetric unit they are
       examined to see if they are special reflections, ones whose phase must
       be either real (0 or PI) or imaginary (PI/2 or 3*PI/2). */
    CheckSpec(*ih, *ik, *il, spec, iptest);

}
/*----------------------------------------------------------------------------*/
void CheckSpec(int ih, int ik, int il, int *spec, int *iptest)

/*
   Checks if the current reflection (H K L) is a special reflection, i.e.,
   whose phase must be either real (0 or PI) or imaginary (PI/2 or 3*PI/2).

   It returns spec=1 if the reflection is special. Otherwise, 0.
   It returns in iptest 0 if real and 90 is imaginary.

   Conditions checked are:

      - H=0 special
      - K=0 special
      - Z=0 special
      - H=K special
      - If for H=0 or K=0 K+H is odd, it indicates an imaginary value for the
        reflection.

      Except the last condition, all other special reflections are real.

      Summary of special reflections.

                 Real               Imaginary
                ------             -----------
                (2n,0,L)            (2n+1,0,L)
                (0,2n,L)            (0,2n+1,L)
                (H,K,0)
                (H,H,L)
                                                                              */

{

    /* Initialization of values to return. */
    *spec = 0;
    *iptest = 0;

    /* Checking the conditions. */
    if (ih == 0 || ik == 0)
    {
        int i2, i = ih + ik;

        i2 = 2 * (i / 2);
        if (i > i2) *iptest = 90;                   /* Imaginary. Otherwise, Real. */
        (*spec)++;
        return;
    }
    if (il == 0 || ih == ik)(*spec)++;                                   /* Real. */

    return;
}
/*----------------------------------------------------------------------------*/
void MatrixMult(int A[], int *ih, int *ik, int *il, int *ip1, int *ip2)

/*
    Does matrix multiplication to bring reflections into the asymmetric unit.

        (H' K' Z' AMP' PHS') = (H K Z AMP PHS) <A>

    where <A> has form:

           A[0]  A[2]   0     0    A[5]
           A[1]  A[3]   0     0    A[6]
            0     0    A[4]   0     0
            0     0     0     1     0
            0     0     0     0    A[7]

     for all cases.

                                                                              */
{
    int ih1;

    /* New indexes H,K,L within asymmetric unit. */
    ih1 = *ih * A[0] + *ik * A[1];                          /* Temporal variable. */
    *ik = *ih * A[2] + *ik * A[3];
    *ih = ih1;
    *il *= A[4];

    /* Computation of the operations over the reflection phase. */
    *ip1 *= A[7];
    *ip2 += *ih * A[5] + *ik * A[6];
}
