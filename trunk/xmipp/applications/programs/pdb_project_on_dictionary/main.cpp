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

#include <data/args.h>
#include <data/matrix2d.h>
#include <data/volume.h>
#include <classification/kSVD.h>
#include <fstream>

double projectOntoDictionary(const Matrix2D<double> &D,
                             double lambda, int S, double lambdaDictionary,
                             bool restore, int patchSize,
                             Matrix3D<double> &Vin, Matrix3D<double> &Vout)
{
    // Get the size of the different vectors
    int N1=patchSize;
    int N0=patchSize;
    int N0_3=N0*N0*N0;
    int N1_3=N1*N1*N1;
    int L1=(patchSize-1)/2;
    int L0=(patchSize-1)/2;

    // Some initialization
    Matrix3D<double> &V0=Vin; // An alias
    Vout.initZeros(Vin);
    Matrix3D<double> Weight;
    Weight=Vout;

    // Separate the two dictionaries
    Matrix2D<double> D1, D2;
    int dim1=2*patchSize*patchSize*patchSize;
    int dim2=  patchSize*patchSize*patchSize;
    D2=D1=D;
    D1.window(0,0,dim1-1,XSIZE(D)-1);
    D2.window(dim1,0,YSIZE(D)-1,XSIZE(D)-1);

    // Readjust the columns of each dictionary so that they are
    // again unitary
    Matrix1D<double> normD1;
    normD1.initZeros(XSIZE(D1));
    Matrix1D<double> normD10;
    normD10.initZeros(XSIZE(D1));
    for (int j=0; j<XSIZE(D1); j++)
    {
        // Compute the norm of the high resolution part of the vector
        for (int i=0; i<N0; i++)
            DIRECT_VEC_ELEM(normD10,j)+=
                DIRECT_MAT_ELEM(D1,i,j)*DIRECT_MAT_ELEM(D1,i,j);
        DIRECT_VEC_ELEM(normD1,j)=DIRECT_VEC_ELEM(normD10,j);

        // Compute the norm of the rest of the vector
        for (int i=N0; i<YSIZE(D1); i++)
            DIRECT_VEC_ELEM(normD1,j)+=
                DIRECT_MAT_ELEM(D1,i,j)*DIRECT_MAT_ELEM(D1,i,j);
        DIRECT_VEC_ELEM(normD10,j)=sqrt(DIRECT_VEC_ELEM(normD10,j));
        DIRECT_VEC_ELEM(normD1,j)=sqrt(DIRECT_VEC_ELEM(normD1,j));
        double inormD1j=1/DIRECT_VEC_ELEM(normD1,j);

        for (int i=0; i<YSIZE(D1); i++)
            DIRECT_MAT_ELEM(D1,i,j)*=inormD1j;

        // Multiply D2 by a number so that it can be used in the recovery
        std::cout << normD10(j) << " " << normD1(j) << std::endl;
        double normD2j=DIRECT_VEC_ELEM(normD1,j)/DIRECT_VEC_ELEM(normD10,j);
        for (int i=0; i<YSIZE(D2); i++)
            DIRECT_MAT_ELEM(D2,i,j)*=normD2j;
    }
    std::cout << "Fin\n";

    // Check if LASSO
    Matrix2D<double> DtD, DtDlambda, DtDlambdaInv;
    if (S==0)
    {
        int K=XSIZE(D1);
        DtD.initZeros(K,K);
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(DtD)
        // Compute the dot product between the columns i and j of D
        for (int k=0; k<YSIZE(D1); k++)
            DIRECT_MAT_ELEM(DtD,i,j)+=
                DIRECT_MAT_ELEM(D1,k,i)*DIRECT_MAT_ELEM(D1,k,j);
        DtDlambda=DtD;
        for (int i=0; i<K; i++)
            DIRECT_MAT_ELEM(DtDlambda,i,i)+=lambda;
        DtDlambda.inv(DtDlambdaInv);
    }

    // Build the pyramid
    Matrix3D<double> V1;
    V0.pyramidReduce(V1);
    STARTINGX(V0)=STARTINGY(V0)=STARTINGZ(V0)=0;

    // Extract the training vectors from the volume
    double error=0, Nerror=0;
    Matrix1D<double> v1(N1*N1*N1+N0*N0*N0);
    Matrix1D<double> alpha(XSIZE(D1)), vp1(N1*N1*N1+N0*N0*N0), vp2(N0*N0*N0);
    init_progress_bar(ZSIZE(V0)-8);
    for (int k0=4; k0<ZSIZE(V0)-4; k0++)
    {
        for (int i0=4; i0<YSIZE(V0)-4; i0++)
            for (int j0=4; j0<XSIZE(V0)-4; j0++)
            {
                // Locate this voxel in the downsampled image
                int k1=ROUND(k0/2.0);
                int i1=ROUND(i0/2.0);
                int j1=ROUND(j0/2.0);

                int idx=0;
                double avg0=0;
                // Copy the pixels at level 0
                for (int kk=-L0; kk<=L0; kk++)
                    for (int ii=-L0; ii<=L0; ii++)
                        for (int jj=-L0; jj<=L0; jj++)
                        {
                            DIRECT_VEC_ELEM(v1,idx)=
                                DIRECT_VOL_ELEM(V0,k0+kk,i0+ii,j0+jj);
                            avg0+=DIRECT_VEC_ELEM(v1,idx);
                            idx++;
                        }
                avg0/=N0_3;

                // Copy the pixels at level 1
                double avg1=0;
                for (int kk=-L1; kk<=L1; kk++)
                    for (int ii=-L1; ii<=L1; ii++)
                        for (int jj=-L1; jj<=L1; jj++)
                        {
                            DIRECT_VEC_ELEM(v1,idx)=
                                DIRECT_VOL_ELEM(V1,k1+kk,i1+ii,j1+jj);
                            avg1+=DIRECT_VEC_ELEM(v1,idx);
                            idx++;
                        }
                avg1/=N0_3;

                // Substract the mean
                for (idx=0; idx<N0; idx++)
                    DIRECT_VEC_ELEM(v1,idx)-=avg0;
                for (idx=N0; idx<XSIZE(v1); idx++)
                    DIRECT_VEC_ELEM(v1,idx)-=avg1;

                // Project this vector onto the dictionary
                if (S!=0) orthogonalMatchingPursuit(v1,D1,S,alpha);
                else lasso(v1,D1,DtD,DtDlambdaInv,lambdaDictionary,alpha);
                vp1.initZeros(); // vp1=D1*alpha
                vp2.initZeros(); // vp2=D2*alpha
                bool proceed=false;
                for (int j=0; j<XSIZE(D1); j++)
                    if (DIRECT_VEC_ELEM(alpha,j)!=0)
                    {
                        proceed=true;
                        for (int i=0; i<YSIZE(D1); i++)
                            DIRECT_VEC_ELEM(vp1,i)+=
                                DIRECT_MAT_ELEM(D1,i,j)*
                                DIRECT_VEC_ELEM(alpha,j);
                        if (restore)
                            for (int i=0; i<YSIZE(D2); i++)
                                DIRECT_VEC_ELEM(vp2,i)+=
                                    DIRECT_MAT_ELEM(D2,i,j)*
                                    DIRECT_VEC_ELEM(alpha,j);
                    }

                if (proceed)
                {
                    if (vp1.module()!=0)
                        std::cout << vp1.module() << " " << vp2.module() << std::endl;

                    // Add the mean
                    for (idx=0; idx<N0; idx++)
                    {
                        DIRECT_VEC_ELEM(vp1,idx)+=avg0;
                        if (restore) DIRECT_VEC_ELEM(vp2,idx)+=avg0;
                    }
                    for (idx=N0; idx<XSIZE(vp1); idx++)
                        DIRECT_VEC_ELEM(vp1,idx)+=avg1;

                    // Measure the projection error
                    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(vp1)
                    error+=ABS(DIRECT_VEC_ELEM(v1,i)-DIRECT_VEC_ELEM(vp1,i));
                    Nerror+=XSIZE(vp1);

                    // Put it back in the volume and update weights
                    idx=0;

                    // Copy the pixels at level 0
                    for (int kk=-L0; kk<=L0; kk++)
                        for (int ii=-L0; ii<=L0; ii++)
                            for (int jj=-L0; jj<=L0; jj++)
                            {
                                if (restore)
                                    DIRECT_VOL_ELEM(Vout,k0+kk,i0+ii,j0+jj)+=
                                        DIRECT_VEC_ELEM(vp2,idx++);
                                else
                                    DIRECT_VOL_ELEM(Vout,k0+kk,i0+ii,j0+jj)+=
                                        DIRECT_VEC_ELEM(vp1,idx++);
                                DIRECT_VOL_ELEM(Weight,k0+kk,i0+ii,j0+jj)++;
                            }
                }
            }
        progress_bar(k0-4);
    }
    progress_bar(ZSIZE(V0)-8);
    error/=Nerror;

    // Normalize the output
    FOR_ALL_ELEMENTS_IN_MATRIX3D(Vout)
    {
        double denominator=lambda+Weight(k,i,j);
        if (denominator!=0)
            Vout(k,i,j)=(lambda*Vin(k,i,j)+Vout(k,i,j))/denominator;
    }
    return error;
}

int main(int argc, char *argv[])
{
    FileName fnIn;     // Volume to project
    FileName fnDict;   // Name of the dictionary
    FileName fnOut;    // Name of the output volume
    int patchSize;     // Patchsize
    double lambda;
    bool restore;      // Perform a restoration

    // Read Parameters
    try
    {
        fnIn      = getParameter(argc,argv,"-i");
        fnOut     = getParameter(argc,argv,"-o","");
        fnDict    = getParameter(argc,argv,"-dict");
        lambda    = textToFloat(getParameter(argc,argv,"-l","0"));
        restore   = checkParameter(argc,argv,"-restore");
    }
    catch (Xmipp_error &XE)
    {
        std::cout << XE;
        std::cout << "Usage: xmipp_pdb_dictionary\n"
                  << "    -i <volume>        : Volume to project\n"
                  << "    -dict <dictionary> : Name of the dictionary\n"
                  << "   [-o <volume>]       : Output volume\n"
                  << "   [-l <lambda=0>]     : Regularization parameter\n"
                  << "   [-restore]          : Perform restoration\n"
                  ;
        exit(1);
    }

    // Call main routine
    try
    {
        // Show parameters
        std::cout << "Input:     " << fnIn    << std::endl
                  << "Output:    " << fnOut   << std::endl
                  << "Dictionary:" << fnDict  << std::endl
                  << "Lambda:    " << lambda  << std::endl
                  << "Restore:   " << restore << std::endl
                  ;

        // Read the dictionary
        std::ifstream fhDict;
        fhDict.open(fnDict.c_str());
        if (!fhDict)
            REPORT_ERROR(1,(std::string)"Cannot open "+fnDict+" for output");
        int S, dictSize, N;
        double lambdaDictionary;
        fhDict >> S >> lambdaDictionary >> dictSize >> N >> patchSize;
        Matrix2D<double> D;
        D.resize(N,dictSize);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(D)
        fhDict >> D(i,j);
        fhDict.close();

        // Read the volume
        VolumeXmipp Vin;
        Vin.read(fnIn);

        // Project the volume
        VolumeXmipp Vout;
        double error=projectOntoDictionary(D,lambda,S,lambdaDictionary,
                                           restore, patchSize, Vin(),Vout());
        std::cout << "Projection error= " << error << std::endl;

        // And write results
        if (fnOut!="") Vout.write(fnOut);
        else Vout.write(fnIn);
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        exit(1);
    }
    exit(0);
}
