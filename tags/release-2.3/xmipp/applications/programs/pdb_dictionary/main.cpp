/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <data/args.h>
#include <data/selfile.h>
#include <reconstruction/convert_pdb2vol.h>
#include <reconstruction/fourier_filter.h>
#include <classification/kSVD.h>
#include <fstream>

//#define DEBUG
void extractTrainingPatches(const FileName &fnPDB, int patchSize,
    int step, double Ts, double resolution1, double resolution2,
    std::vector< Matrix1D<double> > &training)
{
    // Convert PDB to volume
    Prog_PDBPhantom_Parameters pdbconverter;
    pdbconverter.fn_pdb=fnPDB;
    pdbconverter.fn_out="";
    pdbconverter.Ts=Ts;
    pdbconverter.run();
    
    // Filter the volume to the final resolution
    Matrix3D<double> V0=pdbconverter.Vlow();
    FourierMask Filter1;
    Filter1.FilterShape=RAISED_COSINE;
    Filter1.FilterBand=LOWPASS;
    Filter1.w1=Ts/resolution1;
    Filter1.raised_w=0.05;
    V0.setXmippOrigin();
    Filter1.apply_mask_Space(V0);
    V0.threshold("below",0,0);

    // Filter the volume to the restoration resolution
    Matrix3D<double> V0R=pdbconverter.Vlow();
    FourierMask Filter2;
    Filter2.FilterShape=RAISED_COSINE;
    Filter2.FilterBand=LOWPASS;
    Filter2.w1=Ts/resolution2;
    Filter2.raised_w=0.05;
    V0R.setXmippOrigin();
    Filter2.apply_mask_Space(V0R);
    V0R.threshold("below",0,0);

    // Build the pyramid and the differences between pyramid approximations
    V0.resize(2*FLOOR(ZSIZE(V0)/2.0),2*FLOOR(YSIZE(V0)/2.0),
        2*FLOOR(XSIZE(V0)/2.0));
    V0R.resize(V0);
    Matrix3D<double> V1;
    V0.pyramidReduce(V1);
    STARTINGX(V0)=STARTINGY(V0)=STARTINGZ(V0)=0;
    STARTINGX(V0R)=STARTINGY(V0R)=STARTINGZ(V0R)=0;

    // Extract the training vectors from the volume
    int N1=patchSize;
    int N0=patchSize;
    int N0_3=N0*N0*N0;
    int N1_3=N1*N1*N1;
    int L1=(patchSize-1)/2;
    int L0=(patchSize-1)/2;
    Matrix1D<double> v(N1*N1*N1+2*N0*N0*N0);
    std::vector< Matrix1D<double> > auxTraining;
    std::vector<double> energy;
    double bestEnergy=0;
    for (int k0=4; k0<ZSIZE(V0)-4; k0+=step)
        for (int i0=4; i0<YSIZE(V0)-4; i0+=step)
            for (int j0=4; j0<XSIZE(V0)-4; j0+=step)
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
                            DIRECT_VEC_ELEM(v,idx)=
                                DIRECT_VOL_ELEM(V0,k0+kk,i0+ii,j0+jj);
                            avg0+=DIRECT_VEC_ELEM(v,idx);
                            idx++;
                        }
                avg0/=N0_3;

                // Copy the pixels at level 1
                double avg1=0;
                for (int kk=-L1; kk<=L1; kk++)
                    for (int ii=-L1; ii<=L1; ii++)
                        for (int jj=-L1; jj<=L1; jj++)
                        {
                            DIRECT_VEC_ELEM(v,idx)=
                                DIRECT_VOL_ELEM(V1,k1+kk,i1+ii,j1+jj);
                            avg1+=DIRECT_VEC_ELEM(v,idx);
                            idx++;
                        }
                avg1/=N1_3;
                
                // Copy the pixels at level 0R
                double avg0R=0;
                for (int kk=-L0; kk<=L0; kk++)
                    for (int ii=-L0; ii<=L0; ii++)
                        for (int jj=-L0; jj<=L0; jj++)
                        {
                            DIRECT_VEC_ELEM(v,idx)=
                                DIRECT_VOL_ELEM(V0R,k0+kk,i0+ii,j0+jj);
                            avg0R+=DIRECT_VEC_ELEM(v,idx);
                            idx++;
                        }
                avg0R/=N0_3;

                // Substract the mean
                for (idx=0; idx<N0; idx++)
                    DIRECT_VEC_ELEM(v,idx)-=avg0;
                for (idx=N0; idx<N0+N1; idx++)
                    DIRECT_VEC_ELEM(v,idx)-=avg1;
                for (idx=N0+N1; idx<XSIZE(v); idx++)
                    DIRECT_VEC_ELEM(v,idx)-=avg0R;

                auxTraining.push_back(v);
                energy.push_back(v.module());
                if (energy[energy.size()-1]>bestEnergy)
                    bestEnergy=energy[energy.size()-1];
            }

    // Remove vectors with very little energy
    double energyThreshold=0.1*bestEnergy; // In fact it is 0.1^2 of the energy
                                           // because bestEnergy is a square root
    int imax=energy.size();
    for (int i=0; i<imax; i++)
        if (energy[i]>energyThreshold)
        {
            training.push_back(auxTraining[i]);
            #ifdef DEBUG
                VolumeXmipp save(5*3,5,5);
                int k0=0, idx=0;
                for (int kk=0; kk<N0; kk++)
                    for (int ii=0; ii<N0; ii++)
                        for (int jj=0; jj<N0; jj++)
                            save(k0+kk,ii,jj)=
                                DIRECT_VEC_ELEM(auxTraining[i],idx++);
                k0=5;
                for (int kk=0; kk<N0; kk++)
                    for (int ii=0; ii<N0; ii++)
                        for (int jj=0; jj<N0; jj++)
                            save(k0+kk,ii,jj)=
                                DIRECT_VEC_ELEM(auxTraining[i],idx++);
                k0=10;
                for (int kk=0; kk<N0; kk++)
                    for (int ii=0; ii<N0; ii++)
                        for (int jj=0; jj<N0; jj++)
                            save(k0+kk,ii,jj)=
                                DIRECT_VEC_ELEM(auxTraining[i],idx++);
                std::cout << save();
                save.write("PPPtrainingVector.vol");
                save()=V0; save.write("PPPLevel0.vol");
                save()=V0R; save.write("PPPLevel0_Restoration.vol");
                save()=V1; save.write("PPPLevel1.vol");
                std::cout << "Press any key\n";
                char c; std::cin >> c;
            #endif
        }
}
#undef DEBUG

int main(int argc, char *argv[])
{
    FileName fnSel;      // Selfile with the PDBs
    FileName fnOut;      // Name of the dictionary
    int S;               // Number of allowed atoms for a patch (OMP)
    double lambda;       // Regularization parameter (LASSO)
    int dictSize;        // Number of atoms in the dictionary
    int patchSize;       // Size of the window centered at a voxel
    int step;            // Step between patches
    double Ts;           // Sampling rate
    double resolution1;  // Final resolution
    double resolution2;  // Restoration resolution
    bool initRandom;     // Initalize randomly the dictionary
    int projectionMethod;// Projection method
    // Read Parameters
    try
    {
        fnSel = getParameter(argc,argv,"-sel");
        fnOut = getParameter(argc,argv,"-o");
        S = textToInteger(getParameter(argc,argv,"-S","3"));
        lambda = textToFloat(getParameter(argc,argv,"-l","0.2"));
        dictSize = textToInteger(getParameter(argc,argv,"-dictSize","400"));
        patchSize = textToInteger(getParameter(argc,argv,"-patchSize","5"));
        step = textToInteger(getParameter(argc,argv,"-step","2"));
        Ts = textToFloat(getParameter(argc,argv,"-sampling","2"));
        resolution1 = textToFloat(getParameter(argc,argv,"-resolution1","10"));
        resolution2 = textToFloat(getParameter(argc,argv,"-resolution2"," 7"));
        initRandom = checkParameter(argc,argv,"-initRandom");
        projectionMethod = textToInteger(getParameter(argc,argv,
            "-projectionMethod","1"));
        if (projectionMethod==LASSO_PROJECTION) S=0;
        else lambda=0;
    }
    catch (Xmipp_error &XE)
    {
        std::cout << XE;
        std::cout << "Usage: xmipp_pdb_dictionary\n"
                  << "    -sel <selfile>          : Set of PDB files\n"
                  << "    -o <dictionary>         : Name of the dictionary\n"
                  << "   [-S <S=3>]               : Number of allowed atoms (OMP)\n"
                  << "   [-l <l=0.2>]             : Sparsity regularization weight (LASSO)\n"
                  << "   [-dictSize <N=400>]      : Size of the dictionary\n"
                  << "   [-patchSize <W=5>]       : Size of the patch around a voxel\n"
                  << "   [-step <s=2>]            : Distance between learning patches\n"
                  << "   [-sampling <Ts=2>]       : Sampling rate in Angstroms/pixel\n"
                  << "   [-resolution1 <A=10>]    : Final resolution in Angstroms\n"
                  << "   [-resolution2 <A=7>]     : Restoration resolution in Angstroms\n"
                  << "   [-initRandom]            : Initialize the dictionary randomly\n"
                  << "   [-projectionMethod <m=1>]: Projection method\n"
                  << "                              1=OMP\n"
                  << "                              2=LASSO\n"
        ;
        exit(1);
    }

    // Call main routine
    try
    {
        // Show parameters
        std::cout << "Selfile:         " << fnSel            << std::endl
                  << "Output:          " << fnOut            << std::endl
                  << "S:               " << S                << std::endl
                  << "lambda:          " << lambda           << std::endl
                  << "DictSize:        " << dictSize         << std::endl
                  << "patchSize:       " << patchSize        << std::endl
                  << "step:            " << step             << std::endl
                  << "sampling:        " << Ts               << std::endl
                  << "resolution1:     " << resolution1      << std::endl
                  << "resolution2:     " << resolution2      << std::endl
                  << "initRandom:      " << initRandom       << std::endl
                  << "projectionMethod:" << projectionMethod << std::endl
        ;

        // Define variables
        SelFile SF;
        SF.read(fnSel);
        std::vector< Matrix1D<double> > training;
        
        // Extract training patches
        while (!SF.eof())
        {
            FileName fnPDB=SF.NextImg();
            extractTrainingPatches(fnPDB, patchSize, step, Ts, resolution1,
                resolution2, training);
        }
        
        // Initialize the dictionary
        int N=XSIZE(training[0]);
        Matrix2D<double> D;
        D.initZeros(N,dictSize);
        
        if (!initRandom)
        {
            // Fill the columns chosing one of the data samples randomly
            Matrix1D<int> used;
            used.initZeros(training.size());
            for (int d=0; d<dictSize; d++)
            {
                int selected;
                do
                {
                    selected=ROUND(rnd_unif(0,XSIZE(used)-1));
                } while (used(selected));
                used(selected)=1;
                double inorm=1.0/training[selected].module();
                for (int j=0; j<N; j++)
                    D(j,d)=training[selected](j)*inorm;
            }
            used.clear();
        }
        else
        {
            D.initRandom(0,1,"gaussian");

            // Normalize the dictionary
            for (int j=0; j<XSIZE(D); j++)
            {
                // Make sure that the atoms are 0 mean
                double avg=0;
                for (int i=0; i<YSIZE(D); i++)
                    avg+=D(i,j);
                avg/=YSIZE(D);
                for (int i=0; i<YSIZE(D); i++)
                    D(i,j)-=avg;

                // Make them unitary
                double norm=0;
                for (int i=0; i<YSIZE(D); i++)
                    norm+=D(i,j)*D(i,j);
                norm=1/sqrt(norm);
                for (int i=0; i<YSIZE(D); i++)
                    D(i,j)*=norm;
            }
        }

        // Now optimize the dictionary
        std::vector< Matrix1D<double> > Alpha;
        std::cout << "Learning dictionary ...\n";
        kSVD(training, S, D, Alpha, false, 20, 0.01, projectionMethod, lambda);

        // Write the dictionary
        std::ofstream fhOut;
        fhOut.open(fnOut.c_str());
        if (!fhOut)
            REPORT_ERROR(1,(std::string)"Cannot open "+fnOut+" for output");
        fhOut << S << " " << " " << lambda << " " << dictSize << " " << N
              << " " << patchSize << std::endl << D << std::endl;
        fhOut.close();
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        exit(1);
    }
    exit(0);
}
