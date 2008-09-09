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

void extractTrainingPatches(const FileName &fnPDB, int patchSize,
    int step, double Ts, double resolution,
    std::vector< Matrix1D<double> > &training)
{
    // Convert PDB to volume
    Prog_PDBPhantom_Parameters pdbconverter;
    pdbconverter.fn_pdb=fnPDB;
    pdbconverter.fn_out="";
    pdbconverter.Ts=Ts;
    pdbconverter.run();
    Matrix3D<double> &V0=pdbconverter.Vlow();
    
    // Filter the volume to the desired resolution
    FourierMask Filter;
    Filter.FilterShape=RAISED_COSINE;
    Filter.FilterBand=LOWPASS;
    Filter.w1=Ts/resolution;
    Filter.raised_w=0.05;
    V0.setXmippOrigin();
    Filter.generate_mask(V0);
    Filter.apply_mask_Space(V0);
    V0.threshold("below",0,0);

    // Build the pyramid and the differences between pyramid approximations
    V0.resize(4*FLOOR(ZSIZE(V0)/4.0),4*FLOOR(YSIZE(V0)/4.0),
        4*FLOOR(XSIZE(V0)/4.0));
    Matrix3D<double> V1;
    V0.pyramidReduce(V1);
    STARTINGX(V0)=STARTINGY(V0)=STARTINGZ(V0)=0;

    // Extract the training vectors from the volume
    const int L1=2;
    const int L0=2;
    const int N1=2*L1+1;
    const int N0=2*L0+1;
    Matrix1D<double> v(N1*N1*N1+N0*N0*N0);
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
                // Copy the pixels at level 0
                for (int kk=-L0; kk<=L0; kk++)
                    for (int ii=-L0; ii<=L0; ii++)
                        for (int jj=-L0; jj<=L0; jj++)
                            DIRECT_VEC_ELEM(v,idx++)=
                                DIRECT_VOL_ELEM(V0,k0+kk,i0+ii,j0+jj);

                // Copy the pixels at level 1
                for (int kk=-L1; kk<=L1; kk++)
                    for (int ii=-L1; ii<=L1; ii++)
                        for (int jj=-L1; jj<=L1; jj++)
                            DIRECT_VEC_ELEM(v,idx++)=
                                DIRECT_VOL_ELEM(V1,k1+kk,i1+ii,j1+jj);
                
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
            training.push_back(auxTraining[i]);
}

int main(int argc, char *argv[])
{
    FileName fnSel;    // Selfile with the PDBs
    FileName fnOut;    // Name of the dictionary
    int S;             // Number of allowed atoms for a patch
    int dictSize;      // Number of atoms in the dictionary
    int patchSize;     // Size of the window centered at a voxel
    int step;          // Step between patches
    double Ts;         // Sampling rate
    double resolution; // Final resolution
    // Read Parameters
    try
    {
        fnSel = getParameter(argc,argv,"-sel");
        fnOut = getParameter(argc,argv,"-o");
        S = textToInteger(getParameter(argc,argv,"-S","3"));
        dictSize = textToInteger(getParameter(argc,argv,"-dictSize","400"));
        patchSize = textToInteger(getParameter(argc,argv,"-patchSize","5"));
        step = textToInteger(getParameter(argc,argv,"-patchSize","2"));
        Ts = textToFloat(getParameter(argc,argv,"-sampling","2"));
        resolution = textToFloat(getParameter(argc,argv,"-resolution","10"));
    }
    catch (Xmipp_error &XE)
    {
        std::cout << XE;
        std::cout << "Usage: xmipp_pdb_dictionary\n"
                  << "    -sel <selfile>     : Set of PDB files\n"
                  << "    -o <dictionary>    : Name of the dictionary\n"
                  << "   [-S <S=3>]          : Number of allowed atoms\n"
                  << "   [-dictSize <N=400>] : Size of the dictionary\n"
                  << "   [-patchSize <W=5>]  : Size of the patch around a voxel\n"
                  << "   [-step <s=2>]       : Distance between learning patches\n"
                  << "   [-sampling <Ts=2>]  : Sampling rate in Angstroms/pixel\n"
                  << "   [-resolution <A=10>]: Final resolution in Angstroms\n"
        ;
        exit(1);
    }

    // Call main routine
    try
    {
        // Show parameters
        std::cout << "Selfile:   " << fnSel      << std::endl
                  << "Output:    " << fnOut      << std::endl
                  << "S:         " << S          << std::endl
                  << "DictSize:  " << dictSize   << std::endl
                  << "patchSize: " << patchSize  << std::endl
                  << "step:      " << step       << std::endl
                  << "sampling:  " << Ts         << std::endl
                  << "resolution:" << resolution << std::endl
        ;

        // Define variables
        SelFile SF;
        SF.read(fnSel);
        std::vector< Matrix1D<double> > training;
        
        // Extract training patches
        while (!SF.eof())
        {
            FileName fnPDB=SF.NextImg();
            extractTrainingPatches(fnPDB, patchSize, step, Ts, resolution,
                training);
        }
        
        // Initialize the dictionary at random
        int N=XSIZE(training[0]);
        Matrix2D<double> D;
        D.initZeros(N,dictSize);
        Matrix1D<int> used;
        used.initZeros(training.size());
        
        // Fill the first column to DC
        D(0,0)=1/sqrt(N);
        for (int j=1; j<N; j++) D(j,0)=D(0,0);
        
        // Fill the rest of columns chosing one of the data samples randomly
        for (int d=1; d<dictSize; d++)
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
        
        // Now optimize the dictionary
        std::vector< Matrix1D<double> > Alpha;
        std::cout << "Learning dictionary ...\n";
        kSVD(training, S, D, Alpha);
        
        // Write the dictionary
        std::ofstream fhOut;
        fhOut.open(fnOut.c_str());
        if (!fhOut)
            REPORT_ERROR(1,(std::string)"Cannot open "+fnOut+" for output");
        fhOut << S << " " << dictSize << " " << N << std::endl << D;
        fhOut.close();
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        exit(1);
    }
    exit(0);
}
