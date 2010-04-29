/***************************************************************************
 *
 * Authors:     Carlos Oscar Sanchez Sorzano (coss.eps@ceu.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 * Lab. de Bioingenieria, Univ. San Pablo CEU
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

#include <mpi.h>
#include <reconstruction/common_lines.h>

int main(int argc, char **argv)
{
    // Initialize MPI
    int rank, NProcessors;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &NProcessors);

    CommonLine_Parameters prm;
    try
    {
        prm.read(argc, argv);
    }
    catch (Xmipp_error XE)
    {
        if (rank == 0)
        {
            std::cout << XE;
            prm.usage();
        }
        MPI_Finalize();
        exit(1);
    }

    try
    {
        // Everyone produces its own side info
        prm.Nmpi=NProcessors;
        prm.produceSideInfo();
        if (rank==0)
            std::cout << prm;

        // Everyone runs on a little piece of the dataset
        if (rank==0)
            std::cerr << "Computing common lines ...\n";
        prm.run(rank);
        MPI_Barrier(MPI_COMM_WORLD);

        // Finally send all the results to node 0
        std::vector<int> i_j;
        std::vector<double> angi_angj_d;
        if (rank == 0)
        {
            std::cerr << "Gathering information ...\n";
            int toGo = NProcessors-1;
            while (toGo > 0)
            {
                MPI_Status status;
                int N;
                MPI_Recv(&N, 1, MPI_INT, MPI_ANY_SOURCE,
                         MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                i_j.resize(2*N);
                angi_angj_d.resize(3*N);
                MPI_Recv(&i_j[0], 2*N, MPI_INT, status.MPI_SOURCE,
                    MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
                MPI_Recv(&angi_angj_d[0], 3*N, MPI_DOUBLE, status.MPI_SOURCE,
                    MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
                for (int i=0; i<N; i++)
                {
                    long int ii=i_j[2*i];
                    long int jj=i_j[2*i+1];
                    long int idx=ii*prm.Nimg+jj;
                    prm.CLmatrix[idx].angi=angi_angj_d[3*i];
                    prm.CLmatrix[idx].angj=angi_angj_d[3*i+1];
                    prm.CLmatrix[idx].distanceij=angi_angj_d[3*i+2];
                }
                toGo--;
            }
        }
        else
        {
            for (int i=0; i<prm.Nimg; i++)
                for (int j=0; j<prm.Nimg; j++)
                {
                    long int idx=i*prm.Nimg+j;
                    double d=prm.CLmatrix[idx].distanceij;
                    if (d>0)
                    {
                        i_j.push_back(i);
                        i_j.push_back(j);
                        angi_angj_d.push_back(prm.CLmatrix[idx].angi);
                        angi_angj_d.push_back(prm.CLmatrix[idx].angj);
                        angi_angj_d.push_back(d);
                    }
                }
            int N=i_j.size()/2;
            MPI_Send(&N,                1, MPI_INT,    0, 0, MPI_COMM_WORLD); 
            MPI_Send(&i_j[0],         2*N, MPI_INT,    0, 0, MPI_COMM_WORLD);
            MPI_Send(&angi_angj_d[0], 3*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        if (rank == 0)
        {
            prm.qualifyCommonLines();
            prm.writeResults();
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE << std::endl;
        MPI_Finalize();
	return 1 ;
    }
    MPI_Finalize();
    return 0 ;
}
