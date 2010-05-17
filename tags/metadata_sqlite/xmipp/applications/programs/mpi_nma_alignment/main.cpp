/***************************************************************************
 *
 * Authors:  Slavica Jonic slavica.jonic@impmc.jussieu.fr  
 *           Carlos Oscar Sanchez Sorzano coss.eps@ceu.es
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <reconstruction/nma_alignment.h>

#include <mpi.h>

int main(int argc, char **argv)
{
    // Initialize MPI
    int rank, NProcessors;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &NProcessors);

    Prog_nma_alignment_prm prm;
    try
    {
        // Read input parameters
	prm.MPIversion=true;
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
        // Prepare side information
        prm.produce_side_info(rank);
        MetaData SF_in(prm.fn_in);

        // Divide the selfile in chunks
        int imgNbr = prm.get_images_to_process();
        int Nchunk = (int)((float)imgNbr / (float)(NProcessors - 1));
        int myFirst = (rank - 1) * Nchunk;
        int myLast = rank * Nchunk - 1;
        if (rank == NProcessors - 1) myLast = imgNbr - 1;

        // Make the alignment, rank=0 receives all the assignments
        // The rest of the ranks compute the angular parameters for their
        // assigned images
	
	int numberofparam = 7 + prm.modeList.size();
	
        double v[numberofparam];
        if (rank == 0)
        {
            int i=0;
FOR_ALL_OBJECTS_IN_METADATA(SF_in){
   FileName fnImg;
   SF_in.getValue(MDL_IMAGE,fnImg);
   prm.img_names.push_back(fnImg);
Matrix1D<double> dummy;
   prm.listAssignments.push_back(dummy);
   i++;
}
            int toGo = imgNbr;
            MPI_Status status;
            //std::cerr << "Assigning modes and angles ...\n";
            init_progress_bar(imgNbr);
            while (toGo > 0)
            {
                MPI_Recv(v, numberofparam, MPI_DOUBLE, MPI_ANY_SOURCE,
                         MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                int i = (int)v[0];
		Matrix1D<double> aux(numberofparam-1);

	        for (int j=1; j<numberofparam; j++)
                    aux(j-1)=v[j];
		prm.listAssignments[i]=aux;
                toGo--;
            }
            progress_bar(imgNbr);
        }
        else
        {
            for (int i = myFirst; i <= myLast; i++)
            {
                // Read image and estimate angular parameters
                ImageXmipp I;
                FileName tempname;
                SF_in.getValue(MDL_IMAGE,tempname,i);
                
                I.read(tempname, false, false, false);
                I().setXmippOrigin();
                prm.assignParameters(I);

                // Send the alignment parameters to the master
                v[0] = i;
		for (int j=1; j<numberofparam; j++)
                {
                    v[j]=prm.parameters(j-1);
                }
                MPI_Send(v, numberofparam, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
        }

        if (rank == 0) prm.finish_processing();
        MPI_Finalize();
        return 0 ;
    }

    catch (Xmipp_error XE)
    {
        std::cout << XE << std::endl;
        MPI_Finalize();
	return 1 ;
    }
}
