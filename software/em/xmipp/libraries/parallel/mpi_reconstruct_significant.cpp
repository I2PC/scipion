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

#include "mpi_reconstruct_significant.h"

MpiProgReconstructSignificant::MpiProgReconstructSignificant()
{
	node=NULL;
}

MpiProgReconstructSignificant::~MpiProgReconstructSignificant()
{
	delete node;
}

void MpiProgReconstructSignificant::read(int argc, char** argv)
{
    node = new MpiNode(argc, argv);
   	rank = node->rank;
   	Nprocessors = node->size;

   	ProgReconstructSignificant::read(argc, (const char **)argv);
}

void MpiProgReconstructSignificant::synchronize()
{
	node->barrierWait();
}

void MpiProgReconstructSignificant::gatherAlignment()
{
	// Share weights and cc volumes
	MultidimArray<double> aux;
	if (rank==0)
		aux.resizeNoCopy(cc);
	xmipp_MPI_Reduce(MULTIDIM_ARRAY(weight), MULTIDIM_ARRAY(aux), MULTIDIM_SIZE(weight), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)
		weight=aux;
	xmipp_MPI_Reduce(MULTIDIM_ARRAY(cc), MULTIDIM_ARRAY(aux), MULTIDIM_SIZE(cc), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank==0)
		cc=aux;

	// Write all metadatas
	if (rank!=0)
		for (size_t n=0; n<mdReconstructionPartial.size(); ++n)
		{
			FileName fnPartial=formatString("%s/partial_node%03d_%03d.xmd",fnDir.c_str(),(int)rank,(int)n);
			FileName fnProjectionMatching=formatString("%s/projmatch_node%03d_%03d.xmd",fnDir.c_str(),(int)rank,(int)n);

			if (mdReconstructionPartial[n].size()>0)
				mdReconstructionPartial[n].write(fnPartial);
			if (mdReconstructionProjectionMatching[n].size()>0)
				mdReconstructionProjectionMatching[n].write(fnProjectionMatching);
		}
	synchronize();

	// Now the master takes all of them
	if (rank==0)
	{
		MetaData MDAux;
		for (size_t otherRank=1; otherRank<Nprocessors; ++otherRank)
		{
			for (size_t n=0; n<mdReconstructionPartial.size(); ++n)
			{
				FileName fnPartial=formatString("%s/partial_node%03d_%03d.xmd",fnDir.c_str(),(int)otherRank,(int)n);
				FileName fnProjectionMatching=formatString("%s/projmatch_node%03d_%03d.xmd",fnDir.c_str(),(int)otherRank,(int)n);

				if (fnPartial.exists())
				{
					MDAux.read(fnPartial);
					mdReconstructionPartial[n].unionAll(MDAux);
					deleteFile(fnPartial);
				}
				if (fnProjectionMatching.exists())
				{
					MDAux.read(fnProjectionMatching);
					mdReconstructionProjectionMatching[n].unionAll(MDAux);
					deleteFile(fnProjectionMatching);
				}
			}
		}
	}

	// std::cout << "synchronize" << std::endl;
	synchronize();

}
