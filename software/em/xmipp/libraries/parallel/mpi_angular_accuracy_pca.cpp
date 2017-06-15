/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
 *
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

#include "mpi_angular_accuracy_pca.h"

MpiProgAngularAccuracyPCA::MpiProgAngularAccuracyPCA()
{
	node=NULL;
}

MpiProgAngularAccuracyPCA::~MpiProgAngularAccuracyPCA()
{
	delete node;
}

void MpiProgAngularAccuracyPCA::read(int argc, char** argv)
{
    node = new MpiNode(argc, argv);
   	rank = node->rank;
   	Nprocessors = node->size;
   	ProgAngularAccuracyPCA::read(argc, (const char **)argv);
}

void MpiProgAngularAccuracyPCA::synchronize()
{
	node->barrierWait();
}

void MpiProgAngularAccuracyPCA::gatherResults()
{
	if (rank!=0)
	{
		FileName fnPartial=formatString("%s/partial_node%03d.xmd",fnOut.getDir().c_str(),(int)rank);
		if (mdPartial.size()>0)
			mdPartial.write(fnPartial);
	}

	synchronize();

	// Now the master takes all of them
	if (rank==0)
	{
		MetaData MDAux;
		for (size_t otherRank=1; otherRank<Nprocessors; ++otherRank)
		{
				FileName fnP = formatString("%s/partial_node%03d.xmd",fnOut.getDir().c_str(),(int)otherRank);

				if (fnP.exists())
				{
					MDAux.read(fnP);
					mdPartial.unionAll(MDAux);
					deleteFile(fnP);
				}
		}
	}

	synchronize();

}

