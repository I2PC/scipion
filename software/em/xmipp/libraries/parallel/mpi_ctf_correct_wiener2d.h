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

#ifndef MPI_CTF_CORRECT_WIENER2D_H_
#define MPI_CTF_CORRECT_WIENER2D_H_

#include <reconstruction/ctf_correct_wiener2d.h>
#include "parallel/xmipp_mpi.h"

/**@defgroup MpiProgCorrectWiener2D performs CTF correction on 2D images by wiener filtering (MPI)
   @ingroup ParallelLibrary */
//@{

class MpiProgCorrectWiener2D: public ProgCorrectWiener2D
{
public:
	MpiNode *node;
public:
	// Empty constructor
	MpiProgCorrectWiener2D();

	// Destructor
	~MpiProgCorrectWiener2D();

	// Redefine how to read the command line
	void read(int argc, char** argv);

	// Redefine how to synchronize
	void synchronize();

	// Redefine how to gather the alignment
    void gatherResults();
};
//@}

#endif /* MPI_CTF_CORRECT_WIENER2D_H_ */
