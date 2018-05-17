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

#ifndef _PROG_PDB_REDUCE_PSEUDOATOMS_H_ /*SOFTWARE_EM_XMIPP_LIBRARIES_RECONSTRUCTION_PDB_REDUCE_PSEUDOATOMS_H_*/
#define _PROG_PDB_REDUCE_PSEUDOATOMS_H_ /*SOFTWARE_EM_XMIPP_LIBRARIES_RECONSTRUCTION_PDB_REDUCE_PSEUDOATOMS_H_*/

#include <data/xmipp_program.h>

class ProgPdbReduce: public XmippProgram
{
	FileName fn_volume;

	FileName fn_out;

	double thresh;

	int num;

public:
	/*Empty constructor*/
	ProgPdbReduce();

	/** Params definitions */
	void defineParams();
	/** Read from a command line.
		An exception might be thrown by any of the internal conversions,
		this would mean that there is an error in the command line and you
		might show a usage message. */
	void readParams();

	/** Show parameters. */
	void show();

	/** Run. */
	void run();

public:

	void reduceNumberPseudoatoms();

};




#endif /* SOFTWARE_EM_XMIPP_LIBRARIES_RECONSTRUCTION_PDB_REDUCE_PSEUDOATOMS_H_ */
