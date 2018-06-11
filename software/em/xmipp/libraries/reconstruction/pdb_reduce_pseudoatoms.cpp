/***************************************************************************
 *
 * Authors: Javier Mota Garcia (jmota@cnb.csic.es)
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

#include "pdb_reduce_pseudoatoms.h"
#include <data/pdb.h>
#include <data/args.h>
#include <fstream>

ProgPdbReduce::ProgPdbReduce()
{
	thresh = 0.0;
	num = -1.0;
}

void ProgPdbReduce::defineParams()
{
	addUsageLine("Reduce the number of pseudoatoms in a volume.");
	addExampleLine("   xmipp_pdb_reduce -i 1o7d.vol -o 1o7dreduced.vol --threshold 0.15");

	addParamsLine("   -i <pdb_file>                          : File to process");
	addParamsLine("  [-o <fn_root>]                          : Root name for output");
	addParamsLine("  [--number <num=-1.0>]                : Sampling rate (Angstroms/pixel)");
	addParamsLine("  [--threshold <thresh=0.0>]                : Sampling rate (Angstroms/pixel)");
}

void ProgPdbReduce::readParams()
{
	fn_volume = getParam("-i");
	fn_out = checkParam("-o") ? getParam("-o") : fn_volume.withoutExtension();
	num = getDoubleParam("--number");
	thresh = getDoubleParam("--threshold");

}

/* Show -------------------------------------------------------------------- */
void ProgPdbReduce::show()
{
    std::cout << "PDB file:           " << fn_volume           << std::endl
    		  << "Number:             " << num                 << std::endl
              << "Threshold:          " << thresh              << std::endl;
}

void ProgPdbReduce::reduceNumberPseudoatoms()
{
	PDBRichPhantom pdb;
	if(thresh != 0.0)
		pdb.read(fn_volume, -1.0, thresh);
	else
	{
		pdb.read(fn_volume, num);
		std::sort(pdb.intensities.rbegin(), pdb.intensities.rend());
		thresh = pdb.intensities.at(num);
		pdb.atomList.clear();
		pdb.remarks.clear();
		pdb.read(fn_volume, -1.0, thresh);
	}
	pdb.write(fn_out);
}

void ProgPdbReduce::run()
{
	reduceNumberPseudoatoms();
}


