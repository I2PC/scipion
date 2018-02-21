/***************************************************************************
 *
 * Authors:
 *
 * Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "pdb_analysis.h"

void ProgPdbAnalysis::defineParams()
{
	addUsageLine("Analyze a PDB file and produce different statistics");
	addParamsLine("-i <file>    : PDB to study");
	addParamsLine("--operation <op>  : Operation to perform");
	addParamsLine("    where <op>");
	addParamsLine("          distance_histogram <fileOut> <Nnearest=3> <MaxDistance=-1>: Compute the distance histogram between");
	addParamsLine("                                                                    : an atom and its N nearest neighbours");
	addParamsLine("                                                                    : The maximum distance of the neighbours may be limited (in Angstroms)");
	addParamsLine("          rmsd <pdb2>: RMSD between pdb1 and pdb2, they are supposed to have the same atoms and in the same order");
	addExampleLine("Compute the histogram of interatomic distances",false);
	addExampleLine("xmipp_pdb_analysis -i mypdb.pdb --operation distance_histogram distance.hist");
}

void ProgPdbAnalysis::readParams()
{
	fn_pdb=getParam("-i");
	op=getParam("--operation");
	if (op=="distance_histogram")
	{
		fn_hist=getParam("--operation",1);
		Nnearest=getIntParam("--operation",2);
		maxDistance=getDoubleParam("--operation",3);
	}
	else if (op=="rmsd")
	{
		fn_pdb2=getParam("--operation",1);
	}
}

void ProgPdbAnalysis::show()
{
	if (verbose==0)
		return;
	std::cout
	<< "PDB:          " << fn_pdb << std::endl
	<< "Operation:    " << op << std::endl;
	if (op=="distance_histogram")
		std::cout << "Output histogram: " << fn_hist << std::endl
		          << "Nnearest:         " << Nnearest << std::endl
				  << "MaxDistance:      " << maxDistance << std::endl;
	else if (op=="rmsd")
	{
		std::cout << "PDB2:             " << fn_pdb2 << std::endl;
	}
}

void ProgPdbAnalysis::run()
{
	PDBPhantom pdb; // It cannot be a PDBRichAtom because it also has to work with pseudoatomic structures
	pdb.read(fn_pdb);
	if (op=="distance_histogram")
	{
		Histogram1D hist;
		distanceHistogramPDB(pdb,Nnearest,maxDistance,200,hist);
		hist.write(fn_hist);
	}
	else if (op=="rmsd")
	{
		PDBPhantom pdb2;
		pdb2.read(fn_pdb2);
		std::cout << "RMSD= " << computeRMSD(pdb,pdb2) << std::endl;
	}
}
