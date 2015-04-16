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

#include "pdb_nma_deform.h"

void ProgPdbNmaDeform::defineParams()
{
	addUsageLine("Deform a PDB according to a list of NMA modes");
	addUsageLine("+It is assumed that there is a deformation value per mode, and that the number of vectors in the mode is the same as the number of atoms in the PDB");
	addParamsLine("--pdb <file>    : PDB to deform");
	addParamsLine("--nma <metadata_file>  : List of modes, the metadata must use the label NMAModefile");
	addParamsLine("--deformations <...> : Deformation amplitudes in each direction");
	addParamsLine("-o <file>       : Deformed PDB");
	addExampleLine("xmipp_pdb_nma_deform --pdb 2tbv.pdb -o 2tbv_deformed.pdb --nma modelist.xmd --deformations 15");
}

void ProgPdbNmaDeform::readParams()
{
	fn_pdb=getParam("--pdb");
	fn_nma=getParam("--nma");
	fn_out=getParam("-o");
	StringVector listDeformations;
	getListParam("--deformations",listDeformations);
	deformations.resizeNoCopy(listDeformations.size());
	FOR_ALL_ELEMENTS_IN_ARRAY1D(deformations)
	deformations(i)=textToFloat(listDeformations[i]);
}

void ProgPdbNmaDeform::show()
{
	if (verbose==0)
		return;
	std::cout
	<< "PDB:          " << fn_pdb << std::endl
	<< "NMA list:     " << fn_nma << std::endl
	<< "Output:       " << fn_out << std::endl
	<< "Deformations: "
	;
	FOR_ALL_ELEMENTS_IN_ARRAY1D(deformations)
		std::cout << deformations(i) << " ";
	std::cout << std::endl;
}

void ProgPdbNmaDeform::run()
{
	PDBRichPhantom pdb;
	MetaData modes;
	pdb.read(fn_pdb);
	modes.read(fn_nma);
	modes.removeDisabled();
	int j=0;
	FileName fnMode;
	MultidimArray<double> mode;
	mode.resizeNoCopy(pdb.getNumberOfAtoms(),3);
	FOR_ALL_OBJECTS_IN_METADATA(modes)
	{
		// Read mode
		modes.getValue(MDL_NMA_MODEFILE,fnMode,__iter.objId);
		std::ifstream fhMode;
		fhMode.open(fnMode.c_str());
		if (!fhMode)
			REPORT_ERROR(ERR_IO_NOREAD,fnMode);
		fhMode >> mode;
		fhMode.close();

		// Apply mode
		double lambda=A1D_ELEM(deformations,j);
		for (size_t i=0; i<YSIZE(mode); ++i)
		{
			RichAtom& atom_i=pdb.atomList[i];
			atom_i.x+=lambda*DIRECT_A2D_ELEM(mode,i,0);
			atom_i.y+=lambda*DIRECT_A2D_ELEM(mode,i,1);
			atom_i.z+=lambda*DIRECT_A2D_ELEM(mode,i,2);
		}
                j++;
	}
	pdb.write(fn_out);
}
