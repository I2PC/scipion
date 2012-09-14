/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2002)
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

#include "mpi_angular_gcar_commonlines.h"

// Empty constructor =======================================================
ProgAngularGCARCommonLines::ProgAngularGCARCommonLines(int argc, char **argv)
{
	node=new MpiNode(argc,argv);
    if (!node->isMaster())
    	verbose=0;
}

// MPI destructor
ProgAngularGCARCommonLines::~ProgAngularGCARCommonLines()
{
	delete node;
}

// Read arguments ==========================================================
void ProgAngularGCARCommonLines::readParams()
{
	fnIn = getParam("-i");
	fnRoot = getParam("--oroot");
    max_shift_change = getDoubleParam("--max_shift_change");
    psi_step = getDoubleParam("--psi_step");
    shift_step = getDoubleParam("--shift_step");
}

// Show ====================================================================
void ProgAngularGCARCommonLines::show()
{
    if (!verbose)
        return;
    std::cout
    << "Input:            " << fnIn << std::endl
    << "Output root:      " << fnRoot << std::endl
    << "Psi step:         " << psi_step << std::endl
    << "Max shift change: " << max_shift_change << " step: " << shift_step << std::endl
    ;
}

// usage ===================================================================
void ProgAngularGCARCommonLines::defineParams()
{
    addUsageLine("Makes a rotational invariant representation of the image collection");
    addParamsLine("    -i <selfile>               : Selfile with experimental images");
    addParamsLine("    --oroot <rootname>          : Rootname for output");
    addParamsLine("  [--max_shift_change <r=0>]   : Maximum change allowed in shift");
    addParamsLine("  [--psi_step <ang=5>]         : Step in psi in degrees");
    addParamsLine("  [--shift_step <r=1>]         : Step in shift in pixels");
    addExampleLine("Typical use:",false);
    addExampleLine("xmipp_angular_project_library -i referenceVolume.vol -o reference.stk --sampling_rate 5");
    addExampleLine("xmipp_angular_discrete_assign -i projections.sel -o discrete_assignment.xmd --ref reference.doc");
}

// Produce side info =====================================================
void ProgAngularGCARCommonLines::produceSideInfo()
{
	if (node->isMaster())
	{
		FileName(fnRoot+"_matrix.raw").createEmptyFileWithGivenLength(100*100*sizeof(double));
	}
	node->barrierWait();
	W.mapToFile(fnRoot+"_matrix.raw",100,100);
}

// Run ====================================================================
void ProgAngularGCARCommonLines::run()
{
    show();
    produceSideInfo();
}
