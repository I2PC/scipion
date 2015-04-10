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
#ifndef _PROG_RESTORE_WITH_DICTIONARY_HH
#define _PROG_RESTORE_WITH_DICTIONARY_HH

#include "pdb_construct_dictionary.h"

/**@defgroup PDBRestoreWithDictionary Restore a volume with a low/high dictionary
   @ingroup ReconsLibrary */
//@{
/** Restore with Low and High resolution dictionary*/
class ProgRestoreWithPDBDictionary: public ProgPDBDictionary
{
public:
	/** Input and output volume names */
	FileName fnIn, fnOut;

    double lambda;

    int iterations;
public:
    Matrix2D<double> Ui, UitUi;
	Matrix1D<double> wi, v1, v2, y, yp;

public:
    void defineParams();
    void readParams();
    void show();
    void run();

    void selectDictionaryPatches(const MultidimArray<double> &lowResolutionPatch, Matrix1D<double> &lowResolutionPatchSignature,
    		std::vector<size_t> &selectedPatchesIdx, std::vector<double> &weight);

    double approximatePatch(const MultidimArray<double> &lowResolutionPatch,
    		std::vector< size_t > &selectedPatchesIdx, std::vector<double> &weight, Matrix1D<double> &alpha);

    void reconstructPatch(size_t idxTransf, std::vector< size_t > &selectedPatchesIdx, Matrix1D<double> &alpha,
    		MultidimArray<double> &highResolutionPatch);
};
//@}
#endif
