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

#include <data/xmipp_program.h>

/**@defgroup PDBRestoreWithDictionary Restore a volume with a low/high dictionary
   @ingroup ReconsLibrary */
//@{
/** Restore with Low and High resolution dictionary*/
class ProgRestoreWithPDBDictionary: public XmippProgram
{
public:
	/** Input and output volume names */
	FileName fnIn, fnOut;

    /** Dictionary rootname */
    FileName fnDictionaryRoot;

    /** A patch is candidate if its standard deviation is at least this factor of the total standard deviation */
    double stdThreshold;

    double angleThreshold;

    double lambda;

    int iterator;
public:
    /** Low resolution and high resolution dictionary */
    Image<double> dictionaryLow, dictionaryHigh;

    Matrix2D<double> Ui, UitUi;
	Matrix1D<double> wi, v1, v2, y, yp;

public:
    void defineParams();
    void readParams();
    void show();
    void run();

    void selectDictionaryPatches(const MultidimArray<double> &lowResolutionPatch, std::vector<size_t> &idx,
    		std::vector<double> &weight);

    double approximatePatch(const MultidimArray<double> &lowResolutionPatch, std::vector<size_t> &idx,
    		std::vector<double> &weight, Matrix1D<double> &alpha);

    void reconstructPatch(std::vector<size_t> &idx, Matrix1D<double> &alpha,
    		   MultidimArray<double> &highResolutionPatch);
};
//@}
#endif
