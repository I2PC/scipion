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
#ifndef _PROG_CONSTRUCT_DICTIONARY_HH
#define _PROG_CONSTRUCT_DICTIONARY_HH

#include <data/xmipp_program.h>

/**@defgroup PDBConstructDictionary Construct a low and high resolution dictionary
   @ingroup ReconsLibrary */
//@{
/** Construct Low and High resolution dictionary*/
class ProgConstructPDBDictionary: public XmippProgram
{
public:
	/** Metadata with the low resolution volumes */
    FileName fnLow;

    /** Metadata with the high resolution volumes */
    FileName fnHigh;

    /** Output rootname */
    FileName fnRoot;

    /** Patch is of size size x size x size */
    int patchSize;

    /** A patch is candidate if its standard deviation is at least this factor of the total standard deviation */
    double stdThreshold;

    double angleThreshold;
public:
    /** Low resolution and high resolution dictionary */
    std::vector< MultidimArray<double> > dictionaryLow, dictionaryHigh;

public:
    void defineParams();
    void readParams();
    void show();
    void run();

    /** True if the patch is not already in the low resolution dictionary */
    bool notInDictionary(const MultidimArray<double> &candidatePatch) const;

    /** Save dictionaries */
    void saveDictionaries() const;
};
//@}
#endif
