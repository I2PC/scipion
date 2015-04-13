/***************************************************************************
 *
 * Authors:    Carlos Oscar Sorzano  coss@cnb.csic.es (2014)
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
#ifndef _PROG_METADATA_SPLIT3D
#define _PROG_METADATA_SPLIT3D

#include <data/xmipp_program.h>
#include <data/metadata.h>

/**@defgroup MetadataSplit3D Split a metadata according to a 3D volume
   @ingroup ReconsLibrary */
//@{
/** MetadataSplit3D parameters. */
class ProgMetadataSplit3D: public XmippProgram
{
public:
	/** Input metadata with images */
	FileName fn_in;
    /** Reference volume */
    FileName fn_vol;
    /** Filename for the symmetry file */
    FileName fn_sym;
    /** Rootname for output */
    FileName fn_oroot;
    /** Angular sampling rate */
    double angularSampling;
    /** Maximum angular distance */
    double maxDist;
public:
    // Metadata with reference images
    MetaData mdRef;

    // Metadata with experimental images
    MetaData mdIn;

    // Correlates_well vector
    Matrix1D<int> correlatesWell;
public:
    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Usage
    void defineParams();

    /// Run
    void run();
};
//@}
#endif
