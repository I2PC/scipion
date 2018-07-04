/***************************************************************************
 *
 * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
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

#ifndef _PROG_COORDINATES_NOISY_ZONES_FILTER_HH
#define _PROG_COORDINATES_NOISY_ZONES_FILTER_HH

#include <data/xmipp_program.h>
#include <data/micrograph.h>
#include <data/morphology.h>


/**@defgroup CoordinatesMicrographFilters Coordinates Micrograph Filter
   @ingroup ReconsLibrary */
//@{
/** Generic class to filtrate the coordinates in noisy zones */
class ProgCoordinatesNoisyZonesFilter: public XmippProgram
{
public:
    /** Input coordinates */
    FileName fnInCoord;

    /** Input Micrograph */
    FileName fnInMic;

    /** Output coordinates*/
    FileName fnOut;

    /** Patch is of size: size x size */
    int patchSize;

public:
    virtual void defineParams();
    virtual void readParams();
    virtual void show();
    virtual void run();
};
//@}
#endif
