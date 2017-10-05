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
#ifndef _PROG_VOLUME_TEXTURE_HH
#define _PROG_VOLUME_TEXTURE_HH

#include <data/xmipp_program.h>
#include <data/mask.h>

/**@defgroup VolumeTexture Volume texture
   @ingroup ReconsLibrary */
//@{
/** Generic class to analyze the texture within a volume*/
class ProgVolumeTexture: public XmippProgram
{
public:
	/** Input volume */
    FileName fnIn;

    /** Reference volume */
    FileName fnRef;

    /** Root path for output */
    FileName fnRoot;

    /** Patch is of size size x size x size */
    int patchSize;

    /** Mask */
    Mask mask;
public:
    virtual void defineParams();
    virtual void readParams();
    virtual void show();
    virtual void run();
};
//@}
#endif
