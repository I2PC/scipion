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

#ifndef _CORRECTPHASE_HH
#define _CORRECTPHASE_HH

#include "fourier_filter.h"

#include <data/selfile.h>

/**@defgroup CorrectPhase ctf_correct_phase (CTF Phase correction)
   @ingroup ReconsLibraryPrograms */
//@{
/// Correct Phase parameters
class CorrectPhaseParams
{
public:
    /// ctfdat
    FileName fnCtfdat;

    /// Side Info: CTF
    FourierMask ctf;
    
    /// Side Info: ctfdat
    CTFDat ctfdat;
public:
    /** Read parameters from command line. */
    void read(int argc, char **argv);

    /** Show. */
    void show();

    /** Usage. */
    void usage();

    /** Produce side information.
        The CTFdat is read. */
    void produceSideInfo();

    /** Correct phase of the images in the ctfdat.*/
    void run();
};
//@}
#endif
