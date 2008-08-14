/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
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

    /// Epsilon=minimum CTF value to correct
    double epsilon;

#define CORRECT_SETTING_SMALL_TO_ZERO 0
#define CORRECT_LEAVING_SMALL         1
#define CORRECT_AMPLIFYING_NOT_SMALL  2
    /** Correcting method. Valid methods are:
       \\CORRECT_SETTING_SMALL_TO_ZERO: Values where the CTF<epsilon are
          set to 0
       \\CORRECT_LEAVING_SMALL: Values where the CTF<epsilon are
          left as they are
       \\CORRECT_AMPLIFYING_NOT_SMALL: Values where the ABS(CTF)>epsilon are
          divided by the CTF
    */
    int method;

    /// Side Info: CTF
    FourierMask ctf;
    
    /// Side Info: ctfdat
    CTFDat ctfdat;
public:
    /** Empty constructor */
    CorrectPhaseParams(): epsilon(0), method(0)
    {}

    /** Read parameters from command line. */
    void read(int argc, char **argv);

    /** Show. */
    void show();

    /** Usage. */
    void usage();

    /** Produce side information.
        The CTFdat is read. */
    void produceSideInfo();

    /** Read a CTF.. */
    void readCTF(const FileName &fnCTF);

    /** Correct phase of the images in the ctfdat.*/
    void run();

    /** Correct phase of an image.
        An exception is thrown if the input image is not of the same size
        as the ctf or if the CTF is not real */
    void correct(Matrix2D< std::complex<double> > &v);
};
//@}
#endif
