/***************************************************************************
 *
 * Authors:     Roberto Marabini
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

#ifndef CTF_PHASE_FLIPPING_H
#define CTF_PHASE_FLIPPING_H
#include <data/args.h>
#include <data/micrograph.h>
#include <data/xmipp_fftw.h>
#include <data/ctf.h>
#include <data/xmipp_fft.h>

/**@defgroup Micrograph_phase_Flipp (phase flipping)
   @ingroup ReconsLibrary
*/
//@{

class ProgCTFPhaseFlipping: public XmippProgram
{
public:
    /** Input micrograph */
    Micrograph M_in;
    /** Input micrograph file name */
    FileName fn_in;
    /** Out flipped micrograph */
    FileName fn_out;
    /** CTF param file name*/
    FileName fnt_ctf;
    /** Downsampling factor */
    double downsampling;

    /** Define parameters */
    void defineParams();
    /** Read parameters */
    void readParams();
    /** Do_all. */
    void run();
    /** Show parameters. */
    void show();
};

void actualPhaseFlip(MultidimArray<double> &I, CTFDescription ctf);

//@}
#endif
