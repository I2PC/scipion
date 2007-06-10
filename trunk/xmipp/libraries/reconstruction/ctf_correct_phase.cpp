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

#include "ctf_correct_phase.h"

#include <data/args.h>

/* Read parameters from command line. -------------------------------------- */
void CorrectPhaseParams::read(int argc, char **argv)
{
    fnCtfdat = getParameter(argc, argv, "-ctfdat");
    epsilon = AtoF(getParameter(argc, argv, "-small", "0"));
    std::string aux;
    aux = getParameter(argc, argv, "-method", "");
    if (aux == "remove")                  method = CORRECT_SETTING_SMALL_TO_ZERO;
    else if (aux == "leave" || aux == "") method = CORRECT_LEAVING_SMALL;
    else if (aux == "divide")             method = CORRECT_AMPLIFYING_NOT_SMALL;
}

/* Show -------------------------------------------------------------------- */
void CorrectPhaseParams::show()
{
    std::cout << "ctfdat: " << fnCtfdat << std::endl
              << "Small is under " << epsilon << std::endl
              << "Correcting method: ";
    switch (method)
    {
    case CORRECT_SETTING_SMALL_TO_ZERO:
        std::cout << "Set small values to 0\n";
        break;
    case CORRECT_LEAVING_SMALL:
        std::cout << "Leave small values as they are\n";
        break;
    case CORRECT_AMPLIFYING_NOT_SMALL:
        std::cout << "Correct amplitude except for the small values\n";
        break;
    }
}

/* Usage ------------------------------------------------------------------- */
void CorrectPhaseParams::usage()
{
    std::cerr << "   -ctfdat <CTF descr file or selfile> : It must not be centered\n"
              << "  [-small <epsilon=0>]                 : Values under epsilon are small\n"
              << "  [-method <mth=leave>]                : Valid methods are: remove, leave\n"
              << "                                         divide\n";
    ;
}

/* Produce Side information ------------------------------------------------ */
void CorrectPhaseParams::produceSideInfo()
{
    ctfdat.read(fnCtfdat);
}

void CorrectPhaseParams::readCTF(const FileName &fnCTF)
{
    ctf.FilterBand = CTF;
    ctf.ctf.enable_CTFnoise = false;
    ctf.ctf.read(fnCTF);
    ctf.ctf.Produce_Side_Info();
}

/* Correct a single image -------------------------------------------------- */
//#define DEBUG
void CorrectPhaseParams::correct(Matrix2D< complex<double> > &v)
{
    ctf.generate_mask(v);
#ifdef DEBUG
    std::cout << "New image ----------------------------\n";
#endif

    FOR_ALL_ELEMENTS_IN_MATRIX2D(v)
    {
        complex<double> m = ctf.mask2D(i, j);
        if (m.imag() != 0)
            REPORT_ERROR(1, "CorrectPhase::correct: CTF is not real\n");
#ifdef DEBUG
        std::cout << "CTF at (" << j << "," << i << ")="
                  << m << " Value there " << v(i, j);
#endif
        switch (method)
        {
        case CORRECT_SETTING_SMALL_TO_ZERO:
            if (m.real() < 0)
                if (v(i, j).real() < -epsilon) v(i, j) *= -1;
                else                           v(i, j) = 0;
            break;
        case CORRECT_LEAVING_SMALL:
            if (m.real() < -epsilon) v(i, j) *= -1;
            break;
        case CORRECT_AMPLIFYING_NOT_SMALL:
            if (ABS(m.real()) > epsilon) v(i, j) /= m.real();
            break;
        }
#ifdef DEBUG
        std::cout << " Final value " << v(i, j) << std::endl;
#endif
    }
}
#undef DEBUG

/* Correct a set of images ------------------------------------------------- */
void CorrectPhaseParams::run()
{
    Matrix2D< complex<double> > fft;
    ctfdat.goFirstLine();
    cerr << "Correcting CTF phase ...\n";
    int istep = CEIL((double)ctfdat.lineNo() / 60.0);
    init_progress_bar(ctfdat.lineNo());
    int i = 0;
    while (!ctfdat.eof())
    {
        FileName fnProjection, fnCTF;
	ctfdat.getCurrentLine(fnProjection,fnCTF);
	if (fnProjection!="") {
            ImageXmipp I;
            I.read(fnProjection);
            FourierTransform(I(), fft);
	    readCTF(fnCTF);
            correct(fft);
            InverseFourierTransform(fft, I());
            I.write();
	}
        if (i++ % istep == 0) progress_bar(i);
	ctfdat.nextLine();
    }
    progress_bar(ctfdat.lineNo());
}
