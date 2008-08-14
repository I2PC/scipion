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

#include "ctf_correct_amplitude2d.h"

#include <data/args.h>

/* Read parameters from command line. -------------------------------------- */
void CorrectAmplitude2DParams::read(int argc, char **argv)
{
    fnCtfdat = getParameter(argc, argv, "-ctfdat");
    beta = textToFloat(getParameter(argc, argv, "-relax","0.1"));
    Niterations = textToInteger(getParameter(argc, argv, "-iterations","100"));
}

/* Show -------------------------------------------------------------------- */
void CorrectAmplitude2DParams::show()
{
    std::cout << "Ctfdat:     " << fnCtfdat << std::endl
              << "Relaxation: " << beta << std::endl
              << "Iterations: " << Niterations << std::endl;
}

/* Usage ------------------------------------------------------------------- */
void CorrectAmplitude2DParams::usage()
{
    std::cerr << "   -ctfdat <CTF descr file or selfile> : list of particles and CTFs\n"
              << "  [-relax <beta=0.1>]                  : relaxation factor\n"
              << "  [-iterations <N=100>]                : number of iterations\n"
    ;
}

/* Produce Side information ------------------------------------------------ */
void CorrectAmplitude2DParams::produceSideInfo()
{
    ctfdat.read(fnCtfdat);
}

void CorrectAmplitude2DParams::readCTF(const FileName &fnCTF)
{
    ctf.FilterBand = CTF;
    ctf.ctf.enable_CTFnoise = false;
    ctf.ctf.read(fnCTF);
    ctf.ctf.Produce_Side_Info();
}

/* Correct a single image -------------------------------------------------- */
void CorrectAmplitude2DParams::correct(Matrix2D< std::complex<double> > &g)
{
    ctf.generate_mask(g);

    Matrix2D< std::complex<double> > f;
    f=g;
    for (int n=0; n<Niterations; n++)
        FOR_ALL_ELEMENTS_IN_MATRIX2D(f)
            f(i,j)=f(i,j)+beta*ctf.mask2D(i, j)*(g(i,j)-ctf.mask2D(i, j)*f(i,j));
    g=f;
}

/* Correct a set of images ------------------------------------------------- */
void CorrectAmplitude2DParams::run()
{
    Matrix2D< std::complex<double> > fft;
    ctfdat.goFirstLine();
    std::cerr << "Correcting CTF amplitude in 2D ...\n";
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
