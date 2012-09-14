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

#include "ctf_correct_phase.h"

#include <data/args.h>
#include <data/fftw.h>
#include <data/metadata.h>

/* Read parameters from command line. -------------------------------------- */
void CorrectPhaseParams::read(int argc, char **argv)
{
    fnCtfdat = getParameter(argc, argv, "-ctfdat");
}

/* Show -------------------------------------------------------------------- */
void CorrectPhaseParams::show()
{
    std::cout << "ctfdat: " << fnCtfdat << std::endl;
}

/* Usage ------------------------------------------------------------------- */
void CorrectPhaseParams::usage()
{
    std::cerr << "   -ctfdat <CTF descr file or selfile> : list of particles and CTFs\n";
}

/* Produce Side information ------------------------------------------------ */
void CorrectPhaseParams::produceSideInfo()
{
    ctfdat.read(fnCtfdat);
}

/* Correct a set of images ------------------------------------------------- */
void CorrectPhaseParams::run()
{
    MultidimArray< std::complex<double> > fft;
    FourierTransformer transformer;
    std::cerr << "Correcting CTF phase ...\n";
    int istep = CEIL((double)ctfdat.size() / 60.0);
    init_progress_bar(ctfdat.size());
    int i = 0;
    FOR_ALL_OBJECTS_IN_METADATA(ctfdat)
    {
        FileName fnProjection, fnCTF;
        if (!ctfdat.getValue(MDL_IMAGE,fnProjection))
        	REPORT_ERROR(ERR_MD_NOOBJ,(std::string)"Cannot find images in "+fnCtfdat);
        if (!ctfdat.getValue(MDL_CTF_MODEL,fnCTF))
        	REPORT_ERROR(ERR_MD_NOOBJ,(std::string)"Cannot find CTF for "+fnProjection);

        // Read input image and compute its Fourier transform
        Image<double> I;
        I.read(fnProjection);
        transformer.FourierTransform(I(), fft, false);

        // Read the CTF
        ctf.FilterBand = CTF;
        ctf.ctf.enable_CTFnoise = false;
        ctf.ctf.read(fnCTF);
        ctf.ctf.Produce_Side_Info();

        // Apply the phase correction
        Matrix1D<double> w(2);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(fft)
        {
            FFT_IDX2DIGFREQ(i,YSIZE(I()),YY(w));
            FFT_IDX2DIGFREQ(j,XSIZE(I()),XX(w));
            if (ctf.maskValue(w)<0)
                fft(i, j) *= -1;
        }

        // Go back to real space
        transformer.inverseFourierTransform();
        I.write();

        if (i++ % istep == 0)
            progress_bar(i);
    }
    progress_bar(ctfdat.size());
}
