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

#include <data/progs.h>
#include <data/args.h>
#include <data/fftw.h>

class FFT_parameters: public Prog_parameters
{
public:
#define COMPLETE_FFT    0
#define ONLY_AMPLITUDES 1
#define ONLY_PHASE      2
    int FFT_mode;
    bool apply_log;
    bool squared;
    bool do_not_center;
    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        if (checkParameter(argc, argv, "-phase"))      FFT_mode = ONLY_PHASE;
        else if (checkParameter(argc, argv, "-amplitudes")) FFT_mode = ONLY_AMPLITUDES;
        else                                           FFT_mode = COMPLETE_FFT;
        apply_log = checkParameter(argc, argv, "-log10");
        squared = checkParameter(argc, argv, "-squared");
        do_not_center = checkParameter(argc, argv, "-do_not_center");
    }

    void show()
    {
        Prog_parameters::show();
        switch (FFT_mode)
        {
        case COMPLETE_FFT:
            std::cout << "Computing the whole FFT\n";
            break;
        case ONLY_AMPLITUDES:
            std::cout << "Computing Amplitude FFT\n";
            break;
        case ONLY_PHASE:
            std::cout << "Computing Phase FFT\n";
            break;
        }
        if (apply_log) std::cout << "Computing Log of result\n";
        if (squared)   std::cout << "Squaring amplitudes\n";
        if (do_not_center) std::cout << "Not centering result\n";
        else               std::cout << "Centering result\n";
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr << "  [-phase]                  : By default, the whole FFT\n"
        << "  [-amplitudes]             : is returned\n"
        << "  [-log10]                  : Return logarithm of result\n"
        << "  [-squared]                : Return the square of the result\n"
        << "  [-do_not_center]          : By default, the result is centered\n";
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    FFT_parameters *eprm = (FFT_parameters *) prm;
    Matrix2D< std::complex<double> > fftI;
    FourierTransform(img(), fftI);
    switch (eprm->FFT_mode)
    {
    case ONLY_PHASE:
        FFT_phase(fftI, img());
        break;
    case ONLY_AMPLITUDES:
        FFT_magnitude(fftI, img());
        break;
    }
    if (eprm->squared) img() *= img();
    if (eprm->apply_log)
        FOR_ALL_ELEMENTS_IN_MATRIX2D(img())
            MAT_ELEM(img(), i, j) = log10(1 + MAT_ELEM(img(), i, j));
    if (!eprm->do_not_center) CenterFFT(img(), true);
    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    FFT_parameters *eprm = (FFT_parameters *) prm;
    Matrix3D< std::complex<double> > fftV;
    FourierTransform(vol(), fftV);
    switch (eprm->FFT_mode)
    {
    case ONLY_PHASE:
        FFT_phase(fftV, vol());
        break;
    case ONLY_AMPLITUDES:
        FFT_magnitude(fftV, vol());
        break;
    }
    if (eprm->squared) vol() *= vol();
    if (eprm->apply_log)
        FOR_ALL_ELEMENTS_IN_MATRIX3D(vol())
            VOL_ELEM(vol(), k, i, j) = log10(1 + VOL_ELEM(vol(), k, i, j));
    if (!eprm->do_not_center) CenterFFT(vol(), true);
    return true;
}

int main(int argc, char **argv)
{
    FFT_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

