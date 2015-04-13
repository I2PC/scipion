/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *              Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "ctf_phase_flip.h"

void ProgCTFPhaseFlipping::defineParams()
{
    addUsageLine("Correct the phase of micrographs");
    addUsageLine("+This program flips the phase of those frequencies that were already ");
    addUsageLine("+flipped by the CTF. Flipping the phase at the level of the micrograph is recommended.");
    addParamsLine(" -i <file>               : Input micrograph");
    addParamsLine(" -o <file>               : Output micrograph");
    addParamsLine(" --ctf <ctfparam_file>   : CTF description");
    addParamsLine(" [--downsampling <D=1>]  : Downsampling factor of the input micrograph with respect to the original");
    addParamsLine("                         : micrograph.");
}

void ProgCTFPhaseFlipping::readParams()
{
    fn_in        = getParam("-i");
    fn_out       = getParam("-o");
    fnt_ctf      = getParam("--ctf");
    downsampling = getDoubleParam("--downsampling");
}

void ProgCTFPhaseFlipping::show()
{
    if (verbose==0)
        return;
    std::cout
    << "input_micrograph:      " << fn_in        << std::endl
    << "output_micrograph:     " << fn_out       << std::endl
    << "ctf_param_file:        " << fnt_ctf      << std::endl
    << "downsampling:          " << downsampling << std::endl
    ;
}

void ProgCTFPhaseFlipping::run()
{
    show();

    // Read the micrograph in an array of doubles
    Image<double> M_in;
    M_in.read(fn_in);

    // Read CTF
    CTFDescription ctf;
    ctf.clear();
    ctf.read(fnt_ctf);
    ctf.changeSamplingRate(ctf.Tm*downsampling);
    ctf.produceSideInfo();

    actualPhaseFlip(M_in(),ctf);

    M_in.write(fn_out);
}

void actualPhaseFlip(MultidimArray<double> &I, CTFDescription ctf)
{
    // Perform the Fourier transform
    FourierTransformer transformer;
    MultidimArray< std::complex<double> > M_inFourier;
    transformer.FourierTransform(I,M_inFourier,false);

    Matrix1D<double> freq(2); // Frequencies for Fourier plane
    int yDim=YSIZE(I);
    int xDim=XSIZE(I);
    double iTm=1.0/ctf.Tm;
    for (size_t i=0; i<YSIZE(M_inFourier); ++i)
    {
    	FFT_IDX2DIGFREQ(i, yDim, YY(freq));
    	YY(freq) *= iTm;
        for (size_t j=0; j<XSIZE(M_inFourier); ++j)
        {
        	FFT_IDX2DIGFREQ(j, xDim, XX(freq));
        	XX(freq) *= iTm;
            ctf.precomputeValues(XX(freq),YY(freq));
            if (ctf.getValuePureWithoutDampingAt()<0)
                DIRECT_A2D_ELEM(M_inFourier,i,j)*=-1;
        }
    }

    // Perform inverse Fourier transform and finish
    transformer.inverseFourierTransform();
}
