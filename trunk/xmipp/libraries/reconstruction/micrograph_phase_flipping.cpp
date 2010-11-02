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

#include "micrograph_phase_flipping.h"

void Prog_micrograph_phase_flipping::show(void)
{
    std::cout
        << "input_micrograph:      " << fn_in    << std::endl
        << "output_micrograph:     " << fn_out << std::endl
        << "ctf_param_file:        " << fnt_ctf << std::endl
    ;
}

void Prog_micrograph_phase_flipping::run(void)
{
    // Read the micrograph in an array of doubles
    Image<double> M_in;
    M_in.read(fn_in);
    
    // Perform the Fourier transform
    FourierTransformer transformer;
    MultidimArray< std::complex<double> > M_inFourier;
    transformer.FourierTransform(M_in(),M_inFourier,false);

    // Read CTF
    CTFDescription ctf;
    ctf.clear();
    ctf.read(fnt_ctf);
    ctf.Produce_Side_Info();

    Matrix1D<int>    idx(2);  // Indexes for Fourier plane
    Matrix1D<double> freq(2); // Frequencies for Fourier plane
    FOR_ALL_ELEMENTS_IN_ARRAY2D(M_inFourier)
    {
        VECTOR_R2(idx,j,i);
        FFT_idx2digfreq(M_in(),idx,freq);
        digfreq2contfreq(freq, freq, ctf.Tm);
        ctf.precomputeValues(XX(freq),YY(freq));
        if (ctf.CTFpure_without_damping_at()<0)
        {
            M_inFourier(i,j)*=-1;
        }
    }
    
    // Perform inverse Fourier transform and finish
    transformer.inverseFourierTransform();
    M_in.write(fn_out);
}
