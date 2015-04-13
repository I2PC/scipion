/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

#include "ctf_show.h"

void ProgCTFShow::defineParams()
{
    addUsageLine("Convert parametric CTF into image");
    addParamsLine(" -i <file>               : parametric_CTF");
    addParamsLine(" -o <file>               : image CTF");
    addParamsLine(" --dim <xdim>               : image size");
    addExampleLine("Create CTF image from parametric CTF file ctf.param and with dimension 256 pixels:", false);
    addExampleLine("xmipp_ctf_show -i ctf.param -o 20.mrc --dim 256");
}

void ProgCTFShow::readParams()
{
    fn_in   = getParam("-i");
    fn_out  = getParam("-o");
    dim     = getIntParam("--dim");
}
#include <iostream>

void ProgCTFShow::show()
{
    if (verbose==0)
        return;
    std::cout
    << "parametric CTF:      " << fn_in    << std::endl
    << "image CTF:           " << fn_out << std::endl
    ;
}

void ProgCTFShow::run()
{
    show();

    // Read CTF
    CTFDescription ctf;
    ctf.clear();
    ctf.read(fn_in);
    ctf.produceSideInfo();

    Matrix1D<double> freq(2); // Frequencies for Fourier plane
    int yDim=dim;
    int xDim=dim;
    double iTm=1.0/ctf.Tm;
    Image<double> img(xDim,yDim);
    MultidimArray<double> &dummy = img.data;
    int x2Dim = xDim / 2 ;
    int y2Dim = yDim / 2 ;

    for (int i=0; i<yDim; ++i)
    {
        FFT_IDX2DIGFREQ(i, yDim, YY(freq));
        YY(freq) *= iTm;
        for (int j=0; j<xDim; ++j)
        {
            FFT_IDX2DIGFREQ(j, xDim, XX(freq));
            XX(freq) *= iTm;
            ctf.precomputeValues(XX(freq),YY(freq));
            DIRECT_A2D_ELEM(dummy,intWRAP(i+x2Dim,0,xDim-1),
            		              intWRAP(j+y2Dim,0,yDim-1))
                              =ctf.getValuePureWithoutDampingAt();
        }
    }

    // Perform inverse Fourier transform and finish
    img.write(fn_out);
}
