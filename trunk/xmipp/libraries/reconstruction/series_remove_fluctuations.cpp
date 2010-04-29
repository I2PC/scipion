/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2009)
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

#include "series_remove_fluctuations.h"
#include <data/args.h>
#include <data/fftw.h>

// Read from command line --------------------------------------------------
void Series_remove_fluctuations_parameters::read(int argc, char **argv)
{
    fn_sel  = getParameter(argc, argv, "-i");
    fn_root = getParameter(argc, argv, "-oroot", "");
    maxFreq = textToFloat(getParameter(argc, argv, "-lpf", "0.25"));
}

// Produce side info -------------------------------------------------------
void Series_remove_fluctuations_parameters::produceSideInfo()
{
    // Read the selfile into a volume
    SF.read(fn_sel);
    int Zdim, Ydim, Xdim;
    Zdim=SF.ImgNo();
    SF.ImgSize(Ydim,Xdim);
    V.initZeros(Zdim,Ydim,Xdim);
    int k=0;
    while (!SF.eof())
    {
        ImageXmipp I;
        I.read(SF.NextImg());
        V.setSlice(k,I());
        k++;
    }
}

// Show --------------------------------------------------------------------
void Series_remove_fluctuations_parameters::show() const
{
    std::cout << "Removing fluctuations from a tilt series\n";
    std::cout << "Input series:   " << fn_sel  << std::endl
              << "Output root:    " << fn_root << std::endl
              << "Max freq (LPF): " << maxFreq << std::endl
              ;
}

// Usage -------------------------------------------------------------------s
void Series_remove_fluctuations_parameters::usage() const
{
    std::cerr << "series_remove_fluctuations\n"
              << "  -i <selfile>               : Input images\n"
              << " [-oroot <rootname>]         : If not given, the input is rewritten\n"
              << " [-lpf <f=0.25>]             : Max. Freq. for the low pass filter\n"
              ;
}

// Denoise image -----------------------------------------------------------
void Series_remove_fluctuations_parameters::run()
{
    // Filter each line in the series
    int maxPixel=CEIL(ZSIZE(V)*maxFreq);
    std::cout << "Maxpixel=" << maxPixel << std::endl;
    XmippFftw transformer;
    init_progress_bar(YSIZE(V));
    for (int i=0; i<YSIZE(V); i++)
    {
        for (int j=0; j<XSIZE(V); j++)
        {
            // Get line
            Matrix1D<double> line;
            line.initZeros(ZSIZE(V));
            for (int k=0; k<ZSIZE(V); k++)
                line(k)=V(k,i,j);

            // Fourier transform
            transformer.setReal(line);
            transformer.FourierTransform();

            // Filter
            Matrix1D< std::complex<double> > lineFourier;
            transformer.getFourierAlias(lineFourier);
            for (int k=maxPixel; k<XSIZE(lineFourier); k++)
                lineFourier(k)=0;

            // Inverse Fourier transform and back to the volume
            transformer.inverseFourierTransform();
            for (int k=0; k<ZSIZE(V); k++)
                V(k,i,j)=line(k);
        }
        progress_bar(i+1);
    }

    // Write the results
    SF.go_first_ACTIVE();
    int k=0;
    while (!SF.eof())
    {
        ImageXmipp I;
        V.getSlice(k,I());

        FileName fnimg;
        if (fn_root=="") fnimg=SF.NextImg();
        else fnimg.compose(fn_root,k,"xmp");
        I.write(fnimg);

        k++;
    }
}
