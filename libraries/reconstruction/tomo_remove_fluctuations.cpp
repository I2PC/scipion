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

#include "tomo_remove_fluctuations.h"
#include <data/args.h>
#include <data/xmipp_fftw.h>
#include <data/metadata_extension.h>
#include <data/xmipp_image.h>

// Read from command line --------------------------------------------------
void ProgTomoRemoveFluctuations::readParams()
{
    fnIn  = getParam("-i");
    fnOut = getParam("-o");
    maxFreq = getDoubleParam("--lpf");
}

// Usage -------------------------------------------------------------------
void ProgTomoRemoveFluctuations::defineParams()
{
	addUsageLine("Remove the flickering in a tilt series. This is a phenomenon rather ");
	addUsageLine("common in X-ray microscopy. For doing so, a low-pass filter of a given ");
	addUsageLine("frequency (normalized to 0.5) is applied across the time series, i.e., a ");
	addUsageLine("line is formed with all the lines at the same position along the tilt series. ");
	addUsageLine("This time series is the filtered and put back into the tilt series. ");
	addParamsLine("  -i <file>                 : Input images (selfile or stack)");
	addParamsLine("  -o <outputStack>          : Output stack");
    addParamsLine(" [--lpf <f=0.25>]           : Low pass filter");
    addParamsLine("                            : Frequency is normalized to 0.5");
}

// Produce side info -------------------------------------------------------
void ProgTomoRemoveFluctuations::produceSideInfo()
{
    // Read the selfile into a volume
    SF.read(fnIn);
    size_t Zdim, dummy, Ydim, Xdim, Ndim;
    Zdim=SF.size();
    getImageSize(SF,Xdim,Ydim,dummy,Ndim);
    V().setMmap(true);
    V().initZeros(Zdim,Ydim,Xdim);
    int k=0;
    Image<double> I;
    FileName fnImg;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        SF.getValue(MDL_IMAGE,fnImg,__iter.objId);
        I.read(fnImg);
        V().setSlice(k,I());
        k++;
    }
}

// Show --------------------------------------------------------------------
void ProgTomoRemoveFluctuations::show() const
{
	if (verbose==0)
		return;
    std::cout << "Removing fluctuations from a tilt series\n";
    std::cout << "Input series:   " << fnIn    << std::endl
              << "Output root:    " << fnOut   << std::endl
              << "Max freq (LPF): " << maxFreq << std::endl
    ;
}

// Denoise image -----------------------------------------------------------
void ProgTomoRemoveFluctuations::run()
{
	produceSideInfo();

	// Filter each line in the series
	MultidimArray<double> &mV=V();
    int maxPixel=CEIL(ZSIZE(mV)*maxFreq);
    FourierTransformer transformer;
    if (verbose>=1)
    	init_progress_bar(YSIZE(mV));
    MultidimArray<double> line;
    MultidimArray< std::complex<double> > lineFourier;
    std::complex<double> zero=0;
    line.initZeros(ZSIZE(mV));
    for (size_t i=0; i<YSIZE(mV); i++)
    {
        for (size_t j=0; j<XSIZE(mV); j++)
        {
            // Get line
            for (size_t k=0; k<ZSIZE(mV); k++)
                DIRECT_A1D_ELEM(line,k)=DIRECT_A3D_ELEM(mV,k,i,j);

            // Fourier transform
            transformer.setReal(line);
            transformer.FourierTransform();

            // Filter
            transformer.getFourierAlias(lineFourier);
            for (size_t k=maxPixel; k<XSIZE(lineFourier); k++)
            	DIRECT_A1D_ELEM(lineFourier,k)=zero;

            // Inverse Fourier transform and back to the volume
            transformer.inverseFourierTransform();
            for (size_t k=0; k<ZSIZE(mV); k++)
            	DIRECT_A3D_ELEM(mV,k,i,j)=DIRECT_A1D_ELEM(line,k);
        }
        if (verbose>=1)
        	progress_bar(i+1);
    }
    if (verbose>=1)
    	progress_bar(YSIZE(mV));

    // Write the results
    mV.setDimensions(XSIZE(mV),YSIZE(mV),1,ZSIZE(mV));
    V.write(fnOut,ALL_IMAGES,true);
}
