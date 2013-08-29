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

#include <data/xmipp_program.h>
#include <data/xmipp_fftw.h>

class ProgVolumeStructureFactor: public XmippProgram
{
protected:
	FileName fnVol, fnStructure;
	double sampling;

    void defineParams()
    {
        //Usage
        addUsageLine("Calculate the structure factor and the Guinier plot of a volume");
        //Parameters
        addParamsLine("-i <volume>           : Volume to analyze");
        addParamsLine("[-o <structure=\"\">] : Metadata with the structure factor and the Guinier plot");
        addParamsLine("                      : If no name is given, then volume_structure.xmd");
        addParamsLine("[--sampling <T=1>]    : Sampling rate in Angstroms/pixel");
    }

    void readParams()
    {
    	fnVol=getParam("-i");
    	fnStructure=getParam("-o");
    	if (fnStructure=="")
    		fnStructure=fnVol.removeAllExtensions()+"_structure.xmd";
    	sampling=getDoubleParam("--sampling");
    }

    void show()
    {
    	std::cout
    	<< "Input volume:     " << fnVol       << std::endl
    	<< "Output structure: " << fnStructure << std::endl
    	<< "Sampling:         " << sampling    << std::endl
    	;
    }

    void run()
    {
    	show();

    	Image<double> V;
    	FourierTransformer transformer;
    	MultidimArray< std::complex<double> > VFourier;
    	MultidimArray<double> VFourierMag, spectrum;

    	V.read(fnVol);
    	transformer.FourierTransform(V(),VFourier,false);
    	FFT_magnitude(VFourier,VFourierMag);
    	getSpectrum(VFourierMag,spectrum,POWER_SPECTRUM);

    	MetaData MDout;
    	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(spectrum)
    	{
    		if (i==0 || A1D_ELEM(spectrum,i)==0)
    			continue;
    		size_t id=MDout.addObject();
    		double f;
    		FFT_IDX2DIGFREQ(i,XSIZE(V()),f);
    		if (f>0.5)
    			continue;
    		f/=sampling;
    		double fA=1/f;
    		double f2=f*f;

    		MDout.setValue(MDL_RESOLUTION_FREQ,f,id);
    		MDout.setValue(MDL_RESOLUTION_FREQREAL,fA,id);
    		MDout.setValue(MDL_RESOLUTION_FREQ2,f2,id);
    		MDout.setValue(MDL_RESOLUTION_STRUCTURE_FACTOR,A1D_ELEM(spectrum,i),id);
    		double aux=0;
    		if (A1D_ELEM(spectrum,i)>0)
    			aux=log(A1D_ELEM(spectrum,i));
    		MDout.setValue(MDL_RESOLUTION_LOG_STRUCTURE_FACTOR,aux,id);
    	}
    	MDout.write(fnStructure);
    }
};

RUN_XMIPP_PROGRAM(ProgVolumeStructureFactor)
