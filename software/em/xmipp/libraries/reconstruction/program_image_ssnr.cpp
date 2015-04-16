/***************************************************************************
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
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


#include "program_image_ssnr.h"
#include <data/mask.h>
#include <data/xmipp_fftw.h>

void ProgImageSSNR::defineParams()
{
	produces_a_metadata = true;
    each_image_produces_an_output = false;
    keep_input_columns = true;
    remove_disabled = false;
    addUsageLine("Analyze image SSNR. For doing so, we estimate the SSNR of each particle defining as signal "\
    		     "the part of the image within a certain radius, and as noise the part of the image outside.");
    XmippMetadataProgram::defineParams();
    addParamsLine(" [-R <r=-1>]: Particle radius, by default, half size of the image");
    addParamsLine(" [--Rwidth <r=3>]: Transition radius the mask is 1 till R-Rwidth, and zero from R+Rwidth");
    addParamsLine(" [--fmin <f=40>]: Minimum frequency (in Angstroms) to measure");
    addParamsLine(" [--fmax <f=3>]:  Maximum frequency (in Angstroms) to measure");
    addParamsLine(" [--sampling <Ts=1>]:  Sampling rate in Angstroms/pixel");
    addParamsLine(" [--ssnrcut <ssnr=-1>]:  Disable images whose SSNR is below this value");
    addParamsLine(" [--ssnrpercent <p=-1>]:  Disable images whose SSNR is below this percentage");
    addParamsLine(" [--normalizessnr]: Normalize the SSNR by dividing by the maximum SSNR");
    addExampleLine("xmipp_image_ssnr -i images.xmd -o imagesOut.xmd ");
}

void ProgImageSSNR::readParams()
{
    XmippMetadataProgram::readParams();
    R=getDoubleParam("-R");
    Rwidth=getDoubleParam("--Rwidth");
    fmin=getDoubleParam("--fmin");
    fmax=getDoubleParam("--fmax");
    sampling=getDoubleParam("--sampling");
    ssnrcut=getDoubleParam("--ssnrcut");
    ssnrpercent=getDoubleParam("--ssnrpercent");
    normalizessnr=checkParam("--normalizessnr");
}

void ProgImageSSNR::preProcess()
{
	size_t Xdim, Ydim, Zdim, Ndim;
	getImageSize(*getInputMd(), Xdim, Ydim, Zdim, Ndim);
	maskS.initZeros(Ydim,Xdim);
	maskS.setXmippOrigin();

	if (R==-1)
		R=0.5*Xdim-Rwidth;
	RaisedCosineMask(maskS,R-Rwidth,R+Rwidth);
	maskN=1-maskS;

	imin=(size_t)std::max(3.0,0.5*Xdim*(sampling/fmin));
	imax=(size_t)std::min(Xdim-3.0,0.5*Xdim*(sampling/fmax));
}

void ProgImageSSNR::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
	// Copy Enabled
	int enabled;
	rowIn.getValue(MDL_ENABLED,enabled);
	rowOut.setValue(MDL_ENABLED,enabled);

    img.read(fnImg);
    img().setXmippOrigin();

    imgN()=imgS()=img();
    imgS()*=maskS;
    imgN()*=maskN;

    getSpectrum(imgS(),spectrumS,POWER_SPECTRUM);
    getSpectrum(imgN(),spectrumN,POWER_SPECTRUM);

    double SSNR=0;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(spectrumS)
    {
    	if (i>=imin && i<=imax)
    	{
    		double S=DIRECT_A1D_ELEM(spectrumS,i);
    		double N=DIRECT_A1D_ELEM(spectrumN,i);
    		if (S>0 && N>0)
    		SSNR+=log10(S)-log10(N);
    	}
    }
    SSNR*=10.0/(imax-imin+1); // decibels
    rowOut.setValue(MDL_CUMULATIVE_SSNR,SSNR);

    //spectrumS.write("PPPS.txt");
    //spectrumN.write("PPPN.txt");
    //std::cout << fnImg << " SSNR=" << SSNR << std::endl;
    //std::cout << "Press any key\n";
    // char c; std::cin>> c;
}

void thresholdSSNR(MetaData &mdOut, double ssnrcut)
{
	FOR_ALL_OBJECTS_IN_METADATA(mdOut)
	{
		double ssnr;
		mdOut.getValue(MDL_CUMULATIVE_SSNR,ssnr,__iter.objId);
		if (ssnr<ssnrcut)
			mdOut.setValue(MDL_ENABLED,-1,__iter.objId);
	}
}

void normalizeSSNR(MetaData &mdOut)
{
	double maxSSNR=-1e38;
	FOR_ALL_OBJECTS_IN_METADATA(mdOut)
	{
		double ssnr;
		mdOut.getValue(MDL_CUMULATIVE_SSNR,ssnr,__iter.objId);
		if (ssnr>maxSSNR)
			maxSSNR=ssnr;
	}

	if (maxSSNR>0)
	{
		double imaxSSNR=1/maxSSNR;
		FOR_ALL_OBJECTS_IN_METADATA(mdOut)
		{
			double ssnr;
			mdOut.getValue(MDL_CUMULATIVE_SSNR,ssnr,__iter.objId);
			mdOut.setValue(MDL_WEIGHT_SSNR,ssnr*imaxSSNR,__iter.objId);
		}
	}
}

void ProgImageSSNR::postProcess()
{
	if (ssnrcut>0)
		thresholdSSNR(*getOutputMd(),ssnrcut);

	if (ssnrpercent>0)
	{
		std::vector<double> ssnr;
		getOutputMd()->getColumnValues(MDL_CUMULATIVE_SSNR,ssnr);
		std::sort(ssnr.begin(),ssnr.end());
		size_t idx=(size_t)(ssnrpercent/100.0*ssnr.size());
		thresholdSSNR(*getOutputMd(),ssnr[idx]);
	}

	if (normalizessnr)
		normalizeSSNR(*getOutputMd());
	if (fn_out!="")
		getOutputMd()->write(fn_out);
	else
		getOutputMd()->write(fn_in);
}
