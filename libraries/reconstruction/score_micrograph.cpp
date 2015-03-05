/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
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

#include "data/args.h"
#include "data/filters.h"
#include "data/xmipp_fft.h"

#include "score_micrograph.h"
#include "fourier_filter.h"
#include "fringe_processing.h"

/* Read parameters --------------------------------------------------------- */
void ProgScoreMicrograph::readParams()
{

	pieceDim = 256;
	overlap = 0.5;
	skipBorders = 2;
	Nsubpiece = 1;
	bootstrapN = -1;
	estimate_ctf = true;
	defocus_range = 1000;
	fn_micrograph = getParam("--micrograph");
	particleSize = getIntParam("--particleSize");

	 prmEstimateCTFFromPSD.readBasicParams(this);

}

/* Usage ------------------------------------------------------------------- */
void ProgScoreMicrograph::defineParams()
{
	addUsageLine("Estimate the score of a given Micrograph attending on different parameters.");
    addParamsLine("   --micrograph <file>         : File with the micrograph");
    addParamsLine("  [--particleSize <p=100>]       : Size of the particle");
    ProgCTFEstimateFromPSD::defineBasicParams(this);
}

/* Apply ------------------------------------------------------------------- */
void ProgScoreMicrograph::run()
{

	std::cout << particleSize << std::endl;
	char bufferDW[200];
	FileName fn_micrographBadPxl;
	fn_micrographBadPxl = fn_micrograph.withoutExtension()+"bp__tmp."+fn_micrograph.getExtension();
	sprintf(bufferDW, "xmipp_transform_filter -i %s -o %s --bad_pixels outliers 2", fn_micrograph.c_str(), fn_micrographBadPxl.c_str());
	if (!system(bufferDW))
		REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot filter");

	FileName fn_micrographDwn;
	fn_micrographDwn = fn_micrograph.withoutExtension()+"_tmp."+fn_micrograph.getExtension();
	sprintf(bufferDW, "xmipp_transform_downsample -i %s -o %s --step 2.0 --method fourier" , 	fn_micrographBadPxl.c_str(), fn_micrographDwn.c_str());
	if (!system(bufferDW))
		REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot downsample");

	prmEstimateCTFFromPSD.defocus_range = defocus_range;
	prmEstimateCTFFromPSD.downsampleFactor = 2.0;
	prmEstimateCTFFromPSD.Tm = prmEstimateCTFFromPSD.Tm*prmEstimateCTFFromPSD.downsampleFactor;
	prmEstimateCTFFromMicrograph.prmEstimateCTFFromPSD = prmEstimateCTFFromPSD;
	prmEstimateCTFFromMicrograph.pieceDim = pieceDim;
	prmEstimateCTFFromMicrograph.overlap = overlap;
	prmEstimateCTFFromMicrograph.skipBorders = skipBorders;
	prmEstimateCTFFromMicrograph.Nsubpiece = 1;
	prmEstimateCTFFromMicrograph.bootstrapN = bootstrapN;
	prmEstimateCTFFromMicrograph.estimate_ctf = estimate_ctf;
	prmEstimateCTFFromMicrograph.fn_micrograph = fn_micrographDwn;
	prmEstimateCTFFromMicrograph.fn_root = "/home/jvargas/Linux/Proyectos/LUMC/set_001_challenge/kk";
	prmEstimateCTFFromMicrograph.run();

	MetaData MD;
	size_t id=MD.addObject();
	MD.setValue(MDL_MICROGRAPH,fn_micrographDwn,id);
	MD.setValue(MDL_PSD,prmEstimateCTFFromMicrograph.fn_root+".psd",id);
	MD.setValue(MDL_PSD_ENHANCED,prmEstimateCTFFromMicrograph.fn_root+"_enhanced_psd.xmp",id);
	MD.setValue(MDL_PSD_ENHANCED,prmEstimateCTFFromMicrograph.fn_root+"_enhanced_psd.xmp",id);
	MD.setValue(MDL_CTF_MODEL,prmEstimateCTFFromMicrograph.fn_root+".ctfparam",id);
	MD.setValue(MDL_IMAGE1,prmEstimateCTFFromMicrograph.fn_root+"_ctfmodel_quadrant.xmp",id);
	MD.setValue(MDL_IMAGE2,prmEstimateCTFFromMicrograph.fn_root+"_ctfmodel_halfplane.xmp",id);
	MD.setValue(MDL_CTF_DOWNSAMPLE_PERFORMED,2.000000,id);
	MD.write("auxiliaryFile.xmd");

	//prmPSDSort.downsampling =
	prmPSDSort.fn_in = "auxiliaryFile.xmd";
	prmPSDSort.fn_out = "kkauxiliaryFile.xmd";
	// COSS prmPSDSort.downsampling = 2.000000; // this parameter must be written in the metadata

	char buffer[200];
	sprintf(buffer, "xmipp_ctf_sort_psds -i %s -o %s --downsampling %f", 	prmPSDSort.fn_in.c_str(), prmPSDSort.fn_out.c_str(),2.000000);
	if (!system(bufferDW))
		REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot sort psds");

	//Spiral Phase Transform
    MultidimArray<double> im, In, Mod, Part, Ice;
    MultidimArray<bool> ROI1, ROI2;
    Image<double> img;
    img.read(fn_micrographDwn);
    im = img();
    im.statisticsAdjust(0,1);
    In.resizeNoCopy(im);
    Part.resizeNoCopy(im);
    Ice.resizeNoCopy(im);
    Mod.resizeNoCopy(im);
    ROI1.resizeNoCopy(im);
    ROI1.setXmippOrigin();
    ROI2.resizeNoCopy(im);
    ROI2.setXmippOrigin();
    Ice.setXmippOrigin();
    Part.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI1)
    	A2D_ELEM(ROI1,i,j)= true;

    //aprox  Xdim/particleSize is the particle frequency
    size_t Xdim, Ydim, Zdim,Ndim;
    im.getDimensions(Xdim, Ydim,Zdim,Ndim);

    double avg, stdv;
    im.computeAvgStdev(avg,stdv);
    double R = Xdim/particleSize;
    double S = Xdim/particleSize;
    FourierTransformer ftrans(FFTW_BACKWARD);
    normalize(ftrans,im,In,Mod,R,S,ROI1);

    double partDensity = (Mod.sum()/(Xdim*Ydim))*100;

    double th=EntropyOtsuSegmentation(Mod,0.005,false);
    //th *= 2;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI1)
    {
    	if (A2D_ELEM(Mod,i,j) < th)
    	{
    		A2D_ELEM(ROI1,i,j) = false;
    		A2D_ELEM(ROI2,i,j) = true;
    	}
    }

    double min_val, max_val, avg1, avg2, stdv1, stdv2;
    computeStats_within_binary_mask(ROI1, im, min_val, max_val,avg1,stdv1);
    computeStats_within_binary_mask(ROI2, im, min_val, max_val,avg2,stdv2);
    double snr = (stdv1/stdv2)*(stdv1/stdv2);


    double maxFreq;
    double fitting;
    double psdCorr90;
    double psdRadialInt;

    std::cout << prmPSDSort.fn_in.c_str() << std::endl;
    std::cout << prmPSDSort.fn_out.c_str() << std::endl;

    MD.read(prmPSDSort.fn_out);
    FOR_ALL_OBJECTS_IN_METADATA(MD)
    {
    	MD.getValue(MDL_CTF_CRIT_MAXFREQ,maxFreq,__iter.objId);
    	MD.getValue(MDL_CTF_CRIT_FITTINGSCORE,fitting,__iter.objId);
    	MD.getValue(MDL_CTF_CRIT_PSDCORRELATION90,psdCorr90,__iter.objId);
    	MD.getValue(MDL_CTF_CRIT_PSDRADIALINTEGRAL,psdRadialInt,__iter.objId);
    }

    std::cout << "partDensity : " << partDensity << std::endl;
    std::cout << "snr :" << snr << std::endl;
    std::cout << "maxFreq :" << maxFreq << std::endl;
    std::cout << "fitting :" << fitting << std::endl;
    std::cout << "psdCorr90 :" << psdCorr90 << std::endl;
    std::cout << "psdRadialInt :" << psdRadialInt << std::endl;

    img()=Part;
    img.write("kk1.mrc");
    img()=Ice;
    img.write("kk2.mrc");

}

#undef DEBUG
