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

#include "phantom_simulate_microscope.h"

#include <data/args.h>

/* Read parameters --------------------------------------------------------- */
void ProgSimulateMicroscope::readParams()
{
    XmippMetadataProgram::readParams();

    fn_ctf = "";//initialize empty, force recalculation of first time
    pmdIn = getInputMd();
    CTFpresent=true;
    if (checkParam("--ctf"))
    {
        //Fill the input metadata with the value of 'fn_ctf'
    	fn_ctf=getParam("--ctf");
        MDConstGenerator generator(fn_ctf);
        generator.label = MDL_CTF_MODEL;
        generator.fill(*pmdIn);
    }
    else
    {
        if (pmdIn->containsLabel(MDL_CTF_MODEL))
        {
            //sort the images according to the ctf to avoid the recaculation of it
            //beeten images of the same ctf group
            MetaData md(*pmdIn);
            pmdIn->sort(md, MDL_CTF_MODEL);
        }
        else
        	CTFpresent=false;
    }
    after_ctf_noise = checkParam("--after_ctf_noise");
    defocus_change = getDoubleParam("--defocus_change");
    if (checkParam("--noise"))
    {
        sigma = getDoubleParam("--noise");
        low_pass_before_CTF = getDoubleParam("--noise",1);
        estimateSNR=false;
    }
    else if (checkParam("--noNoise"))
    {
        estimateSNR=0;
        sigma=0;
        low_pass_before_CTF=0.5;
    }
    else
    {
        targetSNR = getDoubleParam("--targetSNR");
        low_pass_before_CTF=0.5;
        estimateSNR=true;
    }
    downsampling = getDoubleParam("--downsampling");

}

/* Usage ------------------------------------------------------------------- */
void ProgSimulateMicroscope::defineParams()
{
    each_image_produces_an_output = true;
    save_metadata_stack = true;
    XmippMetadataProgram::defineParams();
    addUsageLine("Simulate the effect of the microscope on ideal projections.");
    addParamsLine("==CTF options==");
    addParamsLine(" [--ctf <CTFdescr=\"\">]       : a CTF description, if this param is not supplied it should come in metadata");
    addParamsLine(" [--after_ctf_noise]      : apply noise after CTF");
    addParamsLine(" [--defocus_change <v=0>] : change in the defocus value (percentage)");
    addParamsLine("==Noise options==");
    addParamsLine(" --noise <stddev> <w=0.5> : noise to be added, this noise is filtered at the frequency specified (<0.5).");
    addParamsLine("or --targetSNR <snr>      : the necessary noise power for a specified SNR is estimated");
    addParamsLine("or --noNoise              : do not add any noise, only simulate the CTF");
    addParamsLine(" [--downsampling <D=1>]  : Downsampling factor of the input micrograph with respect to the original");
    addParamsLine("                         : micrograph.");
    addExampleLine("Generate a set of images with the CTF applied without any noise", false);
    addExampleLine("   xmipp_phantom_simulate_microscope -i g0ta.sel --oroot g1ta --ctf untilt_ARMAavg.ctfparam");
    addExampleLine("Generate a set of images with a target SNR", false);
    addExampleLine("   xmipp_phantom_simulate_microscope -i g0ta.sel --oroot g2ta --ctf untilt_ARMAavg.ctfparam --targetSNR 0.2 --after_ctf_noise");
    addExampleLine("Generate a set of images with the CTF applied and noise before and after CTF", false);
    addExampleLine("   xmipp_phantom_simulate_microscope -i g0ta.sel --oroot g2ta --ctf untilt_ARMAavg.ctfparam --noise 4.15773 --after_ctf_noise");
}

/* Show -------------------------------------------------------------------- */
void ProgSimulateMicroscope::show()
{
    if (!verbose)
        return;
    XmippMetadataProgram::show();
    std::cout
    << "Noise: " << sigma << std::endl
    << "Low pass freq: " << low_pass_before_CTF << std::endl
    << "After CTF noise: " << after_ctf_noise << std::endl
    << "Defocus change: " << defocus_change << std::endl
    ;
    if (estimateSNR)
        std::cout
        << "Target SNR: " << targetSNR << std::endl;
    if (CTFpresent)
    	std::cout << "CTF file: " << fn_ctf << std::endl;
    else
    	std::cout << "No CTF provided\n";
}

void ProgSimulateMicroscope::estimateSigma()
{
	if (CTFpresent)
	{
		ctf.FilterBand = CTF;
		ctf.ctf.enable_CTFnoise = false;
		ctf.ctf.read(getParam("--ctf"));
		ctf.ctf.changeSamplingRate(ctf.ctf.Tm*downsampling);
		ctf.ctf.produceSideInfo();
	}

    size_t N_stats = pmdIn->size();
    std::cout << "N_stats=" << N_stats << std::endl;

    MultidimArray<double> proj_power(N_stats);
    MultidimArray<double> proj_area(N_stats);
    double power_avg, power_stddev, area_avg, area_stddev, avg, dummy;
    if (verbose!=0)
    {
        std::cerr << "Estimating noise power for target SNR=" << targetSNR << std::endl;
        init_progress_bar(N_stats);
    }
    FileName fnImg;
    Image<double> proj;
    size_t nImg=0;
    FOR_ALL_OBJECTS_IN_METADATA(*pmdIn)
    {
        pmdIn->getValue(image_label, fnImg,__iter.objId);
        proj.read(fnImg);
        MultidimArray<double> mProj=proj();

        if (CTFpresent)
        {
			if (nImg == 0)
				ctf.generateMask(mProj);
			ctf.applyMaskSpace(mProj);
        }

        // Compute projection area
        DIRECT_A1D_ELEM(proj_area,nImg) = 0;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mProj)
        if (fabs(DIRECT_MULTIDIM_ELEM(mProj,n)) > 1e-6)
            DIRECT_A1D_ELEM(proj_area,nImg)++;

        // Compute projection power
        mProj.computeStats(avg, DIRECT_A1D_ELEM(proj_power,nImg), dummy, dummy);

        if (nImg++ % 30 == 0 && verbose!=0)
            progress_bar(nImg);
    }
    progress_bar(N_stats);
    Histogram1D hist_proj, hist_area;
    compute_hist(proj_power, hist_proj, 300);
    compute_hist(proj_area, hist_area, 300);
    proj_power.computeStats(power_avg, power_stddev, dummy, dummy);
    std::cout << "# Projection power average: " << power_avg << std::endl
    << "# Projection power stddev:  " << power_stddev << std::endl
    << "# Projection percentil  2.5%: " << hist_proj.percentil(2.5) << std::endl
    << "# Projection percentil 97.5%: " << hist_proj.percentil(97.5) << std::endl;
    proj_area.computeStats(area_avg, area_stddev, dummy, dummy);
    std::cout << "# Projection area average: " << area_avg << std::endl
    << "# Projection area stddev:  " << area_stddev << std::endl
    << "# Area percentil  2.5%:    " << hist_area.percentil(2.5) << std::endl
    << "# Area percentil 97.5%:    " << hist_area.percentil(97.5) << std::endl;

    sigma=sqrt(power_avg*power_avg*Xdim*Ydim / (targetSNR*area_avg));
    std::cout << "Estimated sigma=" << sigma << std::endl;
    updateCtfs();
}

void ProgSimulateMicroscope::setupFourierFilter(FourierFilter &filter, bool isBackground,
        double &power)
{
    static int dXdim = 2 * Xdim, dYdim = 2 * Ydim;
    static MultidimArray<double> aux;

    filter.FilterBand = CTF;
    filter.ctf.read(fn_ctf);
    filter.ctf.enable_CTF = !isBackground;
    filter.ctf.enable_CTFnoise = isBackground;
    filter.ctf.produceSideInfo();
    aux.resizeNoCopy(dYdim, dXdim);
    aux.setXmippOrigin();
    filter.do_generate_3dmask=true;
    filter.generateMask(aux);
    power = filter.maskPower();
}

void ProgSimulateMicroscope::updateCtfs()
{
    double before_power = 0, after_power = 0;
    sigma_after_CTF = sigma_before_CTF = 0;

    if (CTFpresent)
    {
    	setupFourierFilter(ctf, false, before_power);
		if (after_ctf_noise)
			setupFourierFilter(after_ctf, true, after_power);
    }

    // Compute noise balance
    if ((after_power != 0 || before_power != 0) && sigma!=0)
    {
        double p = after_power / (after_power + before_power);
        double K = 1 / sqrt(p * after_power + (1 - p) * before_power);
        sigma_after_CTF = sqrt(p) * K * sigma;
        sigma_before_CTF = sqrt(1 - p) * K * sigma;
    }
    else if (sigma != 0)
    {
        sigma_before_CTF = sigma;
        sigma_after_CTF = 0;
    }
    else
        sigma_before_CTF=sigma_after_CTF=0;
}

/* Produce side information ------------------------------------------------ */
void ProgSimulateMicroscope::preProcess()
{
    size_t dum, dum2;
    getImageSize(*pmdIn, Xdim, Ydim, dum, dum2);

    if (low_pass_before_CTF < 0.5)
    {
        lowpass.FilterBand = LOWPASS;
        lowpass.FilterShape = RAISED_COSINE;
        lowpass.w1 = low_pass_before_CTF;
    }

    if (estimateSNR)
        estimateSigma();
}

void ProgSimulateMicroscope::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    static Image<double> img;
    static FileName last_ctf;
    static bool firstImage=true;
    last_ctf = fn_ctf;
    img.readApplyGeo(fnImg, rowIn);

    rowIn.getValue(MDL_CTF_MODEL, fn_ctf);
    rowOut.setValue(MDL_CTF_MODEL, fn_ctf);
    if (fn_ctf != last_ctf || firstImage)
    {
        updateCtfs();
        firstImage=false;
    }

    if (ZSIZE(img())!=1)
        REPORT_ERROR(ERR_MULTIDIM_DIM,"This process is not intended for volumes");

    apply(img());

    img.write(fnImgOut);
}

/* Apply ------------------------------------------------------------------- */
void ProgSimulateMicroscope::apply(MultidimArray<double> &I)
{
    I.setXmippOrigin();
    I.selfWindow(FIRST_XMIPP_INDEX(2*Ydim), FIRST_XMIPP_INDEX(2*Xdim),
                 LAST_XMIPP_INDEX(2*Ydim), LAST_XMIPP_INDEX(2*Xdim));

    // Add noise before CTF
    MultidimArray<double> noisy;
    noisy.resize(I);
    noisy.initRandom(0, sigma_before_CTF, RND_GAUSSIAN);
    if (low_pass_before_CTF < 0.5)
        lowpass.applyMaskSpace(noisy);
    I += noisy;

    if (CTFpresent)
    {
		// Check if the mask is a defocus changing CTF
		// In that case generate a new mask with a random defocus
		if (defocus_change != 0)
		{
			MultidimArray<double> aux;
			ctf.ctf.DeltafU *= rnd_unif(1 - defocus_change / 100, 1 + defocus_change / 100);
			ctf.ctf.DeltafV *= rnd_unif(1 - defocus_change / 100, 1 + defocus_change / 100);
			aux.initZeros(2*Ydim, 2*Xdim);
			ctf.generateMask(aux);
		}

		// Apply CTF
		ctf.applyMaskSpace(I);

		// Add noise after CTF
		noisy.initRandom(0, sigma_after_CTF, RND_GAUSSIAN);
		if (after_ctf_noise)
			after_ctf.applyMaskSpace(noisy);
		I += noisy;
    }
    I.selfWindow(FIRST_XMIPP_INDEX(Ydim), FIRST_XMIPP_INDEX(Xdim),
                 LAST_XMIPP_INDEX(Ydim), LAST_XMIPP_INDEX(Xdim));
}
