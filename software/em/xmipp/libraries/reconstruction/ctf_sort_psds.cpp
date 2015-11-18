/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Slavica Jonic (slavica.jonic@impmc.jussieu.fr)
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

#include <vector>
#include "ctf_sort_psds.h"
#include "ctf_enhance_psd.h"
#include "ctf_estimate_from_micrograph.h"
#include <data/args.h>
#include <data/filters.h>
#include <data/transformations.h>
#include <data/histogram.h>


/* Constructor ------------------------------------------------------------- */
ProgPSDSort::ProgPSDSort()
{
	produces_an_output=true;
	produces_a_metadata=true;
	keep_input_columns=true;
}

void ProgPSDSort::defineLabelParam()
{
	addParamsLine(" [--label+ <image_label=micrograph>]   : Label to be used to read/write images.");
}

/* Read parameters --------------------------------------------------------- */
void ProgPSDSort::readParams()
{
	XmippMetadataProgram::readParams();
    filter_w1 = getDoubleParam("-f1");
    filter_w2 = getDoubleParam("-f2");
    decay_width = getDoubleParam("-decay");
    mask_w1 = getDoubleParam("-m1");
    mask_w2 = getDoubleParam("-m2");
    // COSS downsampling = getDoubleParam("--downsampling");
}

/* Usage ------------------------------------------------------------------- */
void ProgPSDSort::defineParams()
{
    XmippMetadataProgram::defineParams();
    addUsageLine("Evaluate the CTFs and PSDs of a set of micrographs.");
    addUsageLine("This process is strongly coupled to the output produced by the preprocessing micrographs step of the Xmipp protocols. ");
    addUsageLine("For each input PSD, the program writes its enhanced version since it is used in the computation of some of the criteria.");
    addUsageLine("+The different criteria for evaluating the PSDs are:");
    addUsageLine("+ ");
    addUsageLine("+$ *Damping*: this is the envelope value at the border of the PSD. Micrographs ");
    addUsageLine("+with a high envelope value at border are either wrongly estimated strongly undersampled.");
    addUsageLine("+ ");
    addUsageLine("+$ *First zero average*: this is average in Angstroms of the first zero. ");
    addUsageLine("+Normally, this value should be between 4x and 10x the sampling rate in Angstroms.");
    addUsageLine("+ ");
    addUsageLine("+$ *Maximum frequency*: this is the resolution (in Angstroms) at which the envelope drops below 1% of the maximum envelope");
    addUsageLine("+ ");
    addUsageLine("+$ *First zero disagreement*: if the CTF has been estimated by two different methods ");
    addUsageLine("+(normally Xmipp and Ctffind), then this criterion measures the average disagreement ");
    addUsageLine("+in Angstroms between the first zero in the two estimates. Low disagreements are ");
    addUsageLine("+indicative of correct fit.");
    addUsageLine("+ ");
    addUsageLine("+$ *First zero ratio*: this measures the astigmatism of the CTF by computing the ratio ");
    addUsageLine("+between the largest and smallest axes of the first zero ellipse. Ratios close to 1 ");
    addUsageLine("+indicate no astigmatism.");
    addUsageLine("+ ");
    addUsageLine("+$ *Ratio between the standard deviation at 1st zero and 1st minimum*: the variance in the experimental PSD along the");
    addUsageLine("+first zero and the first CTF minimum should be approximately equal (ratio=1).");
    addUsageLine("+ ");
    addUsageLine("+$ *CTF margin*: ratio between the average difference in the experimental PSD between the 1st Thon");
    addUsageLine("+ring and its previous zero, and the variance of the experimental PSD along the first zero");
    addUsageLine("+first zero and the first CTF minimum should be approximately equal (ratio=1).");
    addUsageLine("+ ");
    addUsageLine("+$ *Fitting score*: the CTF is computed by fitting a theoretical model to the experimentally observed PSD. ");
    addUsageLine("+This criterion is the fitting score. Smaller scores correspond to better fits.");
    addUsageLine("+ ");
    addUsageLine("+$ *Fitting correlation between zeros 1 and 3*: the region between the first and third zeroes ");
    addUsageLine("+is particularly important since it is where the Thon rings are most visible. ");
    addUsageLine("+This criterion reports the correlation between the experimental and theoretical PSDs ");
    addUsageLine("+within this region. High correlations indicate good fits.");
    addUsageLine("+ ");
    addUsageLine("+$ *Non-astigmatic validity*: if we consider the CTF cosine part in the direction U and V ");
    addUsageLine("+and add both as if they were waves, this criterion shows the frequency (in Angstroms) at which ");
    addUsageLine("+both waves would interfere completely destructively. Beyond this frequency, it cannot be assumed that ");
    addUsageLine("+a non-astigmatic CTF correction can manage an astigmatic CTF");
    addUsageLine("+ ");
    addUsageLine("+$ *PSD correlation at 90 degrees*: The PSD of non-astigmatic micrographs correlate well ");
    addUsageLine("+with itself after rotating the micrograph 90 degrees. This is so because non-astigmatic ");
    addUsageLine("+PSDs are circularly symmetrical, while astigmatic micrographs are elliptically symmetrical.");
    addUsageLine("+High correlation when rotating 90 degrees is an indicator of non-astigmatism.");
    addUsageLine("+This criterion is computed on the enhanced PSD. See [[ctf_enhance_psd_v3][ctf_enhance_psd]].");
    addUsageLine("+ ");
    addUsageLine("+$ *PSD radial integral*: this criterion reports the integral of the radially symmetrized PSD.");
    addUsageLine("+This criterion can highlight differences among the background noises of micrographs. ");
    addUsageLine("+This criterion is computed on the enhanced PSD. See [[ctf_enhance_psd_v3][ctf_enhance_psd]].");
    addUsageLine("+ ");
    addUsageLine("+$ *PSD variance*: the PSD is estimated by averaging different PSD local estimates in small regions of the micrograph. ");
    addUsageLine("+This criterion measures the variance of the different PSD local estimates. Untilted micrographs ");
    addUsageLine("+have equal defoci all over the micrograph, and therefore, the variance is due only to noise. ");
    addUsageLine("+However, tilted micrographs have an increased PSD variance since different regions of the micrograph ");
    addUsageLine("+have different defoci. Low variance of the PSD are indicative of non-tilted micrographs");
    addUsageLine("+ ");
    addUsageLine("+$ *PSD Principal Component 1 Variance*: when considering the local PSDs previously defined as vectors ");
    addUsageLine("+in a multidimensional space, we can compute the variance of their projection onto the first principal component axis. ");
    addUsageLine("+Low variance of this projection is indicative of a uniformity of local PSDs, i.e., this is another measure ");
    addUsageLine("+of the presence of tilt in the micrograph.");
    addUsageLine("+ ");
    addUsageLine("+$ *PSD !PCA Runs test*: when computing the projections onto the first principal component, as discussed in the previous criterion, ");
    addUsageLine("+one might expect that the sign of the projection is random for untilted micrographs. Micrographs with a marked ");
    addUsageLine("+non-random pattern of projections are indicative of tilted micrographs. The larger the value of this criterion, the less random the pattern is.");
    addParamsLine("==+ Enhancement filter parameters");
    addParamsLine("  [-f1 <freq_low=0.02>]      : Low freq. for band pass filtration, max 0.5");
    addParamsLine("  [-f2 <freq_high=0.2>]      : High freq. for band pass filtration, max 0.5");
    addParamsLine("  [-decay <freq_decay=0.02>] : Decay for the transition bands");
    addParamsLine("  [-m1 <mfreq_low=0.01>]     : Low freq. for mask, max 0.5");
    addParamsLine("  [-m2 <mfreq_high=0.45>]    : High freq. for mask, max 0.5");
    // COSS addParamsLine("  [--downsampling <K=1>]     : Downsampling factor used");
}

/* Show -------------------------------------------------------------------- */
void ProgPSDSort::show()
{
	if (verbose==0)
		return;
    XmippMetadataProgram::show();
    std::cout
    << "Filter w1:    " << filter_w1 << std::endl
    << "Filter w2:    " << filter_w2 << std::endl
    << "Filter decay: " << decay_width << std::endl
    << "Mask w1:      " << mask_w1 << std::endl
    << "Mask w2:      " << mask_w2 << std::endl ;
    // COSS << "Downsampling: " << downsampling << std::endl;
}

/* Compute Correlation ----------------------------------------------------- */
void ProgPSDSort::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    CTFDescription CTF1, CTF2;
	PSDEvaluation evaluation;
    FileName fnMicrograph, fnPSD, fnCTF, fnCTF2;

	evaluation.ctf_envelope_ssnr = fn_in.getDir()+"envelope.xmd";

	int enabled;
	rowIn.getValue(MDL_ENABLED,enabled);
    rowIn.getValue(MDL_MICROGRAPH, fnMicrograph);
    rowIn.getValue(MDL_PSD,fnPSD);
	if (enabled==-1 || fnPSD == "NA")
	{
		rowOut.setValue(MDL_ENABLED,-1);
		return;
	}

    if (rowIn.containsLabel(MDL_CTF_MODEL))
    {
      rowIn.getValue(MDL_CTF_MODEL, fnCTF);
      if (!fnCTF.exists())
      {
    	  rowOut.setValue(MDL_ENABLED,-1);
    	  return;
      }
      CTF1.read(fnCTF);
    }
    else
      CTF1.readFromMdRow(rowIn);

    FileName fnRoot = fnMicrograph.withoutExtension();

    CTF1.produceSideInfo();
    evaluation.defocusU=CTF1.DeltafU;
    evaluation.defocusV=CTF1.DeltafV;

    if (rowIn.containsLabel(MDL_CTF_MODEL2))
    	rowIn.getValue(MDL_CTF_MODEL2,fnCTF2);

    if (!fnCTF2.empty() && fnCTF2 != "NA")
    {
    	CTF2.read(fnCTF2);
    	CTF2.produceSideInfo();
    }

    // Evaluate beating due to astigmatism
    /* The argument of the cosine part of the CTF is:                     wu=cos(K1*deltafu*u.^2+K2*u.^4)
     * If there is no astigmatism, then in the V direction we would have  wv=cos(K1*deltafv*u.^2+K2*u.^4)
     * and we should have (wu+wv)/2 = wu.
     * If there is astigmatism, the sum (wu+wv)/2 will depart from the behavior of wu
     * http://en.wikipedia.org/wiki/Beat_%28acoustics%29
     * calling argu and argv the argument of the two waves, we have
     *
     * (wu+wv)/2=cos((argu+argv)/2)cos((argu-argv)/2)
     *
     * The term cos((argu-argv)/2) acts as an envelope which may even vanish. Let's analyze this envelope
     * cos(0.5*(K1*deltafu*u.^2+K2*u.^4-K1*deltafv*u.^2+K2*u.^4))=cos(0.5*K1*abs(deltafu-deltafv)*u^2)
     *
     * When this envelope is 0, the interference between the two waves is totally destructive (we cannot apply a
     * non-astigmatic correction to an astigmatic CTF). This happens for
     *
     * 0.5*K1*abs(deltafu-deltafv)*u_0^2=pi/2 ---> u_0=sqrt(PI/(K1*abs(deltafu-deltafv)))
     *
     * This is the expression of critBeating
     */
    evaluation.beating=1.0/sqrt(PI/(CTF1.K1*abs(CTF1.DeltafU-CTF1.DeltafV)));

    // Read input PSD data
    Image<double> PSD;
    PSD.read(fnPSD);
    // Enhance the PSD
    ProgCTFEnhancePSD enhancePSD;
    enhancePSD.filter_w1 = filter_w1;
    enhancePSD.filter_w2 = filter_w2;
    enhancePSD.decay_width = decay_width;
    enhancePSD.mask_w1 = mask_w1;
    enhancePSD.mask_w2 = mask_w2;
    enhancePSD.applyFilter(PSD());

    // Evaluate the radial integral
    PSD().setXmippOrigin();
    Matrix1D< int > center_of_rot(2);
    MultidimArray< double > radial_mean;
    MultidimArray<int> radial_count;
    radialAverage(PSD(),center_of_rot,radial_mean,radial_count);
    radial_mean.selfABS();
    radial_mean/=radial_mean.computeMax();
    evaluation.PSDradialIntegral=radial_mean.sum();

    // Rotate 90 degrees and compute correlation
    Image<double> PSDrotated;
    rotate(LINEAR,PSDrotated(),PSD(),90);
    evaluation.PSDcorrelation90=correlationIndex(PSD(), PSDrotated());

    // Get the fitting score and other quality criteria computed by ctf_estimate_from_micrograph
    MetaData MDctf1;
    MDctf1.read(fnCTF);
    size_t objId1 = MDctf1.firstObject();

#define GET_CTF_CRITERION(labelll,xxx) \
    if (rowIn.containsLabel(labelll)) \
    	rowIn.getValue(labelll,xxx); \
    else if (MDctf1.containsLabel(labelll)) \
    	MDctf1.getValue(labelll,xxx,objId1); \
    else \
    	xxx=0;
    GET_CTF_CRITERION(MDL_CTF_CRIT_FITTINGSCORE,evaluation.fittingScore);
    GET_CTF_CRITERION(MDL_CTF_CRIT_FITTINGCORR13,evaluation.fittingCorr13);
    GET_CTF_CRITERION(MDL_CTF_CRIT_PSDVARIANCE,evaluation.PSDVariance);
    GET_CTF_CRITERION(MDL_CTF_CRIT_PSDPCA1VARIANCE,evaluation.PSDPC1Variance);
    GET_CTF_CRITERION(MDL_CTF_CRIT_PSDPCARUNSTEST,evaluation.PSDPCRunsTest);

    // Explore the CTF
    Matrix1D<double> u(2), freqZero1(2), freqZero2(2), freqMin1(2), pixelZero1(2), pixelMin1(2);
    double wmax=0.5/CTF1.Tm;
    double maxModuleZero=0, minModuleZero=1e38;
    double N=0;
    evaluation.maxDampingAtBorder=0;
    evaluation.firstZeroDisagreement=-1;
    evaluation.firstZeroAvg=0;
    double firstZeroAvgPSD=0;
    double firstZeroStddevPSD=0;
    double firstMinAvgPSD=0;
    double firstMinStddevPSD=0;
    double firstZeroMinAvgPSD=0;
    double firstZeroMinStddevPSD=0;
    evaluation.maxFreq=1000;

    CTF1.precomputeValues(0.0,0.0);
	double idamping0=1.0/CTF1.getValueDampingAt();
	double f2pixel=CTF1.Tm*XSIZE(PSD()); // COSS *downsampling
	if (rowIn.containsLabel(MDL_CTF_DOWNSAMPLE_PERFORMED))
	{
		double aux;
		rowIn.getValue(MDL_CTF_DOWNSAMPLE_PERFORMED,aux);
		f2pixel*=aux;
	}

	MetaData mdEnvelope;
	Matrix1D< double > envelope(100);
	envelope.initZeros();

	double Nalpha = 180;
    for (double alpha=0; alpha<=PI; alpha+=PI/Nalpha, N++)
    {
    	VECTOR_R2(u,cos(alpha),sin(alpha));

    	// Get the zero in the direction of u
    	CTF1.lookFor(1, u, freqZero1, 0);
    	double moduleZero=1.0/freqZero1.module();
    	maxModuleZero=XMIPP_MAX(maxModuleZero,moduleZero);
    	minModuleZero=XMIPP_MIN(minModuleZero,moduleZero);
    	evaluation.firstZeroAvg+=moduleZero;

    	// Get the first minimum (it is at higher frequency than the zero)
    	CTF1.lookFor(1, u, freqMin1, -1);

    	pixelZero1=freqZero1*f2pixel;
    	pixelMin1=freqMin1*f2pixel;
    	double psdZero=PSD().interpolatedElement2D(XX(pixelZero1),YY(pixelZero1),0.0);
    	double psdMin=PSD().interpolatedElement2D(XX(pixelMin1),YY(pixelMin1),0.0);
    	firstMinAvgPSD+=psdMin;
    	firstMinStddevPSD+=psdMin*psdMin;
    	firstZeroAvgPSD+=psdZero;
    	firstZeroStddevPSD+=psdZero*psdZero;
    	double zeroMinDiff=psdMin-psdZero;
    	firstZeroMinAvgPSD+=zeroMinDiff;
    	firstZeroMinStddevPSD+=zeroMinDiff*zeroMinDiff;

    	// Evaluate damping
    	double wx=wmax*XX(u);
    	double wy=wmax*YY(u);
    	CTF1.precomputeValues(wx,wy);
    	double damping=CTF1.getValueDampingAt();
    	damping=damping*damping;
    	evaluation.maxDampingAtBorder=XMIPP_MAX(evaluation.maxDampingAtBorder,damping);

        int idx = 0;
    	for (double w=0; w<wmax; w+=wmax/100.0)
    	{
        	wx=w*XX(u);
        	wy=w*YY(u);
        	CTF1.precomputeValues(wx,wy);
        	double normalizedDamping=fabs(CTF1.getValueDampingAt()*idamping0);
        	if (normalizedDamping>0.1)
        		evaluation.maxFreq=std::min(evaluation.maxFreq,1.0/w);

        	VEC_ELEM(envelope,idx) += double(fabs(CTF1.getValueDampingAt()));
    		idx++;
    	}

    	if (fnCTF2!="") {
        	CTF2.lookFor(1, u, freqZero2, 0);
        	double module2=1.0/freqZero2.module();
        	double diff=ABS(moduleZero-module2);
        	evaluation.firstZeroDisagreement=XMIPP_MAX(evaluation.firstZeroDisagreement,diff);
    	}
    }

    size_t objId2 = mdEnvelope.firstObject();
    int idx=0;
	for (double w=0; w<wmax; w+=wmax/100.0)
	{
		mdEnvelope.setValue(MDL_RESOLUTION_FREQ,w,objId2);
		mdEnvelope.setValue(MDL_CTF_ENVELOPE,VEC_ELEM(envelope,idx)/Nalpha,objId2);
		objId2 = mdEnvelope.addObject();
		idx++;
	}

    evaluation.firstZeroAvg/=N;
    evaluation.firstZeroRatio=maxModuleZero/minModuleZero;
    firstZeroAvgPSD/=N;
    firstZeroStddevPSD=sqrt(fabs(firstZeroStddevPSD/N-firstZeroAvgPSD*firstZeroAvgPSD));
    firstMinAvgPSD/=N;
    firstMinStddevPSD=sqrt(fabs(firstMinStddevPSD/N-firstMinAvgPSD*firstMinAvgPSD));
    firstZeroMinAvgPSD/=N;
    firstZeroMinStddevPSD=sqrt(fabs(firstZeroMinStddevPSD/N-firstZeroMinAvgPSD*firstZeroMinAvgPSD));

    evaluation.firstMinimumStddev_ZeroStddev=1000;
    evaluation.firstMinimumDiffStddev_ZeroStddev=1000;
    if (firstZeroStddevPSD>1e-6)
    {
    	evaluation.firstMinimumStddev_ZeroStddev=firstMinStddevPSD/firstZeroStddevPSD;
    	evaluation.firstMinimumDiffStddev_ZeroStddev=firstZeroMinAvgPSD/firstZeroStddevPSD;
    }

    // Evaluate micrograph normality
	ImageGeneric M;
	M.readMapped(fnMicrograph);
	double avg, stddev, minval, maxval;
	M().computeStats(avg, stddev, minval, maxval);
	Histogram1D hist;
	compute_hist(M(), hist, minval, maxval, 400);
	hist += 1;
	hist /= hist.sum();

	Histogram1D histGaussian;
	histGaussian.initZeros(hist);
	histGaussian.hmin=hist.hmin;
	histGaussian.hmax=hist.hmax;
	histGaussian.step_size=hist.step_size;
	histGaussian.istep_size=hist.istep_size;
	FOR_ALL_ELEMENTS_IN_ARRAY1D(histGaussian) {
		double x;
		hist.index2val(i, x);
		A1D_ELEM(histGaussian,i) = gaussian1D(x, stddev, avg);
	}
	evaluation.histogramNormality=0.5*(KLDistance(hist,histGaussian)+
			                           KLDistance(histGaussian,hist));

	// Write criteria
	rowOut.setValue(MDL_CTF_DEFOCUSU,evaluation.defocusU);
	rowOut.setValue(MDL_CTF_DEFOCUSV,evaluation.defocusV);
	rowOut.setValue(MDL_CTF_CRIT_FIRSTZEROAVG,evaluation.firstZeroAvg);
	rowOut.setValue(MDL_CTF_CRIT_MAXFREQ,evaluation.maxFreq);
	rowOut.setValue(MDL_CTF_CRIT_DAMPING,evaluation.maxDampingAtBorder);
    if (evaluation.firstZeroDisagreement>0)
    	rowOut.setValue(MDL_CTF_CRIT_FIRSTZERODISAGREEMENT,evaluation.firstZeroDisagreement);

    rowOut.setValue(MDL_CTF_CRIT_FIRSTZERORATIO,evaluation.firstZeroRatio);
    rowOut.setValue(MDL_CTF_ENVELOPE_PLOT,evaluation.ctf_envelope_ssnr);
    rowOut.setValue(MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_RATIO,evaluation.firstMinimumStddev_ZeroStddev);
    rowOut.setValue(MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_DIFF_RATIO,evaluation.firstMinimumDiffStddev_ZeroStddev);
    rowOut.setValue(MDL_CTF_CRIT_FITTINGSCORE,evaluation.fittingScore);
    rowOut.setValue(MDL_CTF_CRIT_FITTINGCORR13,evaluation.fittingCorr13);
    rowOut.setValue(MDL_CTF_CRIT_NONASTIGMATICVALIDITY,evaluation.beating);
    rowOut.setValue(MDL_CTF_CRIT_PSDCORRELATION90,evaluation.PSDcorrelation90);
    rowOut.setValue(MDL_CTF_CRIT_PSDRADIALINTEGRAL,evaluation.PSDradialIntegral);
    rowOut.setValue(MDL_CTF_CRIT_PSDVARIANCE,evaluation.PSDVariance);
    rowOut.setValue(MDL_CTF_CRIT_PSDPCA1VARIANCE,evaluation.PSDPC1Variance);
    rowOut.setValue(MDL_CTF_CRIT_PSDPCARUNSTEST,evaluation.PSDPCRunsTest);
    rowOut.setValue(MDL_CTF_CRIT_NORMALITY, evaluation.histogramNormality);


    //mdEnvelope.write(evaluation.ctf_envelope_ssnr,MD_OVERWRITE);
    mdEnvelope.write(evaluation.ctf_envelope_ssnr,MD_OVERWRITE);
    //std::cout << mdEnvelope << std::endl;

}

