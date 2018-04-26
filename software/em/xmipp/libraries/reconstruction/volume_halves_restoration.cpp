/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
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

#include "volume_halves_restoration.h"
#include <data/numerical_tools.h>

// Read arguments ==========================================================
void ProgVolumeHalvesRestoration::readParams()
{
    fnV1 = getParam("--i1");
    fnV2 = getParam("--i2");
    fnRoot = getParam("--oroot");
    NiterReal = getIntParam("--denoising");
    NiterFourier = getIntParam("--deconvolution");
    sigma0 = getDoubleParam("--deconvolution",1);
    lambda = getDoubleParam("--deconvolution",2);
    bankStep = getDoubleParam("--filterBank",0);
    bankOverlap = getDoubleParam("--filterBank",1);
    weightFun = getIntParam("--filterBank",2);
    weightPower = getDoubleParam("--filterBank",3);
    NiterDiff = getIntParam("--difference");
    Kdiff = getDoubleParam("--difference",1);
    if (checkParam("--mask"))
        mask.readParams(this);
}

// Show ====================================================================
void ProgVolumeHalvesRestoration::show()
{
    if (!verbose)
        return;
    std::cout
	<< "Volume1:  " << fnV1 << std::endl
	<< "Volume2:  " << fnV2 << std::endl
	<< "Rootname: " << fnRoot << std::endl
	<< "Denoising Iterations: " << NiterReal << std::endl
	<< "Deconvolution Iterations: " << NiterFourier << std::endl
	<< "Sigma0:   " << sigma0 << std::endl
	<< "Lambda:   " << lambda << std::endl
	<< "Bank step:" << bankStep << std::endl
	<< "Bank overlap:" << bankOverlap << std::endl
	<< "Weight fun:" << weightFun << std::endl
	<< "Weight power:" << weightPower << std::endl
	<< "Difference Iterations: " << NiterDiff << std::endl
	<< "Kdiff: " << Kdiff << std::endl
	;
    mask.show();
}

// usage ===================================================================
void ProgVolumeHalvesRestoration::defineParams()
{
    addUsageLine("Given two halves of a volume (and an optional mask), produce a better estimate of the volume underneath");
    addParamsLine("   --i1 <volume1>              : First half");
    addParamsLine("   --i2 <volume2>              : Second half");
    addParamsLine("  [--oroot <root=\"volumeRestored\">] : Output rootname");
    addParamsLine("  [--denoising <N=0>]          : Number of iterations of denoising in real space");
    addParamsLine("  [--deconvolution <N=0> <sigma0=0.2> <lambda=0.001>]   : Number of iterations of deconvolution in Fourier space, initial sigma and lambda");
    addParamsLine("  [--filterBank <step=0> <overlap=0.5> <weightFun=1> <weightPower=3>] : Frequency step for the filter bank (typically, 0.01; between 0 and 0.5)");
    addParamsLine("                                        : filter overlap is between 0 (no overlap) and 1 (full overlap)");
    addParamsLine("                                : Weight function (0=mean, 1=min, 2=mean*diff");
    addParamsLine("  [--difference <N=0> <K=1.5>]  : Number of iterations of difference evaluation in real space");
    Mask::defineParams(this,INT_MASK);
}

void ProgVolumeHalvesRestoration::produceSideInfo()
{
	V1.read(fnV1);
	V2.read(fnV2);
	V1().setXmippOrigin();
	V2().setXmippOrigin();

	if (mask.fn_mask!="")
	{
		mask.generate_mask();
		pMask=&mask.get_binary_mask();
		pMaskSize=(size_t)(pMask->sum());
	}
	else
		pMask = NULL;

	V1r()=V1(); // Copy the restored volumes
	V2r()=V2();

	// Initialize filter and calculate frequencies
	transformer.FourierTransform(V1(),fVol,false);
	R2.resizeNoCopy(fVol);
	R2.initConstant(-1);
	const MultidimArray<double> &mV1=V1();
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(fVol)
	{
        double fz, fy, fx;
		FFT_IDX2DIGFREQ(k,ZSIZE(mV1),fz);
		FFT_IDX2DIGFREQ(i,YSIZE(mV1),fy);
		FFT_IDX2DIGFREQ(j,XSIZE(mV1),fx);
		A3D_ELEM(R2,k,i,j) = fx*fx+fy*fy+fz*fz;
	}
}

void ProgVolumeHalvesRestoration::run()
{
	show();
	produceSideInfo();

	if (NiterReal>0)
		for (int iter=0; iter<NiterReal; iter++)
		{
			std::cout << "Denoising iteration " << iter << std::endl;
			estimateS();
			significanceRealSpace(V1r(),V1r());
			significanceRealSpace(V2r(),V2r());
		}

	if (NiterFourier>0)
	{
		sigmaConv1=sigmaConv2=sigma0;
		for (int iter=0; iter<NiterFourier; iter++)
		{
			std::cout << "Deconvolution iteration " << iter << std::endl;
			estimateS();

			transformer1.FourierTransform(V1r(),fV1r, false);
			transformer2.FourierTransform(V2r(),fV2r, false);
			transformer.FourierTransform(S(),fVol);
			optimizeSigma();

			deconvolveS();
		}

		S.write(fnRoot+"_deconvolved.vol");
		convolveS();
		S.write(fnRoot+"_convolved.vol");
	}

	if (bankStep>0)
		filterBank();

	if (NiterDiff>0)
		for (int iter=0; iter<NiterDiff; iter++)
		{
			std::cout << "Difference iteration " << iter << std::endl;
			evaluateDifference();
		}

	V1r.write(fnRoot+"_restored1.vol");
	V2r.write(fnRoot+"_restored2.vol");
}

void ProgVolumeHalvesRestoration::estimateS()
{
	S().resizeNoCopy(V1r());
	MultidimArray<double> &mS=S();
	const MultidimArray<double> &mV1r=V1r();
	const MultidimArray<double> &mV2r=V2r();

	// Compute average
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mV1r)
	DIRECT_MULTIDIM_ELEM(mS,n)=0.5*(DIRECT_MULTIDIM_ELEM(mV1r,n)+DIRECT_MULTIDIM_ELEM(mV2r,n));

	// Apply mask
	if (pMask!=NULL)
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mS)
		if (DIRECT_MULTIDIM_ELEM(*pMask,n)==0)
			DIRECT_MULTIDIM_ELEM(mS,n)=0;

	// Apply positivity
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mS)
	if (DIRECT_MULTIDIM_ELEM(mS,n)<0)
		DIRECT_MULTIDIM_ELEM(mS,n)=0;

	// Filter S
	transformer.FourierTransform(mS,fVol,false);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fVol)
		if (DIRECT_MULTIDIM_ELEM(R2,n)>0.25)
			DIRECT_MULTIDIM_ELEM(fVol,n)=0.0;
	transformer.inverseFourierTransform();

	// Calculate HS CDF in real space
	MultidimArray<double> aux;
	if (pMask==NULL)
	{
		aux=S();
		aux*=aux;
	}
	else
	{
		aux.resizeNoCopy(pMaskSize);
		size_t idx=0;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mS)
		if (DIRECT_MULTIDIM_ELEM(*pMask,n)!=0)
			DIRECT_MULTIDIM_ELEM(aux,idx++)=DIRECT_MULTIDIM_ELEM(mS,n)*DIRECT_MULTIDIM_ELEM(mS,n);
	}
	cdfS.calculateCDF(aux);
}

void ProgVolumeHalvesRestoration::significanceRealSpace(const MultidimArray<double> &Vi, MultidimArray<double> &Vir)
{
	// Calculate N=Vi-H*S and its energy CDF
	N().resizeNoCopy(Vi);
	MultidimArray<double> &mN=N();
	const MultidimArray<double> &mS=S();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mN)
	{
		double diff=DIRECT_MULTIDIM_ELEM(Vi,n)-DIRECT_MULTIDIM_ELEM(mS,n);
		DIRECT_MULTIDIM_ELEM(mN,n)=diff*diff;
	}
	CDF cdfN;
	cdfN.calculateCDF(mN);

	// Mask the input volume with a mask value that is proportional to the probability of being larger than noise
	Vir.resizeNoCopy(Vi);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mN)
	{
		double e=DIRECT_MULTIDIM_ELEM(Vi,n)*DIRECT_MULTIDIM_ELEM(Vi,n);
		double pN=cdfN.getProbability(e);
		if (pN<1)
		{
			double pS=cdfS.getProbability(e);
			double pp=pS*pN;
			DIRECT_MULTIDIM_ELEM(Vir,n)=pp*DIRECT_MULTIDIM_ELEM(Vi,n);
		}
		else
			DIRECT_MULTIDIM_ELEM(Vir,n)=DIRECT_MULTIDIM_ELEM(Vi,n);
	}
}

void ProgVolumeHalvesRestoration::deconvolveS()
{
	if (verbose>0)
		std::cout << "   Deconvolving with sigma=" << sigmaConv1  << " " << sigmaConv2 << std::endl;
    double K1=-0.5/(sigmaConv1*sigmaConv1);
    double K2=-0.5/(sigmaConv2*sigmaConv2);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fVol)
    {
		double R2n=DIRECT_MULTIDIM_ELEM(R2,n);
		if (R2n<=0.25)
		{
			double H1=exp(K1*R2n);
			double H2=exp(K2*R2n);
			DIRECT_MULTIDIM_ELEM(fVol,n)=(H1*DIRECT_MULTIDIM_ELEM(fV1r,n)+H2*DIRECT_MULTIDIM_ELEM(fV2r,n))/(H1*H1+H2*H2+lambda*R2n);

			DIRECT_MULTIDIM_ELEM(fV1r,n)/=H1;
			DIRECT_MULTIDIM_ELEM(fV2r,n)/=H2;
		}
    }
    transformer1.inverseFourierTransform();
    transformer2.inverseFourierTransform();

//    MultidimArray<double> &mS=S();
//    transformer.inverseFourierTransform();
//    S.write("PPPS.vol");
}

void ProgVolumeHalvesRestoration::convolveS()
{
	double sigmaConv=(sigmaConv1+sigmaConv2)/2;
    double K=-0.5/(sigmaConv*sigmaConv);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fVol)
    {
		double R2n=DIRECT_MULTIDIM_ELEM(R2,n);
		if (R2n<=0.25)
		{
			double H1=exp(K*R2n);
			DIRECT_MULTIDIM_ELEM(fVol,n)*=H1;
		}
    }
    transformer.inverseFourierTransform();
}

double restorationSigmaCost(double *x, void *_prm)
{
	ProgVolumeHalvesRestoration *prm=(ProgVolumeHalvesRestoration *) _prm;
	double sigma1=x[1];
	double sigma2=x[2];
	if (sigma1<0 || sigma2<0 || sigma1>2 || sigma2>2)
		return 1e38;
    double K1=-0.5/(sigma1*sigma1);
    double K2=-0.5/(sigma2*sigma2);
    double error=0;
    double N=0;
    const MultidimArray< std::complex<double> > &fV=prm->fVol;
    const MultidimArray< std::complex<double> > &fV1r=prm->fV1r;
    const MultidimArray< std::complex<double> > &fV2r=prm->fV2r;
    const MultidimArray<double> &R2=prm->R2;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fV)
	{
		double R2n=DIRECT_MULTIDIM_ELEM(R2,n);
		if (R2n<=0.25)
		{
			double H1=exp(K1*R2n);
			double H2=exp(K2*R2n);
			error+=abs(DIRECT_MULTIDIM_ELEM(fV,n)*H1-DIRECT_MULTIDIM_ELEM(fV1r,n))+
				   abs(DIRECT_MULTIDIM_ELEM(fV,n)*H2-DIRECT_MULTIDIM_ELEM(fV2r,n));
			N++;
		}
	}
    return error;
}

void ProgVolumeHalvesRestoration::optimizeSigma()
{

    Matrix1D<double> p(2), steps(2);
    p(0)=sigmaConv1;
    p(1)=sigmaConv2;
    steps.initConstant(1);
    double cost;
    int iter;
	powellOptimizer(p, 1, 2, &restorationSigmaCost, this, 0.01, cost, iter, steps, verbose>=2);
	sigmaConv1=p(0);
	sigmaConv2=p(1);
}

void ProgVolumeHalvesRestoration::filterBand(const MultidimArray< std::complex<double> > &fV, FourierTransformer &transformer, double w)
{
	// Apply band pass filter
	double w2=w*w;
	double w2Step=(w+bankStep)*(w+bankStep);
	MultidimArray< std::complex<double> >&fVout=transformer.fFourier;
	fVout.initZeros();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fV)
	{
		double R2n=DIRECT_MULTIDIM_ELEM(R2,n);
		if (R2n>=w2 && R2n<w2Step)
			DIRECT_MULTIDIM_ELEM(fVout,n)=DIRECT_MULTIDIM_ELEM(fV,n);
	}
	transformer.inverseFourierTransform();
}

//#define DEBUG
void ProgVolumeHalvesRestoration::filterBank()
{
	MultidimArray<double> Vfiltered1, Vfiltered2;
	Vfiltered1.initZeros(V1r());
	Vfiltered2.initZeros(V2r());

	transformer1.FourierTransform(V1r(),fV1r, false);
	transformer2.FourierTransform(V2r(),fV2r, false);

	FourierTransformer bank1, bank2;
	bank1.setReal(Vfiltered1);
	bank2.setReal(Vfiltered2);
	double filterStep = bankStep*(1-bankOverlap);
	MultidimArray<double> &mV1r=V1r();
	MultidimArray<double> &mV2r=V2r();
	MultidimArray<double> &mS=S();
	mV1r.initZeros(Vfiltered2);
	mV2r.initZeros(Vfiltered2);
	mS.initZeros(Vfiltered2);
	int i=0;
	int imax=ceil(0.5/filterStep);
	std::cerr << "Calculating filter bank ..." << std::endl;
	init_progress_bar(imax);
	for (double w=0; w<0.5; w+=filterStep)
	{
		filterBand(fV1r,bank1,w);
		filterBand(fV2r,bank2,w);

		// Compute energy of the noise
		N().resizeNoCopy(Vfiltered1);
		MultidimArray<double> &mN=N();
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mN)
		{
			double diff=DIRECT_MULTIDIM_ELEM(Vfiltered1,n)-DIRECT_MULTIDIM_ELEM(Vfiltered2,n);
			DIRECT_MULTIDIM_ELEM(mN,n)=0.5*diff*diff;
		}

		CDF cdfN;
		cdfN.calculateCDF(mN);

		// Compute weights
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mN)
		{
			double e1=DIRECT_MULTIDIM_ELEM(Vfiltered1,n)*DIRECT_MULTIDIM_ELEM(Vfiltered1,n);
			double w1=cdfN.getProbability(e1);

			double e2=DIRECT_MULTIDIM_ELEM(Vfiltered2,n)*DIRECT_MULTIDIM_ELEM(Vfiltered2,n);
			double w2=cdfN.getProbability(e2);

			double weight;
			switch (weightFun)
			{
			case 0: weight=0.5*(w1+w2); break;
			case 1: weight=std::min(w1,w2); break;
			case 2: weight=0.5*(w1+w2)*(1-fabs(w1-w2)/(w1+w2)); break;
			}
			weight=pow(weight,weightPower);

			double Vf1w=DIRECT_MULTIDIM_ELEM(Vfiltered1,n)*weight;
			double Vf2w=DIRECT_MULTIDIM_ELEM(Vfiltered2,n)*weight;
			DIRECT_MULTIDIM_ELEM(mV1r,n)+=Vf1w;
			DIRECT_MULTIDIM_ELEM(mV2r,n)+=Vf2w;
			if (e1>e2)
				DIRECT_MULTIDIM_ELEM(mS,n)+=Vf1w;
			else
				DIRECT_MULTIDIM_ELEM(mS,n)+=Vf2w;
		}
		progress_bar(++i);

#ifdef DEBUG
		Image<double> save;
		save()=Vfiltered1;
		save.write("PPP1.vol");
		save()=mV1r;
		save.write("PPP1r.vol");
		save()=Vfiltered2;
		save.write("PPP2.vol");
		save()=mS;
		save.write("PPPS.vol");
		save()=Vfiltered2-Vfiltered1;
		save.write("PPPdiff.vol");
		save()=weight;
		save.write("PPPweight.vol");
		save()=sumWeight;
		save.write("PPPsumWeight.vol");
		std::cout << "w=" << w << " Press";
		char c; std::cin >> c;
#endif
	}
	progress_bar(imax);

	S()  *= 1-bankOverlap;
	mV1r *= 1-bankOverlap;
	mV2r *= 1-bankOverlap;
	S.write(fnRoot+"_filterBank.vol");
}

void ProgVolumeHalvesRestoration::evaluateDifference()
{
	// Compute the difference between the two
	N()=V1r();
	N()-=V2r();

	S()=V1r();
	S()+=V2r();
	S()*=0.5;

	// Compute the std within the signal mask
	double mean, stddev;
	if (pMask==NULL)
		N().computeAvgStdev(mean,stddev);
	else
		N().computeAvgStdev_within_binary_mask(*pMask,mean,stddev);
	stddev*=Kdiff;

	MultidimArray<double> &m1=V1r();
	MultidimArray<double> &m2=V2r();
	MultidimArray<double> &mN=N();
	MultidimArray<double> &mS=S();
	double k=-0.5/(stddev*stddev);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mN)
	{
		double w=exp(k*DIRECT_MULTIDIM_ELEM(mN,n)*DIRECT_MULTIDIM_ELEM(mN,n));
		double s=DIRECT_MULTIDIM_ELEM(mS,n);
		double d1=DIRECT_MULTIDIM_ELEM(m1,n)-s;
		double d2=DIRECT_MULTIDIM_ELEM(m2,n)-s;
		DIRECT_MULTIDIM_ELEM(m1,n)=s+d1*w;
		DIRECT_MULTIDIM_ELEM(m2,n)=s+d2*w;
	}

	S()=V1r();
	S()+=V2r();
	S()*=0.5;
	S.write(fnRoot+"_avgDiff.vol");
}
