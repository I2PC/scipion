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
		sigmaConv=sigma0;
		for (int iter=0; iter<NiterFourier; iter++)
		{
			std::cout << "Deconvolution iteration " << iter << std::endl;
			estimateS();

			transformer.FourierTransform(V1r(),fV1r);
			transformer.FourierTransform(V2r(),fV2r);
			transformer.FourierTransform(S(),fVol);
			optimizeSigma();

			deconvolveS();
		}
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
    double K=-0.5/(sigmaConv*sigmaConv);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fVol)
    {
		double R2n=DIRECT_MULTIDIM_ELEM(R2,n);
		if (R2n<=0.25)
		{
			double H1=exp(K*R2n);
			double H2=H1;
			DIRECT_MULTIDIM_ELEM(fVol,n)=(H1*DIRECT_MULTIDIM_ELEM(fV1r,n)+H2*DIRECT_MULTIDIM_ELEM(fV2r,n))/(H1*H1+H2*H2+lambda*R2n);
		}
    }

    MultidimArray<double> &mS=S();
    transformer.inverseFourierTransform(fVol,mS);
    S.write("PPPS.vol");
}

double restorationSigmaCost(double *x, void *_prm)
{
	ProgVolumeHalvesRestoration *prm=(ProgVolumeHalvesRestoration *) _prm;
	double sigma=x[1];
	if (sigma<0)
		return 1e38;
    double K=-0.5/(sigma*sigma);
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
			double H1=exp(K*R2n);
			error+=abs(DIRECT_MULTIDIM_ELEM(fV,n)*H1-DIRECT_MULTIDIM_ELEM(fV1r,n))+
				   abs(DIRECT_MULTIDIM_ELEM(fV,n)*H1-DIRECT_MULTIDIM_ELEM(fV2r,n));
			N++;
		}
	}
    std::cout << sigma << " "  << error << std::endl;
    return error;
}

void ProgVolumeHalvesRestoration::optimizeSigma()
{

    Matrix1D<double> p(1), steps(1);
    p(0)=sigmaConv;
    steps.initConstant(1);
    double cost;
    int iter;
	powellOptimizer(p, 1, 1, &restorationSigmaCost, this, 0.01, cost, iter, steps, verbose>=2);
}
