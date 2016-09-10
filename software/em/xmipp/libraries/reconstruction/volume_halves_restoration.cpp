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

void CDF::calculateCDF(MultidimArray<double> &V, double probStep)
{
	double *ptr=&DIRECT_MULTIDIM_ELEM(V,0);
	size_t N=MULTIDIM_SIZE(V);
	std::sort(ptr,ptr+N);
	minVal = ptr[0];
	maxVal = ptr[N-1];

	int Nsteps=(int)round(1.0/probStep);
	x.resizeNoCopy(Nsteps);
	probXLessThanx.resizeNoCopy(Nsteps);
	int i=0;
	for (double p=probStep/2; p<1; p+=probStep, i++)
	{
		size_t idx=(size_t)round(p*N);
		A1D_ELEM(probXLessThanx,i)=p;
		A1D_ELEM(x,i)=ptr[idx];
//		std::cout << "i=" << i << " x=" << A1D_ELEM(x,i) << " p=" << A1D_ELEM(probXLessThanx,i) << std::endl;
	}
}

#define INTERP(x,x0,y0,xF,yF) (y0+(x-x0)*(yF-y0)/(xF-x0))

double CDF::getProbability(double xi)
{
	if (xi>maxVal)
		return 1;
	else if (xi<minVal)
		return 0;
	else
	{
		size_t N=XSIZE(x);
		if (xi<DIRECT_A1D_ELEM(x,0))
			return INTERP(xi,minVal,0.0,DIRECT_A1D_ELEM(x,0),DIRECT_A1D_ELEM(probXLessThanx,0));
		else if (xi>DIRECT_A1D_ELEM(x,N-1))
			return INTERP(xi,DIRECT_A1D_ELEM(x,N-1),DIRECT_A1D_ELEM(probXLessThanx,N-1),maxVal,1.0);
		else
		{
			int iLeft=0;
			int iRight=N-1;
			while (iLeft<=iRight)
			{
				int iMiddle = iLeft+(iRight-iLeft)/2;
				if (xi>=DIRECT_A1D_ELEM(x,iMiddle) && xi<=DIRECT_A1D_ELEM(x,iMiddle+1))
				{
//					std::cout << "iMiddle = " << iMiddle << " " << DIRECT_A1D_ELEM(x,iMiddle) << " "  << DIRECT_A1D_ELEM(x,iMiddle+1) << std::endl;
					if (DIRECT_A1D_ELEM(x,iMiddle)==DIRECT_A1D_ELEM(x,iMiddle+1))
						return 0.5*(DIRECT_A1D_ELEM(probXLessThanx,iMiddle)+DIRECT_A1D_ELEM(probXLessThanx,iMiddle+1));
					else
						return INTERP(xi,DIRECT_A1D_ELEM(x,iMiddle),  DIRECT_A1D_ELEM(probXLessThanx,iMiddle),
										 DIRECT_A1D_ELEM(x,iMiddle+1),DIRECT_A1D_ELEM(probXLessThanx,iMiddle+1));
				}
				else if (xi<DIRECT_A1D_ELEM(x,iMiddle))
					iRight=iMiddle;
				else
					iLeft=iMiddle;
			}
			std::cout << "It should never reach here" << std::endl;
		}
	}
}

// Read arguments ==========================================================
void ProgVolumeHalvesRestoration::readParams()
{
    fnV1 = getParam("--i1");
    fnV2 = getParam("--i2");
    fnRoot = getParam("--oroot");
    applyPos = checkParam("--applyPositivity");
    Niter = getIntParam("--Niter");
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
	<< "Niter:    " << Niter << std::endl
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
    addParamsLine("  [--applyPositivity]          : Remove negative values");
    addParamsLine("  [--Niter <N=5>]              : Number of iterations");
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
	H.resizeNoCopy(XSIZE(fVol));
	H.initConstant(1.0);
	Ridx.resizeNoCopy(fVol);
	Ridx.initConstant(-1);
	const MultidimArray<double> &mV1=V1();
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(fVol)
	{
        double fz, fy, fx;
		FFT_IDX2DIGFREQ(k,ZSIZE(mV1),fz);
		FFT_IDX2DIGFREQ(i,YSIZE(mV1),fy);
		FFT_IDX2DIGFREQ(j,XSIZE(mV1),fx);
		double R=sqrt(fx*fx+fy*fy+fz*fz);
		if (R<0.5)
			A3D_ELEM(Ridx,k,i,j) = (int)round(R*XSIZE(mV1));
	}
}

void ProgVolumeHalvesRestoration::run()
{
	show();
	produceSideInfo();

	for (int iter=0; iter<Niter; iter++)
	{
		std::cout << "Iteration " << iter << std::endl;
		estimateS();
		estimateHS();
		significanceRealSpace(V1r(),V1r());
		significanceRealSpace(V2r(),V2r());
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
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mV1r)
	DIRECT_MULTIDIM_ELEM(mS,n)=0.5*(DIRECT_MULTIDIM_ELEM(mV1r,n)+DIRECT_MULTIDIM_ELEM(mV2r,n));
	applyMask(mS);
}

void ProgVolumeHalvesRestoration::applyMask(MultidimArray<double> &V)
{
	if (pMask==NULL)
		return;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V)
	if (DIRECT_MULTIDIM_ELEM(*pMask,n)==0)
		DIRECT_MULTIDIM_ELEM(V,n)=0;
	if (applyPos)
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V)
		if (DIRECT_MULTIDIM_ELEM(V,n)<0)
			DIRECT_MULTIDIM_ELEM(V,n)=0;
}

void ProgVolumeHalvesRestoration::estimateHS()
{
	// Filter S
	HS()=S();
	transformer.FourierTransform(HS(),fVol,false);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fVol)
	{
		int fIdx = DIRECT_MULTIDIM_ELEM(Ridx,n);
		if (fIdx>=0)
			DIRECT_MULTIDIM_ELEM(Ridx,n)*=DIRECT_MULTIDIM_ELEM(H,fIdx);
	}
	transformer.inverseFourierTransform();

	// Calculate HS CDF in real space
	MultidimArray<double> aux;
	if (pMask==NULL)
	{
		aux=HS();
		aux*=aux;
	}
	else
	{
		aux.resizeNoCopy(pMaskSize);
		size_t idx=0;
		const MultidimArray<double> mHS=HS();
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mHS)
		if (DIRECT_MULTIDIM_ELEM(*pMask,n)!=0)
			DIRECT_MULTIDIM_ELEM(aux,idx++)=DIRECT_MULTIDIM_ELEM(mHS,n)*DIRECT_MULTIDIM_ELEM(mHS,n);
	}
	cdfHS.calculateCDF(aux);
}

void ProgVolumeHalvesRestoration::significanceRealSpace(const MultidimArray<double> &Vi, MultidimArray<double> &Vir)
{
	// Calculate N=Vi-H*S and its energy CDF
	N().resizeNoCopy(Vi);
	MultidimArray<double> &mN=N();
	const MultidimArray<double> &mHS=HS();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mN)
	{
		double diff=DIRECT_MULTIDIM_ELEM(Vi,n)-DIRECT_MULTIDIM_ELEM(mHS,n);
		DIRECT_MULTIDIM_ELEM(mN,n)=diff*diff;
	}
	CDF cdfN;
	cdfN.calculateCDF(mN);

	// Mask the input volume with a mask value that is proportional to the probability of being larger than noise
	Vir.resizeNoCopy(Vi);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mN)
	{
		double e=DIRECT_MULTIDIM_ELEM(Vi,n)*DIRECT_MULTIDIM_ELEM(Vi,n);
//		std::cout << "n=" << n << " e=" << e;
		double pN=cdfN.getProbability(e);
//		std::cout << " pN=" << pN;
		if (pN<1)
		{
			double pHS=cdfHS.getProbability(e);
//			std::cout << " pHS=" << pHS << std::endl;
			double pp=pHS*pN;
			DIRECT_MULTIDIM_ELEM(Vir,n)=pp*DIRECT_MULTIDIM_ELEM(Vi,n);
		}
		else
			DIRECT_MULTIDIM_ELEM(Vir,n)=DIRECT_MULTIDIM_ELEM(Vi,n);
	}
}
