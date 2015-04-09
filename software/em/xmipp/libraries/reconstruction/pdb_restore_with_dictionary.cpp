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

#include "pdb_restore_with_dictionary.h"

void ProgRestoreWithPDBDictionary::defineParams()
{
    // Usage
    addUsageLine("This program takes a set of PDB files at low and high resolution and constructs a dictionary for them.");

    // Parameters
    addParamsLine(" -i <volume>      : Volume to restore");
    addParamsLine(" [-o <volume=\"\">] : Restored volume");
    addParamsLine(" [--root <root=dictionary>] : Rootname of the dictionary");
    addParamsLine(" [--stdThreshold <s=0.15>] : Threshold on the standard deviation to include a patch");
    addParamsLine(" [--angleThreshold <s=5>] : Threshold in degrees on the angle to include a patch");
    addParamsLine(" [--regularization <lambda=0.01>] : Regularization factor for the approximation phase");
    addParamsLine(" [--iterator <s=50>] : Number of iterations");
}

void ProgRestoreWithPDBDictionary::readParams()
{
	fnIn=getParam("-i");
	fnOut=getParam("-o");
	if (fnOut=="")
		fnOut=fnIn;
	fnDictionaryRoot=getParam("--root");
	stdThreshold=getDoubleParam("--stdThreshold");
	angleThreshold=cos(DEG2RAD(getDoubleParam("--angleThreshold")));
	lambda=getDoubleParam("--regularization");
	iterator=getIntParam("--iterator");
}

/* Show -------------------------------------------------------------------- */
void ProgRestoreWithPDBDictionary::show()
{
    if (verbose)
		std::cout
		<< "Input volume:        " << fnIn              << std::endl
		<< "Output volume:       " << fnOut             << std::endl
		<< "Dictionary rootname: " << fnDictionaryRoot  << std::endl
		<< "Standard deviation:  " << stdThreshold      << std::endl
		<< "Angle:               " << angleThreshold    << std::endl
		<< "Iterator:            " << iterator          << std::endl
		;
}

//#define DEBUG
void ProgRestoreWithPDBDictionary::run()
{
    show();

    dictionaryLow.read(fnDictionaryRoot+"_low.mrcs");
    dictionaryHigh.read(fnDictionaryRoot+"_high.mrcs");
	const MultidimArray<double> &mDictionaryLow=dictionaryLow();

    Image<double> V, Vhigh;
    MultidimArray<double> weightHigh;
    V.read(fnIn);
	const MultidimArray<double> &mV=V();
	double min, max, mean, std=0;
	mV.computeStats(mean,std,min,max);

    MultidimArray<double> patchLow, patchLowNormalized, patchHigh;
    int patchSize=XSIZE(mDictionaryLow);
    int patchSize_2=patchSize/2;
	size_t Npatches=0, NcandidatePatches=0;

	std::vector<size_t> idx;
	std::vector<double> weight;
	Matrix1D<double> alpha;
	init_progress_bar(ZSIZE(mV));
	patchHigh.resizeNoCopy(patchSize,patchSize,patchSize);
	Vhigh().initZeros(mV);
	weightHigh.initZeros(mV);
	const MultidimArray<double> &mVhigh=Vhigh();
	for (int k=patchSize_2; k<(int)ZSIZE(mV)-patchSize_2; ++k)
	{
		for (int i=patchSize_2; i<(int)ZSIZE(mV)-patchSize_2; ++i)
			for (int j=patchSize_2; j<(int)ZSIZE(mV)-patchSize_2; ++j)
			{
				 ++Npatches;
				 mV.window(patchLow,k,i,j,k+patchSize-1,i+patchSize-1,j+patchSize-1);
				 STARTINGX(patchLow)=STARTINGY(patchLow)=STARTINGZ(patchLow)=0;

				double minPatchLow, maxPatchLow, meanPatchLow, stdPatchLow=0;
				patchLow.computeStats(meanPatchLow,stdPatchLow,minPatchLow,maxPatchLow);
				double R2=1;
				patchHigh.initZeros(patchLow);

				if (stdPatchLow > stdThreshold*std)
				{
					++NcandidatePatches;
					patchLowNormalized=patchLow/sqrt(patchLow.sum2());
					selectDictionaryPatches(patchLowNormalized,idx,weight);
					if (idx.size()>0)
					{
						R2=approximatePatch(patchLow,idx,weight,alpha);
					    reconstructPatch(idx,alpha,patchHigh);
//						Image<double> save;
//						save()=patchLow;
//						save.write("PPPlow.vol");
//						save()=patchHigh;
//						save.write("PPPhigh.vol");
//						std::cout << "Press any key" << std::endl;
//						char c; std::cin >> c;
					}
				}

				// Insert patchHigh in Vhigh
				for (int kk=0; kk<patchSize; ++kk)
					for (int ii=0; ii<patchSize; ++ii)
						for (int jj=0; jj<patchSize; ++jj)
						{
							A3D_ELEM(mVhigh,k+kk-patchSize_2,i+ii-patchSize_2,j+jj-patchSize_2)+=R2*DIRECT_A3D_ELEM(patchHigh,kk,ii,jj);
							A3D_ELEM(weightHigh,k+kk-patchSize_2,i+ii-patchSize_2,j+jj-patchSize_2)+=R2;
						}
			}
		progress_bar(k);
	}
	progress_bar(ZSIZE(mV));

	// Correct by the Vhigh weights
	FOR_ALL_ELEMENTS_IN_ARRAY3D(mVhigh)
	if (A3D_ELEM(weightHigh,k,i,j)>0)
		A3D_ELEM(mVhigh,k,i,j)/=A3D_ELEM(weightHigh,k,i,j);

	Vhigh.write(fnOut);
}

//#define DEBUG
void ProgRestoreWithPDBDictionary::selectDictionaryPatches(const MultidimArray<double> &lowResolutionPatch,
		std::vector<size_t> &idx, std::vector<double> &weight)
{
	idx.clear();
	weight.clear();
	const MultidimArray<double> &mDictionaryLow=dictionaryLow();
	for (size_t i=0; i<NSIZE(mDictionaryLow); ++i)
	{
		const double *ptrDictionaryPatch=&DIRECT_NZYX_ELEM(mDictionaryLow,i,0,0,0);
		double dotProduct=0;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(lowResolutionPatch)
			dotProduct+=(*ptrDictionaryPatch++)*DIRECT_MULTIDIM_ELEM(lowResolutionPatch,n);
		if (dotProduct>angleThreshold)
		{
			idx.push_back(i);
			weight.push_back(1/dotProduct);
#ifdef DEBUG
			Image<double> save;
			save()=lowResolutionPatch;
			save.write("PPPexperimentalPatch.vol");
			memcpy(&DIRECT_MULTIDIM_ELEM(save(),0),&DIRECT_NZYX_ELEM(mDictionaryLow,i,0,0,0),
					MULTIDIM_SIZE(lowResolutionPatch)*sizeof(double));
			save.write("PPPdictionaryPatch.vol");
			std::cout << RAD2DEG(acos(dotProduct)) << ". Press any key" << std::endl;
			char c; std::cin >> c;
#endif
		}
	}
//	std::cout << idx.size() << std::endl;
}
#undef DEBUG

/* See Algorithm 1 in Trinh, Luong, Dibos, Rocchisani, Pham, Nguyen. Novel Exanple-based method for super-resolution
 *    and denoising of medical images. IEEE Trans. Image Processing, 23: 1882-1895 (2014)
 */
double ProgRestoreWithPDBDictionary::approximatePatch(const MultidimArray<double> &lowResolutionPatch,
		std::vector<size_t> &idx, std::vector<double> &weight, Matrix1D<double> &alpha)
{
	const MultidimArray<double> &mDictionaryLow=dictionaryLow();

	alpha.initZeros(idx.size());
    wi=alpha;
	Ui.initZeros(MULTIDIM_SIZE(lowResolutionPatch),VEC_XSIZE(alpha));

	// Prepare wi
	FOR_ALL_ELEMENTS_IN_MATRIX1D(wi)
	VEC_ELEM(wi,i)=lambda*(1+weight[i]);

	// Prepare Ui
	for (size_t j=0; j<MAT_XSIZE(Ui); ++j)
	{
		const double *ptr=&DIRECT_NZYX_ELEM(mDictionaryLow,idx[j],0,0,0);
		for (size_t i=0; i<MAT_YSIZE(Ui); ++i, ++ptr)
			MAT_ELEM(Ui,i,j)=*ptr;
	}
	matrixOperation_AtA(Ui,UitUi);

	// Prepare y
	y.resizeNoCopy(MULTIDIM_SIZE(lowResolutionPatch));
	memcpy(&VEC_ELEM(y,0),&DIRECT_MULTIDIM_ELEM(lowResolutionPatch,0),MULTIDIM_SIZE(lowResolutionPatch)*sizeof(double));

	// Prepare alpha
	alpha.initConstant(1);

	// Iterate
	for (int iter=0; iter<iterator; ++iter)
	{
		matrixOperation_Atx(Ui,y,v1);
		matrixOperation_Ax(UitUi,alpha,v2);
		FOR_ALL_ELEMENTS_IN_MATRIX1D(alpha)
		VEC_ELEM(alpha,i)*=VEC_ELEM(v1,i)/(VEC_ELEM(v2,i)+VEC_ELEM(wi,i));
		// std::cout << "iter=" << iter << " alpha=" << alpha << std::endl;
	}

	// Calculate R2
	matrixOperation_Ax(Ui,alpha,yp);
	double diff2=0.0, norm2=0.0;
	FOR_ALL_ELEMENTS_IN_MATRIX1D(y)
	{
		double diff=VEC_ELEM(y,i)-VEC_ELEM(yp,i);
		diff2+=diff*diff;
		norm2+=VEC_ELEM(y,i)*VEC_ELEM(y,i);
	}
	return 1-diff2/norm2;
}

void ProgRestoreWithPDBDictionary::reconstructPatch(std::vector<size_t> &idx, Matrix1D<double> &alpha,
		   MultidimArray<double> &highResolutionPatch)
{
	const MultidimArray<double> &mDictionaryHigh=dictionaryHigh();
	for (size_t j=0; j<VEC_XSIZE(alpha); ++j)
	{
		const double *ptr=&DIRECT_NZYX_ELEM(mDictionaryHigh,idx[j],0,0,0);
		for (size_t n=0; n<MULTIDIM_SIZE(highResolutionPatch); ++n, ++ptr)
			DIRECT_MULTIDIM_ELEM(highResolutionPatch,n)+=VEC_ELEM(alpha,j)*(*ptr);
	}
}
