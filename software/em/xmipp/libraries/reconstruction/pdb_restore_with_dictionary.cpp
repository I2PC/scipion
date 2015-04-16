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
    addParamsLine(" [--regularization <lambda=0.01>] : Regularization factor for the approximation phase");
    addParamsLine(" [--iter <s=50>] : Number of iterations");
	ProgPDBDictionary::defineParams();
}

void ProgRestoreWithPDBDictionary::readParams()
{
	fnIn=getParam("-i");
	fnOut=getParam("-o");
	if (fnOut=="")
		fnOut=fnIn;
	fnRoot=getParam("--root");
	lambda=getDoubleParam("--regularization");
	iterations=getIntParam("--iter");
	ProgPDBDictionary::readParams();
}

/* Show -------------------------------------------------------------------- */
void ProgRestoreWithPDBDictionary::show()
{
    if (verbose)
    {
		std::cout
		<< "Input volume:        " << fnIn              << std::endl
		<< "Output volume:       " << fnOut             << std::endl
		<< "Iterations:          " << iterations        << std::endl
		;
		ProgPDBDictionary::show();
    }
}

//#define DEBUG
void ProgRestoreWithPDBDictionary::run()
{
    show();
    loadDictionaries();
    constructRotationGroup();

    Image<double> V, Vhigh;
    MultidimArray<double> weightHigh;
    V.read(fnIn);
	const MultidimArray<double> &mV=V();
	double min, max, mean, std=0;
	mV.computeStats(mean,std,min,max);

    MultidimArray<double> patchLow, patchLowNormalized, canonicalPatch, patchHigh;
    int patchSize_2=patchSize/2;
	size_t Npatches=0, NcandidatePatches=0;

	std::vector< size_t > selectedPatchesIdx;
	std::vector<double> weight;
	Matrix1D<double> alpha, canonicalSignature;
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
				 mV.window(patchLow,k-patchSize_2,i-patchSize_2,j-patchSize_2,k+patchSize_2,i+patchSize_2,j+patchSize_2);
				 STARTINGX(patchLow)=STARTINGY(patchLow)=STARTINGZ(patchLow)=0;

				double minPatchLow, maxPatchLow, meanPatchLow, stdPatchLow=0;
				patchLow.computeStats(meanPatchLow,stdPatchLow,minPatchLow,maxPatchLow);
				double R2=0;
				patchHigh.initZeros(patchLow);

				if (stdPatchLow > stdThreshold*std)
				{
					++NcandidatePatches;
					patchLowNormalized=patchLow;
					double norm=sqrt(patchLow.sum2());
					patchLowNormalized*=1.0/norm;
					size_t idxTransf=canonicalOrientation(patchLowNormalized,canonicalPatch,canonicalSignature);
					selectDictionaryPatches(canonicalPatch, canonicalSignature, selectedPatchesIdx, weight);
					if (selectedPatchesIdx.size()>0)
					{
						R2=approximatePatch(canonicalPatch,selectedPatchesIdx,weight,alpha);
					    reconstructPatch(idxTransf,selectedPatchesIdx,alpha,patchHigh);
					    patchHigh*=norm;
					    R2*=norm;
//						Image<double> save;
//						save()=patchLow;
//						save.write("PPPlow.vol");
//						save()=patchHigh;
//						save.write("PPPhigh.vol");
//						std::cout << "R2=" << R2 << " alpha=" << alpha << std::endl;
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
	Image<double> save;
	save()=weightHigh;
	save.write("PPPweightHigh.vol");
	save()=mVhigh;
	save.write("PPPVHigh.vol");
	FOR_ALL_ELEMENTS_IN_ARRAY3D(mVhigh)
	if (A3D_ELEM(weightHigh,k,i,j)>0)
		A3D_ELEM(mVhigh,k,i,j)/=A3D_ELEM(weightHigh,k,i,j);

	Vhigh.write(fnOut);
}

//#define DEBUG
void ProgRestoreWithPDBDictionary::selectDictionaryPatches(const MultidimArray<double> &lowResolutionPatch,
		Matrix1D<double> &lowResolutionPatchSignature, std::vector<size_t> &selectedPatchesIdx, std::vector<double> &weight)
{
	selectedPatchesIdx.clear();
	weight.clear();
	size_t imax=dictionaryLow.size();
	for (size_t i=0; i<imax; ++i)
	{
		const Matrix1D<double> &dictSignature=dictionarySignature[i];
		double dotProduct=XX(lowResolutionPatchSignature)*XX(dictSignature)+
                          YY(lowResolutionPatchSignature)*YY(dictSignature)+
                          ZZ(lowResolutionPatchSignature)*ZZ(dictSignature);
		if (dotProduct>angleThreshold)
		{
			const double *ptrDictionaryPatch=MULTIDIM_ARRAY(dictionaryLow[i]);
			dotProduct=0;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(lowResolutionPatch)
				dotProduct+=(*ptrDictionaryPatch++)*DIRECT_MULTIDIM_ELEM(lowResolutionPatch,n);
			if (dotProduct>angleThreshold)
			{
				if (dotProduct>1)
					dotProduct=1;
				selectedPatchesIdx.push_back(i);
				weight.push_back(1/dotProduct);
	#ifdef DEBUG
				Image<double> save;
				save()=lowResolutionPatch;
				save.write("PPPexperimentalPatch.vol");
				memcpy(&DIRECT_MULTIDIM_ELEM(save(),0),MULTIDIM_ARRAY(dictionaryLow[i]),
						MULTIDIM_SIZE(lowResolutionPatch)*sizeof(double));
				save.write("PPPdictionaryPatch.vol");
				std::cout << RAD2DEG(acos(dotProduct)) << ". Press any key" << std::endl;
				char c; std::cin >> c;
	#endif
			}
		}
	}
#ifdef DEBUG
	std::cout << "Size idx=" << selectedPatchesIdx.size() << std::endl;
	std::cout << "Weights=";
	for (size_t i=0; i<weight.size(); ++i)
		std::cout << weight[i] << " ";
	std::cout << std::endl;
#endif
}
#undef DEBUG

/* See Algorithm 1 in Trinh, Luong, Dibos, Rocchisani, Pham, Nguyen. Novel Exanple-based method for super-resolution
 *    and denoising of medical images. IEEE Trans. Image Processing, 23: 1882-1895 (2014)
 */
double ProgRestoreWithPDBDictionary::approximatePatch(const MultidimArray<double> &lowResolutionPatch,
		std::vector< size_t > &selectedPatchesIdx, std::vector<double> &weight, Matrix1D<double> &alpha)
{
	alpha.initZeros(selectedPatchesIdx.size());
    wi=alpha;
	Ui.initZeros(MULTIDIM_SIZE(lowResolutionPatch),VEC_XSIZE(alpha));

	// Prepare wi
	FOR_ALL_ELEMENTS_IN_MATRIX1D(wi)
	VEC_ELEM(wi,i)=lambda*(1+weight[i]);

	// Prepare Ui
	for (size_t j=0; j<MAT_XSIZE(Ui); ++j)
	{
		const double *ptr=MULTIDIM_ARRAY(dictionaryLow[selectedPatchesIdx[j]]);
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
	for (int iter=0; iter<iterations; ++iter)
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

void ProgRestoreWithPDBDictionary::reconstructPatch(size_t idxTransf, std::vector<size_t> &selectedPatchesIdx,
		Matrix1D<double> &alpha, MultidimArray<double> &highResolutionPatch)
{
	const Matrix2D<double> &A=rotationGroup[idxTransf];
	for (size_t j=0; j<VEC_XSIZE(alpha); ++j)
	{
		applyGeometry(LINEAR,auxPatch,dictionaryHigh[selectedPatchesIdx[j]],A,IS_NOT_INV,DONT_WRAP);
		const double *ptr=MULTIDIM_ARRAY(auxPatch);
		for (size_t n=0; n<MULTIDIM_SIZE(highResolutionPatch); ++n, ++ptr)
			DIRECT_MULTIDIM_ELEM(highResolutionPatch,n)+=VEC_ELEM(alpha,j)*(*ptr);
	}
}
