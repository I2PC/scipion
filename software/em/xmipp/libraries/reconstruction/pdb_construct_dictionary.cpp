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

#include "pdb_construct_dictionary.h"
#include <data/xmipp_image_extension.h>

void ProgPDBDictionary::defineParams()
{
    // Parameters
    addParamsLine(" [--stdThreshold <s=0.25>] : Threshold on the standard deviation to include a patch");
    addParamsLine(" [--angleThreshold <s=20>] : Threshold in degrees on the angle to include a patch");
    addParamsLine(" [--regularization <lambda=0.01>] : Regularization factor for the approximation phase");
    addParamsLine(" [--iter <s=50>] : Number of iterations");
}

void ProgPDBDictionary::readParams()
{
	stdThreshold=getDoubleParam("--stdThreshold");
	angleThreshold=cos(DEG2RAD(getDoubleParam("--angleThreshold")));
	lambda=getDoubleParam("--regularization");
	iterations=getIntParam("--iter");
}

void ProgPDBDictionary::show()
{
    if (verbose)
		std::cout
		<< "Std threshold:       " << stdThreshold   << std::endl
		<< "Angle threshold:     " << angleThreshold << std::endl
		<< "Iterations:          " << iterations     << std::endl
		<< "Regularization:      " << lambda         << std::endl
		;
}

bool ProgPDBDictionary::notInDictionary(const MultidimArray<double> &candidatePatch, MultidimArray<double> &canonicalPatch,
		Matrix1D<double> &canonicalSignature, size_t &canonicalIdx)
{
	canonicalIdx=canonicalOrientation(candidatePatch,canonicalPatch,canonicalSignature);
	size_t imax=dictionaryLow.size();

	for (size_t i=0; i<imax; ++i)
	{
		const Matrix1D<double> &dictSignature=dictionarySignature[i];
		double dotProduct=XX(canonicalSignature)*XX(dictSignature)+
                          YY(canonicalSignature)*YY(dictSignature)+
                          ZZ(canonicalSignature)*ZZ(dictSignature);
		if (dotProduct>angleThreshold)
		{
			const MultidimArray<double> &dictionaryPatch=dictionaryLow[i];
			dotProduct=0;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(dictionaryPatch)
				dotProduct+=DIRECT_MULTIDIM_ELEM(dictionaryPatch,n)*DIRECT_MULTIDIM_ELEM(canonicalPatch,n);
			if (dotProduct>angleThreshold)
				return false;
		}
	}
	return true;
}

//#define DEBUG
void ProgPDBDictionary::selectDictionaryPatches(const MultidimArray<double> &lowResolutionPatch,
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
double ProgPDBDictionary::approximatePatch(const MultidimArray<double> &lowResolutionPatch,
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
	double ymean=y.computeMean();
	FOR_ALL_ELEMENTS_IN_MATRIX1D(y)
	{
		double diff=VEC_ELEM(y,i)-VEC_ELEM(yp,i);
		diff2+=diff*diff;
		diff=VEC_ELEM(y,i)-ymean;
		norm2+=diff*diff;
	}
	return 1-diff2/norm2;
}

void ProgPDBDictionary::reconstructPatch(size_t idxTransf, std::vector<size_t> &selectedPatchesIdx,
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


void ProgPDBDictionary::loadDictionaries()
{
	std::cout << "Loading dictionaries .." << std::endl;
	size_t Xdim, Ydim, Zdim, Ndim;
	FileName fnLow=fnRoot+"_low.mrcs";
	FileName fnHigh=fnRoot+"_high.mrcs";
	FileName fnSignature=fnRoot+"_signature.raw";
	getImageSize(fnLow, Xdim, Ydim, Zdim, Ndim);
	patchSize=Xdim;
	dictionaryLow.reserve(Ndim);
	dictionaryHigh.reserve(Ndim);
	FILE *fhSignature=fopen(fnSignature.c_str(),"r");
	Image<double> aux;
	Matrix1D<double> auxSignature(3);
	for (size_t n=0; n<Ndim; ++n)
	{
		aux.read(fnLow,DATA,n+1);
		dictionaryLow.push_back(aux());
		aux.read(fnHigh,DATA,n+1);
		dictionaryHigh.push_back(aux());
		size_t sizeRead=fread(MATRIX1D_ARRAY(auxSignature),sizeof(double),3,fhSignature);
		if (sizeRead==3)
			dictionarySignature.push_back(auxSignature);
	}
	fclose(fhSignature);
	if (dictionarySignature.size()!=dictionaryLow.size())
		REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,"The dictionary is corrupted: the number of signatures does not correspond to the number of patches");
}

void ProgPDBDictionary::saveDictionaries() const
{
	std::cout << "Saving dictionaries ..." << std::endl;
	size_t imax=dictionaryLow.size();
	FileName fnLow=fnRoot+"_low.mrcs";
	FileName fnHigh=fnRoot+"_high.mrcs";
	FileName fnSignature=fnRoot+"_signature.raw";
	createEmptyFile(fnLow,patchSize,patchSize,patchSize,imax,true);
	createEmptyFile(fnHigh,patchSize,patchSize,patchSize,imax,true);
	Image<double> aux;
	FILE *fhSignature=fopen(fnSignature.c_str(),"wb");
	for (size_t i=0; i<imax; ++i)
	{
		aux()=dictionaryLow[i];
		aux.write(fnLow,i+1,true,WRITE_REPLACE);
		aux()=dictionaryHigh[i];
		aux.write(fnHigh,i+1,true,WRITE_REPLACE);
		fwrite(MATRIX1D_ARRAY(dictionarySignature[i]),sizeof(double),3,fhSignature);
	}
	fclose(fhSignature);
}

void ProgPDBDictionary::constructRotationGroup()
{
	Matrix2D<double> A(4,4), B;
	A.initIdentity();
	rotationGroup.push_back(A);

	A.initZeros();
	A(0,1)=A(2,2)=1; A(1,0)=-1; A(3,3)=1;
	rotationGroup.push_back(A);

	A.initZeros();
	A(2,0)=1; A(1,1)=1; A(0,2)=1; A(3,3)=1;
	rotationGroup.push_back(A);

	A.initZeros();
	A(0,0)=1; A(1,2)=-1; A(2,1)=1; A(3,3)=1;
	rotationGroup.push_back(A);

	A.initZeros();
	A(0,0)=-1; A(2,2)=1; A(1,1)=1; A(3,3)=1;
	rotationGroup.push_back(A);

	A.initZeros();
	A(0,0)=1; A(1,1)=-1; A(2,2)=1; A(3,3)=1;
	rotationGroup.push_back(A);

	bool finished;
	do
	{
		finished=true;
		for (size_t i=0; i<rotationGroup.size(); ++i)
			for (size_t j=0; j<rotationGroup.size(); ++j)
			{
				B=rotationGroup[i]*rotationGroup[j];
				bool found=false;
				for (size_t k=0; k<rotationGroup.size(); ++k)
					if (rotationGroup[k].equal(B))
					{
						found=true;
						break;
					}
				if (!found)
				{
					rotationGroup.push_back(B);
					finished=false;
				}
			}
	} while (!finished);
}

size_t ProgPDBDictionary::canonicalOrientation(const MultidimArray<double> &patch, MultidimArray<double> &canonicalPatch,
		Matrix1D<double> &patchSignature)
{
	patchSignature.resizeNoCopy(3);
	size_t imax=rotationGroup.size();
	double bestMoment=-1e38;
	size_t bestIdx=0;
	for (size_t i=0; i<imax; ++i)
	{
		applyGeometry(LINEAR,auxPatch,patch,rotationGroup[i],IS_INV,DONT_WRAP);
		// Calculate gradients
		double momentX=0, momentY=0, momentZ=0;
		auxPatch.setXmippOrigin();
		FOR_ALL_ELEMENTS_IN_ARRAY3D(auxPatch)
		{
			double val=A3D_ELEM(auxPatch,k,i,j);
			momentX+=j*val;
			momentY+=i*val;
			momentZ+=k*val;
		}
		double moment=momentX+momentY+momentZ;
		if (moment>bestMoment)
		{
			bestMoment=moment;
			bestIdx=i;
			canonicalPatch=auxPatch;
			XX(patchSignature)=momentX;
			YY(patchSignature)=momentY;
			ZZ(patchSignature)=momentZ;
		}
	}
	patchSignature.selfNormalize();
	return bestIdx;
}

void ProgConstructPDBDictionary::defineParams()
{
    // Usage
    addUsageLine("This program takes a set of PDB files at low and high resolution and constructs a dictionary for them.");

    // Parameters
    addParamsLine(" --low <metadata>      : Metadata with the low resolution volumes");
    addParamsLine(" --high <metadata>     : Metadata with the high resolution volumes");
    addParamsLine(" [--oroot <root=dictionary>] : Rootname for the output files");
    addParamsLine(" [--patchSize <n=5>] : Patch size for the dictionary");
    ProgPDBDictionary::defineParams();
}

void ProgConstructPDBDictionary::readParams()
{
	fnLow=getParam("--low");
	fnHigh=getParam("--high");
	fnRoot=getParam("--oroot");
	patchSize=getIntParam("--patchSize");
	ProgPDBDictionary::readParams();
}

void ProgConstructPDBDictionary::show()
{
    if (verbose)
    {
		std::cout
		<< "Input low volumes:   " << fnLow     << std::endl
		<< "Input high volumes:  " << fnHigh    << std::endl
		<< "Output rootname:     " << fnRoot    << std::endl
		<< "Patch size:          " << patchSize << std::endl
		;
		ProgPDBDictionary::show();
    }
}

//#define DEBUG
void ProgConstructPDBDictionary::run()
{
    show();
    if (fileExists(fnRoot+"_low.mrcs"))
    	loadDictionaries();

    MetaData mdlow, mdhigh;
    mdlow.read(fnLow);
    mdhigh.read(fnHigh);
    constructRotationGroup();

    FileName fnVol;
    Image<double> Vlow, Vhigh;
#ifdef DEBUG
    Image<double>Vmask;
#endif
    MultidimArray<double> patchLow, patchHigh, canonicalPatch;
    Matrix1D<double> alpha, canonicalSignature;
    int patchSize_2=patchSize/2;
    size_t canonicalIdx;

    patchLow.resize(patchSize_2,patchSize_2,patchSize_2);
    patchLow.setXmippOrigin();

    FOR_ALL_OBJECTS_IN_METADATA2(mdlow,mdhigh)
    {
    	// Read the low and high resolution volumes
    	mdlow.getValue(MDL_IMAGE,fnVol,__iter.objId);
    	Vlow.read(fnVol);
    	std::cout << "Processing " << fnVol << " and ";
    	mdhigh.getValue(MDL_IMAGE,fnVol,__iter2.objId);
    	std::cout << fnVol << std::endl;
    	Vhigh.read(fnVol);

    	// Go through the volumes and decide whether to move into the dictionary
    	const MultidimArray<double> &mVlow=Vlow();
    	const MultidimArray<double> &mVhigh=Vhigh();

    	double minLow, maxLow, meanLow, stdLow=0;
    	double minHigh, maxHigh, meanHigh, stdHigh=0;
    	mVlow.computeStats(meanLow,stdLow,minLow,maxLow);
    	mVhigh.computeStats(meanHigh,stdHigh,minHigh,maxHigh);
#ifdef DEBUG
    	Vmask().initZeros(Vhigh());
#endif

    	std::vector< size_t > selectedPatchesIdx;
    	std::vector<double> weight;
    	size_t Npatches=0, NcandidatePatches=0, NsuccessfulPatches=0;
    	init_progress_bar(ZSIZE(mVlow));
        for (int k=patchSize_2; k<(int)ZSIZE(mVlow)-patchSize_2; ++k)
        {
            for (int i=patchSize_2; i<(int)ZSIZE(mVlow)-patchSize_2; ++i)
                for (int j=patchSize_2; j<(int)ZSIZE(mVlow)-patchSize_2; ++j)
                {
                	 ++Npatches;
                     mVlow.window(patchLow,k,i,j,k+patchSize-1,i+patchSize-1,j+patchSize-1);
                     mVhigh.window(patchHigh,k,i,j,k+patchSize-1,i+patchSize-1,j+patchSize-1);

                 	double minPatchLow, maxPatchLow, meanPatchLow, stdPatchLow=0;
                 	double minPatchHigh, maxPatchHigh, meanPatchHigh, stdPatchHigh=0;
                 	patchLow.computeStats(meanPatchLow,stdPatchLow,minPatchLow,maxPatchLow);
                 	patchHigh.computeStats(meanPatchHigh,stdPatchHigh,minPatchHigh,maxPatchHigh);

                 	if (stdPatchLow > stdThreshold*stdLow && stdPatchHigh > stdThreshold*stdHigh)
                 	{
                 		++NcandidatePatches;

                 		// Candidate patch
                 		double inormPatchLow=1.0/sqrt(patchLow.sum2());
                 		patchLow*=inormPatchLow;
    					size_t canonicalIdx=canonicalOrientation(patchLow,canonicalPatch,canonicalSignature);
    					selectDictionaryPatches(canonicalPatch, canonicalSignature, selectedPatchesIdx, weight);
    					bool introduceInDictionary=(selectedPatchesIdx.size()==0);
    					if (!introduceInDictionary)
    					{
    						double R2=approximatePatch(canonicalPatch,selectedPatchesIdx,weight,alpha);
    						introduceInDictionary=(R2<angleThreshold);
    					}
                 		if (introduceInDictionary)
                 		{
                 			++NsuccessfulPatches;
                 			selfApplyGeometry(LINEAR,patchHigh,rotationGroup[canonicalIdx],IS_INV,DONT_WRAP);

                 			dictionaryLow.push_back(canonicalPatch);
                     		patchHigh*=inormPatchLow;
                 			dictionaryHigh.push_back(patchHigh);
                 			dictionarySignature.push_back(canonicalSignature);
                 		}

#ifdef DEBUG
                 		Vmask(k,i,j)=1;
#endif
                 	}
                }
            progress_bar(k);
        }
        progress_bar(ZSIZE(mVlow));
        std::cout << "Candidate patches =" << NcandidatePatches << "(" << ((double)NcandidatePatches)/Npatches*100 << "%)"
        		  << "  successful=" << NsuccessfulPatches << "(" << ((double)NsuccessfulPatches)/NcandidatePatches*100 << "%)" << std::endl;
        saveDictionaries();
#ifdef DEBUG
        Vmask.write("PPPmask.vol");
		std::cout << "Press any key\n";
		char c; std::cin >> c;
#endif
    }
}
#undef DEBUG
