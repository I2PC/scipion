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
#include <data/numerical_tools.h>

void ProgPDBDictionary::defineParams()
{
    // Parameters
    addParamsLine(" [--stdThreshold <s=0.25>] : Threshold on the standard deviation to include a patch");
    addParamsLine(" [--angleThreshold <s=20>] : Threshold in degrees on the angle to consider a patch");
    addParamsLine(" [--regularization <lambda=0.01>] : Regularization factor for the approximation phase");
    addParamsLine(" [--iter <s=50>] : Number of iterations");
    addParamsLine(" [--mode <m=\"3D\">] : 3D, 2Dx, 2Dy, 2Dz mode");
}

void ProgPDBDictionary::readParams()
{
	stdThreshold=getDoubleParam("--stdThreshold");
	angleThreshold=cos(DEG2RAD(getDoubleParam("--angleThreshold")));
	lambda=getDoubleParam("--regularization");
	iterations=getIntParam("--iter");
	String modeStr=getParam("--mode");
	if (modeStr=="3D")
		mode=0;
	else if (modeStr=="2Dx")
		mode=1;
	else if (modeStr=="2Dy")
		mode=2;
	else if (modeStr=="2Dz")
		mode=3;
}

void ProgPDBDictionary::show()
{
    if (verbose)
		std::cout
		<< "Std threshold:       " << stdThreshold   << std::endl
		<< "Angle threshold:     " << angleThreshold << std::endl
		<< "Iterations:          " << iterations     << std::endl
		<< "Regularization:      " << lambda         << std::endl
		<< "Mode:                " << mode           << std::endl
		;
}

bool ProgPDBDictionary::notInDictionary(const MultidimArray<double> &candidatePatch, MultidimArray<double> &canonicalPatch,
		Matrix1D<double> &canonicalSignature, size_t &canonicalIdx)
{
	if (mode==0)
		canonicalIdx=canonicalOrientation3D(candidatePatch,canonicalPatch,canonicalSignature);
	else
		canonicalIdx=canonicalOrientation2D(candidatePatch,canonicalPatch,canonicalSignature);
	size_t imax=dictionaryLow.size();

	for (size_t i=0; i<imax; ++i)
	{
		const Matrix1D<double> &dictSignature=dictionarySignature[i];
		double dotProduct=0;
		for (size_t j=0; j<VEC_XSIZE(dictSignature); ++j)
		   dotProduct+=VEC_ELEM(canonicalSignature,j)*VEC_ELEM(dictSignature,j);
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
		double dotProduct=0;
		for (size_t j=0; j<VEC_XSIZE(dictSignature); ++j)
		   dotProduct+=VEC_ELEM(lowResolutionPatchSignature,j)*VEC_ELEM(dictSignature,j);
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
#ifdef NEVERDEFINED
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
	double R2=1-diff2/norm2;
	if (selectedPatchesIdx.size()>30 && R2<0.975)
	{
		Ui.write("PPPUi.txt");
		y.write("PPPy.txt");
		yp.write("PPPyp.txt");
		std::cout << "A strange case detected" << std::endl;
		char c; std::cin >> c;
	}
	return R2;
}
#endif

double ProgPDBDictionary::approximatePatch(const MultidimArray<double> &lowResolutionPatch,
		std::vector< size_t > &selectedPatchesIdx, std::vector<double> &weight, Matrix1D<double> &alpha)
{
	Ui.initZeros(MULTIDIM_SIZE(lowResolutionPatch),selectedPatchesIdx.size());

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
	matrixOperation_Atx(Ui,y,Uity);

	// Solve Least Squares problem
	solveBySVD(UitUi,Uity,alpha,1e-6);

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
	double R2=1-diff2/norm2;
	return R2;
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
	if (Zdim==1)
		auxSignature.resizeNoCopy(4);
	for (size_t n=0; n<Ndim; ++n)
	{
		aux.read(fnLow,DATA,n+1);
		dictionaryLow.push_back(aux());
		aux.read(fnHigh,DATA,n+1);
		dictionaryHigh.push_back(aux());
		size_t sizeRead=fread(MATRIX1D_ARRAY(auxSignature),sizeof(double),VEC_XSIZE(auxSignature),fhSignature);
		if (sizeRead==VEC_XSIZE(auxSignature))
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
	int Zdim;
	if (mode==0)
		Zdim=patchSize;
	else
		Zdim=1;
	createEmptyFile(fnLow,patchSize,patchSize,Zdim,imax,true);
	createEmptyFile(fnHigh,patchSize,patchSize,Zdim,imax,true);
	Image<double> aux;
	FILE *fhSignature=fopen(fnSignature.c_str(),"wb");
	for (size_t i=0; i<imax; ++i)
	{
		aux()=dictionaryLow[i];
		aux.write(fnLow,i+1,true,WRITE_REPLACE);
		aux()=dictionaryHigh[i];
		aux.write(fnHigh,i+1,true,WRITE_REPLACE);
		fwrite(MATRIX1D_ARRAY(dictionarySignature[i]),sizeof(double),VEC_XSIZE(dictionarySignature[i]),fhSignature);
	}
	fclose(fhSignature);
}

void ProgPDBDictionary::extractPatch(const MultidimArray<double> &V, MultidimArray<double> &patch, int k, int i, int j)
{
	if (mode==0)
		V.window(patch,k-patchSize_2,i-patchSize_2,j-patchSize_2,k+patchSize_2,i+patchSize_2,j+patchSize_2);
	else if (mode==1)
	{
		for (int kk=-patchSize_2; kk<=patchSize_2; ++kk)
			for (int ii=-patchSize_2; ii<=patchSize_2; ++ii)
				A2D_ELEM(patch,kk,ii)=A3D_ELEM(V,k+kk,i+ii,j);
	}
	else if (mode==2)
	{
		for (int kk=-patchSize_2; kk<=patchSize_2; ++kk)
			for (int jj=-patchSize_2; jj<=patchSize_2; ++jj)
				A2D_ELEM(patch,kk,jj)=A3D_ELEM(V,k+kk,i,j+jj);
	}
	else if (mode==3)
	{
		for (int ii=-patchSize_2; ii<=patchSize_2; ++ii)
			for (int jj=-patchSize_2; jj<=patchSize_2; ++jj)
				A2D_ELEM(patch,ii,jj)=A3D_ELEM(V,k,i+ii,j+jj);
	}
}

void ProgPDBDictionary::insertPatch(MultidimArray<double> &Vhigh, MultidimArray<double> &weightHigh, const MultidimArray<double> &patchHigh,
		int k, int i, int j, double R2)
{
	if (mode==0)
	{
		for (int kk=-patchSize_2; kk<=patchSize_2; ++kk)
			for (int ii=-patchSize_2; ii<=patchSize_2; ++ii)
				for (int jj=-patchSize_2; jj<=patchSize_2; ++jj)
				{
					A3D_ELEM(Vhigh,k+kk,i+ii,j+jj)+=R2*A3D_ELEM(patchHigh,kk,ii,jj);
					A3D_ELEM(weightHigh,k+kk,i+ii,j+jj)+=R2;
				}
	}
	else if (mode==1)
	{
		for (int kk=-patchSize_2; kk<=patchSize_2; ++kk)
			for (int ii=-patchSize_2; ii<=patchSize_2; ++ii)
			{
				A3D_ELEM(Vhigh,k+kk,i+ii,j)+=R2*A2D_ELEM(patchHigh,kk,ii);
				A3D_ELEM(weightHigh,k+kk,i+ii,j)+=R2;
			}
	}
	else if (mode==2)
	{
		for (int kk=-patchSize_2; kk<=patchSize_2; ++kk)
			for (int jj=-patchSize_2; jj<=patchSize_2; ++jj)
			{
				A3D_ELEM(Vhigh,k+kk,i,j+jj)+=R2*A2D_ELEM(patchHigh,kk,jj);
				A3D_ELEM(weightHigh,k+kk,i,j+jj)+=R2;
			}
	}
	else if (mode==3)
	{
		for (int ii=-patchSize_2; ii<=patchSize_2; ++ii)
			for (int jj=-patchSize_2; jj<=patchSize_2; ++jj)
			{
				A3D_ELEM(Vhigh,k,i+ii,j+jj)+=R2*A2D_ELEM(patchHigh,ii,jj);
				A3D_ELEM(weightHigh,k,i+ii,j+jj)+=R2;
			}
	}
}

void ProgPDBDictionary::constructRotationGroup2D()
{
	Matrix2D<double> A(3,3);
	A.initIdentity();
	rotationGroup.push_back(A);

	A.initZeros();
	A(0,0)=-1; A(1,1)=1; A(2,2)=1;
	rotationGroup.push_back(A);

	A.initZeros();
	A(0,0)=1; A(1,1)=-1; A(2,2)=1;
	rotationGroup.push_back(A);

	A.initZeros();
	A(0,1)=1; A(1,0)=-1; A(2,2)=1;
	rotationGroup.push_back(A);

	constructRotationGroup();
}

void ProgPDBDictionary::constructRotationGroup3D()
{
	Matrix2D<double> A(4,4);
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

	constructRotationGroup();
}

void ProgPDBDictionary::constructRotationGroup()
{
	Matrix2D<double> B;
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
	std::cout << "Size of rotation group=" << rotationGroup.size() << std::endl;
}

size_t ProgPDBDictionary::canonicalOrientation2D(const MultidimArray<double> &patch, MultidimArray<double> &canonicalPatch,
		Matrix1D<double> &patchSignature)
{
	patchSignature.resizeNoCopy(4);
	size_t nmax=rotationGroup.size();
	double bestMoment=-1e38;
	size_t bestIdx=0;
	for (size_t n=0; n<nmax; ++n)
	{
		applyGeometry(LINEAR,auxPatch,patch,rotationGroup[n],IS_INV,DONT_WRAP);
		// Calculate gradients
		double momentX=0, momentY=0, momentXmY=0, momentXY=0;
		auxPatch.setXmippOrigin();
		FOR_ALL_ELEMENTS_IN_ARRAY2D(auxPatch)
		{
			double val=A2D_ELEM(auxPatch,i,j);
			double jval=j*val;
			double ival=i*val;
			momentX+=jval;
			momentY+=ival;
			momentXmY+=jval-ival;
			momentXY+=jval+ival;
		}
		double moment=momentX+momentY;
		if (moment>bestMoment)
		{
			bestMoment=moment;
			bestIdx=n;
			canonicalPatch=auxPatch;
			VEC_ELEM(patchSignature,0)=momentX;
			VEC_ELEM(patchSignature,1)=momentY;
			VEC_ELEM(patchSignature,2)=momentXmY;
			VEC_ELEM(patchSignature,3)=momentXY;
		}
	}
	patchSignature.selfNormalize();
	return bestIdx;
}

size_t ProgPDBDictionary::canonicalOrientation3D(const MultidimArray<double> &patch, MultidimArray<double> &canonicalPatch,
		Matrix1D<double> &patchSignature)
{
	patchSignature.resizeNoCopy(3);
	size_t nmax=rotationGroup.size();
	double bestMoment=-1e38;
	size_t bestIdx=0;
	for (size_t n=0; n<nmax; ++n)
	{
		applyGeometry(LINEAR,auxPatch,patch,rotationGroup[n],IS_INV,DONT_WRAP);
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
			bestIdx=n;
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
    addParamsLine(" [--R2Threshold <s=0.99>] : Threshold in R2 to include a patch");
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
	R2threshold=getDoubleParam("--R2Threshold");
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
		<< "R2 threshold:        " << R2threshold << std::endl
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
    if (mode==0)
    	constructRotationGroup3D();
    else
    	constructRotationGroup2D();

    FileName fnVol;
    Image<double> Vlow, Vhigh;
#ifdef DEBUG
    Image<double>Vmask;
#endif
    MultidimArray<double> patchLow, patchHigh, canonicalPatch;
    Matrix1D<double> alpha, canonicalSignature;
    patchSize_2=patchSize/2;
    size_t canonicalIdx;

    if (mode==0)
    	patchLow.resize(patchSize,patchSize,patchSize);
    else
    	patchLow.resize(patchSize,patchSize);
    patchLow.setXmippOrigin();
    patchHigh=patchLow;

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
                	 extractPatch(mVlow,patchLow,k,i,j);
                	 extractPatch(mVhigh,patchHigh,k,i,j);

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
    					size_t canonicalIdx=0;
    					if (mode==0)
    						canonicalIdx=canonicalOrientation3D(patchLow,canonicalPatch,canonicalSignature);
    					else
    						canonicalIdx=canonicalOrientation2D(patchLow,canonicalPatch,canonicalSignature);
    					selectDictionaryPatches(canonicalPatch, canonicalSignature, selectedPatchesIdx, weight);
    					bool introduceInDictionary=(selectedPatchesIdx.size()==0);
#ifdef DEBUG
    					std::cout << "Evaluating " << k << "," << i << "," << j << " -> selectedPatches=" << selectedPatchesIdx.size() << std::endl;
#endif
    					double R2=-1;
    					if (!introduceInDictionary)
    					{
    						R2=approximatePatch(canonicalPatch,selectedPatchesIdx,weight,alpha);
#ifdef DEBUG
    					std::cout << "   R2=" << R2 << std::endl;
#endif
    						introduceInDictionary=(R2<R2threshold);
    					}
                 		if (introduceInDictionary)
                 		{
                 			++NsuccessfulPatches;
                 			selfApplyGeometry(LINEAR,patchHigh,rotationGroup[canonicalIdx],IS_INV,DONT_WRAP);

                 			dictionaryLow.push_back(canonicalPatch);
                     		patchHigh*=inormPatchLow;
                 			dictionaryHigh.push_back(patchHigh);
                 			dictionarySignature.push_back(canonicalSignature);
#ifdef DEBUG
    					std::cout << "   Introducing it into dictionary" << std::endl;
#endif
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
