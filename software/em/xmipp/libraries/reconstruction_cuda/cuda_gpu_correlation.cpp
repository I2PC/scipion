
//Host includes
#include "cuda_gpu_correlation.h"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <algorithm>


//CUDA includes
#include "cuda_basic_math.h"
#include "cuda_utils.h"

#define PI 3.14159265

__global__ void pointwiseMultiplicationComplexOneManyKernel(cufftDoubleComplex *M, cufftDoubleComplex *manyF, cufftDoubleComplex *MmanyF,
		int nzyxdim, int yxdim)
{
	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned long int idxLow = idx % yxdim;

	if (idx>=nzyxdim)
		return;

	float normFactor = (1.0/yxdim);

	cuDoubleComplex mulOut = cuCmul(manyF[idx], M[idxLow]);

	MmanyF[idx] = make_cuDoubleComplex( cuCreal(mulOut)*normFactor ,  cuCimag(mulOut)*normFactor ) ;

}

__global__ void calculateDenomFunctionKernel(double *MFrealSpace, double *MF2realSpace, double *maskAutocorrelation, double *out,
		int nzyxdim, int yxdim)
{
	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned long int idxLow = idx%yxdim;

	if (idx>=nzyxdim)
		return;

	out[idx] = sqrt(MF2realSpace[idx] - (MFrealSpace[idx]*MFrealSpace[idx]/maskAutocorrelation[idxLow]));
}


__global__ void calculateNccKernel(double *RefExpRealSpace, double *MFrealSpaceRef, double *MFrealSpaceExp, double *denomRef, double *denomExp,
		double *mask, double *NCC, size_t nzyxdim, size_t yxdim, size_t numExp, size_t maskCount)
{

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned long int idxLow = idx % yxdim;

	if(idx>=nzyxdim)
		return;

	int idxDimProj = 0;
	int outIdx = 0;
	int inIdx = 0;

	for(int n=0; n<numExp; n++)
	{
		double num = (RefExpRealSpace[idx + idxDimProj] - ((MFrealSpaceRef[idx]*MFrealSpaceExp[idxLow + inIdx])/(mask[idxLow])) );
		double den1 = denomRef[idx];
		double den2 = denomExp[idxLow + inIdx];

		if(den1!=0.0 && den2!=0.0 && !isnan(den1) && !isnan(den2) && mask[idxLow]>maskCount*0.1)
			NCC[idx + idxDimProj] = num/(den1*den2);
		else
			NCC[idx + idxDimProj] = -1;

		outIdx+=idxDimProj;
		inIdx+=yxdim;
		idxDimProj+=nzyxdim;
	}

}



__global__ void calculateNccRotationKernel(double *RefExpRealSpace, cufftDoubleComplex *polarFFTRef, cufftDoubleComplex *polarFFTExp,
		cufftDoubleComplex *polarSquaredFFTRef, cufftDoubleComplex *polarSquaredFFTExp,	cufftDoubleComplex *maskFFTPolar, double *NCC,
		size_t yxdimFFT, size_t numExp, size_t nzyxdim, size_t yxdim)
{

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=nzyxdim)
		return;

	unsigned long int idxRef = (int)(idx/yxdim)*yxdimFFT;
	unsigned long int idxDimProj = 0;
	unsigned long int idxExp = 0;

	double normValue = 1.0/yxdimFFT;

	double M1M2Polar = cuCreal(cuCmul(maskFFTPolar[0],maskFFTPolar[0]))*normValue;
	double polarValRef = cuCreal(cuCmul(polarFFTRef[idxRef],maskFFTPolar[0]))*normValue;
	double polarSqValRef = cuCreal(cuCmul(polarSquaredFFTRef[idxRef],maskFFTPolar[0]))*normValue;

	for(int n=0; n<numExp; n++)
	{
		double polarValExp = cuCreal(cuCmul(polarFFTExp[idxExp],maskFFTPolar[0]))*normValue;
		double polarSqValExp = cuCreal(cuCmul(polarSquaredFFTExp[idxExp],maskFFTPolar[0]))*normValue;

		double num = (RefExpRealSpace[idx + idxDimProj] - (polarValRef*polarValExp/M1M2Polar) );
		double den1 = sqrt(polarSqValRef - (polarValRef*polarValRef/M1M2Polar) );
		double den2 = sqrt(polarSqValExp - (polarValExp*polarValExp/M1M2Polar) );

		if(den1!=0.0 && den2!=0.0 && !isnan(den1) && !isnan(den2))
			NCC[idx + idxDimProj] = num/(den1*den2);
		else
			NCC[idx + idxDimProj] = -1;

		idxExp+=yxdimFFT;
		idxDimProj+=nzyxdim;
	}

}




__global__ void pointwiseMultiplicationComplexKernel(cufftDoubleComplex *reference, cufftDoubleComplex *experimental, cufftDoubleComplex *RefExpFourier,
		size_t nzyxdim, size_t yxdim, size_t numExp)
{
	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned long int idxLow = idx % yxdim;

	if(idx>=nzyxdim)
		return;

	double normFactor = (1.0/yxdim);
	int outIdx = 0;
	int inIdx = 0;

	for(int n=0; n<numExp; n++){

		cuDoubleComplex mulOut = cuCmul(reference[idx], experimental[idxLow + inIdx]);
		RefExpFourier[idx + outIdx] = make_cuDoubleComplex( cuCreal(mulOut)*normFactor ,  cuCimag(mulOut)*normFactor ) ;

		outIdx+=nzyxdim;
		inIdx+=yxdim;
	}

}


void pointwiseMultiplicationFourier(const GpuMultidimArrayAtGpu< std::complex<double> > &M, const GpuMultidimArrayAtGpu < std::complex<double> >& manyF,
		GpuMultidimArrayAtGpu< std::complex<double> > &MmanyF)
{
    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    manyF.calculateGridSizeVectorized(blockSize, gridSize);

    pointwiseMultiplicationComplexOneManyKernel <<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			((cufftDoubleComplex*)M.d_data, (cufftDoubleComplex*)manyF.d_data, (cufftDoubleComplex*) MmanyF.d_data, manyF.nzyxdim, manyF.yxdim);

}

void calculateDenomFunction(const GpuMultidimArrayAtGpu< double > &MFrealSpace, const GpuMultidimArrayAtGpu < double >& MF2realSpace,
		const GpuMultidimArrayAtGpu < double >& maskAutocorrelation, GpuMultidimArrayAtGpu< double > &out)
{
    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    MFrealSpace.calculateGridSizeVectorized(blockSize, gridSize);

    calculateDenomFunctionKernel <<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			(MFrealSpace.d_data, MF2realSpace.d_data, maskAutocorrelation.d_data, out.d_data, (int)MFrealSpace.nzyxdim, (int)MFrealSpace.yxdim);

}

void GpuCorrelationAux::produceSideInfo()
{
	GpuMultidimArrayAtGpu< std::complex<double> > MF, MF2;
	MF.resize(d_projFFT);
	MF2.resize(d_projSquaredFFT);
	pointwiseMultiplicationFourier(d_maskFFT, d_projFFT, MF);
	pointwiseMultiplicationFourier(d_maskFFT, d_projSquaredFFT, MF2);
	d_projSquaredFFT.clear();

	GpuMultidimArrayAtGpu<double> MF2realSpace(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	MFrealSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	MF.ifft(MFrealSpace);
	MF2.ifft(MF2realSpace);
	MF.clear();
	MF2.clear();

	GpuMultidimArrayAtGpu< std::complex<double> > maskAux(d_projFFT.Xdim, d_projFFT.Ydim);
	pointwiseMultiplicationFourier(d_maskFFT, d_maskFFT, maskAux);
	maskAutocorrelation.resize(Xdim, Ydim);
	maskAux.ifft(maskAutocorrelation);
	maskAux.clear();
	d_denom.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	calculateDenomFunction(MFrealSpace, MF2realSpace, maskAutocorrelation, d_denom);
	MF2realSpace.clear();
}

double** cuda_calculate_correlation_rotation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux)
{
	GpuMultidimArrayAtGpu< std::complex<double> > RefExpFourier(referenceAux.d_projPolarFFT.Xdim, referenceAux.d_projPolarFFT.Ydim,
			referenceAux.d_projPolarFFT.Zdim, referenceAux.d_projPolarFFT.Ndim*experimentalAux.d_projPolarFFT.Ndim);

    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    referenceAux.d_projPolarFFT.calculateGridSizeVectorized(blockSize, gridSize);

    pointwiseMultiplicationComplexKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			((cufftDoubleComplex*)referenceAux.d_projPolarFFT.d_data, (cufftDoubleComplex*)experimentalAux.d_projPolarFFT.d_data,
					(cufftDoubleComplex*)RefExpFourier.d_data, referenceAux.d_projPolarFFT.nzyxdim, referenceAux.d_projPolarFFT.yxdim,
					experimentalAux.d_projPolarFFT.Ndim);

    GpuMultidimArrayAtGpu< double > RefExpRealSpace(referenceAux.XdimPolar, referenceAux.YdimPolar, referenceAux.d_projPolarFFT.Zdim,
    		referenceAux.d_projPolarFFT.Ndim*experimentalAux.d_projPolarFFT.Ndim);
    RefExpFourier.ifft(RefExpRealSpace);
    RefExpFourier.clear();

	unsigned long int tamTotal = referenceAux.XdimPolar*referenceAux.YdimPolar*referenceAux.d_projPolarFFT.Ndim;

 	//XmippDim3 blockSize2(numTh, 1, 1), gridSize2;
 	//referenceAux.MFrealSpacePolar.calculateGridSizeVectorized(blockSize2, gridSize2); //AJ esto es un lio

	const dim3 blockSize2(numTh, 1, 1);
    int numBlkx = (int)(tamTotal)/numTh;
    if((int)(tamTotal)%numTh>0){
    	numBlkx++;
    }
    const dim3 gridSize2(numBlkx, 1, 1);

 	GpuMultidimArrayAtGpu<double> d_NCC(referenceAux.XdimPolar, referenceAux.YdimPolar, referenceAux.d_projPolarFFT.Zdim,
				referenceAux.d_projPolarFFT.Ndim*experimentalAux.d_projPolarFFT.Ndim);
	double **NCC = new double* [experimentalAux.d_projPolarFFT.Ndim];
	for(int n=0; n<experimentalAux.d_projPolarFFT.Ndim; n++)
		NCC[n] = new double [tamTotal];

	calculateNccRotationKernel<<< gridSize2, blockSize2 >>>
			(RefExpRealSpace.d_data, (cufftDoubleComplex*)referenceAux.d_projPolarFFT.d_data, (cufftDoubleComplex*)experimentalAux.d_projPolarFFT.d_data,
					(cufftDoubleComplex*)referenceAux.d_projPolarSquaredFFT.d_data, (cufftDoubleComplex*)experimentalAux.d_projPolarSquaredFFT.d_data,
					(cufftDoubleComplex*)referenceAux.d_maskFFTPolar.d_data, d_NCC.d_data, referenceAux.d_projPolarFFT.yxdim, experimentalAux.d_projPolarFFT.Ndim,
					referenceAux.XdimPolar*referenceAux.YdimPolar*referenceAux.d_projPolarFFT.Ndim, referenceAux.XdimPolar*referenceAux.YdimPolar);

	int pointer=0;
	for(int n=0; n<experimentalAux.d_projPolarFFT.Ndim; n++){
		gpuErrchk(cudaMemcpy(NCC[n], &d_NCC.d_data[pointer], sizeof(double)*(tamTotal), cudaMemcpyDeviceToHost));
		pointer+=tamTotal;
	}

	return NCC;

}


double** cuda_calculate_correlation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux)
{

	GpuMultidimArrayAtGpu< std::complex<double> > RefExpFourier(referenceAux.d_projFFT.Xdim, referenceAux.d_projFFT.Ydim,
			referenceAux.d_projFFT.Zdim, referenceAux.d_projFFT.Ndim*experimentalAux.d_projFFT.Ndim);

    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    referenceAux.d_projFFT.calculateGridSizeVectorized(blockSize, gridSize);

    pointwiseMultiplicationComplexKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			((cufftDoubleComplex*)referenceAux.d_projFFT.d_data, (cufftDoubleComplex*)experimentalAux.d_projFFT.d_data, (cufftDoubleComplex*)RefExpFourier.d_data,
					referenceAux.d_projFFT.nzyxdim, referenceAux.d_projFFT.yxdim, experimentalAux.d_projFFT.Ndim);

    GpuMultidimArrayAtGpu< double > RefExpRealSpace(referenceAux.Xdim, referenceAux.Ydim, referenceAux.d_projFFT.Zdim,
    		referenceAux.d_projFFT.Ndim*experimentalAux.d_projFFT.Ndim);
    RefExpFourier.ifft(RefExpRealSpace);
    RefExpFourier.clear();

 	XmippDim3 blockSize2(numTh, 1, 1), gridSize2;
	referenceAux.MFrealSpace.calculateGridSizeVectorized(blockSize2, gridSize2); //AJ esto es un lio

	GpuMultidimArrayAtGpu<double> d_NCC(referenceAux.Xdim, referenceAux.Ydim, referenceAux.d_projFFT.Zdim,
			referenceAux.d_projFFT.Ndim*experimentalAux.d_projFFT.Ndim);
	double **NCC = new double* [experimentalAux.d_projFFT.Ndim];
	for(int n=0; n<experimentalAux.d_projFFT.Ndim; n++)
		NCC[n] = new double [referenceAux.Xdim*referenceAux.Ydim*referenceAux.d_projFFT.Ndim];


	calculateNccKernel<<< CONVERT2DIM3(gridSize2), CONVERT2DIM3(blockSize2) >>>
			(RefExpRealSpace.d_data, referenceAux.MFrealSpace.d_data, experimentalAux.MFrealSpace.d_data, referenceAux.d_denom.d_data,
					experimentalAux.d_denom.d_data, referenceAux.maskAutocorrelation.d_data, d_NCC.d_data,
					referenceAux.MFrealSpace.nzyxdim, referenceAux.MFrealSpace.yxdim, experimentalAux.MFrealSpace.Ndim, referenceAux.maskCount);

	int pointer=0;
	for(int n=0; n<experimentalAux.d_projFFT.Ndim; n++){
		gpuErrchk(cudaMemcpy(NCC[n], &d_NCC.d_data[pointer], sizeof(double)*(referenceAux.MFrealSpace.nzyxdim), cudaMemcpyDeviceToHost));
		pointer+=referenceAux.MFrealSpace.nzyxdim;
	}

    return NCC;
}


__global__ void cart2polar(double *image, double *polar, double *polar2, int maxRadius, int maxAng, int Nimgs, int Ydim, int Xdim, bool rotate)
{
	int angle = blockDim.x * blockIdx.x + threadIdx.x;
	int radius = blockDim.y * blockIdx.y + threadIdx.y;
	//int batch = blockDim.z * blockIdx.z + threadIdx.z;

	if (radius>=maxRadius || angle>=maxAng)
		return;

	double x = (double)(radius*cosf((float)(angle*PI/180))) + Xdim/2;
	double y = (double)(radius*sinf((float)(angle*PI/180))) + Ydim/2;

	double dx_low = floor(x);
	double dy_low = floor(y);
	int x_low = (int)dx_low;
	int y_low = (int)dy_low;
	double x_x_low=x-dx_low;
	double y_y_low=y-dy_low;
	double one_x=1.0-x_x_low;
	double one_y=1.0-y_y_low;
	double w00=one_y*one_x;
	double w01=one_y*x_x_low;
	double w10=y_y_low*one_x;
	double w11=y_y_low*x_x_low;

	int NXY=Xdim*Ydim;
	int NXYpolar=maxAng*maxRadius;
	int imgIdx00=y_low * Xdim + x_low;
	int imgIdx01=imgIdx00+1;
	int imgIdx10=imgIdx00+Xdim;
	int imgIdx11=imgIdx10+1;
	int imgOffset=0;
	int polarOffset=0;
	int polarIdx;
	if(!rotate)
		polarIdx=angle+(radius*maxAng);
	else
		polarIdx = (maxAng-angle-1)+((maxRadius-radius-1)*maxAng);

	for (int n=0; n<Nimgs; n++)
	{
		double I00 = image[imgIdx00+imgOffset];
		double I01 = image[imgIdx01+imgOffset];
		double I10 = image[imgIdx10+imgOffset];
		double I11 = image[imgIdx11+imgOffset];
		double imVal = I00*w00 + I01*w01 + I10*w10 + I11*w11;
		int finalPolarIndex=polarIdx+polarOffset;
		polar[finalPolarIndex] = imVal;
		polar2[finalPolarIndex] = imVal*imVal;

		imgOffset+=NXY;
		polarOffset+=NXYpolar;
	}

}

void cuda_cart2polar(GpuMultidimArrayAtGpu<double> &image, GpuMultidimArrayAtGpu<double> &polar_image, GpuMultidimArrayAtGpu<double> &polar2_image, bool rotate)
{
    int numTh = 32;
    XmippDim3 blockSize(numTh, numTh, 1), gridSize;
    polar_image.calculateGridSize(blockSize, gridSize);
    cart2polar <<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
    		(image.d_data, polar_image.d_data, polar2_image.d_data, polar_image.Ydim, polar_image.Xdim, polar_image.Ndim, image.Ydim, image.Xdim, rotate);
}
