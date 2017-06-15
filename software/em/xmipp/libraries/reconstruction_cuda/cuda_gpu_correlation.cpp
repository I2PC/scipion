
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


/*
__global__ void compute_max(double *d_in, double *d_out, double *d_out_pos, int size)
{
	extern __shared__ double sdata[];

	unsigned int myId = threadIdx.x + (blockDim.x*blockIdx.x);
	unsigned int tid = threadIdx.x;

	sdata[tid]= d_in[myId];
	sdata[tid+blockDim.x] = myId;

	__syncthreads();

	if(myId >= size)
		sdata[tid] = -1;
	__syncthreads();

	for (unsigned int s= blockDim.x/2; s>0; s>>=1){
		if(tid<s){
			sdata[tid]=max(sdata[tid], sdata[tid+s]);
			sdata[tid+blockDim.x]=(sdata[tid], sdata[tid+s])?myId:myId+s;
		}
		__syncthreads();
	}

	if(tid==0){
		d_out[blockIdx.x]=sdata[0];
		d_out_pos[blockIdx.x]=sdata[blockDim.x];
	}
}
*/

__global__ void matrixMultiplication (float* newMat, float* lastMat, float* result, size_t n){

	unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if(idx>=n)
		return;

	int idx9 = idx*9;
	result[idx9] = newMat[idx9]*lastMat[idx9] + newMat[idx9+1]*lastMat[idx9+3];
	result[idx9+1] = newMat[idx9]*lastMat[idx9+1] + newMat[idx9+1]*lastMat[idx9+4];
	result[idx9+2] = newMat[idx9]*lastMat[idx9+2] + newMat[idx9+1]*lastMat[idx9+5] + newMat[idx9+2];

	result[idx9+3] = newMat[idx9+3]*lastMat[idx9] + newMat[idx9+4]*lastMat[idx9+3];
	result[idx9+4] = newMat[idx9+3]*lastMat[idx9+1] + newMat[idx9+4]*lastMat[idx9+4];
	result[idx9+5] = newMat[idx9+3]*lastMat[idx9+2] + newMat[idx9+4]*lastMat[idx9+5] + newMat[idx9+5];

}

__global__ void pointwiseMultiplicationComplexOneManyKernel(cufftDoubleComplex *M, cufftDoubleComplex *manyF, cufftDoubleComplex *MmanyF,
		int nzyxdim, int yxdim)
{
	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned long int idxLow = idx%(int)yxdim;

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
		double *mask, double *NCC, size_t nzyxdim, size_t yxdim, size_t maskCount)
{

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned long int idxLow = idx % yxdim;

	if(idx>=nzyxdim)
		return;

	double num = (RefExpRealSpace[idx] - ((MFrealSpaceRef[idx]*MFrealSpaceExp[idx])/(mask[idxLow])) );
	double den1 = denomRef[idx];
	double den2 = denomExp[idx];

	if(den1!=0.0 && den2!=0.0 && !isnan(den1) && !isnan(den2) && mask[idxLow]>maskCount*0.1)
		NCC[idx] = num/(den1*den2);
	else
		NCC[idx] = -1;

}

__device__ void wrapping (int &x, int &y, size_t xdim, size_t ydim, int wrap){

	if(wrap==0){ //image copies
		if(x<0)
			x=xdim+x;
		else if(x>=xdim)
			x=x-xdim;
		if(y<0)
			y=ydim+y;
		else if(y>=ydim)
			y=y-ydim;
	}else if(wrap==1){ //mirror
		if(x<0)
			x=-x;
		else if(x>=xdim)
			x=xdim-(x-xdim)-1;
		if(y<0)
			y=-y;
		else if(y>=ydim)
			y=ydim-(y-ydim)-1;
	}else if(wrap==2){ //last pixel copies
		if(x<0)
			x=0;
		else if(x>=xdim)
			x=xdim-1;
		if(y<0)
			y=0;
		else if(y>=ydim)
			y=ydim-1;
	}

}

__global__ void applyTransformKernel(double *d_in, double *d_out, float *transMat, size_t nzyxdim, size_t yxdim,
		size_t xdim, size_t ydim)
{

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned int idxIm = idx%yxdim;
	unsigned int numIm = idx/yxdim;

	//AJ wrapping type: 0 image copies, 1 mirror, 2 last pixel copies
	int wrap = 2;

	if(idx>=nzyxdim)
		return;

	int x = idxIm%xdim;
	int y = idxIm/xdim;
	float x_orig = 0;
	float y_orig = 0;
	if(transMat[(numIm*9)]!=0 || transMat[1+(numIm*9)]!=0 || transMat[3+(numIm*9)]!=0 ||transMat[4+(numIm*9)]!=0){
		x = x - xdim/2;
		y = y - ydim/2;
		x_orig += transMat[(numIm*9)]*x + transMat[1+(numIm*9)]*y + xdim/2;
		y_orig += transMat[3+(numIm*9)]*x + transMat[4+(numIm*9)]*y + ydim/2;
	}
	if(transMat[2+(numIm*9)]!=0 || transMat[5+(numIm*9)]!=0){
		x_orig += transMat[2+(numIm*9)];
		y_orig += transMat[5+(numIm*9)];
	}

	int x_orig00 = (int)floorf(x_orig);
	int y_orig00 = (int)floorf(y_orig);
	int x_orig01 = x_orig00+1;
	int y_orig01 = y_orig00;
	int x_orig10 = x_orig00;
	int y_orig10 = y_orig00+1;
	int x_orig11 = x_orig00+1;
	int y_orig11 = y_orig00+1;

	float x_x_low=x_orig-x_orig00;
	float y_y_low=y_orig-y_orig00;
	float one_x=1.0-x_x_low;
	float one_y=1.0-y_y_low;
	float w00=one_y*one_x;
	float w01=one_y*x_x_low;
	float w10=y_y_low*one_x;
	float w11=y_y_low*x_x_low;

	wrapping (x_orig00, y_orig00, xdim, ydim, wrap);
	wrapping (x_orig01, y_orig01, xdim, ydim, wrap);
	wrapping (x_orig10, y_orig10, xdim, ydim, wrap);
	wrapping (x_orig11, y_orig11, xdim, ydim, wrap);

	int imgIdx00=y_orig00 * xdim + x_orig00;
	int imgIdx01=y_orig01 * xdim + x_orig01;
	int imgIdx10=y_orig10 * xdim + x_orig10;
	int imgIdx11=y_orig11 * xdim + x_orig11;

	int imgOffset = numIm*yxdim;
	double I00 = d_in[imgIdx00+imgOffset];
	double I01 = d_in[imgIdx01+imgOffset];
	double I10 = d_in[imgIdx10+imgOffset];
	double I11 = d_in[imgIdx11+imgOffset];
	double imVal = I00*w00 + I01*w01 + I10*w10 + I11*w11;

	d_out[idx] = imVal;

}



__global__ void calculateNccRotationKernel(double *RefExpRealSpace, cufftDoubleComplex *polarFFTRef, cufftDoubleComplex *polarFFTExp,
		cufftDoubleComplex *polarSquaredFFTRef, cufftDoubleComplex *polarSquaredFFTExp,	double maskFFTPolarReal, double *NCC,
		size_t yxdimFFT, size_t nzyxdim, size_t yxdim)
{

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned long int idxLow = (idx/(int)yxdim)*(int)yxdimFFT;

	if(idx>=nzyxdim)
		return;

	double normValue = 1.0/yxdimFFT;

	cufftDoubleComplex maskFFTPolar = make_cuDoubleComplex(maskFFTPolarReal, 0.0);

	double M1M2Polar = cuCreal(cuCmul(maskFFTPolar,maskFFTPolar))*normValue;
	double polarValRef = cuCreal(cuCmul(polarFFTRef[idxLow],maskFFTPolar))*normValue;
	double polarSqValRef = cuCreal(cuCmul(polarSquaredFFTRef[idxLow],maskFFTPolar))*normValue;

	double polarValExp = cuCreal(cuCmul(polarFFTExp[idxLow],maskFFTPolar))*normValue;
	double polarSqValExp = cuCreal(cuCmul(polarSquaredFFTExp[idxLow],maskFFTPolar))*normValue;

	double num = (RefExpRealSpace[idx] - (polarValRef*polarValExp/M1M2Polar) );
	double den1 = sqrt(polarSqValRef - (polarValRef*polarValRef/M1M2Polar) );
	double den2 = sqrt(polarSqValExp - (polarValExp*polarValExp/M1M2Polar) );

	if(den1!=0.0 && den2!=0.0 && !isnan(den1) && !isnan(den2))
		NCC[idx] = num/(den1*den2);
	else
		NCC[idx] = -1;

}


__global__ void pointwiseMultiplicationComplexKernel(cufftDoubleComplex *reference, cufftDoubleComplex *experimental,
		cufftDoubleComplex *RefExpFourier, size_t nzyxdim, size_t yxdim)
{
	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=nzyxdim)
		return;

	double normFactor = (1.0/yxdim);

	cuDoubleComplex mulOut = cuCmul(reference[idx], experimental[idx]);
	RefExpFourier[idx] = make_cuDoubleComplex( cuCreal(mulOut)*normFactor ,  cuCimag(mulOut)*normFactor );
}


__global__ void maskingKernel(double *d_in, double *d_out, double *d_out2, double *mask,
		size_t xdim, size_t ydim, size_t yxdim, size_t numImag, bool rotation, bool experimental){

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=yxdim)
		return;

	unsigned int x_idx = idx%(int)xdim;
	unsigned int y_idx = idx/(int)xdim;
	unsigned int idxWrite;
	if(!rotation && experimental)
		idxWrite = (ydim-1 - y_idx)*xdim + (xdim-1 - x_idx);
	else
		idxWrite = y_idx*xdim + x_idx;

	int offset=0;
	for(int i=0; i<numImag; i++){
		d_out[idxWrite+offset] = d_in[idx+offset]*mask[idx];
		if(!rotation || !experimental)
			d_out2[idxWrite+offset] = d_out[idxWrite+offset]*d_out[idxWrite+offset];

		offset += yxdim;
	}
}

__global__ void paddingKernel(double *d_orig_image, double *image2_gpu, double *mask, double *padded_image_gpu,
		double *padded_image2_gpu, double *padded_mask_gpu, size_t xdim, size_t ydim, size_t yxdim,
		size_t numImag, size_t pad_xdim, size_t pad_ydim, size_t pad_yxdim, bool rotation, bool experimental){

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=yxdim)
		return;

	int xdim2Im = (int)floorf((pad_xdim-xdim)/2);
	int ydim2Im = (int)floorf((pad_ydim-ydim)/2);
	int xdim2Mask = xdim2Im;
	int ydim2Mask = ydim2Im;
	if(!rotation && experimental && xdim%2==0){
		xdim2Im+=1;
		ydim2Im+=1;
	}
	unsigned int x_idx = idx%(int)xdim;
	unsigned int y_idx = idx/(int)xdim;
	unsigned int idxWrite;
	unsigned int idxWriteMask;

	int offset=0;
	for(int j=0; j<numImag; j++){
		idxWrite = (pad_yxdim*j) + (ydim2Im*pad_xdim) + (y_idx*pad_xdim) + xdim2Im + x_idx;
		if(xdim%2==0)
			idxWriteMask = (pad_yxdim*j) + (ydim2Mask*pad_xdim) + (y_idx*pad_xdim) + xdim2Mask + x_idx;
		else
			idxWriteMask = idxWrite;
		padded_image_gpu[idxWrite] = d_orig_image[idx+offset];
		padded_image2_gpu[idxWrite] = image2_gpu[idx+offset];
		if(j==0)
			padded_mask_gpu[idxWriteMask] = mask[idx];
		offset += yxdim;
	}
}


void padding_masking(GpuMultidimArrayAtGpu<double> &d_orig_image, GpuMultidimArrayAtGpu<double> &mask, GpuMultidimArrayAtGpu<double> &padded_image_gpu,
		GpuMultidimArrayAtGpu<double> &padded_image2_gpu, GpuMultidimArrayAtGpu<double> &padded_mask_gpu, bool rotation, bool experimental){

    int numTh = 1024;
	int numBlk = d_orig_image.yxdim/numTh;
	if(d_orig_image.yxdim%numTh > 0)
		numBlk++;

	GpuMultidimArrayAtGpu<double> image_gpu, image2_gpu;
	image_gpu.resize(d_orig_image);
	image2_gpu.resize(d_orig_image);

	maskingKernel<<< numBlk, numTh >>> (d_orig_image.d_data, image_gpu.d_data, image2_gpu.d_data, mask.d_data,
			d_orig_image.Xdim, d_orig_image.Ydim, d_orig_image.yxdim, d_orig_image.Ndim, rotation, experimental);

	image_gpu.copyGpuToGpu(d_orig_image);

	if(!rotation || !experimental){
	/*	int xdim2 = (int)floor((padded_image_gpu.Xdim-d_orig_image.Xdim)/2);
		int ydim2 = (int)floor((padded_image_gpu.Ydim-d_orig_image.Ydim)/2);

		for(int j=0; j<d_orig_image.Ndim; j++){

			int offsetSmall=(int)(d_orig_image.yxdim*j);
			int offsetLarge=(int)((padded_image_gpu.yxdim*j) + (ydim2*padded_image_gpu.Xdim) + xdim2);

			for(int i=0; i<d_orig_image.Ydim; i++){
				gpuErrchk(cudaMemcpy((void*)&padded_image_gpu.d_data[offsetLarge], (const void*)&d_orig_image.d_data[offsetSmall], sizeof(double)*d_orig_image.Xdim, cudaMemcpyDeviceToDevice));
				gpuErrchk(cudaMemcpy((void*)&padded_image2_gpu.d_data[offsetLarge], (const void*)&image2_gpu.d_data[offsetSmall], sizeof(double)*d_orig_image.Xdim, cudaMemcpyDeviceToDevice));
				if(j==0)
					gpuErrchk(cudaMemcpy((void*)&padded_mask_gpu.d_data[offsetLarge], (const void*)&mask.d_data[offsetSmall], sizeof(double)*d_orig_image.Xdim, cudaMemcpyDeviceToDevice));
				offsetSmall += d_orig_image.Xdim;
				offsetLarge += padded_image_gpu.Xdim;
			}
		}*/
		gpuErrchk(cudaMemset(padded_image_gpu.d_data, 0, padded_image_gpu.nzyxdim*sizeof(double)));
		gpuErrchk(cudaMemset(padded_image2_gpu.d_data, 0, padded_image2_gpu.nzyxdim*sizeof(double)));
		gpuErrchk(cudaMemset(padded_mask_gpu.d_data, 0, padded_mask_gpu.nzyxdim*sizeof(double)));


		paddingKernel<<< numBlk, numTh >>>(d_orig_image.d_data, image2_gpu.d_data, mask.d_data,
				padded_image_gpu.d_data, padded_image2_gpu.d_data, padded_mask_gpu.d_data,
				d_orig_image.Xdim, d_orig_image.Ydim, d_orig_image.yxdim, d_orig_image.Ndim,
				padded_image_gpu.Xdim, padded_image_gpu.Ydim, padded_image_gpu.yxdim, rotation, experimental);

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
	//MF2realSpace.copyGpuToGpu(d_denom);
	calculateDenomFunction(MFrealSpace, MF2realSpace, maskAutocorrelation, d_denom);
	MF2realSpace.clear();

}

void cuda_calculate_correlation_rotation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat)
{
	GpuMultidimArrayAtGpu< std::complex<double> > RefExpFourier(referenceAux.d_projPolarFFT.Xdim, referenceAux.d_projPolarFFT.Ydim,
			referenceAux.d_projPolarFFT.Zdim, referenceAux.d_projPolarFFT.Ndim);

    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    referenceAux.d_projPolarFFT.calculateGridSizeVectorized(blockSize, gridSize);

    pointwiseMultiplicationComplexKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			((cufftDoubleComplex*)referenceAux.d_projPolarFFT.d_data, (cufftDoubleComplex*)experimentalAux.d_projPolarFFT.d_data,
					(cufftDoubleComplex*)RefExpFourier.d_data, referenceAux.d_projPolarFFT.nzyxdim, referenceAux.d_projPolarFFT.yxdim);

    GpuMultidimArrayAtGpu< double > RefExpRealSpace(referenceAux.XdimPolar, referenceAux.YdimPolar, referenceAux.d_projPolarFFT.Zdim,
    		referenceAux.d_projPolarFFT.Ndim);
    RefExpFourier.ifft(RefExpRealSpace);
    RefExpFourier.clear();

    XmippDim3 blockSize2(numTh, 1, 1), gridSize2;
    RefExpRealSpace.calculateGridSizeVectorized(blockSize2, gridSize2);

 	GpuMultidimArrayAtGpu<double> d_NCC(referenceAux.XdimPolar, referenceAux.YdimPolar, referenceAux.d_projPolarFFT.Zdim,
				referenceAux.d_projPolarFFT.Ndim);

	double maskFFTPolar = (referenceAux.XdimPolar*referenceAux.YdimPolar);
	calculateNccRotationKernel<<< CONVERT2DIM3(gridSize2), CONVERT2DIM3(blockSize2) >>>
			(RefExpRealSpace.d_data, (cufftDoubleComplex*)referenceAux.d_projPolarFFT.d_data, (cufftDoubleComplex*)experimentalAux.d_projPolarFFT.d_data,
					(cufftDoubleComplex*)referenceAux.d_projPolarSquaredFFT.d_data, (cufftDoubleComplex*)experimentalAux.d_projPolarSquaredFFT.d_data,
					maskFFTPolar, d_NCC.d_data, referenceAux.d_projPolarFFT.yxdim, RefExpRealSpace.nzyxdim, RefExpRealSpace.yxdim);

	//AJ para calcular el maximo
	/*numTh = 1024;
	int numBlk = d_NCC.yxdim/numTh;
	if(d_NCC.yxdim%numTh > 0)
		numBlk++;

	double *d_NCC_aux;
	gpuErrchk(cudaMalloc((void**)&d_NCC_aux, numBlk * numTh * sizeof(double)));
	gpuErrchk(cudaMemcpy(d_NCC_aux, d_NCC.d_data, sizeof(double)*(d_NCC.yxdim), cudaMemcpyDeviceToDevice));

	double *d_out_max, *d_in_aux, *d_out_max_pos;
	gpuErrchk(cudaMalloc((void**)&d_out_max, numBlk * sizeof(double)));
	gpuErrchk(cudaMalloc((void**)&d_out_max_pos, numBlk * sizeof(double)));

	compute_max<<<numBlk, numTh, 2 * numTh * sizeof(double)>>> (d_NCC_aux, d_out_max, d_out_max_pos, (int)d_NCC.yxdim);

	int numBlk3=numBlk;
	int size_aux3=numBlk;
	while(1){
		if(numBlk3>numTh){
			numBlk3=(int)(numBlk3/numTh);
			if((numBlk3%numTh)>0)
				numBlk3++;
		}else{
			numBlk3=1;
		}
		gpuErrchk(cudaMalloc((void**)&d_in_aux, numBlk3 * numTh * sizeof(double)));
		gpuErrchk(cudaMemcpy(d_in_aux, d_out_max, sizeof(double)*size_aux3, cudaMemcpyDeviceToDevice));

		compute_max<<<numBlk3, numTh, 2 * numTh * sizeof(double)>>> (d_in_aux, d_out_max, d_out_max_pos, size_aux3);
		size_aux3=numBlk3;
		if(numBlk3==1)
			break;
	}
	double max_NCC;
	double position;
	gpuErrchk(cudaMemcpy(&max_NCC, &d_out_max[0], sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(&position, &d_out_max_pos[0], sizeof(double), cudaMemcpyDeviceToHost));
	printf("max_NCC %lf \n", max_NCC);
	printf("position %lf \n", position);*/
	//FIN AJ

	double *max_values = new double[d_NCC.Ndim];
	float *posX = new float[d_NCC.Ndim];
	float *posY = new float[d_NCC.Ndim];
	d_NCC.calculateMax(max_values, posX, posY);

	TransformMatrix<float> result(transMat.Ndim);
	TransformMatrix<float> newMat(transMat.Ndim);
	newMat.setRotation(posX);

	numTh = 1024;
	int numBlk = transMat.Ndim/numTh;
	if(transMat.Ndim%numTh > 0)
		numBlk++;
	matrixMultiplication<<<numBlk, numTh>>> (newMat.d_data, transMat.d_data, result.d_data, transMat.Ndim);

	result.copyMatrix(transMat);

	delete[] max_values;
	delete[] posX;
	delete[] posY;

}


void cuda_calculate_correlation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat)
{

	GpuMultidimArrayAtGpu< std::complex<double> > RefExpFourier(referenceAux.d_projFFT.Xdim, referenceAux.d_projFFT.Ydim,
			referenceAux.d_projFFT.Zdim, referenceAux.d_projFFT.Ndim);

    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    referenceAux.d_projFFT.calculateGridSizeVectorized(blockSize, gridSize);

    pointwiseMultiplicationComplexKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			((cufftDoubleComplex*)referenceAux.d_projFFT.d_data, (cufftDoubleComplex*)experimentalAux.d_projFFT.d_data, (cufftDoubleComplex*)RefExpFourier.d_data,
					referenceAux.d_projFFT.nzyxdim, referenceAux.d_projFFT.yxdim);

    GpuMultidimArrayAtGpu< double > RefExpRealSpace(referenceAux.Xdim, referenceAux.Ydim, referenceAux.d_projFFT.Zdim,
    		referenceAux.d_projFFT.Ndim);
    RefExpFourier.ifft(RefExpRealSpace);
    RefExpFourier.clear();

 	XmippDim3 blockSize2(numTh, 1, 1), gridSize2;
 	RefExpRealSpace.calculateGridSizeVectorized(blockSize2, gridSize2);

	GpuMultidimArrayAtGpu<double> d_NCC(referenceAux.Xdim, referenceAux.Ydim, referenceAux.d_projFFT.Zdim,
			referenceAux.d_projFFT.Ndim);

	calculateNccKernel<<< CONVERT2DIM3(gridSize2), CONVERT2DIM3(blockSize2) >>>
			(RefExpRealSpace.d_data, referenceAux.MFrealSpace.d_data, experimentalAux.MFrealSpace.d_data, referenceAux.d_denom.d_data,
					experimentalAux.d_denom.d_data, referenceAux.maskAutocorrelation.d_data, d_NCC.d_data,
					referenceAux.MFrealSpace.nzyxdim, referenceAux.MFrealSpace.yxdim, referenceAux.maskCount);


	//d_NCC.copyGpuToGpu(experimentalAux.MFrealSpace);

	double *max_values = new double[d_NCC.Ndim];
	float *posX = new float[d_NCC.Ndim];
	float *posY = new float[d_NCC.Ndim];
	d_NCC.calculateMax(max_values, posX, posY);

	TransformMatrix<float> result(transMat.Ndim);
	TransformMatrix<float> newMat(transMat.Ndim);
	newMat.setTranslation(posX, posY);

	numTh = 1024;
	int numBlk = transMat.Ndim/numTh;
	if(transMat.Ndim%numTh > 0)
		numBlk++;
	matrixMultiplication<<<numBlk, numTh>>> (newMat.d_data, transMat.d_data, result.d_data, transMat.Ndim);
	result.copyMatrix(transMat);

	delete[] max_values;
	delete[] posX;
	delete[] posY;

}

void apply_transform(GpuMultidimArrayAtGpu<double> &d_original_image, GpuMultidimArrayAtGpu<double> &d_transform_image, TransformMatrix<float> &transMat){

	int numTh = 1024;
	XmippDim3 blockSize(numTh, 1, 1), gridSize;
	d_original_image.calculateGridSizeVectorized(blockSize, gridSize);

	applyTransformKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			(d_original_image.d_data, d_transform_image.d_data, transMat.d_data,
					d_original_image.nzyxdim, d_original_image.yxdim, d_original_image.Xdim, d_original_image.Ydim);

}


__global__ void cart2polar(double *image, double *polar, double *polar2, int maxRadius, int maxAng,
		int Nimgs, int Ydim, int Xdim, bool rotate)
{
	int angle = blockDim.x * blockIdx.x + threadIdx.x;
	int radius = blockDim.y * blockIdx.y + threadIdx.y;

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
