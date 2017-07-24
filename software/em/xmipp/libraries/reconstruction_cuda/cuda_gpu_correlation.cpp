
//Host includes
#include "cuda_gpu_correlation.h"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <algorithm>


//CUDA includes
#include "cuda_basic_math.h"
#include "cuda_utils.h"
#include <time.h>
#include <sys/time.h>
#include <vector>

#define PI 3.14159265




__global__ void matrixMultiplication (float* newMat, float* lastMat, float* result, size_t n, double maxShift, bool *accepted){

	unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if(idx>=n)
		return;

	int idx9 = idx*9;
	float shiftx = newMat[idx9]*lastMat[idx9+2] + newMat[idx9+1]*lastMat[idx9+5] + newMat[idx9+2];
	float shifty = newMat[idx9+3]*lastMat[idx9+2] + newMat[idx9+4]*lastMat[idx9+5] + newMat[idx9+5];
	if(abs(shiftx)>maxShift || abs(shifty)>maxShift){
		accepted[idx] = false;
		result[idx9] = lastMat[idx9];
		result[idx9+1] = lastMat[idx9+1];
		result[idx9+2] = lastMat[idx9+2];
		result[idx9+3] = lastMat[idx9+3];
		result[idx9+4] = lastMat[idx9+4];
		result[idx9+5] = lastMat[idx9+5];
	}else{
		accepted[idx] = true;
		result[idx9] = newMat[idx9]*lastMat[idx9] + newMat[idx9+1]*lastMat[idx9+3];
		result[idx9+2] = shiftx;
		result[idx9+1] = newMat[idx9]*lastMat[idx9+1] + newMat[idx9+1]*lastMat[idx9+4];
		result[idx9+3] = newMat[idx9+3]*lastMat[idx9] + newMat[idx9+4]*lastMat[idx9+3];
		result[idx9+4] = newMat[idx9+3]*lastMat[idx9+1] + newMat[idx9+4]*lastMat[idx9+4];
		result[idx9+5] = shifty;
	}

}

__global__ void pointwiseMultiplicationComplexOneManyKernel(cufftComplex *M, cufftComplex *manyF, cufftComplex *MmanyF,
		size_t nzyxdim, size_t yxdim)
{
	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned long int idxLow = idx%yxdim;

	if (idx>=nzyxdim)
		return;

	float normFactor = (1.0/yxdim);

	cuComplex mulOut = cuCmulf(manyF[idx], M[idxLow]);

	MmanyF[idx] = make_cuFloatComplex( cuCrealf(mulOut)*normFactor ,  cuCimagf(mulOut)*normFactor ) ;
}

__global__ void calculateDenomFunctionKernel(float *MFrealSpace, float *MF2realSpace, float *maskAutocorrelation, float *out,
		size_t nzyxdim, size_t yxdim)
{
	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned long int idxLow = idx%yxdim;

	if (idx>=nzyxdim)
		return;

	out[idx] = sqrt(MF2realSpace[idx] - (MFrealSpace[idx]*MFrealSpace[idx]/maskAutocorrelation[idxLow]));

}


__global__ void calculateNccKernel(float *RefExpRealSpace, float *MFrealSpaceRef, float *MFrealSpaceExp, float *denomRef, float *denomExp,
		float *mask, float *NCC, size_t nzyxdim, size_t yxdim, size_t xdim, size_t ydim, size_t maskCount, int max_shift)
{

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned long int idxLow = idx % yxdim;

	if(idx>=nzyxdim)
		return;

	int idx_x = idxLow%xdim;
	int idx_y=idxLow/xdim;
	if(idx_x>max_shift && idx_x<xdim-max_shift){
		NCC[idx] = -1;
		return;
	}
	if(idx_y>max_shift && idx_y<ydim-max_shift){
		NCC[idx] = -1;
		return;
	}

	float den1 = denomRef[idx];
	float den2 = denomExp[idx];

	if(den1!=0.0 && den2!=0.0 && !isnan(den1) && !isnan(den2) && mask[idxLow]>maskCount*0.9){
		float num = (RefExpRealSpace[idx] - ((MFrealSpaceRef[idx]*MFrealSpaceExp[idx])/(mask[idxLow])) );
		NCC[idx] = num/(den1*den2);
	}else
		NCC[idx] = -1;



}

__device__ void wrapping (int &x, int &y, size_t xdim, size_t ydim){

	/*if(wrap==0){ //image copies
		if(x<0)
			x=xdim+x;
		else if(x>=xdim)
			x=x-xdim;
		if(y<0)
			y=ydim+y;
		else if(y>=ydim)
			y=y-ydim;*/
	//}else if(wrap==1){ //mirror
		if(x<0)
			x=-x;
		else if(x>=xdim)
			x=xdim-(x-xdim)-1;
		if(y<0)
			y=-y;
		else if(y>=ydim)
			y=ydim-(y-ydim)-1;
	/*}else if(wrap==2){ //last pixel copies
		if(x<0)
			x=0;
		else if(x>=xdim)
			x=xdim-1;
		if(y<0)
			y=0;
		else if(y>=ydim)
			y=ydim-1;
	}*/

}

__global__ void applyTransformKernel(float *d_in, float *d_out, float *transMat, size_t nzyxdim, size_t yxdim,
		size_t xdim, size_t ydim)
{

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned int idxIm = idx%yxdim;
	unsigned int numIm = idx/yxdim;

	//AJ wrapping type: 0 image copies, 1 mirror, 2 last pixel copies
	//int wrap = 0;

	if(idx>=nzyxdim)
		return;

	float x = idxIm%xdim;
	float y = idxIm/xdim;
	float x_orig = 0;
	float y_orig = 0;

	x -= transMat[2+(numIm*9)];
	y -= transMat[5+(numIm*9)];

	x = x - xdim/2;
	y = y - ydim/2;

	x_orig += transMat[(numIm*9)]*x - transMat[1+(numIm*9)]*y + xdim/2;
	y_orig += -transMat[3+(numIm*9)]*x + transMat[4+(numIm*9)]*y + ydim/2;

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

	wrapping (x_orig00, y_orig00, xdim, ydim);
	wrapping (x_orig01, y_orig01, xdim, ydim);
	wrapping (x_orig10, y_orig10, xdim, ydim);
	wrapping (x_orig11, y_orig11, xdim, ydim);

	int imgIdx00=y_orig00 * xdim + x_orig00;
	int imgIdx01=y_orig01 * xdim + x_orig01;
	int imgIdx10=y_orig10 * xdim + x_orig10;
	int imgIdx11=y_orig11 * xdim + x_orig11;

	int imgOffset = numIm*yxdim;
	float I00 = d_in[imgIdx00+imgOffset];
	float I01 = d_in[imgIdx01+imgOffset];
	float I10 = d_in[imgIdx10+imgOffset];
	float I11 = d_in[imgIdx11+imgOffset];
	float imVal = I00*w00 + I01*w01 + I10*w10 + I11*w11;

	d_out[idx] = imVal;

}



__global__ void calculateNccRotationKernel(float *RefExpRealSpace, cufftComplex *polarFFTRef, cufftComplex *polarFFTExp,
		cufftComplex *polarSquaredFFTRef, cufftComplex *polarSquaredFFTExp,	float maskFFTPolarReal, float *NCC,
		size_t yxdimFFT, size_t nzyxdim, size_t yxdim)
{

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned long int idxLow = (idx/(int)yxdim)*(int)yxdimFFT;

	if(idx>=nzyxdim)
		return;

	double normValue = 1.0/yxdimFFT;

	cufftComplex maskFFTPolar = make_cuFloatComplex(maskFFTPolarReal, 0.0);

	float M1M2Polar = cuCrealf(cuCmulf(maskFFTPolar,maskFFTPolar))*normValue;
	float polarValRef = cuCrealf(cuCmulf(polarFFTRef[idxLow],maskFFTPolar))*normValue;
	float polarSqValRef = cuCrealf(cuCmulf(polarSquaredFFTRef[idxLow],maskFFTPolar))*normValue;

	float polarValExp = cuCrealf(cuCmulf(polarFFTExp[idxLow],maskFFTPolar))*normValue;
	float polarSqValExp = cuCrealf(cuCmulf(polarSquaredFFTExp[idxLow],maskFFTPolar))*normValue;

	float num = (RefExpRealSpace[idx] - (polarValRef*polarValExp/M1M2Polar) );
	float den1 = sqrt(polarSqValRef - (polarValRef*polarValRef/M1M2Polar) );
	float den2 = sqrt(polarSqValExp - (polarValExp*polarValExp/M1M2Polar) );

	if(den1!=0.0 && den2!=0.0 && !isnan(den1) && !isnan(den2))
		NCC[idx] = num/(den1*den2);
	else
		NCC[idx] = -1;

}


__global__ void pointwiseMultiplicationComplexKernel(cufftComplex *reference, cufftComplex *experimental,
		cufftComplex *RefExpFourier, size_t nzyxdim, size_t yxdim)
{
	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=nzyxdim)
		return;

	float normFactor = (1.0/yxdim);

	cuComplex mulOut = cuCmulf(reference[idx], experimental[idx]);
	RefExpFourier[idx] = make_cuFloatComplex( cuCrealf(mulOut)*normFactor ,  cuCimagf(mulOut)*normFactor );
}


__global__ void maskingKernel(float *d_in, float *d_out, float *d_out2, float *mask,
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

__global__ void paddingKernel(float *d_orig_image, float *image2_gpu, float *mask, float *padded_image_gpu,
		float *padded_image2_gpu, float *padded_mask_gpu, size_t xdim, size_t ydim, size_t yxdim,
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




/*
void ifft_v2(size_t Xdim, size_t Ydim, size_t Ndim, size_t Zdim, std::complex<double> *in, double *out, mycufftHandle &myhandle){


	int Xfdim=(Xdim/2)+1;


	int nr1[] = {Xdim};   // --- Size of the image in real space
	int nr2[] = {Ydim, Xdim};   // --- Size of the image in real space
	int nr3[] = {Zdim, Ydim, Xdim};   // --- Size of the image in real space

	int nf1[] = {Xfdim};   // --- Size of the Fourier transform
	int nf2[] = {Ydim, Xfdim};   // --- Size of the Fourier transform
	int nf3[] = {Zdim, Ydim, Xfdim};   // --- Size of the Fourier transform
	int *nr=NULL, *nf=NULL;
	int NRANK; // 1D, 2D or 3D FFTs
	if (Ydim==1 && Zdim==1)
	{
		NRANK=1;
		nr=nr1;
		nf=nf1;
	}
	else if (Zdim==1)
	{
		NRANK=2;
		nr=nr2;
		nf=nf2;
	}
	else
	{
		NRANK=3;
		nr=nr3;
		nf=nf3;
	}

	int rstride = 1;				// --- Distance between two successive input/output elements
	int fstride = 1;
	int rdist = Xdim*Ydim*Zdim;	    // --- Distance between batches
	int fdist = Xfdim*Ydim*Zdim;


	//cufftHandle planF;
	cufftHandle *planBptr = new cufftHandle;
	if(myhandle.ptr == NULL){
		printf("creo el plan \n");
		cufftPlanMany(planBptr, NRANK, nr, nf, fstride, fdist, nr, rstride, rdist, CUFFT_Z2D, Ndim);
		myhandle.ptr = (void *)planBptr;
		//createPlanFFT(Xdim, Ydim, Ndim, Zdim, myhandle, *workSize);
		calculateFFTPlanSize(myhandle);
		planBptr=(cufftHandle *)myhandle.ptr;
	}else{
		printf("copio el plan \n");
		calculateFFTPlanSize(myhandle);
		planBptr=(cufftHandle *)myhandle.ptr;
	}

	//cufftHandle planB;
	//gpuErrchkFFT(cufftPlanMany(&planB, NRANK, nr, nf, fstride, fdist, nr, rstride, rdist, CUFFT_Z2D, Ndim));

	gpuErrchkFFT(cufftExecZ2D(*planBptr, (cufftDoubleComplex*)in, (cufftDoubleReal*)out));

	gpuErrchk(cudaDeviceSynchronize());
	//cufftDestroy(planB);

}
*/
/*
void fft_v2(size_t Xdim, size_t Ydim, size_t Ndim, size_t Zdim, double *in, std::complex<double> *out, mycufftHandle &myhandle){


	int Xfdim=(Xdim/2)+1;


	int nr1[] = {Xdim};   // --- Size of the image in real space
	int nr2[] = {Ydim, Xdim};   // --- Size of the image in real space
	int nr3[] = {Zdim, Ydim, Xdim};   // --- Size of the image in real space

	int nf1[] = {Xfdim};   // --- Size of the Fourier transform
	int nf2[] = {Ydim, Xfdim};   // --- Size of the Fourier transform
	int nf3[] = {Zdim, Ydim, Xfdim};   // --- Size of the Fourier transform
	int *nr=NULL, *nf=NULL;
	int NRANK; // 1D, 2D or 3D FFTs
	if (Ydim==1 && Zdim==1)
	{
		NRANK=1;
		nr=nr1;
		nf=nf1;
	}
	else if (Zdim==1)
	{
		NRANK=2;
		nr=nr2;
		nf=nf2;
	}
	else
	{
		NRANK=3;
		nr=nr3;
		nf=nf3;
	}

	int rstride = 1;				// --- Distance between two successive input/output elements
	int fstride = 1;
	int rdist = Xdim*Ydim*Zdim;	    // --- Distance between batches
	int fdist = Xfdim*Ydim*Zdim;


	//cufftHandle planF;
	cufftHandle *planFptr = new cufftHandle;
	if(myhandle.ptr == NULL){
		printf("creo el plan \n");
		cufftPlanMany(planFptr, NRANK, nr, nr, rstride, rdist, nf, fstride, fdist, CUFFT_D2Z, Ndim);
		myhandle.ptr = (void *)planFptr;
		//createPlanFFT(Xdim, Ydim, Ndim, Zdim, myhandle, *workSize);
		calculateFFTPlanSize(myhandle);
		planFptr=(cufftHandle *)myhandle.ptr;
	}else{
		printf("copio el plan \n");
		calculateFFTPlanSize(myhandle);
		planFptr=(cufftHandle *)myhandle.ptr;
	}

	//cufftHandle planF;

	//cufftPlanMany(&planF, NRANK, nr, nr, rstride, rdist, nf, fstride, fdist, CUFFT_D2Z, Ndim);

	gpuErrchkFFT(cufftExecD2Z(*planFptr, (cufftDoubleReal*)in, (cufftDoubleComplex*)out));

	gpuErrchk(cudaDeviceSynchronize());
	//cufftDestroy(planF);

}
*/

/*
void fft_v3(double *in, std::complex<double> *out, cufftHandle planF){

	gpuErrchkFFT(cufftExecD2Z(planF, (cufftDoubleReal*)in, (cufftDoubleComplex*)out));
	gpuErrchk(cudaDeviceSynchronize());

}

void ifft_v3(std::complex<double> *in, double *out, cufftHandle planB){

	gpuErrchkFFT(cufftExecZ2D(planB, (cufftDoubleComplex*)in, (cufftDoubleReal*)out));
	gpuErrchk(cudaDeviceSynchronize());

}
*/


void padding_masking(GpuMultidimArrayAtGpu<float> &d_orig_image, GpuMultidimArrayAtGpu<float> &mask, GpuMultidimArrayAtGpu<float> &padded_image_gpu,
		GpuMultidimArrayAtGpu<float> &padded_image2_gpu, GpuMultidimArrayAtGpu<float> &padded_mask_gpu, bool rotation, bool experimental){

    int numTh = 1024;
	int numBlk = d_orig_image.yxdim/numTh;
	if(d_orig_image.yxdim%numTh > 0)
		numBlk++;

	GpuMultidimArrayAtGpu<float> image_gpu, image2_gpu;
	image_gpu.resize(d_orig_image);
	image2_gpu.resize(d_orig_image);

	maskingKernel<<< numBlk, numTh >>> (d_orig_image.d_data, image_gpu.d_data, image2_gpu.d_data, mask.d_data,
			d_orig_image.Xdim, d_orig_image.Ydim, d_orig_image.yxdim, d_orig_image.Ndim, rotation, experimental);

	image_gpu.copyGpuToGpu(d_orig_image);

	if(!rotation || !experimental){
		gpuErrchk(cudaMemset(padded_image_gpu.d_data, 0, padded_image_gpu.nzyxdim*sizeof(float)));
		gpuErrchk(cudaMemset(padded_image2_gpu.d_data, 0, padded_image2_gpu.nzyxdim*sizeof(float)));
		gpuErrchk(cudaMemset(padded_mask_gpu.d_data, 0, padded_mask_gpu.nzyxdim*sizeof(float)));

		paddingKernel<<< numBlk, numTh >>>(d_orig_image.d_data, image2_gpu.d_data, mask.d_data,
				padded_image_gpu.d_data, padded_image2_gpu.d_data, padded_mask_gpu.d_data,
				d_orig_image.Xdim, d_orig_image.Ydim, d_orig_image.yxdim, d_orig_image.Ndim,
				padded_image_gpu.Xdim, padded_image_gpu.Ydim, padded_image_gpu.yxdim, rotation, experimental);

	}

}


void pointwiseMultiplicationFourier(const GpuMultidimArrayAtGpu< std::complex<float> > &M, const GpuMultidimArrayAtGpu < std::complex<float> >& manyF,
		GpuMultidimArrayAtGpu< std::complex<float> > &MmanyF)
{
    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    manyF.calculateGridSizeVectorized(blockSize, gridSize);

    pointwiseMultiplicationComplexOneManyKernel <<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			((cufftComplex*)M.d_data, (cufftComplex*)manyF.d_data, (cufftComplex*) MmanyF.d_data, manyF.nzyxdim, manyF.yxdim);

}

void calculateDenomFunction(const GpuMultidimArrayAtGpu< float > &MFrealSpace, const GpuMultidimArrayAtGpu < float >& MF2realSpace,
		const GpuMultidimArrayAtGpu < float >& maskAutocorrelation, GpuMultidimArrayAtGpu< float > &out)
{
    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    MFrealSpace.calculateGridSizeVectorized(blockSize, gridSize);

    calculateDenomFunctionKernel <<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			(MFrealSpace.d_data, MF2realSpace.d_data, maskAutocorrelation.d_data, out.d_data, MFrealSpace.nzyxdim, MFrealSpace.yxdim);

}



void GpuCorrelationAux::produceSideInfo(mycufftHandle &myhandlePaddedB, mycufftHandle &myhandleMaskB)
{
	GpuMultidimArrayAtGpu< std::complex<float> > MF, MF2;
	MF.resize(d_projFFT);
	MF2.resize(d_projSquaredFFT);

	pointwiseMultiplicationFourier(d_maskFFT, d_projFFT, MF);
	pointwiseMultiplicationFourier(d_maskFFT, d_projSquaredFFT, MF2);
	d_projSquaredFFT.clear();

	GpuMultidimArrayAtGpu<float> MF2realSpace(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	MFrealSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);

	MF.ifft(MFrealSpace, myhandlePaddedB);
	MF2.ifft(MF2realSpace, myhandlePaddedB);
	MF.clear();
	MF2.clear();

	GpuMultidimArrayAtGpu< std::complex<float> > maskAux(d_projFFT.Xdim, d_projFFT.Ydim);
	pointwiseMultiplicationFourier(d_maskFFT, d_maskFFT, maskAux);
	maskAutocorrelation.resize(Xdim, Ydim);
	maskAux.ifft(maskAutocorrelation, myhandleMaskB);
	maskAux.clear();

	d_denom.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	calculateDenomFunction(MFrealSpace, MF2realSpace, maskAutocorrelation, d_denom);
	MF2realSpace.clear();

}



void GpuCorrelationAux::produceSideInfo(mycufftHandle &myhandlePaddedB, mycufftHandle &myhandleMaskB,
		GpuMultidimArrayAtGpu<float> &maskAutocorr)
{
	GpuMultidimArrayAtGpu< std::complex<float> > MF, MF2;
	MF.resize(d_projFFT);
	MF2.resize(d_projSquaredFFT);

	pointwiseMultiplicationFourier(d_maskFFT, d_projFFT, MF);
	pointwiseMultiplicationFourier(d_maskFFT, d_projSquaredFFT, MF2);
	d_projSquaredFFT.clear();

	GpuMultidimArrayAtGpu<float> MF2realSpace(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	MFrealSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);

	MF.ifft(MFrealSpace, myhandlePaddedB);
	MF2.ifft(MF2realSpace, myhandlePaddedB);
	MF.clear();
	MF2.clear();

	/*GpuMultidimArrayAtGpu< std::complex<float> > maskAux(d_projFFT.Xdim, d_projFFT.Ydim);
	pointwiseMultiplicationFourier(d_maskFFT, d_maskFFT, maskAux);
	maskAutocorrelation.resize(Xdim, Ydim);
	maskAux.ifft(maskAutocorrelation, myhandleMaskB);
	maskAux.clear();
	maskAutocorrelation.resize(Xdim, Ydim);*/
	maskAutocorr.copyGpuToGpu(maskAutocorrelation);

	d_denom.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	calculateDenomFunction(MFrealSpace, MF2realSpace, maskAutocorrelation, d_denom);
	MF2realSpace.clear();

}

void cuda_calculate_correlation_rotation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat,
		float *max_vector, int maxShift, mycufftHandle &myhandlePadded)
{
	GpuMultidimArrayAtGpu< std::complex<float> > RefExpFourier(referenceAux.d_projPolarFFT.Xdim, referenceAux.d_projPolarFFT.Ydim,
			referenceAux.d_projPolarFFT.Zdim, referenceAux.d_projPolarFFT.Ndim);

    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    referenceAux.d_projPolarFFT.calculateGridSizeVectorized(blockSize, gridSize);

    pointwiseMultiplicationComplexKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			((cufftComplex*)referenceAux.d_projPolarFFT.d_data, (cufftComplex*)experimentalAux.d_projPolarFFT.d_data,
					(cufftComplex*)RefExpFourier.d_data, referenceAux.d_projPolarFFT.nzyxdim,
					referenceAux.d_projPolarFFT.yxdim);

    GpuMultidimArrayAtGpu< float > RefExpRealSpace(referenceAux.XdimPolar, referenceAux.YdimPolar, referenceAux.d_projPolarFFT.Zdim,
    		referenceAux.d_projPolarFFT.Ndim);
    RefExpFourier.ifft(RefExpRealSpace, myhandlePadded);
    RefExpFourier.clear();

    XmippDim3 blockSize2(numTh, 1, 1), gridSize2;
    RefExpRealSpace.calculateGridSizeVectorized(blockSize2, gridSize2);

 	GpuMultidimArrayAtGpu<float> d_NCC(referenceAux.XdimPolar, referenceAux.YdimPolar, referenceAux.d_projPolarFFT.Zdim,
				referenceAux.d_projPolarFFT.Ndim);

	double maskFFTPolar = (referenceAux.XdimPolar*referenceAux.YdimPolar);
	calculateNccRotationKernel<<< CONVERT2DIM3(gridSize2), CONVERT2DIM3(blockSize2) >>>
			(RefExpRealSpace.d_data, (cufftComplex*)referenceAux.d_projPolarFFT.d_data, (cufftComplex*)experimentalAux.d_projPolarFFT.d_data,
					(cufftComplex*)referenceAux.d_projPolarSquaredFFT.d_data, (cufftComplex*)experimentalAux.d_projPolarSquaredFFT.d_data,
					maskFFTPolar, d_NCC.d_data, referenceAux.d_projPolarFFT.yxdim, RefExpRealSpace.nzyxdim, RefExpRealSpace.yxdim);

	float *max_values = new float[d_NCC.Ndim];
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

	GpuMultidimArrayAtGpu<bool> accepted(transMat.Ndim);
	bool *acceptedCpu = new bool[transMat.Ndim];
	matrixMultiplication<<<numBlk, numTh>>> (newMat.d_data, transMat.d_data, result.d_data, transMat.Ndim, maxShift, accepted.d_data);
	result.copyMatrix(transMat);
	gpuErrchk(cudaMemcpy(acceptedCpu, accepted.d_data, transMat.Ndim*sizeof(bool), cudaMemcpyDeviceToHost));
	for(int i=0; i<d_NCC.Ndim; i++)
	{
		if(acceptedCpu[i])
			max_vector[i]=max_values[i];
		else
			gpuErrchk(cudaMemcpy(&max_vector[i], &d_NCC.d_data[i*d_NCC.yxdim], sizeof(float), cudaMemcpyDeviceToHost));
	}

	delete[] max_values;
	delete[] posX;
	delete[] posY;

}


void cuda_calculate_correlation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat,
		float *max_vector, int maxShift, mycufftHandle &myhandlePadded)
{

	GpuMultidimArrayAtGpu< std::complex<float> > RefExpFourier(referenceAux.d_projFFT.Xdim, referenceAux.d_projFFT.Ydim,
			referenceAux.d_projFFT.Zdim, referenceAux.d_projFFT.Ndim);

    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    referenceAux.d_projFFT.calculateGridSizeVectorized(blockSize, gridSize);

    pointwiseMultiplicationComplexKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			((cufftComplex*)referenceAux.d_projFFT.d_data, (cufftComplex*)experimentalAux.d_projFFT.d_data, (cufftComplex*)RefExpFourier.d_data,
					referenceAux.d_projFFT.nzyxdim, referenceAux.d_projFFT.yxdim);

    GpuMultidimArrayAtGpu< float > RefExpRealSpace(referenceAux.Xdim, referenceAux.Ydim, referenceAux.d_projFFT.Zdim,
    		referenceAux.d_projFFT.Ndim);
    RefExpFourier.ifft(RefExpRealSpace, myhandlePadded);
    RefExpFourier.clear();

 	XmippDim3 blockSize2(numTh, 1, 1), gridSize2;
 	RefExpRealSpace.calculateGridSizeVectorized(blockSize2, gridSize2);

	GpuMultidimArrayAtGpu<float> d_NCC(referenceAux.Xdim, referenceAux.Ydim, referenceAux.d_projFFT.Zdim,
			referenceAux.d_projFFT.Ndim);

	calculateNccKernel<<< CONVERT2DIM3(gridSize2), CONVERT2DIM3(blockSize2) >>>
			(RefExpRealSpace.d_data, referenceAux.MFrealSpace.d_data, experimentalAux.MFrealSpace.d_data, referenceAux.d_denom.d_data,
					experimentalAux.d_denom.d_data, referenceAux.maskAutocorrelation.d_data, d_NCC.d_data,
					referenceAux.MFrealSpace.nzyxdim, referenceAux.MFrealSpace.yxdim, referenceAux.MFrealSpace.Xdim, referenceAux.MFrealSpace.Ydim,
					referenceAux.maskCount, maxShift);

	//experimentalAux.debug.resize(d_NCC);
	//d_NCC.copyGpuToGpu(experimentalAux.debug);

	float *max_values = new float[d_NCC.Ndim];
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

	GpuMultidimArrayAtGpu<bool> accepted(transMat.Ndim);
	bool *acceptedCpu = new bool[transMat.Ndim];
	matrixMultiplication<<<numBlk, numTh>>> (newMat.d_data, transMat.d_data, result.d_data, transMat.Ndim, maxShift, accepted.d_data);
	result.copyMatrix(transMat);
	gpuErrchk(cudaMemcpy(acceptedCpu, accepted.d_data, transMat.Ndim*sizeof(bool), cudaMemcpyDeviceToHost));
	for(int i=0; i<d_NCC.Ndim; i++)
	{
		if(acceptedCpu[i])
			max_vector[i]=max_values[i];
		else
			gpuErrchk(cudaMemcpy(&max_vector[i], &d_NCC.d_data[i*d_NCC.yxdim], sizeof(float), cudaMemcpyDeviceToHost));
	}


	delete[] max_values;
	delete[] posX;
	delete[] posY;

}

void apply_transform(GpuMultidimArrayAtGpu<float> &d_original_image, GpuMultidimArrayAtGpu<float> &d_transform_image, TransformMatrix<float> &transMat){

	int numTh = 1024;
	XmippDim3 blockSize(numTh, 1, 1), gridSize;
	d_original_image.calculateGridSizeVectorized(blockSize, gridSize);

	applyTransformKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			(d_original_image.d_data, d_transform_image.d_data, transMat.d_data,
					d_original_image.nzyxdim, d_original_image.yxdim, d_original_image.Xdim, d_original_image.Ydim);

}


__global__ void cart2polar(float *image, float *polar, float *polar2, int maxRadius, int maxAng,
		int Nimgs, int Ydim, int Xdim, bool rotate)
{
	int angle = blockDim.x * blockIdx.x + threadIdx.x;
	int radius = blockDim.y * blockIdx.y + threadIdx.y;

	if (radius>=maxRadius || angle>=maxAng)
		return;

	float x = (float)(radius*cosf((float)(angle*PI/180))) + Xdim/2;
	float y = (float)(radius*sinf((float)(angle*PI/180))) + Ydim/2;

	float dx_low = floor(x);
	float dy_low = floor(y);
	int x_low = (int)dx_low;
	int y_low = (int)dy_low;
	float x_x_low=x-dx_low;
	float y_y_low=y-dy_low;
	float one_x=1.0-x_x_low;
	float one_y=1.0-y_y_low;
	float w00=one_y*one_x;
	float w01=one_y*x_x_low;
	float w10=y_y_low*one_x;
	float w11=y_y_low*x_x_low;

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
		float I00 = image[imgIdx00+imgOffset];
		float I01 = image[imgIdx01+imgOffset];
		float I10 = image[imgIdx10+imgOffset];
		float I11 = image[imgIdx11+imgOffset];
		float imVal = I00*w00 + I01*w01 + I10*w10 + I11*w11;
		int finalPolarIndex=polarIdx+polarOffset;
		polar[finalPolarIndex] = imVal;
		polar2[finalPolarIndex] = imVal*imVal;

		imgOffset+=NXY;
		polarOffset+=NXYpolar;
	}

}

void cuda_cart2polar(GpuMultidimArrayAtGpu<float> &image, GpuMultidimArrayAtGpu<float> &polar_image, GpuMultidimArrayAtGpu<float> &polar2_image, bool rotate)
{
    int numTh = 32;
    XmippDim3 blockSize(numTh, numTh, 1), gridSize;
    polar_image.calculateGridSize(blockSize, gridSize);
    cart2polar <<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
    		(image.d_data, polar_image.d_data, polar2_image.d_data, polar_image.Ydim, polar_image.Xdim, polar_image.Ndim, image.Ydim, image.Xdim, rotate);
}

/*
void produceSideInfoCuda(GpuCorrelationAux &aux, mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask,
		mycufftHandle &myhandlePolar, size_t memoryUsed){

	//AJ cambiar la mascara para crear un plan menos, es decir, la mascara se repite tantas veces como imagen
	cufftHandle planFPadded, planBPadded, planFPolar, planBPolar, planFMask, planBMask;
	createPlanFFT(aux.d_padded_image.Xdim, aux.d_padded_image.Ydim, aux.d_padded_image.Ndim, aux.d_padded_image.Zdim, planFPadded, planBPadded);
	createPlanFFT(aux.d_polar_image.Xdim, aux.d_polar_image.Ydim, aux.d_polar_image.Ndim, aux.d_polar_image.Zdim, planFPolar, planBPolar);
	createPlanFFT(aux.d_mask.Xdim, aux.d_mask.Ydim, aux.d_mask.Ndim, aux.d_mask.Zdim, planFMask, planBMask);

	myhandlePadded.ptr = &planFPadded;
	myhandlePolar.ptr = &planFPolar;
	myhandleMask.ptr = &planFMask;
	printf("En produceSideInfoCuda\n");
	calculateFFTPlanSize(myhandlePadded);
	calculateFFTPlanSize(myhandlePolar);
	calculateFFTPlanSize(myhandleMask);

	int Xfdim=(aux.d_padded_image.Xdim/2)+1;
	aux.d_projFFT.resize(Xfdim, aux.d_padded_image.Ydim, 1, aux.d_padded_image.Ndim);
	aux.d_projSquaredFFT.resize(Xfdim, aux.d_padded_image.Ydim, 1, aux.d_padded_image.Ndim);
	fft_v3(aux.d_padded_image.d_data, aux.d_projFFT.d_data, planFPadded);
	fft_v3(aux.d_padded2_image.d_data, aux.d_projSquaredFFT.d_data, planFPadded);
	aux.d_padded2_image.clear();

	aux.d_maskFFT.resize(Xfdim, aux.d_mask.Ydim, 1, 1);
	fft_v3(aux.d_mask.d_data, aux.d_maskFFT.d_data, planFMask);

	Xfdim=(aux.d_polar_image.Xdim/2)+1;
	aux.d_projPolarFFT.resize(Xfdim, aux.d_polar_image.Ydim, 1, aux.d_polar_image.Ndim);
	aux.d_projPolarSquaredFFT.resize(Xfdim, aux.d_polar_image.Ydim, 1, aux.d_polar_image.Ndim);
	fft_v3(aux.d_polar_image.d_data, aux.d_projPolarFFT.d_data, planFPolar);
	aux.d_polar_image.clear();
	fft_v3(aux.d_polar2_image.d_data, aux.d_projPolarSquaredFFT.d_data, planFPolar);
	aux.d_polar2_image.clear();



	GpuMultidimArrayAtGpu< std::complex<double> > MF, MF2;
	MF.resize(aux.d_projFFT);
	MF2.resize(aux.d_projSquaredFFT);

	pointwiseMultiplicationFourier(aux.d_maskFFT, aux.d_projFFT, MF);
	pointwiseMultiplicationFourier(aux.d_maskFFT, aux.d_projSquaredFFT, MF2);
	aux.d_projSquaredFFT.clear();

	GpuMultidimArrayAtGpu<double> MF2realSpace(aux.d_padded_image.Xdim, aux.d_padded_image.Ydim,
			aux.d_padded_image.Zdim, aux.d_padded_image.Ndim);
	aux.MFrealSpace.resize(aux.d_padded_image.Xdim, aux.d_padded_image.Ydim,
			aux.d_padded_image.Zdim, aux.d_padded_image.Ndim);
	ifft_v3(MF.d_data, aux.MFrealSpace.d_data, planBPadded);
	ifft_v3(MF2.d_data, MF2realSpace.d_data, planBPadded);

	aux.d_padded_image.clear();
	MF.clear();
	MF2.clear();

	GpuMultidimArrayAtGpu< std::complex<double> > maskAux(aux.d_maskFFT.Xdim, aux.d_maskFFT.Ydim);
	pointwiseMultiplicationFourier(aux.d_maskFFT, aux.d_maskFFT, maskAux);
	aux.maskAutocorrelation.resize(aux.d_mask.Xdim, aux.d_mask.Ydim);
	ifft_v3(maskAux.d_data, aux.maskAutocorrelation.d_data, planBMask);
	maskAux.clear();

	aux.d_denom.resize(aux.MFrealSpace.Xdim, aux.MFrealSpace.Ydim, aux.MFrealSpace.Zdim, aux.MFrealSpace.Ndim);
	calculateDenomFunction(aux.MFrealSpace, MF2realSpace, aux.maskAutocorrelation, aux.d_denom);
	MF2realSpace.clear();

	float memory[3]={0, 0, 0};
    cuda_check_gpu_memory(memory);
    printf("After Total %f, free %f, used %f\n\n",memory[0],memory[1], memory[2]);
    memoryUsed = memory[2];

    cufftDestroy(planFPadded);
    cufftDestroy(planBPadded);
    cufftDestroy(planFPolar);
    cufftDestroy(planBPolar);
    cufftDestroy(planFMask);
    cufftDestroy(planBMask);

}
*/

/*
void cuda_preprocess_experimental_image (double *expImages, GpuCorrelationAux &d_referenceAux,
		GpuCorrelationAux &d_experimentalAux, GpuMultidimArrayAtGpu<double> &mask,
		bool rotation, size_t Xdim, size_t Ydim, size_t pad_Xdim, size_t pad_Ydim, size_t angles, size_t radius,
		int numImagesRef, cufftHandle planFPadded, cufftHandle planFMask, cufftHandle planFPolar, int expIdx){

	d_experimentalAux.d_original_image.resize(Xdim,Ydim,1,numImagesRef);

	for(size_t i=0; i<numImagesRef; i++)
		d_experimentalAux.d_original_image.fillImageToGpu(&expImages[expIdx*Xdim*Ydim],i);

	GpuMultidimArrayAtGpu<double> image_stack_gpu(Xdim,Ydim,1,numImagesRef);
	d_experimentalAux.d_original_image.copyGpuToGpu(image_stack_gpu);

	if(!rotation){
		d_experimentalAux.d_padded_image.resize(pad_Xdim, pad_Ydim, 1, numImagesRef);
		d_experimentalAux.d_padded2_image.resize(pad_Xdim, pad_Ydim, 1, numImagesRef);
		d_experimentalAux.d_mask.resize(pad_Xdim, pad_Ydim, 1, 1);

		padding_masking(d_experimentalAux.d_original_image, mask, d_experimentalAux.d_padded_image, d_experimentalAux.d_padded2_image,
				d_experimentalAux.d_mask, rotation, true);

		int Xfdim=(d_experimentalAux.d_padded_image.Xdim/2)+1;
		d_experimentalAux.d_projFFT.resize(Xfdim, d_experimentalAux.d_padded_image.Ydim, 1, d_experimentalAux.d_padded_image.Ndim);
		d_experimentalAux.d_projSquaredFFT.resize(Xfdim, d_experimentalAux.d_padded_image.Ydim, 1, d_experimentalAux.d_padded_image.Ndim);
		fft_v3(d_experimentalAux.d_padded_image.d_data, d_experimentalAux.d_projFFT.d_data, planFPadded);
		fft_v3(d_experimentalAux.d_padded2_image.d_data, d_experimentalAux.d_projSquaredFFT.d_data, planFPadded);
		d_experimentalAux.d_padded2_image.clear();

		//d_experimentalAux.d_maskFFT.resize(Xfdim, d_experimentalAux.d_mask.Ydim, 1, 1);
		//fft_v3(d_experimentalAux.d_mask.d_data, d_experimentalAux.d_maskFFT.d_data, planFMask);

	}else{

		d_experimentalAux.d_polar_image.resize(360,radius,1,numImagesRef);
		d_experimentalAux.d_polar2_image.resize(360,radius,1,numImagesRef);
		cuda_cart2polar(d_experimentalAux.d_original_image, d_experimentalAux.d_polar_image,
				d_experimentalAux.d_polar2_image, true);

		int Xfdim=(d_experimentalAux.d_polar_image.Xdim/2)+1;
		d_experimentalAux.d_projPolarFFT.resize(Xfdim, d_experimentalAux.d_polar_image.Ydim, 1, d_experimentalAux.d_polar_image.Ndim);
		d_experimentalAux.d_projPolarSquaredFFT.resize(Xfdim, d_experimentalAux.d_polar_image.Ydim, 1, d_experimentalAux.d_polar_image.Ndim);
		fft_v3(d_experimentalAux.d_polar_image.d_data, d_experimentalAux.d_projPolarFFT.d_data, planFPolar);
		d_experimentalAux.d_polar_image.clear();
		fft_v3(d_experimentalAux.d_polar2_image.d_data, d_experimentalAux.d_projPolarSquaredFFT.d_data, planFPolar);
		d_experimentalAux.d_polar2_image.clear();

	}

}
*/


/*
void produceSideInfoCudaExp (GpuCorrelationAux &d_referenceAux, GpuCorrelationAux &d_experimentalAux, cufftHandle planBPadded){

	GpuMultidimArrayAtGpu< std::complex<double> > MF, MF2;
	MF.resize(d_experimentalAux.d_projFFT);
	MF2.resize(d_experimentalAux.d_projSquaredFFT);

	pointwiseMultiplicationFourier(d_referenceAux.d_maskFFT, d_experimentalAux.d_projFFT, MF);
	pointwiseMultiplicationFourier(d_referenceAux.d_maskFFT, d_experimentalAux.d_projSquaredFFT, MF2);
	d_experimentalAux.d_projSquaredFFT.clear();

	GpuMultidimArrayAtGpu<double> MF2realSpace(d_experimentalAux.d_padded_image.Xdim, d_experimentalAux.d_padded_image.Ydim,
			d_experimentalAux.d_padded_image.Zdim, d_experimentalAux.d_padded_image.Ndim);
	d_experimentalAux.MFrealSpace.resize(d_experimentalAux.d_padded_image.Xdim, d_experimentalAux.d_padded_image.Ydim,
			d_experimentalAux.d_padded_image.Zdim, d_experimentalAux.d_padded_image.Ndim);
	ifft_v3(MF.d_data, d_experimentalAux.MFrealSpace.d_data, planBPadded);
	ifft_v3(MF2.d_data, MF2realSpace.d_data, planBPadded);

	d_experimentalAux.d_padded_image.clear();
	MF.clear();
	MF2.clear();

	d_referenceAux.maskAutocorrelation.copyGpuToGpu(d_experimentalAux.maskAutocorrelation);

	d_experimentalAux.d_denom.resize(d_experimentalAux.MFrealSpace.Xdim, d_experimentalAux.MFrealSpace.Ydim,
			d_experimentalAux.MFrealSpace.Zdim, d_experimentalAux.MFrealSpace.Ndim);
	calculateDenomFunction(d_experimentalAux.MFrealSpace, MF2realSpace,
			d_experimentalAux.maskAutocorrelation, d_experimentalAux.d_denom);
	MF2realSpace.clear();

}
*/

/*
void cuda_align_experimental_image(double *expImages, GpuCorrelationAux &d_referenceAux, GpuCorrelationAux &d_experimentalAux,
		TransformMatrix<float> &transMat_tr, TransformMatrix<float> &transMat_rt, double *max_vector_tr, double *max_vector_rt,
		int numImagesRef, GpuMultidimArrayAtGpu<double> &mask, double maxShift,
		mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar, size_t mdExpSize)
{
	size_t Xdim, Ydim, pad_Xdim, pad_Ydim, angles, radius;
	pad_Xdim = d_referenceAux.Xdim;
	pad_Ydim = d_referenceAux.Ydim;
	Xdim = (pad_Xdim+1)/2;
	Ydim = (pad_Ydim+1)/2;
	angles=d_referenceAux.XdimPolar;
	radius=d_referenceAux.YdimPolar;

	//AJ creating FFT plans (forward and backward)
	cufftHandle planFPadded, planBPadded, planFPolar, planBPolar, planFMask, planBMask;
	createPlanFFT(pad_Xdim, pad_Ydim, numImagesRef, 1, planFPadded, planBPadded);
	createPlanFFT(angles, radius, numImagesRef, 1, planFPolar, planBPolar);
	//createPlanFFT(pad_Xdim, pad_Ydim, 1, 1, planFMask, planBMask);

	for(int expIdx=0; expIdx<mdExpSize; expIdx++){ //loop for every experimental image read bu the cpu (in cpu memory)

		//preprocess
		//BUCLEEEEEEEEEEEE
		bool rotation=false;
		cuda_preprocess_experimental_image (expImages, d_referenceAux, d_experimentalAux, mask,
				rotation, Xdim, Ydim, pad_Xdim, pad_Ydim, angles, radius, numImagesRef,
				planFPadded, planFMask, planFPolar, expIdx);

		//produceSideInfo
		if(!rotation){
			d_experimentalAux.maskCount=d_referenceAux.maskCount;
			produceSideInfoCudaExp (d_referenceAux, d_experimentalAux, planBPadded);

			float memory[3]={0, 0, 0};
		    cuda_check_gpu_memory(memory);
		    printf("After Total %f, free %f, used %f\n\n",memory[0],memory[1], memory[2]);

		}

	}//AJ END FOR mdExpSize

}
*/
