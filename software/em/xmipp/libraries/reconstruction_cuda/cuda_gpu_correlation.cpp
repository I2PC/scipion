
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



__global__ void sumRadiusKernel(float *d_in, float *d_out, float *d_out_max, float *d_out_zero, size_t dim, size_t radius, size_t ndim){

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned int numIm = floorf(idx/360);
	unsigned int angle = idx%360;

	if(idx>=dim)
		return;

	d_out[idx]=0.0;
	d_out_max[idx]=-100000;
	int idxRead=360*radius*numIm;
	for(int i=0; i<radius; i++){
		if(d_in[idxRead+(360*i)+angle]==-1.0){
			continue;
		}
		d_out[idx] += d_in[idxRead+(360*i)+angle];
		if(d_in[idxRead+(360*i)+angle]>d_out_max[idx])
			d_out_max[idx] = d_in[idxRead+(360*i)+angle];

		if(i==0)
			d_out_zero[idx] = d_in[idxRead+angle];
	}

}


__global__ void calculateMax(float *d_in, float *d_out, float *position, size_t yxdim, int Ndim, bool firstCall){

	extern __shared__ float sdata[];

	unsigned int idx = threadIdx.x;
	unsigned int blockSize = blockDim.x;


	/*//Version 3
	unsigned int i = blockIdx.x * blockSize + idx;
	int index = 0;
	for(int imN=0; imN<Ndim; imN++){

		sdata[idx]=d_in[i+index];
		if(firstCall)
			sdata[idx+blockSize] = (float)i; //AJ position
		else
			sdata[idx+blockSize] = position[i+index];

		__syncthreads();

		if(i>=yxdim)
			sdata[idx]=-1.0;
		__syncthreads();

		for(unsigned int s=ceilf((float)blockSize/2); s>0; s>>=1){
			if(idx<s){
				sdata[idx]=fmaxf(sdata[idx], sdata[idx+s]);
				sdata[idx+blockSize] = (sdata[idx]==sdata[idx+s]) ? sdata[idx+blockSize+s] : sdata[idx+blockSize];
			}
			__syncthreads();
		}

		if(idx==0){
			d_out[blockIdx.x+(gridDim.x*imN)] = sdata[0];
			position[blockIdx.x+(gridDim.x*imN)] = sdata[blockSize];
			if(imprime){
				//printf("%i d_out %lf, position %lf \n", blockIdx.x+(gridDim.x*imN), d_out[blockIdx.x+(gridDim.x*imN)], position[blockIdx.x+(gridDim.x*imN)]);
			}
		}

		__syncthreads();

		index = index+(int)yxdim;

	}*/

	//Version 6
	unsigned int i = blockIdx.x * blockSize*2 + idx;
	int index = 0;
	for(int imN=0; imN<Ndim; imN++){

		if(i<yxdim){
			if(i+blockSize < yxdim)
				sdata[idx]=fmaxf(d_in[i+index], d_in[i+blockSize+index]);
			else
				sdata[idx]=d_in[i+index];

			if(firstCall)
				sdata[idx+blockSize] = (sdata[idx]==d_in[i+index]) ? (float)i : (float)(i+blockSize); //AJ position
			else
				sdata[idx+blockSize] = (sdata[idx]==d_in[i+index]) ? position[i+index] : position[i+blockSize+index];
		}

		__syncthreads();

		if(i>=yxdim)
			sdata[idx]=-1.0;
		__syncthreads();

		//for(unsigned int s=ceil((float)blockSize/2); s>0; s>>=1){
			//if(idx<s){
				//sdata[idx]=fmaxf(sdata[idx], sdata[idx+s]);
				//sdata[idx+blockSize] = (sdata[idx]==sdata[idx+s]) ? sdata[idx+blockSize+s] : sdata[idx+blockSize];
			//}
			//__syncthreads();
		//}
		if(blockSize >= 1024){
			if(idx<512){
				sdata[idx] = fmaxf(sdata[idx], sdata[idx+512]);
				sdata[idx+blockSize] = (sdata[idx]==sdata[idx+512]) ? sdata[idx+blockSize+512] : sdata[idx+blockSize];
			}
			__syncthreads();
		}
		if(blockSize >= 512){
			if(idx<256){
				sdata[idx] = fmaxf(sdata[idx], sdata[idx+256]);
				sdata[idx+blockSize] = (sdata[idx]==sdata[idx+256]) ? sdata[idx+blockSize+256] : sdata[idx+blockSize];
			}
			__syncthreads();
		}
		if(blockSize >= 256){
			if(idx<128){
				sdata[idx] = fmaxf(sdata[idx], sdata[idx+128]);
				sdata[idx+blockSize] = (sdata[idx]==sdata[idx+128]) ? sdata[idx+blockSize+128] : sdata[idx+blockSize];
			}
			__syncthreads();
		}
		if(blockSize >= 128){
			if(idx<64){
				sdata[idx] = fmaxf(sdata[idx], sdata[idx+64]);
				sdata[idx+blockSize] = (sdata[idx]==sdata[idx+64]) ? sdata[idx+blockSize+64] : sdata[idx+blockSize];
			}
			__syncthreads();
		}
		if(idx<32){
			if(blockSize>=64){
				if(idx<32){
					sdata[idx] = fmaxf(sdata[idx], sdata[idx+32]);
					sdata[idx+blockSize] = (sdata[idx]==sdata[idx+32]) ? sdata[idx+blockSize+32] : sdata[idx+blockSize];
				}
			}
			if(blockSize>=32){
				if(idx<16){
					sdata[idx] = fmaxf(sdata[idx], sdata[idx+16]);
					sdata[idx+blockSize] = (sdata[idx]==sdata[idx+16]) ? sdata[idx+blockSize+16] : sdata[idx+blockSize];
				}
			}
			if(blockSize>=16){
				if(idx<8){
					sdata[idx] = fmaxf(sdata[idx], sdata[idx+8]);
					sdata[idx+blockSize] = (sdata[idx]==sdata[idx+8]) ? sdata[idx+blockSize+8] : sdata[idx+blockSize];
				}
			}
			if(blockSize>=8){
				if(idx<4){
					sdata[idx] = fmaxf(sdata[idx], sdata[idx+4]);
					sdata[idx+blockSize] = (sdata[idx]==sdata[idx+4]) ? sdata[idx+blockSize+4] : sdata[idx+blockSize];
				}
			}
			if(blockSize>=4){
				if(idx<2){
					sdata[idx] = fmaxf(sdata[idx], sdata[idx+2]);
					sdata[idx+blockSize] = (sdata[idx]==sdata[idx+2]) ? sdata[idx+blockSize+2] : sdata[idx+blockSize];
				}
			}
			if(blockSize>=2){
				if(idx<1){
					sdata[idx] = fmaxf(sdata[idx], sdata[idx+1]);
					sdata[idx+blockSize] = (sdata[idx]==sdata[idx+1]) ? sdata[idx+blockSize+1] : sdata[idx+blockSize];
				}
			}
		}

		if(idx==0){
			d_out[blockIdx.x+(gridDim.x*imN)] = sdata[0];
			position[blockIdx.x+(gridDim.x*imN)] = sdata[blockSize];
		}

		__syncthreads();

		index = index+(int)yxdim;

	}

}

__global__ void matrixMultiplication (float* newMat, float* lastMat, float* result, size_t n, double maxShift,
		float *maxGpu, float *NCC, size_t NCC_yxdim){

	unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
	if(idx>=n)
		return;

	double maxShift2 = maxShift*maxShift;

	int idx9 = idx*9;
	float shiftx = newMat[idx9]*lastMat[idx9+2] + newMat[idx9+1]*lastMat[idx9+5] + newMat[idx9+2]*lastMat[idx9+8];
	float shifty = newMat[idx9+3]*lastMat[idx9+2] + newMat[idx9+4]*lastMat[idx9+5] + newMat[idx9+5]*lastMat[idx9+8];
	float radShift = shiftx*shiftx + shifty*shifty;
	//if(abs(shiftx)>maxShift || abs(shifty)>maxShift)
	if(radShift > maxShift2){
		//printf("NO\n");
		result[idx9] = lastMat[idx9];
		result[idx9+1] = lastMat[idx9+1];
		result[idx9+2] = lastMat[idx9+2];
		result[idx9+3] = lastMat[idx9+3];
		result[idx9+4] = lastMat[idx9+4];
		result[idx9+5] = lastMat[idx9+5];
		maxGpu[idx] = NCC[idx*NCC_yxdim];
	}else{
		//printf("SI\n");
		result[idx9] = newMat[idx9]*lastMat[idx9] + newMat[idx9+1]*lastMat[idx9+3] + newMat[idx9+2]*lastMat[idx9+6];
		result[idx9+2] = shiftx;
		result[idx9+1] = newMat[idx9]*lastMat[idx9+1] + newMat[idx9+1]*lastMat[idx9+4] + newMat[idx9+2]*lastMat[idx9+7];
		result[idx9+3] = newMat[idx9+3]*lastMat[idx9] + newMat[idx9+4]*lastMat[idx9+3] + newMat[idx9+5]*lastMat[idx9+6];
		result[idx9+4] = newMat[idx9+3]*lastMat[idx9+1] + newMat[idx9+4]*lastMat[idx9+4] + newMat[idx9+5]*lastMat[idx9+7];
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

	MmanyF[idx] = mulOut*normFactor; //make_cuFloatComplex( cuCrealf(mulOut)*normFactor ,  cuCimagf(mulOut)*normFactor ) ;
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


__global__ void calculateNccKernel(float *RefExpRealSpace, float *MFrealSpaceRef, float *MFrealSpaceExp,
		float *MF2realSpaceRef, float *MF2realSpaceExp, float *mask, float *NCC,
		size_t nzyxdim, size_t yxdim, size_t xdim, size_t ydim, size_t maskCount, int max_shift)
{

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned long int idxLow = idx % yxdim;

	if(idx>=nzyxdim)
		return;

	int idx_x = idxLow%xdim;
	int idx_y=idxLow/xdim;
	if(idx_x>=max_shift && idx_x<xdim-max_shift){
		NCC[idx] = -1;
		return;
	}
	if(idx_y>=max_shift && idx_y<ydim-max_shift){
		NCC[idx] = -1;
		return;
	}

	//float den1 = denomRef[idx];
	//float den2 = denomExp[idx];
	float den1 = sqrt(MF2realSpaceRef[idx] - (MFrealSpaceRef[idx]*MFrealSpaceRef[idx]/mask[idxLow]));
	float den2 = sqrt(MF2realSpaceExp[idx] - (MFrealSpaceExp[idx]*MFrealSpaceExp[idx]/mask[idxLow]));

	if(den1!=0.0 && den2!=0.0 && !isnan(den1) && !isnan(den2) && mask[idxLow]>maskCount*0.9){
		float num = (RefExpRealSpace[idx] - ((MFrealSpaceRef[idx]*MFrealSpaceExp[idx])/(mask[idxLow])) );
		NCC[idx] = num/(den1*den2);
	}else
		NCC[idx] = -1;


/*
	//AJ to store in a vector only the NCC values inside the square allowed by the max_shift
	int numIm = idx/yxdim;
	int a = idx_x/(int)(xdim/2);
	int b = idx_y/(int)(ydim/2);
	int xAux, yAux, idxAux;
	int tamSq = max_shift*max_shift;
	if(a==0 && b==0){
		idxAux = idx_y*max_shift + idx_x;
	}else if(a!=0 && b==0){
		xAux = idx_x-(xdim-max_shift);
		idxAux = idx_y*max_shift + xAux + tamSq;
	}else if(a==0 && b!=0){
		yAux = idx_y-(ydim-max_shift);
		idxAux = yAux*max_shift + idx_x + 2*tamSq;
	}else{
		xAux = idx_x-(xdim-max_shift);
		yAux = idx_y-(ydim-max_shift);
		idxAux = yAux*max_shift + xAux + 3*tamSq;
	}
	NCC_aux[idxAux + numIm*4*tamSq] = NCC[idx];
*/

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
/*		if(x<0)
			x=-x;
		else if(x>=xdim)
			x=xdim-(x-xdim)-1;
		if(y<0)
			y=-y;
		else if(y>=ydim)
			y=ydim-(y-ydim)-1;
*/	/*}else if(wrap==2){ //last pixel copies
		if(x<0)
			x=0;
		else if(x>=xdim)
			x=xdim-1;
		if(y<0)
			y=0;
		else if(y>=ydim)
			y=ydim-1;
	}*/

	if(x<0)
		x=0;
	else if(x>=xdim)
		x=0;
	if(y<0)
		y=0;
	else if(y>=ydim)
		y=0;

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

	//******
	int xaux=x_orig00;
	int yaux=y_orig00;

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

	/*/******
	if(xaux<0 || xaux>=xdim || yaux<0 || yaux>=ydim){
		imVal = 0;
	}*/

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

	float normValue = 1.0/yxdimFFT;
	float maskNorm = maskFFTPolarReal*normValue;

	//cufftComplex maskFFTPolar = make_cuFloatComplex(maskFFTPolarReal, 0.0);

	float M1M2Polar = maskFFTPolarReal*maskNorm; //cuCrealf(cuCmulf(maskFFTPolar,maskFFTPolar))*normValue;
	float polarValRef = cuCrealf(polarFFTRef[idxLow])*maskNorm; //cuCrealf(cuCmulf(polarFFTRef[idxLow],maskFFTPolar))*normValue;
	float polarSqValRef = cuCrealf(polarSquaredFFTRef[idxLow])*maskNorm; //cuCrealf(cuCmulf(polarSquaredFFTRef[idxLow],maskFFTPolar))*normValue;

	float polarValExp = cuCrealf(polarFFTExp[idxLow])*maskNorm; //cuCrealf(cuCmulf(polarFFTExp[idxLow],maskFFTPolar))*normValue;
	float polarSqValExp = cuCrealf(polarSquaredFFTExp[idxLow])*maskNorm; //cuCrealf(cuCmulf(polarSquaredFFTExp[idxLow],maskFFTPolar))*normValue;

	float num = (RefExpRealSpace[idx] - (polarValRef*polarValExp/M1M2Polar) );
	float den1 = sqrt(polarSqValRef - (polarValRef*polarValRef/M1M2Polar) );
	float den2 = sqrt(polarSqValExp - (polarValExp*polarValExp/M1M2Polar) );

	if(den1!=0.0 && den2!=0.0 && !isnan(den1) && !isnan(den2))
		NCC[idx] = num/(den1*den2);
	else
		NCC[idx] = -1.0;

}


__global__ void pointwiseMultiplicationComplexKernel(cufftComplex *reference, cufftComplex *experimental,
		cufftComplex *RefExpFourier, size_t nzyxdim, size_t yxdim)
{
	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=nzyxdim)
		return;

	float normFactor = (1.0/yxdim);

	cuComplex mulOut = cuCmulf(reference[idx], experimental[idx]);
	RefExpFourier[idx] = mulOut*normFactor; //make_cuFloatComplex( cuCrealf(mulOut)*normFactor ,  cuCimagf(mulOut)*normFactor );
}

/*
__global__ void maskingKernel(float *d_in, float *d_out, float *d_out2, float *mask,
		size_t xdim, size_t ydim, size_t yxdim, size_t numImag, bool experimental){

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=yxdim)
		return;

	unsigned int x_idx = idx%(int)xdim;
	unsigned int y_idx = idx/(int)xdim;
	unsigned int idxWrite;
	if(experimental)
		idxWrite = (ydim-1 - y_idx)*xdim + (xdim-1 - x_idx);
	else
		idxWrite = y_idx*xdim + x_idx;

	int offset=0;
	for(int i=0; i<numImag; i++){
		d_out[idxWrite+offset] = d_in[idx+offset]*mask[idx];
		d_out2[idxWrite+offset] = d_out[idxWrite+offset]*d_out[idxWrite+offset];

		offset += yxdim;
	}
}
*/

__global__ void maskingPaddingKernel(float *d_in, float *mask, float *padded_image_gpu,
		float *padded_image2_gpu, float *padded_mask_gpu, size_t xdim, size_t ydim, size_t yxdim,
		size_t numImag, size_t pad_xdim, size_t pad_ydim, size_t pad_yxdim, bool experimental){

	unsigned long int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=yxdim)
		return;

	unsigned int x_idx1 = idx%(int)xdim;
	unsigned int y_idx1 = idx/(int)xdim;
	unsigned int idxWriteToMask;
	if(experimental)
		idxWriteToMask = (ydim-1 - y_idx1)*xdim + (xdim-1 - x_idx1);
	else
		idxWriteToMask = y_idx1*xdim + x_idx1;

	int xdim2Im = (int)floorf((pad_xdim-xdim)/2);
	int ydim2Im = (int)floorf((pad_ydim-ydim)/2);
	int xdim2Mask = xdim2Im;
	int ydim2Mask = ydim2Im;
	if(experimental && xdim%2==0){
		xdim2Im+=1;
		ydim2Im+=1;
	}

	unsigned int x_idx = idxWriteToMask%(int)xdim;
	unsigned int y_idx = idxWriteToMask/(int)xdim;
	unsigned int idxWrite;
	unsigned int idxWriteMask;
	float d_out, d_out2;

	int offset=0;
	for(int j=0; j<numImag; j++){

		d_out = d_in[idx+offset]*mask[idx];
		d_out2 = d_out*d_out;

		idxWrite = (pad_yxdim*j) + (ydim2Im*pad_xdim) + (y_idx*pad_xdim) + xdim2Im + x_idx;
		if(xdim%2==0)
			idxWriteMask = (pad_yxdim*j) + (ydim2Mask*pad_xdim) + (y_idx*pad_xdim) + xdim2Mask + x_idx;
		else
			idxWriteMask = idxWrite;
		padded_image_gpu[idxWrite] = d_out;
		padded_image2_gpu[idxWrite] = d_out2;
		if(j==0 && padded_mask_gpu!=NULL)
			padded_mask_gpu[idxWriteMask] = mask[idx];
		offset += yxdim;
	}
}


void padding_masking(GpuMultidimArrayAtGpu<float> &d_orig_image, GpuMultidimArrayAtGpu<float> &mask, GpuMultidimArrayAtGpu<float> &padded_image_gpu,
		GpuMultidimArrayAtGpu<float> &padded_image2_gpu, GpuMultidimArrayAtGpu<float> &padded_mask_gpu, bool experimental){

    int numTh = 1024;
	int numBlk = d_orig_image.yxdim/numTh;
	if(d_orig_image.yxdim%numTh > 0)
		numBlk++;

	gpuErrchk(cudaMemset(padded_image_gpu.d_data, 0, padded_image_gpu.nzyxdim*sizeof(float)));
	gpuErrchk(cudaMemset(padded_image2_gpu.d_data, 0, padded_image2_gpu.nzyxdim*sizeof(float)));
	if(padded_mask_gpu.d_data!=NULL)
		gpuErrchk(cudaMemset(padded_mask_gpu.d_data, 0, padded_mask_gpu.nzyxdim*sizeof(float)));

	maskingPaddingKernel<<< numBlk, numTh >>>(d_orig_image.d_data, mask.d_data,
			padded_image_gpu.d_data, padded_image2_gpu.d_data, padded_mask_gpu.d_data,
			d_orig_image.Xdim, d_orig_image.Ydim, d_orig_image.yxdim, d_orig_image.Ndim,
			padded_image_gpu.Xdim, padded_image_gpu.Ydim, padded_image_gpu.yxdim, experimental);

/*
	maskingKernel<<< numBlk, numTh >>>(d_orig_image.d_data, image_gpu.d_data, image2_gpu.d_data, mask.d_data,
			d_orig_image.Xdim, d_orig_image.Ydim, d_orig_image.yxdim, d_orig_image.Ndim, rotation, experimental);



	paddingKernel<<< numBlk, numTh >>>(image_gpu.d_data, image2_gpu.d_data, mask.d_data,
			padded_image_gpu.d_data, padded_image2_gpu.d_data, padded_mask_gpu.d_data,
			d_orig_image.Xdim, d_orig_image.Ydim, d_orig_image.yxdim, d_orig_image.Ndim,
			padded_image_gpu.Xdim, padded_image_gpu.Ydim, padded_image_gpu.yxdim, experimental);
*/

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



void GpuCorrelationAux::produceSideInfo(mycufftHandle &myhandlePaddedB, mycufftHandle &myhandleMaskB, StructuresAux &myStructureAux)
{
	//GpuMultidimArrayAtGpu< std::complex<float> > MF, MF2;
	myStructureAux.MF.resize(d_projFFT);
	myStructureAux.MF2.resize(d_projSquaredFFT);

	pointwiseMultiplicationFourier(d_maskFFT, d_projFFT, myStructureAux.MF);
	pointwiseMultiplicationFourier(d_maskFFT, d_projSquaredFFT, myStructureAux.MF2);
	//d_projSquaredFFT.clear();

	MF2realSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	MFrealSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);

	myStructureAux.MF.ifft(MFrealSpace, myhandlePaddedB);
	myStructureAux.MF2.ifft(MF2realSpace, myhandlePaddedB);
	//MF.clear();
	//MF2.clear();

	GpuMultidimArrayAtGpu< std::complex<float> > maskAux(d_projFFT.Xdim, d_projFFT.Ydim);
	pointwiseMultiplicationFourier(d_maskFFT, d_maskFFT, maskAux);
	maskAutocorrelation.resize(Xdim, Ydim);
	maskAux.ifft(maskAutocorrelation, myhandleMaskB);
	maskAux.clear();

	//d_denom.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	//calculateDenomFunction(MFrealSpace, MF2realSpace, maskAutocorrelation, d_denom);
	//MF2realSpace.clear();

}



void GpuCorrelationAux::produceSideInfo(mycufftHandle &myhandlePaddedB, mycufftHandle &myhandleMaskB, StructuresAux &myStructureAux,
		GpuMultidimArrayAtGpu<float> &maskAutocorr)
{
	//GpuMultidimArrayAtGpu< std::complex<float> > MF, MF2;
	myStructureAux.MF.resize(d_projFFT);
	myStructureAux.MF2.resize(d_projSquaredFFT);

	pointwiseMultiplicationFourier(d_maskFFT, d_projFFT, myStructureAux.MF);
	pointwiseMultiplicationFourier(d_maskFFT, d_projSquaredFFT, myStructureAux.MF2);
	//d_projSquaredFFT.clear();

	MF2realSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	MFrealSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);

	myStructureAux.MF.ifft(MFrealSpace, myhandlePaddedB);
	myStructureAux.MF2.ifft(MF2realSpace, myhandlePaddedB);
	//MF.clear();
	//MF2.clear();

	maskAutocorr.copyGpuToGpu(maskAutocorrelation);

	//d_denom.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	//calculateDenomFunction(MFrealSpace, MF2realSpace, maskAutocorrelation, d_denom);
	//MF2realSpace.clear();

}



void calculateMaxNew1D(float *max_values, float *posX, int fixPadding, int xdim, int Ndim, float *d_data, float *auxMax,
		GpuMultidimArrayAtGpu<float> d_out, GpuMultidimArrayAtGpu<float> d_pos){

    int numTh = 1024;
    int numBlk = xdim/1024;
    if(xdim%1024!=0)
    	numBlk++;
    numBlk=ceil((float)numBlk/2); //V4+
    int numBlk2, size_aux2;

    d_out.resize(numBlk*Ndim);
    d_pos.resize(numBlk*Ndim);

	//AJ TIME
	timeval start1, end1;
	double secs1;
	gettimeofday(&start1, NULL);

	//printf("1. numTh %i, numBlk %i yxdim %i \n", numTh, numBlk, yxdim);
	calculateMax<<<numBlk, numTh, 2*numTh * sizeof(float)>>>(d_data, d_out.d_data, d_pos.d_data, xdim, Ndim, true);
	/*float *posi, *max;
	posi = new float[numBlk*Ndim];
	max = new float[numBlk*Ndim];
	cudaMemcpy(max, d_out.d_data, numBlk*Ndim*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(posi, d_pos.d_data, numBlk*Ndim*sizeof(float), cudaMemcpyDeviceToHost);
	for(int j=0; j<numBlk*Ndim; j++)
		//printf("max %lf posi %lf ", max[j], posi[j]);
	//printf("\n");*/

	numBlk2=numBlk;
	size_aux2=numBlk;
	while(1){
		if(numBlk2>numTh){
		   numBlk2=numBlk2/numTh;
		   if(numBlk2%numTh!=0)
			   numBlk2++;
		   numBlk=ceil((float)numBlk/2); //V4+
		}else{
			numTh=ceil((float)size_aux2/2); //V4+
			float aux1 = log((float)numTh)/log(2.0); //V4+
			int aux2 = (int)aux1; //V4+
			float error = aux1-(float)aux2; //V4+
			if(error>0.001) //V4+
				aux2++;  //V4+
			numTh=pow(2,aux2); //V4+
			numBlk2=1;
		}
		//printf("2. numTh %i, numBlk %i yxdim %i \n", numTh, numBlk2, size_aux2);
		calculateMax<<<numBlk2, numTh, 2*numTh * sizeof(float)>>> (d_out.d_data, d_out.d_data, d_pos.d_data, size_aux2, Ndim, false);
		size_aux2=numBlk2;
		if(numBlk2==1)
			break;
   }

	//AJ TIME
	gettimeofday(&end1, NULL);
	secs1 = timeval_diff(&end1, &start1);
	//printf("MAX 2: %.16g miliseconds\n", secs1 * 1000.0);

	//AJ TIME
	timeval start2, end2;
	double secs2;
	gettimeofday(&start2, NULL);

	float h_pos;
	for(int i=0; i<Ndim; i++){

		cudaMemcpy(&h_pos, &d_pos.d_data[i], sizeof(float), cudaMemcpyDeviceToHost);
		int position = (int)h_pos;
		cudaMemcpy(&max_values[i], &auxMax[position+(i*360)], sizeof(float), cudaMemcpyDeviceToHost);

		float posX_aux = (float)(position%xdim);
		float Xdim2 = (float)(xdim/2);

		/*if(posX_aux>=Xdim2){
			posX[i] = xdim-1-posX_aux;
		}else if(posX_aux<Xdim2){
			posX[i] = -(posX_aux+1);*/
		if(posX_aux<Xdim2){
			posX[i] = -(posX_aux+1);
		}else if(posX_aux>=Xdim2){
			posX[i] = xdim-1-posX_aux;
		}

		//Fixing padding problem ¿?¿?
		posX[i]+=fixPadding;

	    //AJ TIME
	    gettimeofday(&end2, NULL);
	    secs2 = timeval_diff(&end2, &start2);
	    //printf("MAX 3: %.16g miliseconds\n", secs2 * 1000.0);

	}

}


void calculateMaxNew2D(float *max_values, float *posX, float *posY, int fixPadding, int yxdim, int Xdim, int Ydim, int Ndim, float *d_data,
		GpuMultidimArrayAtGpu<float> d_out, GpuMultidimArrayAtGpu<float> d_pos){

    int numTh = 1024;
    int numBlk = yxdim/1024;
    if(yxdim%1024!=0)
    	numBlk++;
    numBlk=ceil((float)numBlk/2); //V4+
    int numBlk2, size_aux2;

    d_out.resize(numBlk*Ndim);
    d_pos.resize(numBlk*Ndim);

	//AJ TIME
	timeval start1, end1;
	double secs1;
	gettimeofday(&start1, NULL);

	//printf("1. numTh %i, numBlk %i yxdim %i \n", numTh, numBlk, yxdim);
	calculateMax<<<numBlk, numTh, 2*numTh * sizeof(float)>>>(d_data, d_out.d_data, d_pos.d_data, yxdim, Ndim, true);
	/*float *posi, *max;
	posi = new float[numBlk*Ndim];
	max = new float[numBlk*Ndim];
	cudaMemcpy(max, d_out.d_data, numBlk*Ndim*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(posi, d_pos.d_data, numBlk*Ndim*sizeof(float), cudaMemcpyDeviceToHost);
	for(int j=0; j<numBlk*Ndim; j++)
		//printf("max %lf posi %lf ", max[j], posi[j]);
	//printf("\n");*/

	numBlk2=numBlk;
	size_aux2=numBlk;
	while(1){
		if(numBlk2>numTh){
		   numBlk2=numBlk2/numTh;
		   if(numBlk2%numTh!=0)
			   numBlk2++;
		   numBlk=ceil((float)numBlk/2); //V4+
		}else{
			numTh=ceil((float)size_aux2/2); //V4+
			float aux1 = log((float)numTh)/log(2.0); //V4+
			int aux2 = (int)aux1; //V4+
			float error = aux1-(float)aux2; //V4+
			if(error>0.001) //V4+
				aux2++;  //V4+
			numTh=pow(2,aux2); //V4+
			numBlk2=1;
		}
		//printf("2. numTh %i, numBlk %i yxdim %i \n", numTh, numBlk2, size_aux2);
		calculateMax<<<numBlk2, numTh, 2*numTh * sizeof(float)>>> (d_out.d_data, d_out.d_data, d_pos.d_data, size_aux2, Ndim, false);
		size_aux2=numBlk2;
		if(numBlk2==1)
			break;
   }

	//AJ TIME
	gettimeofday(&end1, NULL);
	secs1 = timeval_diff(&end1, &start1);
	//printf("MAX 2: %.16g miliseconds\n", secs1 * 1000.0);

	//AJ TIME
	timeval start2, end2;
	double secs2;
	gettimeofday(&start2, NULL);

	float h_pos;
	for(int i=0; i<Ndim; i++){

		cudaMemcpy(&h_pos, &d_pos.d_data[i], sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(&max_values[i], &d_out.d_data[i], sizeof(float), cudaMemcpyDeviceToHost);

		int position = (int)h_pos;

		float posX_aux = (float)(position%Xdim);
		float posY_aux = (float)(position/Xdim);
		float Xdim2 = (float)(Xdim/2);
		float Ydim2 = (float)(Ydim/2);

		if(posX_aux>=Xdim2 && posY_aux>=Ydim2){
			posX[i] = Xdim-1-posX_aux;
			posY[i] = Ydim-1-posY_aux;
		}else if(posX_aux<Xdim2 && posY_aux>=Ydim2){
			posX[i] = -(posX_aux+1);
			posY[i] = Ydim-1-posY_aux;
		}else if(posX_aux<Xdim2 && posY_aux<Ydim2){
			posX[i] = -(posX_aux+1);
			posY[i] = -(posY_aux+1);
		}else if(posX_aux>=Xdim2 && posY_aux<Ydim2){
			posX[i] = Xdim-1-posX_aux;
			posY[i] = -(posY_aux+1);
		}

		//Fixing padding problem ¿?¿?
		posX[i]+=fixPadding;
		posY[i]+=fixPadding;

	    //AJ TIME
	    gettimeofday(&end2, NULL);
	    secs2 = timeval_diff(&end2, &start2);
	    //printf("MAX 3: %.16g miliseconds\n", secs2 * 1000.0);

	}

}


void cuda_calculate_correlation_rotation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat,
		float *max_vector, int maxShift, mycufftHandle &myhandlePadded, bool mirror, StructuresAux &myStructureAux)
{
	myStructureAux.RefExpFourierPolar.resize(referenceAux.d_projPolarFFT.Xdim, referenceAux.d_projPolarFFT.Ydim,
			referenceAux.d_projPolarFFT.Zdim, referenceAux.d_projPolarFFT.Ndim);

    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    referenceAux.d_projPolarFFT.calculateGridSizeVectorized(blockSize, gridSize);

    pointwiseMultiplicationComplexKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			((cufftComplex*)referenceAux.d_projPolarFFT.d_data, (cufftComplex*)experimentalAux.d_projPolarFFT.d_data,
					(cufftComplex*)myStructureAux.RefExpFourierPolar.d_data, referenceAux.d_projPolarFFT.nzyxdim,
					referenceAux.d_projPolarFFT.yxdim);

    myStructureAux.RefExpRealSpacePolar.resize(referenceAux.XdimPolar, referenceAux.YdimPolar, referenceAux.d_projPolarFFT.Zdim,
    		referenceAux.d_projPolarFFT.Ndim);
    myStructureAux.RefExpFourierPolar.ifft(myStructureAux.RefExpRealSpacePolar, myhandlePadded);
    //RefExpFourier.clear();

    XmippDim3 blockSize2(numTh, 1, 1), gridSize2;
    myStructureAux.RefExpRealSpacePolar.calculateGridSizeVectorized(blockSize2, gridSize2);

    myStructureAux.d_NCCPolar.resize(referenceAux.XdimPolar, referenceAux.YdimPolar, referenceAux.d_projPolarFFT.Zdim,
				referenceAux.d_projPolarFFT.Ndim);

	double maskFFTPolar = (referenceAux.XdimPolar*referenceAux.YdimPolar);
	calculateNccRotationKernel<<< CONVERT2DIM3(gridSize2), CONVERT2DIM3(blockSize2) >>>
			(myStructureAux.RefExpRealSpacePolar.d_data, (cufftComplex*)referenceAux.d_projPolarFFT.d_data, (cufftComplex*)experimentalAux.d_projPolarFFT.d_data,
					(cufftComplex*)referenceAux.d_projPolarSquaredFFT.d_data, (cufftComplex*)experimentalAux.d_projPolarSquaredFFT.d_data,
					maskFFTPolar, myStructureAux.d_NCCPolar.d_data, referenceAux.d_projPolarFFT.yxdim, myStructureAux.RefExpRealSpacePolar.nzyxdim,
					myStructureAux.RefExpRealSpacePolar.yxdim);

	/*float *max_values = new float[myStructureAux.d_NCCPolar.Ndim];
	float *posX = new float[myStructureAux.d_NCCPolar.Ndim];
	float *posY = new float[myStructureAux.d_NCCPolar.Ndim];
	myStructureAux.d_NCCPolar.calculateMax(max_values, posX, posY, 0);*/


	//AJ sum along the radius
    numTh = 1024;
    int numBlk = (myStructureAux.d_NCCPolar.Xdim*myStructureAux.d_NCCPolar.Ndim)/numTh;
    if((myStructureAux.d_NCCPolar.Xdim*myStructureAux.d_NCCPolar.Ndim)%numTh!=0)
    	numBlk++;

    myStructureAux.d_NCCPolar1D.resize(myStructureAux.d_NCCPolar.Xdim,1,1,myStructureAux.d_NCCPolar.Ndim);
    myStructureAux.auxMax.resize(myStructureAux.d_NCCPolar.Xdim,1,1,myStructureAux.d_NCCPolar.Ndim);
    myStructureAux.auxZero.resize(myStructureAux.d_NCCPolar.Xdim,1,1,myStructureAux.d_NCCPolar.Ndim);
    sumRadiusKernel<<< numBlk, numTh >>>(myStructureAux.d_NCCPolar.d_data, myStructureAux.d_NCCPolar1D.d_data, myStructureAux.auxMax.d_data,
    		myStructureAux.auxZero.d_data, myStructureAux.d_NCCPolar.Xdim*myStructureAux.d_NCCPolar.Ndim, myStructureAux.d_NCCPolar.Ydim,
			myStructureAux.d_NCCPolar.Ndim);


    /*float *test = new float[myStructureAux.d_NCCPolar.Xdim*myStructureAux.d_NCCPolar.Ndim];
    cudaMemcpy((void*)test, (void*)myStructureAux.d_NCCPolar1D.d_data, sizeof(float)*myStructureAux.d_NCCPolar.Xdim*myStructureAux.d_NCCPolar.Ndim, cudaMemcpyDeviceToHost);
    for(int h=0; h<myStructureAux.d_NCCPolar.Xdim*myStructureAux.d_NCCPolar.Ndim; h++)
    	printf("%f ", test[h]);
    printf("\n");*/


	float *max_values = new float[myStructureAux.d_NCCPolar.Ndim];
	float *posX = new float[myStructureAux.d_NCCPolar.Ndim];
	float *posY = new float[myStructureAux.d_NCCPolar.Ndim];

	//float *max_values2 = new float[myStructureAux.d_NCCPolar.Ndim];
	//float *posX2 = new float[myStructureAux.d_NCCPolar.Ndim];
	//float *posY2 = new float[myStructureAux.d_NCCPolar.Ndim];

	calculateMaxNew1D(max_values, posX, 0, myStructureAux.d_NCCPolar1D.Xdim, myStructureAux.d_NCCPolar1D.Ndim, myStructureAux.d_NCCPolar1D.d_data,
			myStructureAux.auxMax.d_data, myStructureAux.d_out_polar_max, myStructureAux.d_pos_polar_max);


	//for(int h=0; h<myStructureAux.d_NCCPolar.Ndim; h++){
		 //printf("1D posX[%i] %f \n", h, posX[h]);
		 //printf("1D max_values[%i] %f \n", h, max_values[h]);
	//}

	//calculateMaxNew2D(max_values, posX, posY, 0, myStructureAux.d_NCCPolar.yxdim, myStructureAux.d_NCCPolar.Xdim, myStructureAux.d_NCCPolar.Ydim,
    		//myStructureAux.d_NCCPolar.Ndim, myStructureAux.d_NCCPolar.d_data, myStructureAux.d_out_polar_max, myStructureAux.d_pos_polar_max);


    //for(int h=0; h<myStructureAux.d_NCCPolar.Ndim; h++)
    		 //printf("2D posX[%i] %f \n", h,posX[h]);

	//printf("posX %f \n", posX[0]);
	//printf("posY %f \n", posY[0]);

	//***
	//posX[0] = 21;

    TransformMatrix<float> result(transMat.Ndim);
	TransformMatrix<float> newMat(transMat.Ndim);
	newMat.setRotation(posX);

	numTh = 1024;
	numBlk = transMat.Ndim/numTh;
	if(transMat.Ndim%numTh > 0)
		numBlk++;
/*
	float *matrixCpu = new float[9];
	printf("ANTES max %f  x %f y %f \n", max_values[0], posX[0], posY[0]);
	transMat.copyMatrixToCpu(matrixCpu);
    for(int h=0; h<9; h++){ //myStructureAux.d_NCC.Ndim
    	printf("ANTES transMat[%i] %f ", h, matrixCpu[h]);
    }
    printf("\n");
*/
	GpuMultidimArrayAtGpu<float> maxGpu(myStructureAux.d_NCCPolar1D.Ndim);
	gpuErrchk(cudaMemcpy(maxGpu.d_data, max_values, myStructureAux.d_NCCPolar1D.Ndim*sizeof(float), cudaMemcpyHostToDevice));
	matrixMultiplication<<<numBlk, numTh>>> (newMat.d_data, transMat.d_data, result.d_data, transMat.Ndim, 2*maxShift,
			maxGpu.d_data, myStructureAux.auxZero.d_data, myStructureAux.d_NCCPolar1D.yxdim);
	result.copyMatrix(transMat);

	gpuErrchk(cudaMemcpy(max_vector, maxGpu.d_data, myStructureAux.d_NCCPolar1D.Ndim*sizeof(float), cudaMemcpyDeviceToHost));

/*	for(int h=0; h<myStructureAux.d_NCCPolar.Ndim; h++){
		 printf("R1D posX[%i] %f ", h, posX[h]);
		 printf("R1D max_vector[%i] %f \n", h, max_vector[h]);
	}
*/
/*	float *matrixCpu = new float[9];
	printf("max %f  x %f y %f \n", max_values[0], posX[0], posY[0]);
	transMat.copyMatrixToCpu(matrixCpu);
    for(int h=0; h<9; h++){ //myStructureAux.d_NCC.Ndim
    	printf("transMat[%i] %f ", h, matrixCpu[h]);
    }
    printf("\n");
*/
	delete[] max_values;
	delete[] posX;
	delete[] posY;

}


void cuda_calculate_correlation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat,
		float *max_vector, int maxShift, mycufftHandle &myhandlePadded, bool mirror, StructuresAux &myStructureAux)
{

	//AJ TIME
	timeval start0, end0;
	double secs0;
    gettimeofday(&start0, NULL);

    myStructureAux.RefExpFourier.resize(referenceAux.d_projFFT.Xdim, referenceAux.d_projFFT.Ydim,
			referenceAux.d_projFFT.Zdim, referenceAux.d_projFFT.Ndim);

	//AJ TIME
	gettimeofday(&end0, NULL);
	secs0 = timeval_diff(&end0, &start0);
	//printf("CORR resize1: %.16g miliseconds\n", secs0 * 1000.0);

    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    referenceAux.d_projFFT.calculateGridSizeVectorized(blockSize, gridSize);

	//AJ TIME
	timeval start1, end1;
	double secs1;
    gettimeofday(&start1, NULL);

    pointwiseMultiplicationComplexKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			((cufftComplex*)referenceAux.d_projFFT.d_data, (cufftComplex*)experimentalAux.d_projFFT.d_data, (cufftComplex*)myStructureAux.RefExpFourier.d_data,
					referenceAux.d_projFFT.nzyxdim, referenceAux.d_projFFT.yxdim);

	//AJ TIME
	gettimeofday(&end1, NULL);
	secs1 = timeval_diff(&end1, &start1);
	//printf("CORR pointwiseMultiplicationComplexKernel: %.16g miliseconds\n", secs1 * 1000.0);
	//AJ TIME
	timeval start2, end2;
	double secs2;
    gettimeofday(&start2, NULL);

    myStructureAux.RefExpRealSpace.resize(referenceAux.Xdim, referenceAux.Ydim, referenceAux.d_projFFT.Zdim,
    		referenceAux.d_projFFT.Ndim);

	//AJ TIME
	gettimeofday(&end2, NULL);
	secs2 = timeval_diff(&end2, &start2);
	//printf("CORR resize2: %.16g miliseconds\n", secs2 * 1000.0);

	//AJ TIME
	timeval start3, end3;
	double secs3;
    gettimeofday(&start3, NULL);

    myStructureAux.RefExpFourier.ifft(myStructureAux.RefExpRealSpace, myhandlePadded);
    //RefExpFourier.clear();

	//AJ TIME
	gettimeofday(&end3, NULL);
	secs3 = timeval_diff(&end3, &start3);
	//printf("CORR ifft: %.16g miliseconds\n", secs3 * 1000.0);

 	XmippDim3 blockSize2(numTh, 1, 1), gridSize2;
 	myStructureAux.RefExpRealSpace.calculateGridSizeVectorized(blockSize2, gridSize2);

	//AJ TIME
	timeval start4, end4;
	double secs4;
    gettimeofday(&start4, NULL);

    myStructureAux.d_NCC.resize(referenceAux.Xdim, referenceAux.Ydim, referenceAux.d_projFFT.Zdim,
			referenceAux.d_projFFT.Ndim);


	//AJ TIME
	gettimeofday(&end4, NULL);
	secs4 = timeval_diff(&end4, &start4);
	//printf("CORR resize3: %.16g miliseconds\n", secs4 * 1000.0);


	//AJ TIME
	timeval start5, end5;
	double secs5;
    gettimeofday(&start5, NULL);

	calculateNccKernel<<< CONVERT2DIM3(gridSize2), CONVERT2DIM3(blockSize2) >>>
			(myStructureAux.RefExpRealSpace.d_data, referenceAux.MFrealSpace.d_data, experimentalAux.MFrealSpace.d_data, referenceAux.MF2realSpace.d_data,
					experimentalAux.MF2realSpace.d_data, referenceAux.maskAutocorrelation.d_data, myStructureAux.d_NCC.d_data, referenceAux.MFrealSpace.nzyxdim,
					referenceAux.MFrealSpace.yxdim, referenceAux.MFrealSpace.Xdim, referenceAux.MFrealSpace.Ydim, referenceAux.maskCount, maxShift);

	//AJ TIME
	gettimeofday(&end5, NULL);
	secs5 = timeval_diff(&end5, &start5);
	//printf("CORR calculateNccKernel: %.16g miliseconds\n", secs5 * 1000.0);

	//experimentalAux.debug.resize(d_NCC);
	//d_NCC.copyGpuToGpu(experimentalAux.debug);

	int fixPadding=0;
	if(referenceAux.XdimOrig%2==0)
		fixPadding=1;
	if(referenceAux.XdimOrig%2!=0)
		fixPadding=-1;

	//AJ TIME
	timeval start6, end6;
	double secs6;
    gettimeofday(&start6, NULL);

	/*float *max_values = new float[myStructureAux.d_NCC.Ndim];
	float *posX = new float[myStructureAux.d_NCC.Ndim];
	float *posY = new float[myStructureAux.d_NCC.Ndim];
    myStructureAux.d_NCC.calculateMax(max_values, posX, posY, fixPadding);*/


    //AJ test new maximum calculation
	float *max_values = new float[myStructureAux.d_NCC.Ndim];
	float *posX = new float[myStructureAux.d_NCC.Ndim];
	float *posY = new float[myStructureAux.d_NCC.Ndim];

    calculateMaxNew2D(max_values, posX, posY, fixPadding, myStructureAux.d_NCC.yxdim, myStructureAux.d_NCC.Xdim, myStructureAux.d_NCC.Ydim,
    		myStructureAux.d_NCC.Ndim, myStructureAux.d_NCC.d_data, myStructureAux.d_out_max, myStructureAux.d_pos_max);

	//AJ TIME
	gettimeofday(&end6, NULL);
	secs6 = timeval_diff(&end6, &start6);
	//printf("CORR calculateMax: %.16g miliseconds\n", secs6 * 1000.0);

	//***
	//posX[0] = 10;
	//posY[0] = 0;

	TransformMatrix<float> result(transMat.Ndim);
	TransformMatrix<float> newMat(transMat.Ndim);
	newMat.setTranslation(posX, posY);

	numTh = 1024;
	int numBlk = transMat.Ndim/numTh;
	if(transMat.Ndim%numTh > 0)
		numBlk++;

	//AJ TIME
	timeval start7, end7;
	double secs7;
    gettimeofday(&start7, NULL);

	/*float *matrixCpu = new float[9];
	printf("ANTES max %f  x %f y %f \n", max_values[0], posX[0], posY[0]);
	transMat.copyMatrixToCpu(matrixCpu);
    for(int h=0; h<9; h++){ //myStructureAux.d_NCC.Ndim
    	printf("ANTES transMat[%i] %f ", h, matrixCpu[h]);
    }
    printf("\n");
*/
	GpuMultidimArrayAtGpu<float> maxGpu(myStructureAux.d_NCC.Ndim);
	gpuErrchk(cudaMemcpy(maxGpu.d_data, max_values, myStructureAux.d_NCC.Ndim*sizeof(float), cudaMemcpyHostToDevice));
	matrixMultiplication<<<numBlk, numTh>>> (newMat.d_data, transMat.d_data, result.d_data, transMat.Ndim, 2*maxShift,
			maxGpu.d_data, myStructureAux.d_NCC.d_data, myStructureAux.d_NCC.yxdim);
	result.copyMatrix(transMat);

	gpuErrchk(cudaMemcpy(max_vector, maxGpu.d_data, myStructureAux.d_NCC.Ndim*sizeof(float), cudaMemcpyDeviceToHost));

	//AJ TIME
	gettimeofday(&end7, NULL);
	secs7 = timeval_diff(&end7, &start7);
	//printf("CORR final: %.16g miliseconds\n", secs7 * 1000.0);

/*	float *matrixCpu = new float[9];
	printf("max %f  x %f y %f \n", max_values[0], posX[0], posY[0]);
	transMat.copyMatrixToCpu(matrixCpu);
    for(int h=0; h<9; h++){ //myStructureAux.d_NCC.Ndim
    	printf("transMat[%i] %f ", h, matrixCpu[h]);
    }
    printf("\n");
*/

/*	for(int h=0; h<myStructureAux.d_NCC.Ndim; h++){
		 printf("T posX[%i] %f ", h, posX[h]);
		 printf("posY[%i] %f ", h, posY[h]);
		 printf("max_vector[%i] %f \n", h, max_vector[h]);
	}
*/
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



