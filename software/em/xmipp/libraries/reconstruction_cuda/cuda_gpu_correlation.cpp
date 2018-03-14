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


__global__ void calcAbsKernel(cufftComplex *d_in, float *d_out, int dim){

	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=dim)
		return;

	d_out[idx]=d_in[idx].x;

}

__global__ void sumRadiusKernel(float *d_in, float *d_out, float *d_out_max, float *d_out_zero,
		int dim, int radius, int ndim){

	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int numIm = floorf(idx/360.0f);
	int angle = idx%360;

	if(idx>=dim)
		return;

	float out = 0.0;
	float out_max = -100000;
	int idxRead=360*radius*numIm;
	for(int i=0; i<radius; i++){
		if(d_in[idxRead+(360*i)+angle]==-1.0){
			continue;
		}
		out += d_in[idxRead+(360*i)+angle];
		if(d_in[idxRead+(360*i)+angle]>d_out_max[idx]){
			out_max = d_in[idxRead+(360*i)+angle];
		}

		if(i==0)
			d_out_zero[idx] = d_in[idxRead+angle];
	}
	d_out[idx] = out;
	d_out_max[idx] = out_max;

}


__global__ void calculateMax2(float *d_in, float *d_out, float *position, int yxdim, int Ndim, bool firstCall){

	extern __shared__ float sdata[];

	int idx = threadIdx.x;
	int blockSize = blockDim.x;


	//Version 6
	int i = blockIdx.x * blockSize + idx;

	//printf("d_in[%i] %f \n", i, d_in[i]);

	int index = 0;
	for(int imN=0; imN<Ndim; imN++){

		if(i<yxdim*gridDim.x){
			sdata[idx]=d_in[i+index];

			if(firstCall)
				sdata[idx+blockSize] = (float)idx; //AJ position
			else
				sdata[idx+blockSize] = position[i+index];
		}
		//if (idx==0)
			//printf("i %i, sdata %f \n", i, sdata[idx]);

		__syncthreads();

		if(i>=yxdim*gridDim.x)
			sdata[idx]=-1.0;
		__syncthreads();


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
			//printf("idx %i sdata[0] %f sdata[blockSize] %f \n", blockIdx.x+(gridDim.x*imN), sdata[0], sdata[blockSize]);
			d_out[blockIdx.x+(gridDim.x*imN)] = sdata[0];
			position[blockIdx.x+(gridDim.x*imN)] = sdata[blockSize];
		}

		__syncthreads();

		index = index+(int)yxdim;

	}

}


__global__ void pointwiseMultiplicationComplexOneManyKernel_three(cufftComplex *M,
		cufftComplex *manyF, cufftComplex *MmanyF, cufftComplex *manyF_sq, cufftComplex *MmanyF_sq,
		cufftComplex *maskAux, int nzyxdim, int yxdim, int ndim, bool power2)
{

	int idx = threadIdx.x;
	int nIm = blockIdx.x;
	int blockSize = blockDim.x;

	int posTh = nIm*yxdim + idx;
	int n = ceilf(((float)yxdim/(float)blockSize));

	if (idx>=yxdim)
		return;

	float normFactor = 1.0f/yxdim;

	int myIdx = posTh;
	int myIdxMask = idx;

	cufftComplex myM = M[myIdxMask];
	cuComplex mulOut = cuCmulf(manyF[myIdx], myM);
	cuComplex mulOut_sq = cuCmulf(manyF_sq[myIdx], myM);

	MmanyF[myIdx] = mulOut*normFactor;
	MmanyF_sq[myIdx] = mulOut_sq*normFactor;
	if(nIm==0)
		maskAux[myIdx] = cuCmulf(myM, myM)*normFactor;

	myIdx+=blockSize;
	myIdxMask+=blockSize;

	for (int i=1; i<n; i++){

		if (posTh+i*blockSize < yxdim*(nIm+1)){

			myM = M[myIdxMask];
			mulOut = cuCmulf(manyF[myIdx], myM);
			mulOut_sq = cuCmulf(manyF_sq[myIdx], myM);

			MmanyF[myIdx] = mulOut*normFactor;
			MmanyF_sq[myIdx] = mulOut_sq*normFactor;
			if(nIm==0)
				maskAux[myIdx] = cuCmulf(myM, myM)*normFactor;

			myIdx+=blockSize;
			myIdxMask+=blockSize;

		}
	}
}


__global__ void pointwiseMultiplicationComplexOneManyKernel_two(cufftComplex *M,
		cufftComplex *manyF, cufftComplex *MmanyF, cufftComplex *manyF_sq, cufftComplex *MmanyF_sq,
		int nzyxdim, int yxdim, int ndim, bool power2)
{

	int idx = threadIdx.x;
	int nIm = blockIdx.x;
	int blockSize = blockDim.x;

	int posTh = nIm*yxdim + idx;
	int n = ceilf(((float)yxdim/(float)blockSize));

	if (idx>=yxdim)
		return;

	float normFactor = 1.0f/yxdim;

	int myIdx = posTh;
	int myIdxMask = idx;

	cufftComplex myM = M[myIdxMask];
	cuComplex mulOut = cuCmulf(manyF[myIdx], myM);
	cuComplex mulOut_sq = cuCmulf(manyF_sq[myIdx], myM);

	MmanyF[myIdx] = mulOut*normFactor;
	MmanyF_sq[myIdx] = mulOut_sq*normFactor;

	myIdx+=blockSize;
	myIdxMask+=blockSize;

	for (int i=1; i<n; i++){

		if (posTh+i*blockSize < yxdim*(nIm+1)){

			myM = M[myIdxMask];
			mulOut = cuCmulf(manyF[myIdx], myM);
			mulOut_sq = cuCmulf(manyF_sq[myIdx], myM);

			MmanyF[myIdx] = mulOut*normFactor;
			MmanyF_sq[myIdx] = mulOut_sq*normFactor;

			myIdx+=blockSize;
			myIdxMask+=blockSize;

		}
	}
}



__global__ void pointwiseMultiplicationComplexOneManyKernel(cufftComplex *M, cufftComplex *manyF, cufftComplex *MmanyF,
		int nzyxdim, int yxdim, bool power2)
{

	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int idxLow;
	if (power2==true)
		idxLow = idx & (yxdim-1);
	else
		idxLow = idx%yxdim;

	if (idx>=nzyxdim)
		return;

	float normFactor = (1.0f/yxdim);

	cuComplex mulOut = cuCmulf(manyF[idx], M[idxLow]);
	MmanyF[idx] = mulOut*normFactor;

}

//NOW WE DONT USE THIS KERNEL
__global__ void calculateDenomFunctionKernel(float *MFrealSpace, float *MF2realSpace, float *maskAutocorrelation,
		float *out, int nzyxdim, int yxdim, bool power2)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int idxLow;
	if (power2==true)
		idxLow = idx & (yxdim-1);
	else
		idxLow = idx%yxdim;

	if (idx>=nzyxdim)
		return;

	out[idx] = sqrt(MF2realSpace[idx] - (MFrealSpace[idx]*MFrealSpace[idx]/maskAutocorrelation[idxLow]));

}


__global__ void calculateNccKernel(float *RefExpRealSpace, float *MFrealSpaceRef, float *MFrealSpaceExp,
		float *MF2realSpaceRef, float *MF2realSpaceExp, float *mask, float *NCC,
		int nzyxdim, int yxdim, int xdim, int ydim, int maskCount, int max_shift, bool power2yx, bool power2x)
{

	//change unsigned... and size_t
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int idxLow;
	if (power2yx==true)
		idxLow = idx & (yxdim-1);
	else
		idxLow = idx%yxdim;

	if(idx>=nzyxdim)
		return;

	int idx_x;
	if (power2x==true)
		idx_x = idxLow & (xdim-1);
	else
		idx_x = idxLow%xdim;

	int idx_y=idxLow/xdim;
	if(idx_x>=max_shift && idx_x<xdim-max_shift){
		NCC[idx] = -1;
		return;
	}
	if(idx_y>=max_shift && idx_y<ydim-max_shift){
		NCC[idx] = -1;
		return;
	}

	float maskTh = mask[idxLow];
	float MF2realSpaceRefTh = MF2realSpaceRef[idx];
	float MFrealSpaceRefTh = MFrealSpaceRef[idx];
	float MF2realSpaceExpTh = MF2realSpaceExp[idx];
	float MFrealSpaceExpTh = MFrealSpaceExp[idx];
	float RefExpRealSpaceTh = RefExpRealSpace[idx];

	float den1 = sqrt(MF2realSpaceRefTh - (MFrealSpaceRefTh*MFrealSpaceRefTh/maskTh));
	float den2 = sqrt(MF2realSpaceExpTh - (MFrealSpaceExpTh*MFrealSpaceExpTh/maskTh));

	float num;
	if(den1!=0.0 && den2!=0.0 && !isnan(den1) && !isnan(den2) && maskTh>maskCount*0.9){
		num = (RefExpRealSpaceTh - ((MFrealSpaceRefTh*MFrealSpaceExpTh)/(maskTh)) );
		NCC[idx] = num/(den1*den2);
	}else
		NCC[idx] = -1;

}

__device__ void wrapping (int &x, int &y, int xdim, int ydim){
		//mirror wrapping
		if(x<0)
			x=-x;
		else if(x>=xdim)
			x=xdim-(x-xdim)-1;
		if(y<0)
			y=-y;
		else if(y>=ydim)
			y=ydim-(y-ydim)-1;
}

__global__ void applyTransformKernel(float *d_in, float *d_out, float *transMat, int nzyxdim, int yxdim,
		int xdim, int ydim, bool power2yx, bool power2x)
{

	extern __shared__ float transMatSD[];

	int idxIm = blockDim.x * blockIdx.x + threadIdx.x;
	int idxTh = threadIdx.x;
	int numIm = blockIdx.y;


	if(idxIm>=yxdim)
		return;

	if(idxTh<9)
		transMatSD[idxTh] = transMat[idxTh+(numIm*9)];
	__syncthreads();

	float x;
	if (power2x==true)
		x = idxIm & (xdim-1);
	else
		x = idxIm%xdim;

	float y = idxIm/xdim;
	float x_orig = 0;
	float y_orig = 0;

	//transMat into shared memory
	x -= transMatSD[2]; //transMat[2+(numIm*9)];
	y -= transMatSD[5]; //transMat[5+(numIm*9)];

	x = x - xdim/2.0f;
	y = y - ydim/2.0f;

	x_orig += transMatSD[0]*x - transMatSD[1]*y + xdim/2.0f;
	y_orig += -transMatSD[3]*x + transMatSD[4]*y + xdim/2.0f;
	// x_orig += transMat[(numIm*9)]*x - transMat[1+(numIm*9)]*y + xdim/2.0f;
	//y_orig += -transMat[3+(numIm*9)]*x + transMat[4+(numIm*9)]*y + ydim/2.0f;

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

	d_out[idxIm+(numIm*yxdim)] = imVal;

}



__global__ void calculateNccRotationKernel(float *RefExpRealSpace, cufftComplex *polarFFTRef, cufftComplex *polarFFTExp,
		cufftComplex *polarSquaredFFTRef, cufftComplex *polarSquaredFFTExp,	float maskFFTPolarReal, float *NCC,
		int yxdimFFT, int nzyxdim, int yxdim)
{

	//change unsigned lont int for int and all the size_t
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int idxLow = (idx/yxdim)*yxdimFFT;

	if(idx>=nzyxdim)
		return;

	float normValue = 1.0f/yxdimFFT;
	float maskNorm = maskFFTPolarReal*normValue;

	float M1M2Polar = maskFFTPolarReal*maskNorm;
	float polarValRef = cuCrealf(polarFFTRef[idxLow])*maskNorm;
	float polarSqValRef = cuCrealf(polarSquaredFFTRef[idxLow])*maskNorm;

	float polarValExp = cuCrealf(polarFFTExp[idxLow])*maskNorm;
	float polarSqValExp = cuCrealf(polarSquaredFFTExp[idxLow])*maskNorm;

	float num = (RefExpRealSpace[idx] - (polarValRef*polarValExp/M1M2Polar) );
	float den1 = sqrt(polarSqValRef - (polarValRef*polarValRef/M1M2Polar) );
	float den2 = sqrt(polarSqValExp - (polarValExp*polarValExp/M1M2Polar) );

	if(den1!=0.0 && den2!=0.0 && !isnan(den1) && !isnan(den2))
		NCC[idx] = num/(den1*den2);
	else
		NCC[idx] = -1.0;

}


__global__ void pointwiseMultiplicationComplexKernel(cufftComplex *reference, cufftComplex *experimental,
		cufftComplex *RefExpFourier, int nzyxdim, int yxdim)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=nzyxdim)
		return;

	float normFactor = (1.0f/yxdim);

	cuComplex mulOut = cuCmulf(reference[idx], experimental[idx]);
	RefExpFourier[idx] = mulOut*normFactor;
}

//TODO: try this instead with blockIdx.y being the number of images maybe with a subset,
//i.e., every thread processing a bunch of images, instead of just one
__global__ void maskingPaddingKernel(float *d_in, float *mask, float *padded_image_gpu,
		float *padded_image2_gpu, float *padded_mask_gpu, int xdim, int ydim, int yxdim,
		int numImag, int pad_xdim, int pad_ydim, int pad_yxdim, bool experimental, bool power2x){

	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int numIm = blockIdx.y;

	if(idx>=yxdim)
		return;

	int x_idx1;
	if (power2x==true)
		x_idx1 = idx & (xdim-1);
	else
		x_idx1 = idx%xdim;

	int y_idx1 = idx/xdim;
	int idxWriteToMask;
	if(experimental)
		idxWriteToMask = (ydim-1 - y_idx1)*xdim + (xdim-1 - x_idx1);
	else
		idxWriteToMask = y_idx1*xdim + x_idx1;

	int xdim2Im = (int)floorf((pad_xdim-xdim)/2.0f);
	int ydim2Im = (int)floorf((pad_ydim-ydim)/2.0f);
	int xdim2Mask = xdim2Im;
	int ydim2Mask = ydim2Im;
	if(experimental && (xdim&1)==0){
		xdim2Im+=1;
		ydim2Im+=1;
	}

	int x_idx;
	if (power2x==true)
		x_idx = idxWriteToMask & (xdim-1);
	else
		x_idx = idxWriteToMask%xdim;

	int y_idx = idxWriteToMask/xdim;
	int idxWrite;
	int idxWriteMask;
	float d_out, d_out2;

	int offset=0;
	//for(int j=0; j<numImag; j++){

	int j = numIm;

		offset=j*yxdim;

		d_out = d_in[idx+offset]*mask[idx];
		d_out2 = d_out*d_out;

		idxWrite = (pad_yxdim*j) + (ydim2Im*pad_xdim) + (y_idx*pad_xdim) + xdim2Im + x_idx;
		if((xdim&1)==0)
			idxWriteMask = (pad_yxdim*j) + (ydim2Mask*pad_xdim) + (y_idx*pad_xdim) + xdim2Mask + x_idx;
		else
			idxWriteMask = idxWrite;
		padded_image_gpu[idxWrite] = d_out;
		padded_image2_gpu[idxWrite] = d_out2;
		if(j==0 && padded_mask_gpu!=NULL)
			padded_mask_gpu[idxWriteMask] = mask[idx];

	//}
}


__global__ void buildTranslationMatrix (float *d_pos, float* lastMat, float* result,
		float *maxGpu, float *NCC, int Xdim, int Ydim, int Ndim, int NCC_yxdim,
		int fixPadding, double maxShift2, bool power2x){

	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=Ndim)
		return;

	int position = (int)d_pos[idx];

	float posX_aux;
	if (power2x==true)
		posX_aux = position & (Xdim-1);
	else
		posX_aux = position%Xdim;

	float posY_aux = (float)(position/Xdim);
	float Xdim2 = (Xdim/2.0f);
	float Ydim2 = (Ydim/2.0f);
	float posX, posY;

	if(posX_aux>=Xdim2 && posY_aux>=Ydim2){
		posX = Xdim-1-posX_aux;
		posY = Ydim-1-posY_aux;
	}else if(posX_aux<Xdim2 && posY_aux>=Ydim2){
		posX = -(posX_aux+1);
		posY = Ydim-1-posY_aux;
	}else if(posX_aux<Xdim2 && posY_aux<Ydim2){
		posX = -(posX_aux+1);
		posY = -(posY_aux+1);
	}else if(posX_aux>=Xdim2 && posY_aux<Ydim2){
		posX = Xdim-1-posX_aux;
		posY = -(posY_aux+1);
	}

	//Fixing padding problem
	posX+=fixPadding;
	posY+=fixPadding;

	int idx9 = idx*9;

	float new0 = 1.0f;
	float new1 = 0.0f;
	float new2 = -posX;
	float new3 = 0.0f;
	float new4 = 1.0f;
	float new5 = -posY;
	//float new6 = 0.0f;
	//float new7 = 0.0f;
	//float new8 = 1.0f;

	float last0 = lastMat[idx9];
	float last1 = lastMat[idx9+1];
	float last2 = lastMat[idx9+2];
	float last3 = lastMat[idx9+3];
	float last4 = lastMat[idx9+4];
	float last5 = lastMat[idx9+5];
	float last6 = lastMat[idx9+6];
	float last7 = lastMat[idx9+7];
	float last8 = lastMat[idx9+8];

	float shiftx = new0*last2 + new1*last5 + new2*last8;
	float shifty = new3*last2 + new4*last5 + new5*last8;
	float radShift = shiftx*shiftx + shifty*shifty;
	if(radShift > maxShift2){
		result[idx9] = last0;
		result[idx9+1] = last1;
		result[idx9+2] = last2;
		result[idx9+3] = last3;
		result[idx9+4] = last4;
		result[idx9+5] = last5;
		maxGpu[idx] = NCC[idx*NCC_yxdim];
	}else{
		result[idx9] = new0*last0 + new1*last3 + new2*last6;
		result[idx9+2] = shiftx;
		result[idx9+1] = new0*last1 + new1*last4 + new2*last7;
		result[idx9+3] = new3*last0 + new4*last3 + new5*last6;
		result[idx9+4] = new3*last1 + new4*last4 + new5*last7;
		result[idx9+5] = shifty;
	}

}



__global__ void buildRotationMatrix (float *d_pos, float* lastMat, float* result,
		float *maxGpu, float *auxMax, float *NCC, int Xdim, int Ndim, int NCC_yxdim,
		int fixPadding, double maxShift2, bool power2x){

	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=Ndim)
		return;

	float posX;
	int position = (int)d_pos[idx];
	maxGpu[idx] = auxMax[position+(idx*360)];

	float posX_aux;
	if (power2x==true)
		posX_aux = position & (Xdim-1);
	else
		posX_aux = position%Xdim;

	float Xdim2 = (Xdim/2.0f);


	if(posX_aux<Xdim2){
		posX = -(posX_aux+1);
	}else if(posX_aux>=Xdim2){
		posX = Xdim-1-posX_aux;
	}

	//Fixing padding problem
	posX+=fixPadding;

	float rad = (float)(-posX*PI/180);

	int idx9 = idx*9;

	float new0 = cosf(rad);
	float new1 = -sinf(rad);
	float new2 = 0.0f;
	float new3 = sinf(rad);
	float new4 = cosf(rad);
	float new5 = 0.0f;
	//float new6 = 0.0f;
	//float new7 = 0.0f;
	//float new8 = 1.0f;

	float last0 = lastMat[idx9];
	float last1 = lastMat[idx9+1];
	float last2 = lastMat[idx9+2];
	float last3 = lastMat[idx9+3];
	float last4 = lastMat[idx9+4];
	float last5 = lastMat[idx9+5];
	float last6 = lastMat[idx9+6];
	float last7 = lastMat[idx9+7];
	float last8 = lastMat[idx9+8];

	float shiftx = new0*last2 + new1*last5 + new2*last8;
	float shifty = new3*last2 + new4*last5 + new5*last8;
	float radShift = shiftx*shiftx + shifty*shifty;
	if(radShift > maxShift2){
		result[idx9] = last0;
		result[idx9+1] = last1;
		result[idx9+2] = last2;
		result[idx9+3] = last3;
		result[idx9+4] = last4;
		result[idx9+5] = last5;
		maxGpu[idx] = NCC[idx*NCC_yxdim];
	}else{
		result[idx9] = new0*last0 + new1*last3 + new2*last6;
		result[idx9+2] = shiftx;
		result[idx9+1] = new0*last1 + new1*last4 + new2*last7;
		result[idx9+3] = new3*last0 + new4*last3 + new5*last6;
		result[idx9+4] = new3*last1 + new4*last4 + new5*last7;
		result[idx9+5] = shifty;
	}

}



__global__ void cart2polar(float *image, float *polar, float *polar2, int maxRadius, int maxAng,
		int Nimgs, int Ydim, int Xdim, bool rotate)
{
	int angle = blockDim.x * blockIdx.x + threadIdx.x;
	int radius = blockDim.y * blockIdx.y + threadIdx.y;

	int numIm = blockIdx.z;

	if (radius>=maxRadius || angle>=maxAng)
		return;

	float x = (float)(radius)*cosf((float)angle*(float)PI/180.0f) + (float)Xdim/2.0f;
	float y = (float)(radius)*sinf((float)angle*(float)PI/180.0f) + (float)Ydim/2.0f;

	float dx_low = floorf(x);
	float dy_low = floorf(y);
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

		int n=numIm;

		imgOffset = n*NXY;
		polarOffset = n*NXYpolar;
		float I00 = image[imgIdx00+imgOffset];
		float I01 = image[imgIdx01+imgOffset];
		float I10 = image[imgIdx10+imgOffset];
		float I11 = image[imgIdx11+imgOffset];
		float imVal = I00*w00 + I01*w01 + I10*w10 + I11*w11;
		int finalPolarIndex=polarIdx+polarOffset;
		polar[finalPolarIndex] = imVal;
		polar2[finalPolarIndex] = imVal*imVal;

}



__global__ void calculateMaxThreads (float *d_in, float *d_out, float *position,
		int yxdim, int Ndim, int yxdim2, int Ndim2){

	int idx = threadIdx.x;
	int nIm = blockIdx.x;
	int blockSize = blockDim.x;

	int posTh = nIm*yxdim + idx;
	int n = ceilf(((float)yxdim/(float)blockSize));

	float tmp, tmpPos;

	if(idx>=yxdim){
		tmp = -1.0f;
		tmpPos = -1.0f;
	}else{
		tmp = d_in[posTh];
		tmpPos = idx;
		for (int i=1; i<n; i++){
			if (posTh+i*blockSize < yxdim*(nIm+1)){
				tmp = fmaxf(tmp, d_in[posTh+i*blockSize]);
				tmpPos = (tmp==d_in[posTh+i*blockSize]) ? idx+i*blockSize : tmpPos;
			}
		}
	}
	int posOut = nIm*blockSize + idx;
	d_out[posOut] = tmp;
	position[posOut] = tmpPos;

	__syncthreads();


	extern __shared__ float sdata[];


	//Version 6
	int i = blockIdx.x * blockSize + idx;


	int index = 0;
	for(int imN=0; imN<Ndim2; imN++){

		if(i<yxdim2*gridDim.x){
			sdata[idx]=d_out[i+index];
			sdata[idx+blockSize] = position[i+index];
		}

		__syncthreads();

		if(i>=yxdim2*gridDim.x)
			sdata[idx]=-1.0;
		__syncthreads();


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

		index = index+(int)yxdim2;

	}

}


void padding_masking(GpuMultidimArrayAtGpu<float> &d_orig_image, GpuMultidimArrayAtGpu<float> &mask,
		GpuMultidimArrayAtGpu<float> &padded_image_gpu, GpuMultidimArrayAtGpu<float> &padded_image2_gpu,
		GpuMultidimArrayAtGpu<float> &padded_mask_gpu, bool experimental, myStreamHandle &myStream){

    int numTh = 1024;
	int numBlk = d_orig_image.yxdim/numTh;
	if(d_orig_image.yxdim%numTh > 0)
		numBlk++;

	dim3 blockSize(numTh,1,1);
	dim3 gridSize(numBlk, d_orig_image.Ndim, 1);

	cudaStream_t *stream = (cudaStream_t*) myStream.ptr;
	gpuErrchk(cudaMemsetAsync(padded_image_gpu.d_data, 0, padded_image_gpu.nzyxdim*sizeof(float), *stream));
	gpuErrchk(cudaMemsetAsync(padded_image2_gpu.d_data, 0, padded_image2_gpu.nzyxdim*sizeof(float), *stream));
	if(padded_mask_gpu.d_data!=NULL)
		gpuErrchk(cudaMemsetAsync(padded_mask_gpu.d_data, 0, padded_mask_gpu.nzyxdim*sizeof(float), *stream));

	bool power2;
	if (d_orig_image.Xdim & (d_orig_image.Xdim-1))
		power2 = false;
	else
		power2 = true;
	maskingPaddingKernel<<< gridSize, blockSize, 0, *stream>>>(d_orig_image.d_data, mask.d_data,
			padded_image_gpu.d_data, padded_image2_gpu.d_data, padded_mask_gpu.d_data,
			d_orig_image.Xdim, d_orig_image.Ydim, d_orig_image.yxdim, d_orig_image.Ndim,
			padded_image_gpu.Xdim, padded_image_gpu.Ydim, padded_image_gpu.yxdim, experimental, power2);

}


void pointwiseMultiplicationFourier(const GpuMultidimArrayAtGpu< std::complex<float> > &M, const GpuMultidimArrayAtGpu < std::complex<float> >& manyF,
		GpuMultidimArrayAtGpu< std::complex<float> > &MmanyF, myStreamHandle &myStream)
{
	cudaStream_t *stream = (cudaStream_t*) myStream.ptr;
    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    manyF.calculateGridSizeVectorized(blockSize, gridSize);

    bool power2;
    if (manyF.yxdim & (manyF.yxdim-1))
    	power2 = false;
    else
    	power2 = true;
    pointwiseMultiplicationComplexOneManyKernel <<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize), 0, *stream >>>
			((cufftComplex*)M.d_data, (cufftComplex*)manyF.d_data, (cufftComplex*) MmanyF.d_data, manyF.nzyxdim, manyF.yxdim, power2);

}

void pointwiseMultiplicationFourier_two(const GpuMultidimArrayAtGpu< std::complex<float> > &M,
		const GpuMultidimArrayAtGpu < std::complex<float> >& manyF,
		GpuMultidimArrayAtGpu< std::complex<float> > &MmanyF,
		const GpuMultidimArrayAtGpu < std::complex<float> >& manyF_sq,
		GpuMultidimArrayAtGpu< std::complex<float> > &MmanyF_sq,
		myStreamHandle &myStream)
{
	cudaStream_t *stream = (cudaStream_t*) myStream.ptr;
    int numTh = 1024;
    int numBlk = manyF.Ndim;

    bool power2;
    if (manyF.yxdim & (manyF.yxdim-1))
    	power2 = false;
    else
    	power2 = true;

    pointwiseMultiplicationComplexOneManyKernel_two <<< numBlk, numTh, 0, *stream >>>
			((cufftComplex*)M.d_data, (cufftComplex*)manyF.d_data, (cufftComplex*) MmanyF.d_data,
					(cufftComplex*)manyF_sq.d_data, (cufftComplex*) MmanyF_sq.d_data,
					manyF.nzyxdim, manyF.yxdim, manyF.Ndim, power2);

}

void pointwiseMultiplicationFourier_three(const GpuMultidimArrayAtGpu< std::complex<float> > &M,
		const GpuMultidimArrayAtGpu < std::complex<float> >& manyF,
		GpuMultidimArrayAtGpu< std::complex<float> > &MmanyF,
		const GpuMultidimArrayAtGpu < std::complex<float> >& manyF_sq,
		GpuMultidimArrayAtGpu< std::complex<float> > &MmanyF_sq,
		myStreamHandle &myStream,
		GpuMultidimArrayAtGpu< std::complex<float> > &maskAux)
{
	cudaStream_t *stream = (cudaStream_t*) myStream.ptr;
    int numTh = 1024;
    int numBlk = manyF.Ndim;

    bool power2;
    if (manyF.yxdim & (manyF.yxdim-1))
    	power2 = false;
    else
    	power2 = true;

    pointwiseMultiplicationComplexOneManyKernel_three <<< numBlk, numTh, 0, *stream >>>
			((cufftComplex*)M.d_data, (cufftComplex*)manyF.d_data, (cufftComplex*) MmanyF.d_data,
					(cufftComplex*)manyF_sq.d_data, (cufftComplex*) MmanyF_sq.d_data,
					(cufftComplex*) maskAux.d_data,
					manyF.nzyxdim, manyF.yxdim, manyF.Ndim, power2);

}


void GpuCorrelationAux::produceSideInfo(mycufftHandle &myhandlePaddedB, mycufftHandle &myhandleMaskB,
		StructuresAux &myStructureAux, myStreamHandle &myStream)
{
	myStructureAux.MF.resize(d_projFFT);
	myStructureAux.MF2.resize(d_projSquaredFFT);
	GpuMultidimArrayAtGpu< std::complex<float> > maskAux(d_projFFT.Xdim, d_projFFT.Ydim);

	pointwiseMultiplicationFourier_three(d_maskFFT, d_projFFT, myStructureAux.MF, d_projSquaredFFT,
			myStructureAux.MF2, myStream, maskAux);

	MF2realSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	MFrealSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);

	GpuMultidimArrayAtGpu< std::complex<float> > dull;

	myStructureAux.MF.ifftStream(MFrealSpace, myhandlePaddedB, myStream, false, dull);
	myStructureAux.MF2.ifftStream(MF2realSpace, myhandlePaddedB, myStream, false, dull);

	maskAutocorrelation.resize(Xdim, Ydim);
	maskAux.ifftStream(maskAutocorrelation, myhandleMaskB, myStream, false, dull);
	maskAux.clear();

}



void GpuCorrelationAux::produceSideInfo(mycufftHandle &myhandlePaddedB, mycufftHandle &myhandleMaskB,
		StructuresAux &myStructureAux, GpuMultidimArrayAtGpu<float> &maskAutocorr, myStreamHandle &myStream)
{

	myStructureAux.MF.resize(d_projFFT);
	myStructureAux.MF2.resize(d_projSquaredFFT);


	pointwiseMultiplicationFourier_two(d_maskFFT, d_projFFT, myStructureAux.MF, d_projSquaredFFT,
			myStructureAux.MF2, myStream);

	GpuMultidimArrayAtGpu< std::complex<float> > dull;

	MF2realSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	MFrealSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);

	myStructureAux.MF.ifftStream(MFrealSpace, myhandlePaddedB, myStream, false, dull);
	myStructureAux.MF2.ifftStream(MF2realSpace, myhandlePaddedB, myStream, false, dull);

	maskAutocorr.copyGpuToGpuStream(maskAutocorrelation, myStream);

}


void GpuCorrelationAux::produceSideInfo(mycufftHandle &myhandlePaddedB, mycufftHandle &myhandleMaskB,
		StructuresAux &myStructureAux, GpuMultidimArrayAtGpu<float> &maskAutocorr, myStreamHandle &myStream,
		bool skip, mycufftHandle &ifftcb)
{

	GpuMultidimArrayAtGpu< std::complex<float> > dull;

	MF2realSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);
	MFrealSpace.resize(Xdim, Ydim, d_projFFT.Zdim, d_projFFT.Ndim);

	d_projFFT.ifftStream(MFrealSpace, ifftcb, myStream, true, d_maskFFT);
	d_projSquaredFFT.ifftStream(MF2realSpace, ifftcb, myStream, true, d_maskFFT);

	maskAutocorr.copyGpuToGpuStream(maskAutocorrelation, myStream);


}


void calculateMaxNew2DNew(int yxdim, int Ndim, float *d_data,
		GpuMultidimArrayAtGpu<float> &d_out, GpuMultidimArrayAtGpu<float> &d_pos, myStreamHandle &myStream){

	cudaStream_t *stream = (cudaStream_t*) myStream.ptr;

    int numTh = 1024;
    int numBlk = Ndim;

    d_out.resize(numTh * numBlk);
    d_pos.resize(numTh * numBlk);

    calculateMaxThreads<<<numBlk, numTh, 2*numTh * sizeof(float), *stream>>>(d_data, d_out.d_data, d_pos.d_data, yxdim, Ndim, numTh, 1);

}


void cuda_calculate_correlation_rotation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat,
		float *max_vector, int maxShift, mycufftHandle &myhandlePadded, bool mirror,
		StructuresAux &myStructureAux, myStreamHandle &myStream, TransformMatrix<float> &resultRT)
{


	cudaStream_t *stream = (cudaStream_t*) myStream.ptr;

	myStructureAux.RefExpFourierPolar.resize(referenceAux.d_projPolarFFT.Xdim, referenceAux.d_projPolarFFT.Ydim,
			referenceAux.d_projPolarFFT.Zdim, referenceAux.d_projPolarFFT.Ndim);

    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    referenceAux.d_projPolarFFT.calculateGridSizeVectorized(blockSize, gridSize);

    pointwiseMultiplicationComplexKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize), 0, *stream >>>
			((cufftComplex*)referenceAux.d_projPolarFFT.d_data, (cufftComplex*)experimentalAux.d_projPolarFFT.d_data,
					(cufftComplex*)myStructureAux.RefExpFourierPolar.d_data, referenceAux.d_projPolarFFT.nzyxdim,
					referenceAux.d_projPolarFFT.yxdim);

    GpuMultidimArrayAtGpu< std::complex<float> > dull;
    myStructureAux.RefExpRealSpacePolar.resize(referenceAux.XdimPolar, referenceAux.YdimPolar, referenceAux.d_projPolarFFT.Zdim,
    		referenceAux.d_projPolarFFT.Ndim);
    myStructureAux.RefExpFourierPolar.ifftStream(myStructureAux.RefExpRealSpacePolar, myhandlePadded, myStream, false, dull);

    XmippDim3 blockSize2(numTh, 1, 1), gridSize2;
    myStructureAux.RefExpRealSpacePolar.calculateGridSizeVectorized(blockSize2, gridSize2);

    myStructureAux.d_NCCPolar.resize(referenceAux.XdimPolar, referenceAux.YdimPolar, referenceAux.d_projPolarFFT.Zdim,
				referenceAux.d_projPolarFFT.Ndim);

	double maskFFTPolar = (referenceAux.XdimPolar*referenceAux.YdimPolar);
	calculateNccRotationKernel<<< CONVERT2DIM3(gridSize2), CONVERT2DIM3(blockSize2), 0, *stream >>>
			(myStructureAux.RefExpRealSpacePolar.d_data, (cufftComplex*)referenceAux.d_projPolarFFT.d_data, (cufftComplex*)experimentalAux.d_projPolarFFT.d_data,
					(cufftComplex*)referenceAux.d_projPolarSquaredFFT.d_data, (cufftComplex*)experimentalAux.d_projPolarSquaredFFT.d_data,
					maskFFTPolar, myStructureAux.d_NCCPolar.d_data, referenceAux.d_projPolarFFT.yxdim, myStructureAux.RefExpRealSpacePolar.nzyxdim,
					myStructureAux.RefExpRealSpacePolar.yxdim);

	//AJ sum along the radius
    numTh = 1024;
    int numBlk = (myStructureAux.d_NCCPolar.Xdim*myStructureAux.d_NCCPolar.Ndim)/numTh;
    if((myStructureAux.d_NCCPolar.Xdim*myStructureAux.d_NCCPolar.Ndim)%numTh!=0)
    	numBlk++;

    myStructureAux.d_NCCPolar1D.resize(myStructureAux.d_NCCPolar.Xdim,1,1,myStructureAux.d_NCCPolar.Ndim);
    myStructureAux.auxMax.resize(myStructureAux.d_NCCPolar.Xdim,1,1,myStructureAux.d_NCCPolar.Ndim);
    myStructureAux.auxZero.resize(myStructureAux.d_NCCPolar.Xdim,1,1,myStructureAux.d_NCCPolar.Ndim);
    sumRadiusKernel<<< numBlk, numTh, 0, *stream >>>(myStructureAux.d_NCCPolar.d_data, myStructureAux.d_NCCPolar1D.d_data, myStructureAux.auxMax.d_data,
    		myStructureAux.auxZero.d_data, myStructureAux.d_NCCPolar.Xdim*myStructureAux.d_NCCPolar.Ndim, myStructureAux.d_NCCPolar.Ydim,
			myStructureAux.d_NCCPolar.Ndim);

	calculateMaxNew2DNew(myStructureAux.d_NCCPolar1D.Xdim, myStructureAux.d_NCCPolar1D.Ndim, myStructureAux.d_NCCPolar1D.d_data,
			myStructureAux.d_out_polar_max, myStructureAux.d_pos_polar_max, myStream);

	numTh = 1024;
	numBlk = transMat.Ndim/numTh;
	if(transMat.Ndim%numTh > 0)
		numBlk++;

	bool _power2x;
	if (myStructureAux.d_NCCPolar1D.Xdim & (myStructureAux.d_NCCPolar1D.Xdim-1))
		_power2x = false;
	else
		_power2x = true;
	double maxShift2 = (2*maxShift)*(2*maxShift);
	myStructureAux.maxGpu.resize(myStructureAux.d_NCCPolar1D.Ndim);
	buildRotationMatrix<<<numBlk, numTh, 0, *stream>>> (myStructureAux.d_pos_polar_max.d_data, transMat.d_data,
			resultRT.d_data, myStructureAux.maxGpu.d_data, myStructureAux.auxMax.d_data, myStructureAux.auxZero.d_data,
			myStructureAux.d_NCCPolar1D.Xdim, myStructureAux.d_NCCPolar1D.Ndim,
			myStructureAux.d_NCCPolar1D.yxdim, 0, maxShift2, _power2x);

	resultRT.copyMatrix(transMat, myStream);

	gpuErrchk(cudaMemcpyAsync(max_vector, myStructureAux.maxGpu.d_data, myStructureAux.maxGpu.Ndim*sizeof(float), cudaMemcpyDeviceToHost, *stream));

}


void cuda_calculate_correlation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat,
		float *max_vector, int maxShift, mycufftHandle &myhandlePadded, bool mirror,
		StructuresAux &myStructureAux, myStreamHandle &myStream, TransformMatrix<float> &resultTR,
		bool saveMaxVector)
{

	cudaStream_t *stream = (cudaStream_t*) myStream.ptr;

    myStructureAux.RefExpFourier.resize(referenceAux.d_projFFT.Xdim, referenceAux.d_projFFT.Ydim,
			referenceAux.d_projFFT.Zdim, referenceAux.d_projFFT.Ndim);

    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    referenceAux.d_projFFT.calculateGridSizeVectorized(blockSize, gridSize);

    pointwiseMultiplicationComplexKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize), 0, *stream >>>
			((cufftComplex*)referenceAux.d_projFFT.d_data, (cufftComplex*)experimentalAux.d_projFFT.d_data, (cufftComplex*)myStructureAux.RefExpFourier.d_data,
					referenceAux.d_projFFT.nzyxdim, referenceAux.d_projFFT.yxdim);


    myStructureAux.RefExpRealSpace.resize(referenceAux.Xdim, referenceAux.Ydim, referenceAux.d_projFFT.Zdim,
    		referenceAux.d_projFFT.Ndim);

    GpuMultidimArrayAtGpu< std::complex<float> > dull;
    myStructureAux.RefExpFourier.ifftStream(myStructureAux.RefExpRealSpace, myhandlePadded, myStream, false, dull);


 	XmippDim3 blockSize2(numTh, 1, 1), gridSize2;
 	myStructureAux.RefExpRealSpace.calculateGridSizeVectorized(blockSize2, gridSize2);

    myStructureAux.d_NCC.resize(referenceAux.Xdim, referenceAux.Ydim, referenceAux.d_projFFT.Zdim,
			referenceAux.d_projFFT.Ndim);

    bool power2yx, power2x;
	if (referenceAux.MFrealSpace.yxdim & (referenceAux.MFrealSpace.yxdim-1))
		power2yx = false;
	else
		power2yx = true;
	if (referenceAux.MFrealSpace.Xdim & (referenceAux.MFrealSpace.Xdim-1))
		power2x = false;
	else
		power2x = true;
	calculateNccKernel<<< CONVERT2DIM3(gridSize2), CONVERT2DIM3(blockSize2), 0, *stream >>>
			(myStructureAux.RefExpRealSpace.d_data, referenceAux.MFrealSpace.d_data, experimentalAux.MFrealSpace.d_data, referenceAux.MF2realSpace.d_data,
					experimentalAux.MF2realSpace.d_data, referenceAux.maskAutocorrelation.d_data, myStructureAux.d_NCC.d_data, referenceAux.MFrealSpace.nzyxdim,
					referenceAux.MFrealSpace.yxdim, referenceAux.MFrealSpace.Xdim, referenceAux.MFrealSpace.Ydim, referenceAux.maskCount, maxShift, power2yx, power2x);

	int fixPadding=0;
	if(referenceAux.XdimOrig%2==0 && referenceAux.Xdim%2==0)
		fixPadding=1;
	if(referenceAux.XdimOrig%2==0 && referenceAux.Xdim%2!=0)
		fixPadding=0;
	if(referenceAux.XdimOrig%2!=0 && referenceAux.Xdim%2==0)
		fixPadding=-1;
	if(referenceAux.XdimOrig%2!=0 && referenceAux.Xdim%2!=0)
		fixPadding=0;

	calculateMaxNew2DNew(myStructureAux.d_NCC.yxdim, myStructureAux.d_NCC.Ndim,
			myStructureAux.d_NCC.d_data, myStructureAux.d_out_max, myStructureAux.d_pos_max, myStream);

	numTh = 1024;
	int numBlk = transMat.Ndim/numTh;
	if(transMat.Ndim%numTh > 0)
		numBlk++;

	bool _power2x;
		if (myStructureAux.d_NCC.Xdim & (myStructureAux.d_NCC.Xdim-1))
			_power2x = false;
		else
			_power2x = true;
	double maxShift2 = (2*maxShift)*(2*maxShift);
	buildTranslationMatrix<<<numBlk, numTh, 0, *stream>>> (myStructureAux.d_pos_max.d_data, transMat.d_data, resultTR.d_data,
			myStructureAux.d_out_max.d_data, myStructureAux.d_NCC.d_data, myStructureAux.d_NCC.Xdim, myStructureAux.d_NCC.Ydim,
			myStructureAux.d_NCC.Ndim, myStructureAux.d_NCC.yxdim, fixPadding, maxShift2, _power2x);

	resultTR.copyMatrix(transMat, myStream);

	if(saveMaxVector)
		gpuErrchk(cudaMemcpyAsync(max_vector, myStructureAux.d_out_max.d_data, myStructureAux.d_NCC.Ndim*sizeof(float), cudaMemcpyDeviceToHost, *stream));

}


void cuda_calculate_correlation_two(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAuxTR,
		TransformMatrix<float> &transMatTR, float *max_vectorTR, int maxShift,
		mycufftHandle &myhandlePaddedTR, bool mirror, StructuresAux &myStructureAuxTR,
		myStreamHandle &myStreamTR,
		GpuCorrelationAux &experimentalAuxRT, TransformMatrix<float> &transMatRT,
		float *max_vectorRT, mycufftHandle &myhandlePaddedRT,
		StructuresAux &myStructureAuxRT, myStreamHandle &myStreamRT,
		TransformMatrix<float> &resultTR, TransformMatrix<float> &resultRT,
		mycufftHandle &ifftcb, bool saveMaxVector)
{

	cudaStream_t *streamTR = (cudaStream_t*) myStreamTR.ptr;
	cudaStream_t *streamRT = (cudaStream_t*) myStreamRT.ptr;


    myStructureAuxTR.RefExpFourier.resize(referenceAux.d_projFFT.Xdim, referenceAux.d_projFFT.Ydim,
			referenceAux.d_projFFT.Zdim, referenceAux.d_projFFT.Ndim);
    myStructureAuxTR.RefExpRealSpace.resize(referenceAux.Xdim, referenceAux.Ydim, referenceAux.d_projFFT.Zdim,
    		referenceAux.d_projFFT.Ndim);
    myStructureAuxTR.d_NCC.resize(referenceAux.Xdim, referenceAux.Ydim, referenceAux.d_projFFT.Zdim,
			referenceAux.d_projFFT.Ndim);

	myStructureAuxRT.RefExpFourierPolar.resize(referenceAux.d_projPolarFFT.Xdim, referenceAux.d_projPolarFFT.Ydim,
			referenceAux.d_projPolarFFT.Zdim, referenceAux.d_projPolarFFT.Ndim);
    myStructureAuxRT.RefExpRealSpacePolar.resize(referenceAux.XdimPolar, referenceAux.YdimPolar, referenceAux.d_projPolarFFT.Zdim,
    		referenceAux.d_projPolarFFT.Ndim);
    myStructureAuxRT.d_NCCPolar.resize(referenceAux.XdimPolar, referenceAux.YdimPolar, referenceAux.d_projPolarFFT.Zdim,
				referenceAux.d_projPolarFFT.Ndim);
    myStructureAuxRT.d_NCCPolar1D.resize(myStructureAuxRT.d_NCCPolar.Xdim,1,1,myStructureAuxRT.d_NCCPolar.Ndim);
    myStructureAuxRT.auxMax.resize(myStructureAuxRT.d_NCCPolar.Xdim,1,1,myStructureAuxRT.d_NCCPolar.Ndim);
    myStructureAuxRT.auxZero.resize(myStructureAuxRT.d_NCCPolar.Xdim,1,1,myStructureAuxRT.d_NCCPolar.Ndim);
    myStructureAuxRT.maxGpu.resize(myStructureAuxRT.d_NCCPolar1D.Ndim);



    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    referenceAux.d_projFFT.calculateGridSizeVectorized(blockSize, gridSize);


    pointwiseMultiplicationComplexKernel<<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize), 0, *streamTR >>>
			((cufftComplex*)referenceAux.d_projFFT.d_data, (cufftComplex*)experimentalAuxTR.d_projFFT.d_data,
					(cufftComplex*)myStructureAuxTR.RefExpFourier.d_data,
					referenceAux.d_projFFT.nzyxdim, referenceAux.d_projFFT.yxdim);

    XmippDim3 blockSize3(numTh, 1, 1), gridSize3;
    referenceAux.d_projPolarFFT.calculateGridSizeVectorized(blockSize3, gridSize3);


    pointwiseMultiplicationComplexKernel<<< CONVERT2DIM3(gridSize3), CONVERT2DIM3(blockSize3), 0, *streamRT >>>
		((cufftComplex*)referenceAux.d_projPolarFFT.d_data, (cufftComplex*)experimentalAuxRT.d_projPolarFFT.d_data,
					(cufftComplex*)myStructureAuxRT.RefExpFourierPolar.d_data, referenceAux.d_projPolarFFT.nzyxdim,
					referenceAux.d_projPolarFFT.yxdim);


    GpuMultidimArrayAtGpu< std::complex<float> > dull;
    myStructureAuxTR.RefExpFourier.ifftStream(myStructureAuxTR.RefExpRealSpace, myhandlePaddedTR, myStreamTR, false, dull);

    myStructureAuxRT.RefExpFourierPolar.ifftStream(myStructureAuxRT.RefExpRealSpacePolar, myhandlePaddedRT, myStreamRT, false, dull);

 	XmippDim3 blockSize2(numTh, 1, 1), gridSize2;
 	myStructureAuxTR.RefExpRealSpace.calculateGridSizeVectorized(blockSize2, gridSize2);

 	bool power2yx, power2x;
	if (referenceAux.MFrealSpace.yxdim & (referenceAux.MFrealSpace.yxdim-1))
		power2yx = false;
	else
		power2yx = true;
	if (referenceAux.MFrealSpace.Xdim & (referenceAux.MFrealSpace.Xdim-1))
		power2x = false;
	else
		power2x = true;
	calculateNccKernel<<< CONVERT2DIM3(gridSize2), CONVERT2DIM3(blockSize2), 0, *streamTR >>>
			(myStructureAuxTR.RefExpRealSpace.d_data, referenceAux.MFrealSpace.d_data, experimentalAuxTR.MFrealSpace.d_data, referenceAux.MF2realSpace.d_data,
					experimentalAuxTR.MF2realSpace.d_data, referenceAux.maskAutocorrelation.d_data, myStructureAuxTR.d_NCC.d_data, referenceAux.MFrealSpace.nzyxdim,
					referenceAux.MFrealSpace.yxdim, referenceAux.MFrealSpace.Xdim, referenceAux.MFrealSpace.Ydim, referenceAux.maskCount, maxShift, power2yx, power2x);


	int fixPadding=0;
	if(referenceAux.XdimOrig%2==0 && referenceAux.Xdim%2==0)
		fixPadding=1;
	if(referenceAux.XdimOrig%2==0 && referenceAux.Xdim%2!=0)
		fixPadding=0;
	if(referenceAux.XdimOrig%2!=0 && referenceAux.Xdim%2==0)
		fixPadding=-1;
	if(referenceAux.XdimOrig%2!=0 && referenceAux.Xdim%2!=0)
		fixPadding=0;

	numTh = 1024;
    XmippDim3 blockSize4(numTh, 1, 1), gridSize4;
    myStructureAuxRT.RefExpRealSpacePolar.calculateGridSizeVectorized(blockSize4, gridSize4);

	double maskFFTPolar = (referenceAux.XdimPolar*referenceAux.YdimPolar);
	calculateNccRotationKernel<<< CONVERT2DIM3(gridSize4), CONVERT2DIM3(blockSize4), 0, *streamRT >>>
			(myStructureAuxRT.RefExpRealSpacePolar.d_data, (cufftComplex*)referenceAux.d_projPolarFFT.d_data, (cufftComplex*)experimentalAuxRT.d_projPolarFFT.d_data,
					(cufftComplex*)referenceAux.d_projPolarSquaredFFT.d_data, (cufftComplex*)experimentalAuxRT.d_projPolarSquaredFFT.d_data,
					maskFFTPolar, myStructureAuxRT.d_NCCPolar.d_data, referenceAux.d_projPolarFFT.yxdim, myStructureAuxRT.RefExpRealSpacePolar.nzyxdim,
					myStructureAuxRT.RefExpRealSpacePolar.yxdim);

	//AJ sum along the radius
    numTh = 1024;
    int numBlk = (myStructureAuxRT.d_NCCPolar.Xdim*myStructureAuxRT.d_NCCPolar.Ndim)/numTh;
    if((myStructureAuxRT.d_NCCPolar.Xdim*myStructureAuxRT.d_NCCPolar.Ndim)%numTh!=0)
    	numBlk++;

    sumRadiusKernel<<< numBlk, numTh, 0, *streamRT >>>(myStructureAuxRT.d_NCCPolar.d_data, myStructureAuxRT.d_NCCPolar1D.d_data, myStructureAuxRT.auxMax.d_data,
    		myStructureAuxRT.auxZero.d_data, myStructureAuxRT.d_NCCPolar.Xdim*myStructureAuxRT.d_NCCPolar.Ndim, myStructureAuxRT.d_NCCPolar.Ydim,
			myStructureAuxRT.d_NCCPolar.Ndim);

	calculateMaxNew2DNew(myStructureAuxTR.d_NCC.yxdim, myStructureAuxTR.d_NCC.Ndim,
			myStructureAuxTR.d_NCC.d_data, myStructureAuxTR.d_out_max, myStructureAuxTR.d_pos_max, myStreamTR);

	calculateMaxNew2DNew(myStructureAuxRT.d_NCCPolar1D.Xdim, myStructureAuxRT.d_NCCPolar1D.Ndim, myStructureAuxRT.d_NCCPolar1D.d_data,
			myStructureAuxRT.d_out_polar_max, myStructureAuxRT.d_pos_polar_max, myStreamRT);

	numTh = 1024;
	numBlk = transMatTR.Ndim/numTh;
	if(transMatTR.Ndim%numTh > 0)
		numBlk++;

	bool _power2x;
		if (myStructureAuxTR.d_NCC.Xdim & (myStructureAuxTR.d_NCC.Xdim-1))
			_power2x = false;
		else
			_power2x = true;
	double maxShift2 = (2*maxShift)*(2*maxShift);
	buildTranslationMatrix<<<numBlk, numTh, 0, *streamTR>>> (myStructureAuxTR.d_pos_max.d_data, transMatTR.d_data, resultTR.d_data,
			myStructureAuxTR.d_out_max.d_data, myStructureAuxTR.d_NCC.d_data, myStructureAuxTR.d_NCC.Xdim, myStructureAuxTR.d_NCC.Ydim,
			myStructureAuxTR.d_NCC.Ndim, myStructureAuxTR.d_NCC.yxdim, fixPadding, maxShift2, _power2x);

	numBlk = transMatRT.Ndim/numTh;
	if(transMatRT.Ndim%numTh > 0)
		numBlk++;

	bool __power2x;
		if (myStructureAuxRT.d_NCCPolar1D.Xdim & (myStructureAuxRT.d_NCCPolar1D.Xdim-1))
			__power2x = false;
		else
			__power2x = true;
	buildRotationMatrix<<<numBlk, numTh, 0, *streamRT>>> (myStructureAuxRT.d_pos_polar_max.d_data, transMatRT.d_data,
			resultRT.d_data, myStructureAuxRT.maxGpu.d_data, myStructureAuxRT.auxMax.d_data, myStructureAuxRT.auxZero.d_data,
			myStructureAuxRT.d_NCCPolar1D.Xdim, myStructureAuxRT.d_NCCPolar1D.Ndim,
			myStructureAuxRT.d_NCCPolar1D.yxdim, 0, maxShift2, __power2x);


	resultTR.copyMatrix(transMatTR, myStreamTR);

	resultRT.copyMatrix(transMatRT, myStreamRT);

	if(saveMaxVector){
		gpuErrchk(cudaMemcpyAsync(max_vectorTR, myStructureAuxTR.d_out_max.d_data, myStructureAuxTR.d_NCC.Ndim*sizeof(float), cudaMemcpyDeviceToHost, *streamTR));
		gpuErrchk(cudaMemcpyAsync(max_vectorRT, myStructureAuxRT.maxGpu.d_data, myStructureAuxRT.maxGpu.Ndim*sizeof(float), cudaMemcpyDeviceToHost, *streamRT));
	}

}




void apply_transform(GpuMultidimArrayAtGpu<float> &d_original_image, GpuMultidimArrayAtGpu<float> &d_transform_image,
		TransformMatrix<float> &transMat, myStreamHandle &myStream){

	cudaStream_t *stream = (cudaStream_t*) myStream.ptr;

	int numTh = 1024;

	int numBlk = d_transform_image.yxdim/numTh;
	if(d_transform_image.yxdim%numTh > 0)
		numBlk++;
	dim3 blockSize(numTh, 1, 1);
	dim3 gridSize(numBlk, d_transform_image.Ndim, 1);

	bool power2yx, power2x;
	if (d_original_image.yxdim & (d_original_image.yxdim-1))
		power2yx = false;
	else
		power2yx = true;
	if (d_original_image.Xdim & (d_original_image.Xdim-1))
		power2x = false;
	else
		power2x = true;
	applyTransformKernel<<< gridSize, blockSize, 9*sizeof(float), *stream >>>
			(d_original_image.d_data, d_transform_image.d_data, transMat.d_data,
					d_original_image.nzyxdim, d_original_image.yxdim, d_original_image.Xdim,
					d_original_image.Ydim, power2yx, power2x);

}



void cuda_cart2polar(GpuMultidimArrayAtGpu<float> &image, GpuMultidimArrayAtGpu<float> &polar_image,
		GpuMultidimArrayAtGpu<float> &polar2_image, bool rotate, myStreamHandle &myStream)
{
	cudaStream_t *stream = (cudaStream_t*) myStream.ptr;
    int numTh = 32;
    int numBlkx = polar_image.Xdim/numTh;
    	if(polar_image.Xdim%numTh > 0)
    		numBlkx++;
	int numBlky = polar_image.Ydim/numTh;
		if(polar_image.Ydim%numTh > 0)
			numBlky++;

    dim3 blockSize(numTh, numTh, 1);
    dim3 gridSize(numBlkx, numBlky, polar_image.Ndim);

    cart2polar <<< gridSize, blockSize, 0, *stream>>>
    		(image.d_data, polar_image.d_data, polar2_image.d_data, polar_image.Ydim, polar_image.Xdim, polar_image.Ndim, image.Ydim, image.Xdim, rotate);
}

void waitGpu (myStreamHandle &myStream, bool allStreams){
	if(!allStreams){
		cudaStream_t *stream = (cudaStream_t*) myStream.ptr;
		gpuErrchk(cudaStreamSynchronize(*stream));
	}else
		gpuErrchk(cudaDeviceSynchronize());
}

void calculateAbs (std::complex<float> *data, float *out, int size, myStreamHandle &myStream){

	cudaStream_t *stream = (cudaStream_t*) myStream.ptr;
	int numTh = 1024;
	int	numBlk = size/numTh;
	if(size%numTh > 0)
		numBlk++;
	calcAbsKernel <<< numBlk, numTh, 0, *stream>>> ((cufftComplex*)data, out, size);

}
