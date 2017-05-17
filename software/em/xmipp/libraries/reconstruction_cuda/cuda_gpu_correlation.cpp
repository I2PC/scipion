
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

__global__ void pointwiseMultiplicationKernel_complex(cufftDoubleComplex *M, cufftDoubleComplex *manyF, cufftDoubleComplex *MmanyF,
		size_t nzyxdim, size_t yxdim)
{
	size_t idx = blockDim.x * blockIdx.x + threadIdx.x;
	size_t idxLow = idx % yxdim;

	if (idx>=nzyxdim)
		return;
//
//	int idxDimProj = Xdim*Ydim*numProj;
//	int idxDim = Xdim*Ydim;
//	float normFactor = (1.0/(Xdim*Ydim));
//	int outIdx = 0;
//	int inIdx = 0;
//
//	for(int n=0; n<numExp; n++){
//		d_out_proj_exp[idx + outIdx] = make_cuDoubleComplex( cuCreal(cuCmul(projFFT[idx], expFFT[idxLow + inIdx]))*normFactor ,  cuCimag(cuCmul(projFFT[idx], expFFT[idxLow + inIdx]))*normFactor ) ;
//		d_out_proj_mask[idx + outIdx] = make_cuDoubleComplex( cuCreal(cuCmul(projFFT[idx], maskFFT[idxLow]))*normFactor ,  cuCimag(cuCmul(projFFT[idx], maskFFT[idxLow]))*normFactor ) ;
//		d_out_mask_exp[idx + outIdx] = make_cuDoubleComplex( cuCreal(cuCmul(maskFFT[idxLow], expFFT[idxLow + inIdx]))*normFactor ,  cuCimag(cuCmul(maskFFT[idxLow], expFFT[idxLow + inIdx]))*normFactor ) ;
//		d_out_mask_mask[idx + outIdx] = make_cuDoubleComplex( cuCreal(cuCmul(maskFFT[idxLow], maskFFT[idxLow]))*normFactor ,  cuCimag(cuCmul(maskFFT[idxLow], maskFFT[idxLow]))*normFactor ) ;
//		d_out_projSquared_mask[idx + outIdx] = make_cuDoubleComplex( cuCreal(cuCmul(projSqFFT[idx], maskFFT[idxLow]))*normFactor ,  cuCimag(cuCmul(projSqFFT[idx], maskFFT[idxLow]))*normFactor ) ;
//		d_out_mask_expSquared[idx + outIdx] = make_cuDoubleComplex( cuCreal(cuCmul(maskFFT[idxLow], expSqFFT[idxLow + inIdx]))*normFactor ,  cuCimag(cuCmul(maskFFT[idxLow], expSqFFT[idxLow + inIdx]))*normFactor ) ;
//		//d_out_proj_exp[idx + outIdx] = expFFT[idxLow + n*Xdim*Ydim];
//		outIdx+=idxDimProj;
//		inIdx+=idxDim;
//	}

}

void pointwiseMultiplicationFourier(const GpuMultidimArrayAtGpu< std::complex<double> > &M, const GpuMultidimArrayAtGpu < std::complex<double> >& manyF,
		GpuMultidimArrayAtGpu< std::complex<double> > &MmanyF)
{
    int numTh = 1024;
    XmippDim3 blockSize(numTh, 1, 1), gridSize;
    manyF.calculateGridSizeVectorized(blockSize, gridSize);

	pointwiseMultiplicationKernel_complex <<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
			((cufftDoubleComplex*)M.d_data, (cufftDoubleComplex*)manyF.d_data, (cufftDoubleComplex*) MmanyF.d_data, manyF.nzyxdim, manyF.yxdim);
}

void GpuCorrelationAux::produceSideInfo()
{
	GpuMultidimArrayAtGpu< std::complex<double> > MF, MF2;
	MF.resize(d_projFFT);
	MF2.resize(d_projFFT);


	d_projSquaredFFT.clear();
}

#ifdef NEVERDEFINED
__global__ void calculate_ncc(double *d_rec_proj_exp, double *d_rec_proj_mask, double *d_rec_mask_exp, double *d_rec_mask_mask,
		double *d_rec_projSquared_mask, double *d_rec_mask_expSquared, double *NCC, int Xdim, int Ydim, int numProj, int counting)
{

	int idx = blockDim.x * blockIdx.x + threadIdx.x;

	if(idx>=(Xdim*Ydim*numProj))
		return;

	double num = (d_rec_proj_exp[idx] - ((d_rec_proj_mask[idx]*d_rec_mask_exp[idx])/(d_rec_mask_mask[idx])) );
	double den1 = sqrt( d_rec_projSquared_mask[idx] - (d_rec_proj_mask[idx]*d_rec_proj_mask[idx]/d_rec_mask_mask[idx]) );
	double den2 = sqrt( d_rec_mask_expSquared[idx] - (d_rec_mask_exp[idx]*d_rec_mask_exp[idx]/d_rec_mask_mask[idx]) );

	if(den1!=0 && den2!=0 && d_rec_mask_mask[idx]>counting*0.1)
		NCC[idx] = num/(den1*den2);
	else
		NCC[idx] = -1;

}

__global__ void multiply_kernel(cufftDoubleComplex *projFFT, cufftDoubleComplex *projSqFFT, cufftDoubleComplex *expFFT,
		cufftDoubleComplex *expSqFFT, cufftDoubleComplex *maskFFT, cufftDoubleComplex *d_out_proj_exp, cufftDoubleComplex *d_out_proj_mask,
		cufftDoubleComplex *d_out_mask_exp, cufftDoubleComplex *d_out_mask_mask, cufftDoubleComplex *d_out_projSquared_mask,
		cufftDoubleComplex *d_out_mask_expSquared, int Xdim, int Ydim, int numProj, int numExp)
{
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	int idxLow = fmod((double)idx, (double)Xdim*Ydim);

	if(idx>=(Xdim*Ydim*numProj))
		return;

	int idxDimProj = Xdim*Ydim*numProj;
	int idxDim = Xdim*Ydim;
	float normFactor = (1.0/(Xdim*Ydim));
	int outIdx = 0;
	int inIdx = 0;

	for(int n=0; n<numExp; n++){
		d_out_proj_exp[idx + outIdx] = make_cuDoubleComplex( cuCreal(cuCmul(projFFT[idx], expFFT[idxLow + inIdx]))*normFactor ,  cuCimag(cuCmul(projFFT[idx], expFFT[idxLow + inIdx]))*normFactor ) ;
		d_out_proj_mask[idx + outIdx] = make_cuDoubleComplex( cuCreal(cuCmul(projFFT[idx], maskFFT[idxLow]))*normFactor ,  cuCimag(cuCmul(projFFT[idx], maskFFT[idxLow]))*normFactor ) ;
		d_out_mask_exp[idx + outIdx] = make_cuDoubleComplex( cuCreal(cuCmul(maskFFT[idxLow], expFFT[idxLow + inIdx]))*normFactor ,  cuCimag(cuCmul(maskFFT[idxLow], expFFT[idxLow + inIdx]))*normFactor ) ;
		d_out_mask_mask[idx + outIdx] = make_cuDoubleComplex( cuCreal(cuCmul(maskFFT[idxLow], maskFFT[idxLow]))*normFactor ,  cuCimag(cuCmul(maskFFT[idxLow], maskFFT[idxLow]))*normFactor ) ;
		d_out_projSquared_mask[idx + outIdx] = make_cuDoubleComplex( cuCreal(cuCmul(projSqFFT[idx], maskFFT[idxLow]))*normFactor ,  cuCimag(cuCmul(projSqFFT[idx], maskFFT[idxLow]))*normFactor ) ;
		d_out_mask_expSquared[idx + outIdx] = make_cuDoubleComplex( cuCreal(cuCmul(maskFFT[idxLow], expSqFFT[idxLow + inIdx]))*normFactor ,  cuCimag(cuCmul(maskFFT[idxLow], expSqFFT[idxLow + inIdx]))*normFactor ) ;
		//d_out_proj_exp[idx + outIdx] = expFFT[idxLow + n*Xdim*Ydim];
		outIdx+=idxDimProj;
		inIdx+=idxDim;
	}

}

double** cuda_calculate_correlation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux){
	//Following the formula (21) in Padfiled's paper, proj=F1 and exp=F2
	cufftDoubleComplex* d_out_proj_exp;
    gpuErrchk(cudaMalloc((void**)&d_out_proj_exp, sizeof(cufftDoubleComplex)*((Xdim/2)+1)*Ydim*numProj*numExp));
	cufftDoubleComplex* d_out_proj_mask;
    gpuErrchk(cudaMalloc((void**)&d_out_proj_mask, sizeof(cufftDoubleComplex)*((Xdim/2)+1)*Ydim*numProj*numExp));
	cufftDoubleComplex* d_out_mask_exp;
    gpuErrchk(cudaMalloc((void**)&d_out_mask_exp, sizeof(cufftDoubleComplex)*((Xdim/2)+1)*Ydim*numProj*numExp));
	cufftDoubleComplex* d_out_mask_mask;
    gpuErrchk(cudaMalloc((void**)&d_out_mask_mask, sizeof(cufftDoubleComplex)*((Xdim/2)+1)*Ydim*numProj*numExp));
	cufftDoubleComplex* d_out_projSquared_mask;
    gpuErrchk(cudaMalloc((void**)&d_out_projSquared_mask, sizeof(cufftDoubleComplex)*((Xdim/2)+1)*Ydim*numProj*numExp));
	cufftDoubleComplex* d_out_mask_expSquared;
    gpuErrchk(cudaMalloc((void**)&d_out_mask_expSquared, sizeof(cufftDoubleComplex)*((Xdim/2)+1)*Ydim*numProj*numExp));

    int numTh = 1024;
    const dim3 blockSize(numTh, 1, 1);
    int numBlkx = (int)(((Xdim/2)+1)*Ydim*numProj)/numTh;
    if((int)(((Xdim/2)+1)*Ydim*numProj)%numTh>0){
    	numBlkx++;
    }
    const dim3 gridSize(numBlkx, 1, 1);
	multiply_kernel<<< gridSize, blockSize >>>
			((cufftDoubleComplex*)d_projFFTPointer, (cufftDoubleComplex*)d_projSquaredFFTPointer, (cufftDoubleComplex*)d_expFFTPointer,
					(cufftDoubleComplex*)d_expSquaredFFTPointer, (cufftDoubleComplex*)d_maskFFTPointer, d_out_proj_exp, d_out_proj_mask,
					d_out_mask_exp, d_out_mask_mask, d_out_projSquared_mask, d_out_mask_expSquared, ((Xdim/2)+1), Ydim, numProj, numExp);

	double *d_rec_proj_exp = cuda_ifft(d_out_proj_exp, Xdim, Ydim, 1, numProj*numExp);
	double *d_rec_proj_mask = cuda_ifft(d_out_proj_mask, Xdim, Ydim, 1, numProj*numExp);
	double *d_rec_mask_exp = cuda_ifft(d_out_mask_exp, Xdim, Ydim, 1, numProj*numExp);
	double *d_rec_mask_mask = cuda_ifft(d_out_mask_mask, Xdim, Ydim, 1, numProj*numExp);
	double *d_rec_projSquared_mask = cuda_ifft(d_out_projSquared_mask, Xdim, Ydim, 1, numProj*numExp);
	double *d_rec_mask_expSquared = cuda_ifft(d_out_mask_expSquared, Xdim, Ydim, 1, numProj*numExp);

	int pointer=0;
	numBlkx = (int)(Xdim*Ydim*numProj)/numTh;
	if((int)(Xdim*Ydim*numProj)%numTh>0){
		numBlkx++;
	}
	const dim3 gridSize2(numBlkx, 1, 1);
	double **NCC = new double* [numExp];
	GpuMultidimArrayAtGpu<double> d_NCC(Xdim, Ydim, 1, numProj);
	for(int n=0;n<numExp;n++){

		NCC[n] = new double [Xdim*Ydim*numProj];

		calculate_ncc<<< gridSize2, blockSize >>>
				(&d_rec_proj_exp[pointer], &d_rec_proj_mask[pointer], &d_rec_mask_exp[pointer], &d_rec_mask_mask[pointer],
						&d_rec_projSquared_mask[pointer], &d_rec_mask_expSquared[pointer], d_NCC.d_data, Xdim, Ydim, numProj, counting);

		gpuErrchk(cudaMemcpy(NCC[n], d_NCC.d_data, sizeof(double)*Ydim*Xdim*numProj, cudaMemcpyDeviceToHost));

		pointer+=Ydim*Xdim*numProj;
	}

	gpuErrchk(cudaFree(d_rec_proj_exp));
	gpuErrchk(cudaFree(d_rec_proj_mask));
	gpuErrchk(cudaFree(d_rec_mask_exp));
	gpuErrchk(cudaFree(d_rec_mask_mask));
	gpuErrchk(cudaFree(d_rec_projSquared_mask));
	gpuErrchk(cudaFree(d_rec_mask_expSquared));

    return NCC;

}
#endif

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
