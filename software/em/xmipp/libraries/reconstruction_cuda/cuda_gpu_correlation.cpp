
//Host includes
#include "cuda_gpu_correlation.h"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <algorithm>


//CUDA includes
#include <cufft.h>
#include <cuComplex.h>

#include "cuda_basic_math.h"
#include "cuda_copy_data.h"
#include "cuda_check_errors.h"
#include "cuda_utils.h"
#include "cuda_memory_check.h"

#define PI 3.14159265


void cuda_check_gpu_memory(float* data){
	float3 memData = getAvailableMem();
	//printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n", memData.z/1024.0/1024.0, memData.y/1024.0/1024.0, memData.x/1024.0/1024.0);
	data[0]=memData.x; //total
	data[1]=memData.y; //free
	data[2]=memData.z; //used
}

__global__ void complextomodule(cufftDoubleComplex *out, double *module, double *angle, int maxSize)
{
	int x = blockDim.x * blockIdx.x + threadIdx.x;

	if(x>=maxSize)
		return;

	module[x] = log(cuCabs(out[x]));
	angle[x] = atan2(cuCimag(out[x]), cuCreal(out[x]));

}

std::complex<double>* cuda_fft(double* &d_in_data, double* recovered_image, double *module, double *angle,
		size_t Xdim, size_t Ydim, size_t Zdim, size_t batch)
{

	int NX = Ydim;
	int NY = Xdim;
	int NRANK = 2;		 			// --- 2D FFTs
	int BATCH = batch;
	int n[] = {NX, NY}; 			// --- Size of the Fourier transform
    int istride = 1;				// --- Distance between two successive input/output elements
	int ostride = 1;   			
    int idist = NX*NY;				// --- Distance between batches
	int odist = NX*((NY/2)+1);
    int inembed[] = {NX, NY};       // --- Input size with pitch
    int onembed[] = {NX, ((NY/2)+1)};

    cufftHandle planF, planB;

    cufftDoubleComplex *d_out_data;
    gpuErrchk(cudaMalloc((void**)&d_out_data, sizeof(cufftDoubleComplex)*NX*((NY/2)+1)*BATCH));
    GpuMultidimArrayAtGpu<double> d_rec_data;
    if(recovered_image!=NULL)
    	d_rec_data.resize(NX, NY, 1, BATCH);

     /* Forward transform */
    gpuErrchkFFT(cufftPlanMany(&planF, NRANK, n,
              inembed, istride, idist,
              onembed, ostride, odist, CUFFT_D2Z, BATCH));
    gpuErrchkFFT(cufftExecD2Z(planF, (cufftDoubleReal*)d_in_data, d_out_data));
    gpuErrchk(cudaDeviceSynchronize());
    cufftDestroy(planF);

    //AJ kernel to calculate the module and angle of the FFT
    if(module!=NULL & angle!=NULL){

        GpuMultidimArrayAtGpu<double> d_module(NX, ((NY/2)+1), 1, BATCH);
        GpuMultidimArrayAtGpu<double> d_angle(NX, ((NY/2)+1), 1, BATCH);

        int numTh = 1024;
        XmippDim3 blockSize(numTh, 1, 1), gridSize;
        d_module.calculateGridSizeVectorized(blockSize, gridSize);
    	complextomodule <<< CONVERT2DIM3(gridSize), CONVERT2DIM3(blockSize) >>>
    			(d_out_data, d_module.d_data, d_angle.d_data, (BATCH*NX*((NY/2)+1)));
    	gpuErrchk(cudaMemcpy(module, d_module.d_data, sizeof(double)*BATCH*NX*((NY/2)+1), cudaMemcpyDeviceToHost));
    	gpuErrchk(cudaMemcpy(angle, d_angle.d_data, sizeof(double)*BATCH*NX*((NY/2)+1), cudaMemcpyDeviceToHost));

    }

    if(recovered_image!=NULL){

    	/* Backward transform */
    	gpuErrchkFFT(cufftPlanMany(&planB, NRANK, n,
    			onembed, ostride, odist,
				inembed, istride, idist, CUFFT_Z2D, BATCH));
    	gpuErrchkFFT(cufftExecZ2D(planB, d_out_data, (cufftDoubleReal*)d_rec_data.d_data));
    	gpuErrchk(cudaDeviceSynchronize());
    	cufftDestroy(planB);

    	gpuErrchk(cudaMemcpy(recovered_image, d_rec_data.d_data, sizeof(double)*NX*NY*BATCH, cudaMemcpyDeviceToHost));
    	//recovered_image.copyFromGpu(d_rec_data);

    }

    return (std::complex<double>*)d_out_data;

}

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

double* cuda_ifft(cufftDoubleComplex* &data, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch)
{

	int NX = Ydim;
	int NY = Xdim;
	int NRANK = 2;		 			// --- 2D FFTs
	int BATCH = batch;
	int n[] = {NX, NY}; 			// --- Size of the Fourier transform
    int istride = 1;				// --- Distance between two successive input/output elements
	int ostride = 1;
    int idist = NX*NY;				// --- Distance between batches
	int odist = NX*((NY/2)+1);
    int inembed[] = {NX, NY};       // --- Input size with pitch
    int onembed[] = {NX, ((NY/2)+1)};

    cufftDoubleReal *d_rec_data;
    cufftHandle planB;

    gpuErrchk(cudaMalloc((void**)&d_rec_data, sizeof(cufftDoubleReal)*NX*NY*BATCH));

	/* Backward transform */
	gpuErrchkFFT(cufftPlanMany(&planB, NRANK, n,
			onembed, ostride, odist,
			inembed, istride, idist, CUFFT_Z2D, BATCH));
	gpuErrchkFFT(cufftExecZ2D(planB, data, d_rec_data));
	gpuErrchk(cudaDeviceSynchronize());
	cufftDestroy(planB);

	gpuErrchk(cudaFree(data));

	return (double*)d_rec_data;
}

double** cuda_calculate_correlation(std::complex<double>* d_projFFTPointer, std::complex<double>* d_projSquaredFFTPointer, std::complex<double>* d_expFFTPointer, std::complex<double>* d_expSquaredFFTPointer, std::complex<double>* d_maskFFTPointer, size_t Xdim, size_t Ydim, size_t Zdim, size_t numProj, size_t numExp, int counting){

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


