
//Host includes
#include "cuda_gpu_correlation.h"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include "cuda_basic_math.h"

//CUDA includes
#include <cuda_runtime.h>
#include <cufft.h>
#include <cuComplex.h>

#include "cuda_copy_data.h"
#include "cuda_check_errors.h"
#include "cuda_utils.h"



//CUDA functions

__global__ void complextomodule(cufftComplex *out, float *module, float *angle, size_t Xdim, size_t Ydim)
{
	int x = blockDim.x * blockIdx.x + threadIdx.x;
	int y = blockDim.y * blockIdx.y + threadIdx.y;

	if(x>=Xdim || y>=Ydim){
		return;
	}

	module[(y * Xdim + x)] = logf(cuCabsf(out[(y * Xdim + x)]));
	angle[(y * Xdim + x)] = atan2(cuCimagf(out[(y * Xdim + x)]), cuCrealf(out[(y * Xdim + x)]));

}

void cuda_fft(float *image, float *rotated_image, float *module, float *angle, size_t Xdim, size_t Ydim, size_t Zdim){

	int NX = Xdim;
	int NY = Ydim;
	int NRANK = 2;
	int BATCH = 1;

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);

    size_t free_byte, total_byte;

    gpuErrchk(cudaMemGetInfo( &free_byte, &total_byte ));

    double free_db = (double)free_byte ;
    double total_db = (double)total_byte ;
    double used_db = total_db - free_db ;

    //printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n", used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);


    cufftHandle planF, planB;
    float *in_data;
    cufftComplex *out_data;
    cufftReal *rec_data;


    gpuErrchk(cudaMalloc((void**)&in_data, sizeof(float)*NX*NY));
    gpuErrchk(cudaMemcpy(in_data, image, sizeof(float)*NX*NY, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMalloc((void**)&out_data, sizeof(cufftComplex)*NX*((int)(floor(NY/2)+1))));
    gpuErrchk(cudaMalloc((void**)&rec_data, sizeof(cufftReal)*NX*NY));


    /* Forward transform */
    gpuErrchkFFT(cufftPlan2d(&planF, NX, NY, CUFFT_R2C));
    gpuErrchkFFT(cufftExecR2C(planF, (cufftReal*)in_data, out_data));
    gpuErrchk(cudaDeviceSynchronize());
    cufftDestroy(planF);



    //AJ kernel to calculate the module and angle of the FFT
    float *d_module;
    float *d_angle;
    gpuErrchk(cudaMalloc((void**)&d_module, sizeof(float)*NX*((int)(floor(NY/2)+1))));
    gpuErrchk(cudaMalloc((void**)&d_angle, sizeof(float)*NX*((int)(floor(NY/2)+1))));

    int numTh = 32;
    const dim3 blockSize(numTh, numTh, 1);
    int numBlkx = (int)(Xdim)/numTh;
    if((Xdim)%numTh>0){
    	numBlkx++;
    }
    int numBlky = (int)(Ydim)/numTh;
    if((Ydim)%numTh>0){
    	numBlky++;
    }
    const dim3 gridSize(numBlkx, numBlky, 1);
    complextomodule <<< gridSize, blockSize >>> (out_data, d_module, d_angle, Xdim, Ydim);
    gpuErrchk(cudaMemcpy(module, d_module, sizeof(float)*NX*((int)(floor(NY/2)+1)), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(angle, d_angle, sizeof(float)*NX*((int)(floor(NY/2)+1)), cudaMemcpyDeviceToHost));


    /* Backward transform */
    gpuErrchkFFT(cufftPlan2d(&planB, NX, NY, CUFFT_C2R));
    gpuErrchkFFT(cufftExecC2R(planB, out_data, rec_data));
    gpuErrchk(cudaDeviceSynchronize());
    cufftDestroy(planB);

    gpuErrchk(cudaMemcpy(rotated_image, rec_data, sizeof(float)*NX*NY, cudaMemcpyDeviceToHost));


    gpuErrchk(cudaFree(in_data));
    gpuErrchk(cudaFree(out_data));
    gpuErrchk(cudaFree(rec_data));
    gpuErrchk(cudaFree(d_module));
    gpuErrchk(cudaFree(d_angle));

    /*

	gpuErrchk(cudaMalloc((void **)&d_output0, Xdim[i]*Ydim[i]*sizeof(float)));
	gpuErrchk(cudaMalloc((void**)&d_angle0, angle_size * sizeof(double)));
	gpuErrchk(cudaMemcpy(d_angle0, ang[i], angle_size * sizeof(double), cudaMemcpyHostToDevice));

	int numTh = 32;
	const dim3 blockSize(numTh, numTh, 1);
	int numBlkx = (int)(maxXdim)/numTh;
	if((maxXdim)%numTh>0){
		numBlkx++;
	}
	int numBlky = (int)(maxYdim)/numTh;
	if((maxYdim)%numTh>0){
		numBlky++;
	}
	const dim3 gridSize(numBlkx, numBlky, 1);


	rotate_kernel_normalized_2D_parallel<<<gridSize, blockSize, 0, stream[0]>>>(d_output0, Xdim[0], Ydim[0], d_angle0, texImage[0]);
	gpuErrchk(cudaMemcpyAsync(rotated_image[0], d_output0, Xdim[0]*Ydim[0]*sizeof(float), cudaMemcpyDeviceToHost, stream[0]));

	gpuErrchk(cudaDeviceSynchronize());
	gpuErrchk(cudaFree(d_output0));
	gpuErrchk(cudaFree(d_angle0));*/

}


