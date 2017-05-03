
//Host includes
#include "cuda_gpu_correlation.h"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include "cuda_basic_math.h"

//CUDA includes
#include <cuda_runtime.h>
#include <cufft.h>
#include <cuComplex.h>

#include "cuda_copy_data.h"
#include "cuda_check_errors.h"
#include "cuda_utils.h"

#define PI 3.14159265

//CUDA functions

__global__ void complextomodule(cufftComplex *out, float *module, float *angle, size_t maxSize)
{
	int x = blockDim.x * blockIdx.x + threadIdx.x;

	if(x>=maxSize){
		return;
	}

	module[x] = logf(cuCabsf(out[x]));
	angle[x] = atan2(cuCimagf(out[x]), cuCrealf(out[x]));

}


__global__ void cart2polar(float *in, float *out, int maxRadio, int maxAng, int maxBatch, int Xdim, int Ydim)
{
	int angle = blockDim.x * blockIdx.x + threadIdx.x;
	int radio = blockDim.y * blockIdx.y + threadIdx.y;
	//int batch = blockDim.z * blockIdx.z + threadIdx.z;

	if(radio>=maxRadio || angle>=maxAng){
		return;	
	}

	int x = (int)llroundf((float)(radio*cosf((float)(angle*PI/180)))) + maxRadio;
	int y = (int)llroundf((float)(radio*sinf((float)(angle*PI/180)))) + maxRadio;
	
	float imVal;
	for (int i=0; i<maxBatch; i++){
		imVal = in[(y * Xdim + x)+(i*Xdim*Ydim)];
		out[(angle+(radio*maxAng))+(i*maxAng*maxRadio)] = imVal;
	}
	
}

#define DEBUG
void cuda_fft(float *image, float *recovered_image, float *module, float *angle, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch)
{

	int NX = Xdim;
	int NY = Ydim;
	int NRANK = 2;		 			// --- 2D FFTs
	int BATCH = batch;
	int n[] = {NX, NY}; 			// --- Size of the Fourier transform
    int istride = 1;				// --- Distance between two successive input/output elements
	int ostride = 1;   			
    int idist = NX*NY;				// --- Distance between batches
	int odist = NX*((NY/2)+1); 		
    int inembed[] = {NX, NY};       // --- Input size with pitch
    int onembed[] = {NX, ((NY/2)+1)};

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);

    /*size_t free_byte, total_byte;

    gpuErrchk(cudaMemGetInfo( &free_byte, &total_byte ));

    double free_db = (double)free_byte ;
    double total_db = (double)total_byte ;
    double used_db = total_db - free_db ;

    printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n", used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);*/


    cufftHandle planF, planB;
    float *d_in_data;
    cufftComplex *d_out_data;
    cufftReal *d_rec_data;

    gpuErrchk(cudaMalloc((void**)&d_in_data, sizeof(float)*NX*NY*BATCH));
    gpuErrchk(cudaMemcpy(d_in_data, image, sizeof(float)*NX*NY*BATCH, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMalloc((void**)&d_out_data, sizeof(cufftComplex)*NX*((NY/2)+1)*BATCH));
    gpuErrchk(cudaMalloc((void**)&d_rec_data, sizeof(cufftReal)*NX*NY*BATCH));


     /* Forward transform */
    gpuErrchkFFT(cufftPlanMany(&planF, NRANK, n,
              inembed, istride, idist,
              onembed, ostride, odist, CUFFT_R2C, BATCH));
    gpuErrchkFFT(cufftExecR2C(planF, (cufftReal*)d_in_data, d_out_data));
    gpuErrchk(cudaDeviceSynchronize());
    cufftDestroy(planF);


    //AJ kernel to calculate the module and angle of the FFT
    float *d_module;
    float *d_angle;
    gpuErrchk(cudaMalloc((void**)&d_module, sizeof(float)*BATCH*NX*((NY/2)+1) ));
    gpuErrchk(cudaMalloc((void**)&d_angle, sizeof(float)*BATCH*NX*((NY/2)+1) ));

    int numTh = 1024;
    const dim3 blockSize(numTh, 1, 1);
    int numBlkx = (int)(BATCH*NX*((NY/2)+1))/numTh;
    if((int)(BATCH*NX*((NY/2)+1))%numTh>0){
    	numBlkx++;
    }
    const dim3 gridSize(numBlkx, 1, 1);
    complextomodule <<< gridSize, blockSize >>> (d_out_data, d_module, d_angle, (size_t)(BATCH*NX*((NY/2)+1)));
    gpuErrchk(cudaMemcpy(module, d_module, sizeof(float)*BATCH*NX*((NY/2)+1), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(angle, d_angle, sizeof(float)*BATCH*NX*((NY/2)+1), cudaMemcpyDeviceToHost));


    /* Backward transform */
    gpuErrchkFFT(cufftPlanMany(&planB, NRANK, n,
              onembed, ostride, odist,
              inembed, istride, idist, CUFFT_C2R, BATCH));
    gpuErrchkFFT(cufftExecC2R(planB, d_out_data, d_rec_data));
    gpuErrchk(cudaDeviceSynchronize());
    cufftDestroy(planB);

    gpuErrchk(cudaMemcpy(recovered_image, d_rec_data, sizeof(float)*NX*NY*BATCH, cudaMemcpyDeviceToHost));


    gpuErrchk(cudaFree(d_in_data));
    gpuErrchk(cudaFree(d_out_data));
    gpuErrchk(cudaFree(d_rec_data));
    gpuErrchk(cudaFree(d_module));
    gpuErrchk(cudaFree(d_angle));


}



void cuda_cart2polar(float *image, float *polar_image, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch){


	//AJ suponemos que todas las imagenes son de igual tamaño
	int NX = Xdim;
	int NY = Ydim;
	int BATCH = batch;

	float *d_in_data;
	float *d_out_data;

    gpuErrchk(cudaMalloc((void**)&d_in_data, sizeof(float)*NX*NY*BATCH));
    gpuErrchk(cudaMemcpy(d_in_data, image, sizeof(float)*NX*NY*BATCH, cudaMemcpyHostToDevice));

	int radios = std::max(NX/2, NY/2);
	int angles = 360;

    gpuErrchk(cudaMalloc((void**)&d_out_data, sizeof(float)*radios*angles*BATCH));


    //AJ kernel to calculate the polar coordinates
    float *d_module;
    float *d_angle;
    gpuErrchk(cudaMalloc((void**)&d_module, sizeof(float)*BATCH*NX*((NY/2)+1) ));
    gpuErrchk(cudaMalloc((void**)&d_angle, sizeof(float)*BATCH*NX*((NY/2)+1) ));

    int numTh = 32;
    const dim3 blockSize(numTh, numTh, 1);
    int numBlkx = (int)(angles)/numTh;
    if((int)(angles)%numTh>0){
    	numBlkx++;
    }
    int numBlky = (int)(radios)/numTh;
    if((int)(radios)%numTh>0){
    	numBlky++;
    }
    const dim3 gridSize(numBlkx, numBlky, 1);
    cart2polar <<< gridSize, blockSize >>> (d_in_data, d_out_data, radios, angles, BATCH, NX, NY);
    gpuErrchk(cudaMemcpy(polar_image, d_out_data, sizeof(float)*radios*angles*BATCH, cudaMemcpyDeviceToHost));

    gpuErrchk(cudaFree(d_in_data));
    gpuErrchk(cudaFree(d_out_data));



}


