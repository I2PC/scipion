
//Host includes
#include "cuda_gpu_rotate_image_linear_nearest_threads.h"

#include <iostream>
#include <stdio.h>
//CUDA includes
#include "cuda_utils.h"
#include "cuda_interpolation2D_basic_rotation_threads.h"

#define NUM_TH 2
cudaTextureObject_t texImageTh[NUM_TH];
float *d_output[NUM_TH];
double *d_angle[NUM_TH];


void CreateTexture(float *image, size_t Xdim, size_t Ydim, size_t Zdim, int interp, int wrap, int i)
{

		float *d_inputIm;//Device Array
		size_t pitch;
		gpuErrchk(cudaMallocPitch(&d_inputIm, &pitch, sizeof(float)*Xdim, Ydim));
		gpuErrchk(cudaMemcpy2D(d_inputIm, pitch, image, sizeof(float)*Xdim, sizeof(float)*Xdim, Ydim, cudaMemcpyHostToDevice));

		cudaResourceDesc resDesc;
		memset(&resDesc, 0, sizeof(resDesc));
		resDesc.resType = cudaResourceTypePitch2D;
		resDesc.res.pitch2D.devPtr = d_inputIm;
		resDesc.res.pitch2D.pitchInBytes =  pitch;
		resDesc.res.pitch2D.width = Xdim;
		resDesc.res.pitch2D.height = Ydim;
		resDesc.res.pitch2D.desc = cudaCreateChannelDesc<float>(); // bits per channel

		cudaTextureDesc texDesc;
		memset(&texDesc, 0, sizeof(texDesc));
		texDesc.readMode = cudaReadModeElementType;
		texDesc.normalizedCoords = true;
		if (interp==0){
			texDesc.filterMode = cudaFilterModePoint;
		}else{
		    texDesc.filterMode = cudaFilterModeLinear;
		}
		texDesc.addressMode[0] = (cudaTextureAddressMode)wrap;
		texDesc.addressMode[1] = (cudaTextureAddressMode)wrap;
		gpuErrchk(cudaCreateTextureObject(&texImageTh[i], &resDesc, &texDesc, NULL));

		std::cerr << "Create texture: " << i << std::endl;

}


void cuda_rotate_image_linear_nearest_threads(float *image, float *rotated_image, size_t Xdim,
		size_t Ydim, size_t Zdim, double* ang, int interp, int wrap, int thIdx){

	cudaSetDevice(0);

	size_t free_byte, total_byte;
	gpuErrchk(cudaMemGetInfo( &free_byte, &total_byte ));
	double free_db = (double)free_byte ;
	double total_db = (double)total_byte ;
	double used_db = total_db - free_db ;
	printf("Thread %i . Before texture - GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n", thIdx, used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);


	int angle_size = 9;


	CreateTexture(image, Xdim, Ydim, Zdim, interp, wrap, thIdx);

	size_t free_byte2, total_byte2;
	gpuErrchk(cudaMemGetInfo( &free_byte2, &total_byte2 ));
	double free_db2 = (double)free_byte2 ;
	double total_db2 = (double)total_byte2 ;
	double used_db2 = total_db2 - free_db2 ;
	printf("Thread %i . After texture - GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n", thIdx, used_db2/1024.0/1024.0, free_db2/1024.0/1024.0, total_db2/1024.0/1024.0);


	gpuErrchk(cudaMalloc((void **)&d_output[thIdx], Xdim*Ydim*sizeof(float)));
	gpuErrchk(cudaMalloc((void**)&d_angle[thIdx], angle_size * sizeof(double)));
	gpuErrchk(cudaMemcpy(d_angle[thIdx], ang, angle_size * sizeof(double), cudaMemcpyHostToDevice));


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

	std::cerr << "Llamada al kernel del hilo " << thIdx << std::endl;
	std::cerr << "con blockSize de " << numTh << " y gridSize de " << numBlkx << " x " << numBlky << std::endl;

	rotate_kernel_normalized_2D_threads<<<gridSize, blockSize>>>(d_output[thIdx], Xdim, Ydim, d_angle[thIdx], texImageTh[thIdx], thIdx);

	std::cerr << "Final del kernel del hilo " << thIdx << std::endl;

	gpuErrchk(cudaDeviceSynchronize());


	gpuErrchk(cudaMemcpy(rotated_image, d_output[thIdx], Xdim*Ydim*sizeof(float), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaFree(d_output[thIdx]));
	gpuErrchk(cudaFree(d_angle[thIdx]));



}
