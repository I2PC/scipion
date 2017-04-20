
//Host includes
#include "cuda_gpu_rotate_image_bspline.h"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include "cuda_basic_math.h"

//CUDA includes
#include <cuda_runtime.h>
#include "cuda_copy_data.h"
#include "cuda_prefiltering_rotation.h"
#include "cuda_interpolation2D_rotation.h"
#include "cuda_interpolation3D_rotation.h"
#include "cuda_check_errors.h"
#include "cuda_utils.h"


//CUDA functions

void cuda_rotate_image_bspline(float *image, float *rotated_image, size_t Xdim, size_t Ydim, size_t Zdim,
		double *ang, int wrap, int first_call, struct ioTime* mytime){


	if(first_call==1){
		gpuErrchk(cudaSetDevice(0));
		gpuErrchk(cudaFree(0));
	}

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_ini_gpu_mem_inS, NULL);
}

	struct cudaPitchedPtr bsplineCoeffs, cudaOutput;
	bsplineCoeffs = CopyVolumeHostToDevice(image, (uint)Xdim, (uint)Ydim, (uint)Zdim);

	// Init texture
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaArray* cuArray;

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_fin_gpu_mem_inS, NULL);
    mytime->secs_gpu_mem_inS = timeval_diff(&mytime->t_fin_gpu_mem_inS, &mytime->t_ini_gpu_mem_inS);
}


    if(Zdim==1){

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_ini_gpu_filtS, NULL);
}

    	//Filtering process (first step)
    	CubicBSplinePrefilter2D((float*)bsplineCoeffs.ptr, (uint)bsplineCoeffs.pitch, (uint)Xdim, (uint)Ydim);

    	gpuErrchk(cudaMallocArray(&cuArray, &channelDesc, Xdim, Ydim));
    	// Copy to device memory some data located at address h_data in host memory
    	gpuErrchk(cudaMemcpy2DToArray(cuArray, 0, 0, bsplineCoeffs.ptr, bsplineCoeffs.pitch, Xdim * sizeof(float), Ydim, cudaMemcpyDeviceToDevice));

    	// Bind the array to the texture reference
    	gpuErrchk(cudaBindTextureToArray(texRef, cuArray, channelDesc));
    	gpuErrchk(cudaFree(bsplineCoeffs.ptr));

    	// Specify texture object parameters
    	texRef.filterMode = cudaFilterModeLinear;
    	texRef.normalized = false;
    	texRef.addressMode[0] = (cudaTextureAddressMode)wrap;
    	texRef.addressMode[1] = (cudaTextureAddressMode)wrap;

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_fin_gpu_filtS, NULL);
    mytime->secs_gpu_filtS = timeval_diff(&mytime->t_fin_gpu_filtS, &mytime->t_ini_gpu_filtS);
}

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_ini_gpu_kernelS, NULL);
}

    	//Interpolation (second step)
    	cudaOutput = interpolate2D(Xdim, Ydim, ang);
    	gpuErrchk(cudaDeviceSynchronize());

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_fin_gpu_kernelS, NULL);
    mytime->secs_gpu_kernelS = timeval_diff(&mytime->t_fin_gpu_kernelS, &mytime->t_ini_gpu_kernelS);
}

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_ini_gpu_mem_outS, NULL);
}

    	CopyVolumeDeviceToHost(rotated_image, cudaOutput, Xdim, Ydim, Zdim);

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_fin_gpu_mem_outS, NULL);
    mytime->secs_gpu_mem_outS = timeval_diff(&mytime->t_fin_gpu_mem_outS, &mytime->t_ini_gpu_mem_outS);
}

    }else if (Zdim>1){

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_ini_gpu_filtS, NULL);
}

    	//Filtering process (first step)
    	CubicBSplinePrefilter3D((float*)bsplineCoeffs.ptr, (uint)bsplineCoeffs.pitch, (uint)Xdim, (uint)Ydim, (uint)Zdim);

    	cudaExtent volumeExtent = make_cudaExtent(Xdim, Ydim, Zdim);
    	gpuErrchk(cudaMalloc3DArray(&cuArray, &channelDesc, volumeExtent));
    	cudaMemcpy3DParms p = {0};
    	p.extent   = volumeExtent;
    	p.srcPtr   = bsplineCoeffs;
    	p.dstArray = cuArray;
    	p.kind     = cudaMemcpyDeviceToDevice;
    	gpuErrchk(cudaMemcpy3D(&p));
    	// bind array to 3D texture
    	gpuErrchk(cudaBindTextureToArray(texRefVol, cuArray, channelDesc));
    	gpuErrchk(cudaFree(bsplineCoeffs.ptr));

    	// Specify texture object parameters
  		texRefVol.filterMode = cudaFilterModeLinear;
    	texRefVol.normalized = false;
    	texRefVol.addressMode[0] = (cudaTextureAddressMode)wrap;
    	texRefVol.addressMode[1] = (cudaTextureAddressMode)wrap;
    	texRefVol.addressMode[2] = (cudaTextureAddressMode)wrap;

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_fin_gpu_filtS, NULL);
    mytime->secs_gpu_filtS = timeval_diff(&mytime->t_fin_gpu_filtS, &mytime->t_ini_gpu_filtS);
}

    	float *d_output;
    	gpuErrchk(cudaMalloc((void **)&d_output, Xdim * Ydim * Zdim * sizeof(float)));



if(mytime->calcTime==true){
    gettimeofday(&mytime->t_ini_gpu_kernelS, NULL);
}

    	//Interpolation (second step)
    	interpolate3D(Xdim, Ydim, Zdim, ang, d_output);
    	gpuErrchk(cudaDeviceSynchronize());

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_fin_gpu_kernelS, NULL);
    mytime->secs_gpu_kernelS = timeval_diff(&mytime->t_fin_gpu_kernelS, &mytime->t_ini_gpu_kernelS);
}

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_ini_gpu_mem_outS, NULL);
}

		gpuErrchk(cudaMemcpy(rotated_image, d_output, Xdim * Ydim * Zdim * sizeof(float), cudaMemcpyDeviceToHost));
		gpuErrchk(cudaFree(d_output));

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_fin_gpu_mem_outS, NULL);
    mytime->secs_gpu_mem_outS = timeval_diff(&mytime->t_fin_gpu_mem_outS, &mytime->t_ini_gpu_mem_outS);
}

    }


    gpuErrchk(cudaFreeArray(cuArray));

}


