
//Host includes
#include <iostream>
#include <stdio.h>
#include <math.h>

//CUDA includes
#include <cuda_runtime.h>

#include "cuda_copy_data.h"
#include "cuda_gpu_rotate_image_v2.h"
#include "cuda_prefiltering_rotation.h"
#include "cuda_interpolation3D_rotation.h"
#include "cuda_interpolation2D_rotation.h"


//CUDA function
void cuda_rotate_image_v2(float *image, float *rotated_image, size_t Xdim, size_t Ydim, size_t Zdim, double *ang){


	//std::cerr  << "Inside CUDA function " << ang << std::endl;

	//CUDA code

	struct cudaPitchedPtr bsplineCoeffs, cudaOutput;
	bsplineCoeffs = CopyVolumeHostToDevice(image, (uint)Xdim, (uint)Ydim, (uint)Zdim);

	// Init texture
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaArray* cuArray;


    if(Zdim==1){

    	//Filtering process (first step)
    	CubicBSplinePrefilter2D((float*)bsplineCoeffs.ptr, (uint)bsplineCoeffs.pitch, (uint)Xdim, (uint)Ydim);

    	cudaMallocArray(&cuArray, &channelDesc, Xdim, Ydim);
    	// Copy to device memory some data located at address h_data in host memory
    	cudaMemcpy2DToArray(cuArray, 0, 0, bsplineCoeffs.ptr, bsplineCoeffs.pitch, Xdim * sizeof(float), Ydim, cudaMemcpyDeviceToDevice);

    	// Bind the array to the texture reference
    	cudaBindTextureToArray(texRef, cuArray, channelDesc);
    	cudaFree(bsplineCoeffs.ptr);

    	// Specify texture object parameters
    	texRef.filterMode = cudaFilterModeLinear;
    	texRef.normalized = false;

    	//Interpolation (second step)
    	cudaOutput = interpolate2D(Xdim, Ydim, ang);
    	cudaDeviceSynchronize();

    	CopyVolumeDeviceToHost(rotated_image, cudaOutput, Xdim, Ydim, Zdim);


    }else if (Zdim>1){

    	//Filtering process (first step)
    	CubicBSplinePrefilter3D((float*)bsplineCoeffs.ptr, (uint)bsplineCoeffs.pitch, (uint)Xdim, (uint)Ydim, (uint)Zdim);

    	cudaExtent volumeExtent = make_cudaExtent(Xdim, Ydim, Zdim);
    	cudaMalloc3DArray(&cuArray, &channelDesc, volumeExtent);
    	cudaMemcpy3DParms p = {0};
    	p.extent   = volumeExtent;
    	p.srcPtr   = bsplineCoeffs;
    	p.dstArray = cuArray;
    	p.kind     = cudaMemcpyDeviceToDevice;
    	cudaMemcpy3D(&p);
    	// bind array to 3D texture
    	cudaBindTextureToArray(texRefVol, cuArray, channelDesc);
    	cudaFree(bsplineCoeffs.ptr);

    	// Specify texture object parameters
  		texRefVol.filterMode = cudaFilterModeLinear;
    	texRefVol.normalized = false;

    	float *d_output;
    	cudaMalloc((void **)&d_output, Xdim * Ydim * Zdim * sizeof(float));

    	//Interpolation (second step)
    	interpolate3D(Xdim, Ydim, Zdim, ang, d_output);
    	cudaDeviceSynchronize();

    	cudaMemcpy(rotated_image, d_output, Xdim * Ydim * Zdim * sizeof(float), cudaMemcpyDeviceToHost);
    	cudaFree(d_output);

    }


    cudaFree(cuArray);

}


