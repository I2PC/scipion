
//Host includes
#include <iostream>
#include <stdio.h>
#include <math.h>
//CUDA includes
#include <cuda_runtime.h>
#include "cuda_copy_data.cpp"
#include "cuda_prefiltering_rotation.cpp"
#include "cuda_interpolation_rotation.cpp"
#include "cuda_gpu_rotate_image_v2.h"



// 2D float texture
texture<float, 2, cudaReadModeElementType> texRef;


//CUDA functions

void cuda_rotate_image_v2(float *image, float *rotated_image, size_t Xdim, size_t Ydim, float ang){

	std::cerr  << "Inside CUDA function " << ang << std::endl;

	//CUDA code
	//size_t matSize=Xdim*Ydim*sizeof(float);

	//Filtering process (first step)
	struct cudaPitchedPtr bsplineCoeffs, cudaOutput;
	bsplineCoeffs = CopyVolumeHostToDevice(image, (uint)Xdim, (uint)Ydim, 1); //AJ el 1 habra que cambiarlo cuando sea volumen
	CubicBSplinePrefilter2D((float*)bsplineCoeffs.ptr, (uint)bsplineCoeffs.pitch, (uint)Xdim, (uint)Ydim);

	// Init texture
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaArray* cuArray;
    cudaMallocArray(&cuArray, &channelDesc, Xdim, Ydim);
    // Copy to device memory some data located at address h_data in host memory
    cudaMemcpy2DToArray(cuArray, 0, 0, bsplineCoeffs.ptr, bsplineCoeffs.pitch, Xdim * sizeof(float), Ydim, cudaMemcpyDeviceToDevice);

    // Bind the array to the texture reference
    cudaBindTextureToArray(texRef, cuArray, channelDesc);

    // Specify texture object parameters
    texRef.filterMode = cudaFilterModeLinear;
    texRef.normalized = false;

    //Interpolation (second step)
    cudaOutput = interpolate(Xdim, Ydim, ang);

    CopyVolumeDeviceToHost(rotated_image, cudaOutput, Xdim, Ydim, 1); //AJ el 1 habra que cambiarlo cuando sea volumen

    cudaFree(cuArray);

}


