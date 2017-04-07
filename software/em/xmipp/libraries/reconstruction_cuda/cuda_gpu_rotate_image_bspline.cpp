
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


//CUDA functions

void cuda_rotate_image_bspline(float *image, float *rotated_image, size_t Xdim, size_t Ydim, size_t Zdim, double *ang, int wrap, int first_call){

//#define TIME

	if(first_call==1){
		gpuErrchk(cudaSetDevice(0));
		gpuErrchk(cudaFree(0));
	}


#ifdef TIME
	clock_t t_ini0, t_fin0;
	double secs0;
	t_ini0 = clock();
#endif


	struct cudaPitchedPtr bsplineCoeffs, cudaOutput;
	bsplineCoeffs = CopyVolumeHostToDevice(image, (uint)Xdim, (uint)Ydim, (uint)Zdim);

	// Init texture
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaArray* cuArray;

#ifdef TIME
    t_fin0 = clock();
    secs0 = (double)(t_fin0 - t_ini0) / CLOCKS_PER_SEC;
    printf("CopyVolumeHostToDevice: %.16g milisegundos\n", (secs0) * 1000.0);
#endif


    if(Zdim==1){

#ifdef TIME
	clock_t t_ini2, t_fin2;
	double secs2;
	t_ini2 = clock();
#endif

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

#ifdef TIME
    t_fin2 = clock();
    secs2 = (double)(t_fin2 - t_ini2) / CLOCKS_PER_SEC;
    printf("Memory host to device and prefiltering: %.16g milisegundos\n", (secs2) * 1000.0);
#endif

#ifdef TIME
	clock_t t_ini3, t_fin3;
	double secs3;
	t_ini3 = clock();
#endif

    	//Interpolation (second step)
    	cudaOutput = interpolate2D(Xdim, Ydim, ang);
    	gpuErrchk(cudaDeviceSynchronize());

#ifdef TIME
    t_fin3 = clock();
    secs3 = (double)(t_fin3 - t_ini3) / CLOCKS_PER_SEC;
    printf("Kernel (interpolation): %.16g milisegundos\n", secs3 * 1000.0);
#endif

#ifdef TIME
	clock_t t_ini4, t_fin4;
	double secs4;
	t_ini4 = clock();
#endif

    	CopyVolumeDeviceToHost(rotated_image, cudaOutput, Xdim, Ydim, Zdim);

#ifdef TIME
    t_fin4 = clock();
    secs4 = (double)(t_fin4 - t_ini4) / CLOCKS_PER_SEC;
    printf("Memory device to host: %.16g milisegundos\n", secs4 * 1000.0);
#endif


    }else if (Zdim>1){


#ifdef TIME
	clock_t t_ini5, t_fin5;
	double secs5;
	t_ini5 = clock();
#endif

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

    	float *d_output;
    	gpuErrchk(cudaMalloc((void **)&d_output, Xdim * Ydim * Zdim * sizeof(float)));

#ifdef TIME
    t_fin5 = clock();
    secs5 = (double)(t_fin5 - t_ini5) / CLOCKS_PER_SEC;
    printf("Memory host to device and prefiltering: %.16g milisegundos\n", (secs5) * 1000.0);
#endif

#ifdef TIME
	clock_t t_ini6, t_fin6;
	double secs6;
	t_ini6 = clock();
#endif

    	//Interpolation (second step)
    	interpolate3D(Xdim, Ydim, Zdim, ang, d_output);
    	gpuErrchk(cudaDeviceSynchronize());

#ifdef TIME
    t_fin6 = clock();
    secs6 = (double)(t_fin6 - t_ini6) / CLOCKS_PER_SEC;
    printf("Kernel (interpolation): %.16g milisegundos\n", secs6 * 1000.0);
#endif

#ifdef TIME
	clock_t t_ini7, t_fin7;
	double secs7;
	t_ini7 = clock();
#endif

		gpuErrchk(cudaMemcpy(rotated_image, d_output, Xdim * Ydim * Zdim * sizeof(float), cudaMemcpyDeviceToHost));
		gpuErrchk(cudaFree(d_output));

#ifdef TIME
    t_fin7 = clock();
    secs7 = (double)(t_fin7 - t_ini7) / CLOCKS_PER_SEC;
    printf("Memory device to host: %.16g milisegundos\n", secs7 * 1000.0);
#endif

    }


    gpuErrchk(cudaFreeArray(cuArray));

}


