
//Host includes
#include "cuda_gpu_rotate_image_linear_nearest.h"

#include <iostream>
#include <stdio.h>
//CUDA includes
#include <cuda_runtime.h>

#include "cuda_copy_data.h"
#include "cuda_interpolation2D_basic_rotation.h"
#include "cuda_interpolation3D_basic_rotation.h"
#include "cuda_check_errors.h"


void cuda_rotate_image_linear_nearest(float *image, float *rotated_image, size_t Xdim,
		size_t Ydim, size_t Zdim, double* ang, int interp, int wrap, int first_call){

		//CUDA code
		size_t matSize=Xdim*Ydim*Zdim*sizeof(float);

		if(first_call==1){
			gpuErrchk(cudaSetDevice(0));
			gpuErrchk(cudaFree(0));
		}

//#define TIME

		/*//AJ prueba con pinned memory
		float *h_image;
		gpuErrchk( cudaMallocHost((void**)&h_image, matSize) );
		memcpy(h_image, image, matSize);*/


#ifdef TIME
	clock_t t_ini8, t_fin8;
	double secs8;
	t_ini8 = clock();
#endif


	struct cudaPitchedPtr bsplineCoeffs;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaArray* cuArray;

	bsplineCoeffs = CopyVolumeHostToDevice(image, (uint)Xdim, (uint)Ydim, (uint)Zdim);

	if(Zdim==1){

		// Init texture
		gpuErrchk(cudaMallocArray(&cuArray, &channelDesc, Xdim, Ydim));
		// Copy to device memory some data located at address h_data in host memory
		gpuErrchk(cudaMemcpy2DToArray(cuArray, 0, 0, bsplineCoeffs.ptr, bsplineCoeffs.pitch, Xdim * sizeof(float), Ydim, cudaMemcpyDeviceToDevice));
		// Bind the array to the texture reference
		gpuErrchk(cudaBindTextureToArray(texRefBasic, cuArray, channelDesc));

		// Specify texture object parameters
		texRefBasic.addressMode[0] = (cudaTextureAddressMode)wrap;
		texRefBasic.addressMode[1] = (cudaTextureAddressMode)wrap;
		if (interp==0){
			texRefBasic.filterMode = cudaFilterModePoint;
		}else{
			texRefBasic.filterMode = cudaFilterModeLinear;
		}
	    texRefBasic.normalized = true;


	}else if (Zdim>1){

		cudaExtent volumeExtent = make_cudaExtent(Xdim, Ydim, Zdim);
		gpuErrchk(cudaMalloc3DArray(&cuArray, &channelDesc, volumeExtent));
		cudaMemcpy3DParms p = {0};
		p.extent   = volumeExtent;
		p.srcPtr   = bsplineCoeffs;
		p.dstArray = cuArray;
    	p.kind     = cudaMemcpyDeviceToDevice;
    	gpuErrchk(cudaMemcpy3D(&p));
    	// bind array to 3D texture
    	gpuErrchk(cudaBindTextureToArray(texRefVolBasic, cuArray, channelDesc));
    	// Specify texture object parameters
    	texRefVolBasic.addressMode[0] = (cudaTextureAddressMode)wrap;
    	texRefVolBasic.addressMode[1] = (cudaTextureAddressMode)wrap;
    	texRefVolBasic.addressMode[2] = (cudaTextureAddressMode)wrap;
    	if (interp==0){
    		texRefVolBasic.filterMode = cudaFilterModePoint;
    	}else{
    		texRefVolBasic.filterMode = cudaFilterModeLinear;
    	}
   		texRefVolBasic.normalized = true;


    }
	gpuErrchk(cudaFree(bsplineCoeffs.ptr));

    // Allocate result of transformation in device memory
    float *d_output;
    gpuErrchk(cudaMalloc((void **)&d_output, matSize));
	double* d_angle;
	int angle_size = (Zdim==1) ? 9 : 16;
	gpuErrchk(cudaMalloc((void**)&d_angle, angle_size * sizeof(double)));
	gpuErrchk(cudaMemcpy(d_angle, ang, angle_size * sizeof(double), cudaMemcpyHostToDevice));


#ifdef TIME
    t_fin8 = clock();
    secs8 = (double)(t_fin8 - t_ini8) / CLOCKS_PER_SEC;
    printf("Memory host to device : %.16g milisegundos\n", secs8 * 1000.0);
#endif

#ifdef TIME
	clock_t t_ini9, t_fin9;
	double secs9;
	t_ini9 = clock();
#endif

	//Kernel
	if(Zdim==1){

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
		rotate_kernel_normalized_2D<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, d_angle);


	}else if(Zdim>1){

		int numTh = 10;
		const dim3 blockSize(numTh, numTh, numTh);
		int numBlkx = (int)(Xdim)/numTh;
		if((Xdim)%numTh>0){
			numBlkx++;
		}
		int numBlky = (int)(Ydim)/numTh;
		if((Ydim)%numTh>0){
			numBlky++;
		}
		int numBlkz = (int)(Zdim)/numTh;
		if((Zdim)%numTh>0){
			numBlkz++;
		}
		const dim3 gridSize(numBlkx, numBlky, numBlkz);
		rotate_kernel_normalized_3D<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, Zdim, d_angle);


	}
	gpuErrchk(cudaDeviceSynchronize());

#ifdef TIME
    t_fin9 = clock();
    secs9 = (double)(t_fin9 - t_ini9) / CLOCKS_PER_SEC;
    printf("Kernel : %.16g milisegundos\n", secs9 * 1000.0);
#endif

#ifdef TIME
	clock_t t_ini10, t_fin10;
	double secs10;
	t_ini10 = clock();
#endif

	gpuErrchk(cudaMemcpy(rotated_image, d_output, matSize, cudaMemcpyDeviceToHost));

	gpuErrchk(cudaFreeArray(cuArray));
	gpuErrchk(cudaFree(d_output));
	gpuErrchk(cudaFree(d_angle));

#ifdef TIME
    t_fin10 = clock();
    secs10 = (double)(t_fin10 - t_ini10) / CLOCKS_PER_SEC;
    printf("Memory device to host : %.16g milisegundos\n", secs10 * 1000.0);
#endif

    //AJ prueba con pinned memory
    //gpuErrchk(cudaFreeHost(h_image));

}
