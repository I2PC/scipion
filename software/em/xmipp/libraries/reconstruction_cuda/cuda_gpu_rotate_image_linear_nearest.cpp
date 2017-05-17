
//Host includes
#include "cuda_gpu_rotate_image_linear_nearest.h"

#include <iostream>
#include <stdio.h>

//CUDA includes
#include "cuda_utils.h"
#include "cuda_xmipp_utils.h"
#include "cuda_interpolation2D_basic_rotation.h"
#include "cuda_interpolation3D_basic_rotation.h"


void cuda_rotate_image_linear_nearest(float *image, float *rotated_image, size_t Xdim, size_t Ydim, size_t Zdim,
		double* ang, int interp, int wrap, int first_call, struct ioTime* mytime){

		//CUDA code
		size_t matSize=Xdim*Ydim*Zdim*sizeof(float);

		if(first_call==1){
			gpuErrchk(cudaSetDevice(0));
			gpuErrchk(cudaFree(0));
		}


		/*//AJ prueba con pinned memory
		float *h_image;
		gpuErrchk( cudaMallocHost((void**)&h_image, matSize) );
		memcpy(h_image, image, matSize);*/

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_ini_gpu_mem_in, NULL);
}

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

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_fin_gpu_mem_in, NULL);
    mytime->secs_gpu_mem_in = timeval_diff(&mytime->t_fin_gpu_mem_in, &mytime->t_ini_gpu_mem_in);
}

	gpuErrchk(cudaFree(bsplineCoeffs.ptr));

    // Allocate result of transformation in device memory
    float *d_output;
    gpuErrchk(cudaMalloc((void **)&d_output, matSize));
	double* d_angle;
	int angle_size = (Zdim==1) ? 9 : 16;
	gpuErrchk(cudaMalloc((void**)&d_angle, angle_size * sizeof(double)));
	gpuErrchk(cudaMemcpy(d_angle, ang, angle_size * sizeof(double), cudaMemcpyHostToDevice));


if(mytime->calcTime==true){
    gettimeofday(&mytime->t_ini_gpu_kernel, NULL);
}

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

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_fin_gpu_kernel, NULL);
    mytime->secs_gpu_kernel = timeval_diff(&mytime->t_fin_gpu_kernel, &mytime->t_ini_gpu_kernel);
}

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_ini_gpu_mem_out, NULL);
}

	gpuErrchk(cudaMemcpy(rotated_image, d_output, matSize, cudaMemcpyDeviceToHost));

if(mytime->calcTime==true){
    gettimeofday(&mytime->t_fin_gpu_mem_out, NULL);
    mytime->secs_gpu_mem_out = timeval_diff(&mytime->t_fin_gpu_mem_out, &mytime->t_ini_gpu_mem_out);
}

	gpuErrchk(cudaFreeArray(cuArray));
	gpuErrchk(cudaFree(d_output));
	gpuErrchk(cudaFree(d_angle));



    //AJ prueba con pinned memory
    //gpuErrchk(cudaFreeHost(h_image));

}
