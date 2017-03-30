
//Host includes
#include "cuda_gpu_rotate_image.h"
#include <iostream>
#include <stdio.h>
//CUDA includes
#include <cuda_runtime.h>
//#include "cuda_copy_data.h"




//CUDA functions

void cuda_rotate_image(float *image, float *rotated_image, size_t Xdim, size_t Ydim, size_t Zdim, double* ang, int interp){

	//std::cerr  << "Inside CUDA function " << ang << std::endl;

	//CUDA code
	size_t matSize=Xdim*Ydim*Zdim*sizeof(float);
	struct cudaPitchedPtr bsplineCoeffs;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaArray* cuArray;

	bsplineCoeffs = CopyVolumeHostToDevice(image, (uint)Xdim, (uint)Ydim, (uint)Zdim);

	if(Zdim==1){

		// Init texture
		cudaMallocArray(&cuArray, &channelDesc, Xdim, Ydim);
		// Copy to device memory some data located at address h_data in host memory
		cudaMemcpy2DToArray(cuArray, 0, 0, bsplineCoeffs.ptr, bsplineCoeffs.pitch, Xdim * sizeof(float), Ydim, cudaMemcpyDeviceToDevice);
		// Bind the array to the texture reference
		cudaBindTextureToArray(texRef, cuArray, channelDesc);

		// Specify texture object parameters
		//texRef.addressMode[0] = cudaAddressModeBorder;
		//texRef.addressMode[1] = cudaAddressModeBorder;
		if (interp==0){
		    texRef.filterMode = cudaFilterModePoint;
		}else{
		    texRef.filterMode = cudaFilterModeLinear;
		}
		if (interp<2){
		    texRef.normalized = true;
		}else{
			texRef.normalized = false;
		}

	}else if (Zdim>1){

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

    	// Specify texture object parameters
    	//texRefVol.addressMode[0] = cudaAddressModeWrap;
    	//texRefVol.addressMode[1] = cudaAddressModeWrap;
    	//texRefVol.addressMode[2] = cudaAddressModeWrap;
    	if (interp==0){
    		texRefVol.filterMode = cudaFilterModePoint;
    	}else{
    		texRefVol.filterMode = cudaFilterModeLinear;
    	}
    	if (interp<2){
    		texRefVol.normalized = true;
    	}else{
    		texRefVol.normalized = false;
    	}

    }
	cudaFree(bsplineCoeffs.ptr);

    // Allocate result of transformation in device memory
    float *d_output;
	cudaMalloc((void **)&d_output, matSize);
	double* d_angle;
	cudaMalloc((void**)&d_angle, 9 * sizeof(double));
	cudaMemcpy(d_angle, ang, 9 * sizeof(double), cudaMemcpyHostToDevice);


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
		if(interp<2){
			rotate_kernel_normalized_2D<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, d_angle);
		}else{
			rotate_kernel_unnormalized_2D<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, d_angle);
		}

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
		if(interp<2){
			rotate_kernel_normalized_3D<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, Zdim, d_angle);
		}else{
			rotate_kernel_unnormalized_3D<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, Zdim, d_angle);
		}

	}
	cudaDeviceSynchronize();

	cudaMemcpy(rotated_image, d_output, matSize, cudaMemcpyDeviceToHost);

	cudaFree(cuArray);
	cudaFree(d_output);
	cudaFree(d_angle);

}
