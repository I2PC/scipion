
//Host includes
#include "cuda_gpu_rotate_image.h"
#include <iostream>
#include <stdio.h>
//CUDA includes
#include <cuda_runtime.h>

// 2D float texture
texture<float, cudaTextureType2D, cudaReadModeElementType> texRef;


//CUDA functions
__global__ void
rotate_kernel(float *output, size_t Xdim, size_t Ydim, float ang)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    float u = x / (float)Xdim;
    float v = y / (float)Ydim;

    // Transform coordinates
    u -= 0.5f;
    v -= 0.5f;
    float tu = u * cosf(ang) - v * sinf(ang) + 0.5f;
    float tv = v * cosf(ang) + u * sinf(ang) + 0.5f;

    // Read from texture and write to global memory
    output[y * Xdim + x] = tex2D(texRef, tu, tv);
}

void cuda_rotate_image(float *image, float *rotated_image, size_t Xdim, size_t Ydim, float ang, int interp){

	std::cerr  << "Inside CUDA function " << ang << std::endl;

	//CUDA code
	size_t matSize=Xdim*Ydim*sizeof(float);

	// Allocate CUDA array in device memory
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    cudaArray* cuArray;
    cudaMallocArray(&cuArray, &channelDesc, num, num);
    // Copy to device memory some data located at address h_data in host memory
    cudaMemcpyToArray(cuArray, 0, 0, image, matSize, cudaMemcpyHostToDevice);

    // Specify texture object parameters
    texRef.addressMode[0] = cudaAddressModeWrap;
    texRef.addressMode[1] = cudaAddressModeWrap;
    if (interp==0){
    	texRef.filterMode = cudaFilterModePoint;
    }else{
    	texRef.filterMode = cudaFilterModeLinear;
    }
    texRef.normalized = true;

    // Bind the array to the texture reference
    cudaBindTextureToArray(texRef, cuArray, channelDesc);

    // Allocate result of transformation in device memory
    float *d_output;
	cudaMalloc((void **)&d_output, matSize);


	//Kernel
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

	rotate_kernel<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, ang);

	cudaDeviceSynchronize();

	cudaMemcpy(rotated_image, d_output, matSize, cudaMemcpyDeviceToHost);


	cudaFree(cuArray);
	cudaFree(d_output);

}
