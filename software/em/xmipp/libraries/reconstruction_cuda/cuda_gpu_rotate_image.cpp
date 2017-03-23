
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
rotate_kernel(float *output, cudaTextureObject_t texObj, int num, int ang)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if(x<num || y<num){
    	return;
    }

    float u = x / (float)num;
    float v = y / (float)num;

    // Transform coordinates
    u -= 0.5f;
    v -= 0.5f;
    float tu = u * cosf(ang) - v * sinf(ang) + 0.5f;
    float tv = v * cosf(ang) + u * sinf(ang) + 0.5f;

    // Read from texture and write to global memory
    output[y * num + x] = tex2D<float>(texRef, tu, tv);
}

void cuda_rotate_image(float *image, float *rotated_image, int ang){

	std::cerr  << "Inside CUDA function " << ang << std::endl;

	//CUDA code
	int num=2;
	size_t matSize=num*num*sizeof(float);

	// Allocate CUDA array in device memory
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    cudaArray* cuArray;
    cudaMallocArray(&cuArray, &channelDesc, num, num);
    // Copy to device memory some data located at address h_data in host memory
    cudaMemcpyToArray(cuArray, 0, 0, image, matSize, cudaMemcpyHostToDevice);

    // Specify texture
    struct cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = cuArray;

    // Specify texture object parameters
    struct cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.addressMode[0]   = cudaAddressModeWrap;
    texDesc.addressMode[1]   = cudaAddressModeWrap;
    texDesc.filterMode       = cudaFilterModePoint;
    texDesc.readMode         = cudaReadModeElementType;
    texDesc.normalizedCoords = 1;

    // Create texture object
    cudaTextureObject_t texObj = 0;
    cudaCreateTextureObject(&texRef, &resDesc, &texDesc, NULL);

    // Allocate result of transformation in device memory
    float *d_output;
	cudaMalloc((void **)&d_output, matSize);


	//Kernel
	int numTh = 32;
	const dim3 blockSize(numTh, numTh, 1);
	int numBlkx = (int)(num)/numTh;
	if((num)%numTh>0){
		numBlkx++;
	}
	int numBlky = (int)(num)/numTh;
	if((num)%numTh>0){
	  numBlky++;
	}
	const dim3 gridSize(numBlkx, numBlky, 1);

	rotate_kernel<<<gridSize, blockSize>>>(d_output, texRef, num, ang);

	cudaMemcpy(rotated_image, d_output, matSize, cudaMemcpyDeviceToHost);

	// Destroy texture object
	cudaDestroyTextureObject(texRef);

	cudaFree(cuArray);
	cudaFree(d_output);

}
