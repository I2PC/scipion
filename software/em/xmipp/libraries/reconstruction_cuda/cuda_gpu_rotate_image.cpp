
//Host includes
#include "cuda_gpu_rotate_image.h"
#include <iostream>
#include <stdio.h>
//CUDA includes
#include <cuda_runtime.h>

// 2D float texture
texture<float, cudaTextureType2D, cudaReadModeElementType> texRef;


//CUDA functions

// Cubic B-spline function
// The 3rd order Maximal Order and Minimum Support function, that it is maximally differentiable.
__device__ float bspline(float t)
{
	t = fabs(t);
	const float a = 2.0f - t;

	if (t < 1.0f) return 2.0f/3.0f - 0.5f*t*t*a;
	else if (t < 2.0f) return a*a*a / 6.0f;
	else return 0.0f;
}


//! Bicubic interpolated texture lookup, using unnormalized coordinates.
//! Straight forward implementation, using 16 nearest neighbour lookups.
//! @param tex  2D texture
//! @param x  unnormalized x texture coordinate
//! @param y  unnormalized y texture coordinate
__device__ float cubicTex2DSimple(texture<float, cudaTextureType2D, cudaReadModeElementType> tex, float x, float y)
{
	// transform the coordinate from [0,extent] to [-0.5, extent-0.5]
	const float2 coord_grid = make_float2(x - 0.5f, y - 0.5f);
	float2 index = make_float2(floor(coord_grid.x), floor(coord_grid.y));
	const float2 fraction = make_float2(coord_grid.x - index.x, coord_grid.y - index.y);
	index.x += 0.5f;  //move from [-0.5, extent-0.5] to [0, extent]
	index.y += 0.5f;  //move from [-0.5, extent-0.5] to [0, extent]

	float result = 0.0f;
	for (float y=-1; y < 2.5f; y++)
	{
		float bsplineY = bspline(y-fraction.y);
		float v = index.y + y;
		for (float x=-1; x < 2.5f; x++)
		{
			float bsplineXY = bspline(x-fraction.x) * bsplineY;
			float u = index.x + x;
			result += bsplineXY * tex2D(tex, u, v);
		}
	}
	return result;
}


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
    cudaMallocArray(&cuArray, &channelDesc, Xdim, Ydim);
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
