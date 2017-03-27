
//Host includes
#include <iostream>
#include <stdio.h>
#include <math.h>
//CUDA includes
#include <cuda_runtime.h>
#include "cuda_copy_data.cpp"
#include "cuda_prefiltering_rotation.cpp"
#include "cuda_gpu_rotate_image_v2.h"



// 2D float texture
texture<float, 2, cudaReadModeElementType> texRef;


//CUDA functions

__device__ void bspline_weights(float2 fraction, float2& w0, float2& w1, float2& w2, float2& w3)
{
	const float2 one_frac = 1.0f - fraction;
	const float2 squared = fraction * fraction;
	const float2 one_sqd = one_frac * one_frac;

	w0 = 1.0f/6.0f * one_sqd * one_frac;
	w1 = 2.0f/3.0f - 0.5f * squared * (2.0f-fraction);
	w2 = 2.0f/3.0f - 0.5f * one_sqd * (2.0f-one_frac);
	w3 = 1.0f/6.0f * squared * fraction;
}

template<class floatN>
__device__ floatN cubicTex2D(texture<float, 2, cudaReadModeElementType> tex, float x, float y)
{
	// transform the coordinate from [0,extent] to [-0.5, extent-0.5]
	const float2 coord_grid = make_float2(x - 0.5f, y - 0.5f);
	const float2 index = floor(coord_grid);
	const float2 fraction = coord_grid - index;
	float2 w0, w1, w2, w3;
	bspline_weights(fraction, w0, w1, w2, w3);

	const float2 g0 = w0 + w1;
	const float2 g1 = w2 + w3;
	const float2 h0 = (w1 / g0) - make_float2(0.5f) + index;  //h0 = w1/g0 - 1, move from [-0.5, extent-0.5] to [0, extent]
	const float2 h1 = (w3 / g1) + make_float2(1.5f) + index;  //h1 = w3/g1 + 1, move from [-0.5, extent-0.5] to [0, extent]

	// fetch the four linear interpolations
	floatN tex00 = tex2D(tex, h0.x, h0.y);
	floatN tex10 = tex2D(tex, h1.x, h0.y);
	floatN tex01 = tex2D(tex, h0.x, h1.y);
	floatN tex11 = tex2D(tex, h1.x, h1.y);

	// weigh along the y-direction
	tex00 = g0.y * tex00 + g1.y * tex01;
	tex10 = g0.y * tex10 + g1.y * tex11;

	// weigh along the x-direction
	return (g0.x * tex00 + g1.x * tex10);
}


__global__ void
interpolate_kernel(float* output, uint width, float2 extent, float2 a, float2 shift)
{
	uint x = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	uint y = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;
	uint i = __umul24(y, width) + x;

	float x0 = (float)x;
	float y0 = (float)y;
	float x1 = a.x * x0 - a.y * y0 + shift.x;
	float y1 = a.x * y0 + a.y * x0 + shift.y;

	output[i] = cubicTex2D<float>(texRef, x1, y1);

}


cudaPitchedPtr interpolate(uint width, uint height, float angle)
{
	// Prepare the geometry
	float2 a = make_float2((float)cos(angle), (float)sin(angle));
	float xOrigin = floor(width/2);
	float yOrigin = floor(height/2);
	float x0 = a.x * (xOrigin) - a.y * (yOrigin);
	float y0 = a.y * (xOrigin) + a.x * (yOrigin);
	float xShift = xOrigin - x0;
	float yShift = yOrigin - y0;

	// Allocate the output image
	float* output;
	cudaMalloc((void**)&output, width * height * sizeof(float));

	// Visit all pixels of the output image and assign their value
	dim3 blockSize(min(PowTwoDivider(width), 16), min(PowTwoDivider(height), 16));
	dim3 gridSize(width / blockSize.x, height / blockSize.y);
	float2 shift = make_float2((float)xShift, (float)yShift);
	float2 extent = make_float2((float)width, (float)height);
	interpolate_kernel<<<gridSize, blockSize>>>(output, width, extent, a, shift);

	return make_cudaPitchedPtr(output, width * sizeof(float), width, height);
}





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
