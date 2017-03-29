
//Host includes
#include <iostream>
#include <stdio.h>
#include <math.h>
//CUDA includes
#include <cuda_runtime.h>

#include "cuda_copy_data.h"
#include "cuda_prefiltering_rotation.cpp"
#include "cuda_gpu_rotate_image_v2.h"



// 2D float texture
texture<float, cudaTextureType2D, cudaReadModeElementType> texRef;

// 3D float texture
texture<float, cudaTextureType3D, cudaReadModeElementType> texRefVol;


//CUDA functions
template<class floatN>
__device__ void bspline_weights(floatN fraction, floatN& w0, floatN& w1, floatN& w2, floatN& w3)
{
	const floatN one_frac = 1.0f - fraction;
	const floatN squared = fraction * fraction;
	const floatN one_sqd = one_frac * one_frac;

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


template<class floatN>
__device__ floatN cubicTex3D(texture<float, 3, cudaReadModeElementType> tex, float3 coord)
{
	// shift the coordinate from [0,extent] to [-0.5, extent-0.5]
	const float3 coord_grid = coord - 0.5f;
	const float3 index = floor(coord_grid);
	const float3 fraction = coord_grid - index;
	float3 w0, w1, w2, w3;
	bspline_weights(fraction, w0, w1, w2, w3);

	const float3 g0 = w0 + w1;
	const float3 g1 = w2 + w3;
	const float3 h0 = (w1 / g0) - 0.5f + index;  //h0 = w1/g0 - 1, move from [-0.5, extent-0.5] to [0, extent]
	const float3 h1 = (w3 / g1) + 1.5f + index;  //h1 = w3/g1 + 1, move from [-0.5, extent-0.5] to [0, extent]

	// fetch the eight linear interpolations
	// weighting and fetching is interleaved for performance and stability reasons
	floatN tex000 = tex3D(tex, h0.x, h0.y, h0.z);
	floatN tex100 = tex3D(tex, h1.x, h0.y, h0.z);
	tex000 = g0.x * tex000 + g1.x * tex100;  //weigh along the x-direction
	floatN tex010 = tex3D(tex, h0.x, h1.y, h0.z);
	floatN tex110 = tex3D(tex, h1.x, h1.y, h0.z);
	tex010 = g0.x * tex010 + g1.x * tex110;  //weigh along the x-direction
	tex000 = g0.y * tex000 + g1.y * tex010;  //weigh along the y-direction
	floatN tex001 = tex3D(tex, h0.x, h0.y, h1.z);
	floatN tex101 = tex3D(tex, h1.x, h0.y, h1.z);
	tex001 = g0.x * tex001 + g1.x * tex101;  //weigh along the x-direction
	floatN tex011 = tex3D(tex, h0.x, h1.y, h1.z);
	floatN tex111 = tex3D(tex, h1.x, h1.y, h1.z);
	tex011 = g0.x * tex011 + g1.x * tex111;  //weigh along the x-direction
	tex001 = g0.y * tex001 + g1.y * tex011;  //weigh along the y-direction

	return (g0.z * tex000 + g1.z * tex001);  //weigh along the z-direction
}



__global__ void
interpolate_kernel2D(float* output, uint width, float2 extent, float2 a, float2 shift)
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


__global__ void
interpolate_kernel3D(float* output, uint width, uint height, float3 extent, float2 a, float3 shift)
{
	uint x = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	uint y = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;
	uint z = __umul24(blockIdx.z, blockDim.z) + threadIdx.z;
	uint i = (width * height * z) + (width * y) + x;

	float x0 = (float)x;
	float y0 = (float)y;
	float z0 = (float)z;
	//AJ aqui se hace una transformacion simple en un solo eje
	float x1 = a.x * x0 - a.y * y0 + shift.x;
	float y1 = a.x * y0 + a.y * x0 + shift.y;
	float z1 = z0 + shift.z;

	output[i] = cubicTex3D<float>(texRefVol, make_float3(x1, y1, z1));

}


cudaPitchedPtr interpolate2D(uint width, uint height, float angle)
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
	interpolate_kernel2D<<<gridSize, blockSize>>>(output, width, extent, a, shift);

	return make_cudaPitchedPtr(output, width * sizeof(float), width, height);
}


void interpolate3D(uint width, uint height, uint depth, float angle, float* output)
{
	// Prepare the geometry
	//AJ mal, esto es una transformacion simple en un solo eje, hay que generalizar
	float2 a = make_float2((float)cos(angle), (float)sin(angle));
	float xOrigin = floor(width/2);
	float yOrigin = floor(height/2);
	float zOrigin = floor(depth/2);

	float x0 = a.x * (xOrigin) - a.y * (yOrigin);
	float y0 = a.y * (xOrigin) + a.x * (yOrigin);
	float z0 = zOrigin;

	float xShift = xOrigin - x0;
	float yShift = yOrigin - y0;
	float zShift = zOrigin - z0;

	// Visit all pixels of the output image and assign their value
	int numTh = 10;
	const dim3 blockSize(numTh, numTh, numTh);
	int numBlkx = (int)(width)/numTh;
	if((width)%numTh>0){
		numBlkx++;
	}
	int numBlky = (int)(height)/numTh;
	if((height)%numTh>0){
		numBlky++;
	}
	int numBlkz = (int)(depth)/numTh;
	if((depth)%numTh>0){
		numBlkz++;
	}
	const dim3 gridSize(numBlkx, numBlky, numBlkz);

	float3 shift = make_float3((float)xShift, (float)yShift, (float)zShift);
	float3 extent = make_float3((float)width, (float)height, (float)depth);
	interpolate_kernel3D<<<gridSize, blockSize>>>(output, width, height, extent, a, shift);

}


void cuda_rotate_image_v2(float *image, float *rotated_image, size_t Xdim, size_t Ydim, size_t Zdim, float ang){

	std::cerr  << "Inside CUDA function " << ang << std::endl;

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


