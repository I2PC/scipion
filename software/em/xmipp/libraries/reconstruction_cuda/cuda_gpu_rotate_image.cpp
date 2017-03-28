
//Host includes
#include "cuda_gpu_rotate_image.h"
#include <iostream>
#include <stdio.h>
//CUDA includes
#include <cuda_runtime.h>
#include "cuda_copy_data.h"


// 2D float texture
texture<float, cudaTextureType2D, cudaReadModeElementType> texRef;

// 3D float texture
texture<float, cudaTextureType3D, cudaReadModeElementType> texRefVol;


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


__device__ float cubicTex3DSimple(texture<float, cudaTextureType3D, cudaReadModeElementType> tex, float3 coord)
{
	// transform the coordinate from [0,extent] to [-0.5, extent-0.5]
	const float3 coord_grid = make_float3(coord.x - 0.5f, coord.y - 0.5f, coord.z - 0.5f);
	float3 index = make_float3(floor(coord_grid.x), floor(coord_grid.y), floor(coord_grid.z));
	const float3 fraction = make_float3(coord_grid.x - index.x, coord_grid.y - index.y, coord_grid.z - index.z);
	index.x += 0.5f;  //move from [-0.5, extent-0.5] to [0, extent]
	index.y += 0.5f;  //move from [-0.5, extent-0.5] to [0, extent]
	index.z += 0.5f;  //move from [-0.5, extent-0.5] to [0, extent]

	float result = 0.0f;
	for (float z=-1; z < 2.5f; z++)  //range [-1, 2]
	{
		float bsplineZ = bspline(z-fraction.z);
		float w = index.z + z;
		for (float y=-1; y < 2.5f; y++)
		{
			float bsplineYZ = bspline(y-fraction.y) * bsplineZ;
			float v = index.y + y;
			for (float x=-1; x < 2.5f; x++)
			{
				float bsplineXYZ = bspline(x-fraction.x) * bsplineYZ;
				float u = index.x + x;
				result += bsplineXYZ * tex3D(tex, u, v, w);
			}
		}
	}
	return result;
}


__global__ void
rotate_kernel_normalized_2D(float *output, size_t Xdim, size_t Ydim, float ang)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    // Transform coordinates
    float u = x / (float)Xdim;
    float v = y / (float)Ydim;
    u -= 0.5f;
    v -= 0.5f;

    float tu = u * cosf(ang) - v * sinf(ang) + 0.5f;
    float tv = v * cosf(ang) + u * sinf(ang) + 0.5f;

    // Read from texture and write to global memory
   	output[y * Xdim + x] = tex2D(texRef, tu, tv);

}


__global__ void
rotate_kernel_unnormalized_2D(float *output, size_t Xdim, size_t Ydim, float ang)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    // Transform coordinates
    float u = x / (float)Xdim;
    float v = y / (float)Ydim;
    u -= 0.5f;
    v -= 0.5f;

    float tu = u * cosf(ang) - v * sinf(ang) + 0.5f;
    float tv = v * cosf(ang) + u * sinf(ang) + 0.5f;

    tu = tu*(float)Xdim;
    tv = tv*(float)Ydim;

    // Read from texture and write to global memory
   	output[y * Xdim + x] = cubicTex2DSimple(texRef, tu, tv);
}



__global__ void
rotate_kernel_normalized_3D(float *output, size_t Xdim, size_t Ydim, size_t Zdim, float ang)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;

    // Transform coordinates
    float u = x / (float)Xdim;
    float v = y / (float)Ydim;
    float w = z / (float)Zdim;
    u -= 0.5f;
    v -= 0.5f;
    w -= 0.5f;

    ang = 0;

    float tu = u * cosf(ang) - v * sinf(ang) + 0.5f;
    float tv = v * cosf(ang) + u * sinf(ang) + 0.5f;
    float tw = w * cosf(ang) + w * sinf(ang) + 0.5f;

    // Read from texture and write to global memory
    printf("Escribe hilo %i", (int)(y * Xdim + x) + (Xdim * Ydim * z));
   	output[(y * Xdim + x) + (Xdim * Ydim * z)] = tex3D(texRefVol, tu, tv, tw);

}


__global__ void
rotate_kernel_unnormalized_3D(float *output, size_t Xdim, size_t Ydim, size_t Zdim, float ang)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;
    int z = blockDim.z * blockIdx.z + threadIdx.z;

    // Transform coordinates
    float u = x / (float)Xdim;
    float v = y / (float)Ydim;
    float w = z / (float)Zdim;
    u -= 0.5f;
    v -= 0.5f;
    w -= 0.5f;

    float tu = u * cosf(ang) - v * sinf(ang) + 0.5f;
    float tv = v * cosf(ang) + u * sinf(ang) + 0.5f;
    float tw = w * cosf(ang) + w * sinf(ang) + 0.5f;

    tu = tu*(float)Xdim;
    tv = tv*(float)Ydim;
    tw = tw*(float)Zdim;

    // Read from texture and write to global memory
   	output[(y * Xdim + x) + (Xdim * Ydim * z)] = cubicTex3DSimple(texRefVol, make_float3(tu, tv, tw));
}


void cuda_rotate_image(float *image, float *rotated_image, size_t Xdim, size_t Ydim, size_t Zdim, float ang, int interp){

	std::cerr  << "Inside CUDA function " << ang << std::endl;

	//CUDA code
	size_t matSize=Xdim*Ydim*Zdim*sizeof(float);
	struct cudaPitchedPtr bsplineCoeffs, cudaOutput;
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

    }

    // Specify texture object parameters
    texRef.addressMode[0] = cudaAddressModeWrap;
    texRef.addressMode[1] = cudaAddressModeWrap;
    if(Zdim>1){
    	texRef.addressMode[2] = cudaAddressModeWrap;
    }
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



    // Allocate result of transformation in device memory
    float *d_output;
	cudaMalloc((void **)&d_output, matSize);


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
			rotate_kernel_normalized_2D<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, ang);
		}else{
			rotate_kernel_unnormalized_2D<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, ang);
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
			rotate_kernel_normalized_3D<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, Zdim, ang);
		}else{
			rotate_kernel_unnormalized_3D<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, Zdim, ang);
		}

	}
	cudaDeviceSynchronize();

	cudaMemcpy(rotated_image, d_output, matSize, cudaMemcpyDeviceToHost);


	cudaFree(cuArray);
	cudaFree(d_output);

}
