
//Host includes
#include "cuda_gpu_rotate_image.h"
#include <iostream>
#include <stdio.h>
//CUDA includes
#include <cuda_runtime.h>
#include "cuda_copy_data.h"
#include "cuda_interpolation2D_rotation.h"


// 2D float texture
texture<float, cudaTextureType2D, cudaReadModeElementType> texRef;

// 3D float texture
texture<float, cudaTextureType3D, cudaReadModeElementType> texRefVol;


//CUDA functions





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
rotate_kernel_normalized_2D(float *output, size_t Xdim, size_t Ydim, double* angle)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    // Transform coordinates
    float u = x / (float)Xdim;
    float v = y / (float)Ydim;
    u -= 0.5f;
    v -= 0.5f;

    float tu = u * (float)angle[0] + v * (float)angle[1] + 0.5f;
    float tv = u * (float)angle[3] + v * (float)angle[4] + 0.5f;

    // Read from texture and write to global memory
   	output[y * Xdim + x] = tex2D(texRef, tu, tv);

}


__global__ void
rotate_kernel_unnormalized_2D(float *output, size_t Xdim, size_t Ydim, double* angle)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    // Transform coordinates
    float u = x / (float)Xdim;
    float v = y / (float)Ydim;
    u -= 0.5f;
    v -= 0.5f;

    float tu = u * (float)angle[0] + v * (float)angle[1] + 0.5f;
    float tv = u * (float)angle[3] + v * (float)angle[4] + 0.5f;

    tu = tu*(float)Xdim;
    tv = tv*(float)Ydim;

    // Read from texture and write to global memory
   	output[y * Xdim + x] = cubicTex2DSimple(texRef, tu, tv);
}



__global__ void
rotate_kernel_normalized_3D(float *output, size_t Xdim, size_t Ydim, size_t Zdim, double* angle)
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

    float tu = u * (float)angle[0] + v * (float)angle[1] + w * (float)angle[2] + 0.5f;
    float tv = u * (float)angle[3] + v * (float)angle[4] + w * (float)angle[5] + 0.5f;
    float tw = u * (float)angle[6] + v * (float)angle[7] + w * (float)angle[8] + 0.5f;

    // Read from texture and write to global memory
   	output[(y * Xdim + x) + (Xdim * Ydim * z)] = tex3D(texRefVol, tu, tv, tw);

}


__global__ void
rotate_kernel_unnormalized_3D(float *output, size_t Xdim, size_t Ydim, size_t Zdim, double* angle)
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

    float tu = u * (float)angle[0] + v * (float)angle[1] + w * (float)angle[2] + 0.5f;
    float tv = u * (float)angle[3] + v * (float)angle[4] + w * (float)angle[5] + 0.5f;
    float tw = u * (float)angle[6] + v * (float)angle[7] + w * (float)angle[8] + 0.5f;

    tu = tu*(float)Xdim;
    tv = tv*(float)Ydim;
    tw = tw*(float)Zdim;

    // Read from texture and write to global memory
   	output[(y * Xdim + x) + (Xdim * Ydim * z)] = cubicTex3DSimple(texRefVol, make_float3(tu, tv, tw));
}


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
