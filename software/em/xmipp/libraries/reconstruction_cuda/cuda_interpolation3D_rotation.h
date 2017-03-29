/***************************************************************************
 * Authors:     Amaya Jimenez (ajimenez@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/



#ifndef CUDABASICMATH
#define CUDABASICMATH

#include "cuda_basic_math.h"

#endif

#ifndef CUDABSPLINE
#define CUDABSPLINE

#include "cuda_bspline_rotation.h"

#endif

// 3D float texture
texture<float, cudaTextureType3D, cudaReadModeElementType> texRefVol;


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
interpolate_kernel3D(float* output, uint width, uint height, float3 extent, double* angle, float3 shift)
{
	uint x = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	uint y = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;
	uint z = __umul24(blockIdx.z, blockDim.z) + threadIdx.z;
	uint i = (width * height * z) + (width * y) + x;

	float x0 = (float)x;
	float y0 = (float)y;
	float z0 = (float)z;

	float x1 = (float)angle[0] * x0 + (float)angle[1] * y0 + (float)angle[2] * z0 + shift.x;
	float y1 = (float)angle[3] * x0 + (float)angle[4] * y0 + (float)angle[5] * z0 + shift.y;
	float z1 = (float)angle[6] * x0 + (float)angle[7] * y0 + (float)angle[8] * z0 + shift.z;

	output[i] = cubicTex3D<float>(texRefVol, make_float3(x1, y1, z1));

}


void interpolate3D(uint width, uint height, uint depth, double* angle, float* output)
{
	// Prepare the geometry
	float xOrigin = floor(width/2);
	float yOrigin = floor(height/2);
	float zOrigin = floor(depth/2);

	float x0 = (float)angle[0] * xOrigin + (float)angle[1] * yOrigin + (float)angle[2] * zOrigin;
	float y0 = (float)angle[3] * xOrigin + (float)angle[4] * yOrigin + (float)angle[5] * zOrigin;
	float z0 = (float)angle[6] * xOrigin + (float)angle[7] * yOrigin + (float)angle[8] * zOrigin;

	float xShift = xOrigin - x0;
	float yShift = yOrigin - y0;
	float zShift = zOrigin - z0;

	double* d_angle;
	cudaMalloc((void**)&d_angle, 9 * sizeof(double));
	cudaMemcpy(d_angle, angle, 9 * sizeof(double), cudaMemcpyHostToDevice);

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
	interpolate_kernel3D<<<gridSize, blockSize>>>(output, width, height, extent, d_angle, shift);

	cudaDeviceSynchronize();
	cudaFree(d_angle);

}
