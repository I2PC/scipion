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

// 2D float texture
texture<float, cudaTextureType2D, cudaReadModeElementType> texRef;

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
interpolate_kernel2D(float* output, uint width, float2 extent, double* angle, float2 shift)
{
	uint x = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	uint y = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;
	uint i = __umul24(y, width) + x;

	float x0 = (float)x;
	float y0 = (float)y;

	float x1 = (float)angle[0] * x0 + (float)angle[1] * y0 + shift.x;
	float y1 = (float)angle[3] * x0 + (float)angle[4] * y0 + shift.y;

	output[i] = cubicTex2D<float>(texRef, x1, y1);

}

cudaPitchedPtr interpolate2D(uint width, uint height, double* angle)
{
	// Prepare the geometry
	float xOrigin = floor(width/2);
	float yOrigin = floor(height/2);
	float x0 = (float)angle[0] * (xOrigin) + (float)angle[1] * (yOrigin);
	float y0 = (float)angle[3] * (xOrigin) + (float)angle[4] * (yOrigin);
	float xShift = xOrigin - x0;
	float yShift = yOrigin - y0;

	// Allocate the output image
	float* output;
	cudaMalloc((void**)&output, width * height * sizeof(float));
	double* d_angle;
	cudaMalloc((void**)&d_angle, 9 * sizeof(double));
	cudaMemcpy(d_angle, angle, 9 * sizeof(double), cudaMemcpyHostToDevice);


	// Visit all pixels of the output image and assign their value
	dim3 blockSize(min(PowTwoDivider(width), 16), min(PowTwoDivider(height), 16));
	dim3 gridSize(width / blockSize.x, height / blockSize.y);
	float2 shift = make_float2((float)xShift, (float)yShift);
	float2 extent = make_float2((float)width, (float)height);
	interpolate_kernel2D<<<gridSize, blockSize>>>(output, width, extent, d_angle, shift);

	cudaDeviceSynchronize();
	cudaFree(d_angle);

	return make_cudaPitchedPtr(output, width * sizeof(float), width, height);
}

