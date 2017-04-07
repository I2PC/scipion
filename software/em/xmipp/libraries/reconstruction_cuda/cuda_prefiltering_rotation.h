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



#ifndef CUDA_PREFILTER_H
#define CUDA_PREFILTER_H

//#include "cuda_basic_math.h"



#include <stdio.h>

template<class floatN>
__device__ floatN InitialCausalCoefficient(
	floatN* c,			// coefficients
	uint DataLength,	// number of coefficients
	int step)			// element interleave in bytes
{
	const uint Horizon = UMIN(12, DataLength);

	// this initialization corresponds to clamping boundaries
	// accelerated loop
	float zn = Pole;
	float Sum = *c;
	for (uint n = 0; n < Horizon; n++) {
		Sum += zn * *c;
		zn *= Pole;
		c = (floatN*)((uchar*)c + step);
	}
	return(Sum);
}

template<class floatN>
__device__ floatN InitialAntiCausalCoefficient(
	floatN* c,			// last coefficient
	uint DataLength,	// number of samples or coefficients
	int step)			// element interleave in bytes
{
	// this initialization corresponds to clamping boundaries
	return((Pole / (Pole - 1.0f)) * *c);
}

template<class floatN>
__device__ void ConvertToInterpolationCoefficients(
	floatN* coeffs,		// input samples --> output coefficients
	uint DataLength,	// number of samples or coefficients
	int step)			// element interleave in bytes
{
	// compute the overall gain
	const float Lambda = (1.0f - Pole) * (1.0f - 1.0f / Pole);

	// causal initialization
	floatN* c = coeffs;
	floatN previous_c;  //cache the previously calculated c rather than look it up again (faster!)
	*c = previous_c = Lambda * InitialCausalCoefficient(c, DataLength, step);
	// causal recursion
	for (uint n = 1; n < DataLength; n++) {
		c = (floatN*)((uchar*)c + step);
		*c = previous_c = Lambda * *c + Pole * previous_c;
	}
	// anticausal initialization
	*c = previous_c = InitialAntiCausalCoefficient(c, DataLength, step);
	// anticausal recursion
	for (int n = DataLength - 2; 0 <= n; n--) {
		c = (floatN*)((uchar*)c - step);
		*c = previous_c = Pole * (previous_c - *c);
	}
}

template<class floatN>
__global__ void SamplesToCoefficients2DX(
	floatN* image,		// in-place processing
	uint pitch,			// width in bytes
	uint width,			// width of the image
	uint height)		// height of the image
{
	// process lines in x-direction
	const uint y = blockIdx.x * blockDim.x + threadIdx.x;
	floatN* line = (floatN*)((uchar*)image + y * pitch);  //direct access

	ConvertToInterpolationCoefficients(line, width, sizeof(floatN));
}

template<class floatN>
__global__ void SamplesToCoefficients2DY(
	floatN* image,		// in-place processing
	uint pitch,			// width in bytes
	uint width,			// width of the image
	uint height)		// height of the image
{
	// process lines in x-direction
	const uint x = blockIdx.x * blockDim.x + threadIdx.x;
	floatN* line = image + x;  //direct access

	ConvertToInterpolationCoefficients(line, height, pitch);
}



template<class floatN>
__global__ void SamplesToCoefficients3DX(
	floatN* volume,		// in-place processing
	uint pitch,			// width in bytes
	uint width,			// width of the volume
	uint height,		// height of the volume
	uint depth)			// depth of the volume
{
	// process lines in x-direction
	const uint y = blockIdx.x * blockDim.x + threadIdx.x;
	const uint z = blockIdx.y * blockDim.y + threadIdx.y;
	const uint startIdx = (z * height + y) * pitch;

	floatN* ptr = (floatN*)((uchar*)volume + startIdx);
	ConvertToInterpolationCoefficients(ptr, width, sizeof(floatN));
}

template<class floatN>
__global__ void SamplesToCoefficients3DY(
	floatN* volume,		// in-place processing
	uint pitch,			// width in bytes
	uint width,			// width of the volume
	uint height,		// height of the volume
	uint depth)			// depth of the volume
{
	// process lines in y-direction
	const uint x = blockIdx.x * blockDim.x + threadIdx.x;
	const uint z = blockIdx.y * blockDim.y + threadIdx.y;
	const uint startIdx = z * height * pitch;

	floatN* ptr = (floatN*)((uchar*)volume + startIdx);
	ConvertToInterpolationCoefficients(ptr + x, height, pitch);
}

template<class floatN>
__global__ void SamplesToCoefficients3DZ(
	floatN* volume,		// in-place processing
	uint pitch,			// width in bytes
	uint width,			// width of the volume
	uint height,		// height of the volume
	uint depth)			// depth of the volume
{
	// process lines in z-direction
	const uint x = blockIdx.x * blockDim.x + threadIdx.x;
	const uint y = blockIdx.y * blockDim.y + threadIdx.y;
	const uint startIdx = y * pitch;
	const uint slice = height * pitch;

	floatN* ptr = (floatN*)((uchar*)volume + startIdx);
	ConvertToInterpolationCoefficients(ptr + x, depth, slice);
}


template<class floatN>
void CubicBSplinePrefilter2D(floatN* image, uint pitch, uint width, uint height)
{

	dim3 dimBlockX(min(PowTwoDivider(height), 64));
	dim3 dimGridX(height / dimBlockX.x);
	SamplesToCoefficients2DX<floatN><<<dimGridX, dimBlockX>>>(image, pitch, width, height);

	dim3 dimBlockY(min(PowTwoDivider(width), 64));
	dim3 dimGridY(width / dimBlockY.x);
	SamplesToCoefficients2DY<floatN><<<dimGridY, dimBlockY>>>(image, pitch, width, height);

}


template<class floatN>
extern void CubicBSplinePrefilter3D(floatN* volume, uint pitch, uint width, uint height, uint depth)
{
	// Try to determine the optimal block dimensions
	uint dimX = min(min(PowTwoDivider(width), PowTwoDivider(height)), 64);
	uint dimY = min(min(PowTwoDivider(depth), PowTwoDivider(height)), 512/dimX);
	dim3 dimBlock(dimX, dimY);

	// Replace the voxel values by the b-spline coefficients
	dim3 dimGridX(height / dimBlock.x, depth / dimBlock.y);
	SamplesToCoefficients3DX<floatN><<<dimGridX, dimBlock>>>(volume, pitch, width, height, depth);

	dim3 dimGridY(width / dimBlock.x, depth / dimBlock.y);
	SamplesToCoefficients3DY<floatN><<<dimGridY, dimBlock>>>(volume, pitch, width, height, depth);

	dim3 dimGridZ(width / dimBlock.x, height / dimBlock.y);
	SamplesToCoefficients3DZ<floatN><<<dimGridZ, dimBlock>>>(volume, pitch, width, height, depth);

}
#endif
