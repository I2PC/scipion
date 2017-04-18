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

#ifndef CUDA_INTERPOLATION2D_H
#define CUDA_INTERPOLATION2D_H

#include "cuda_basic_math.h"
#include "cuda_bspline_rotation.h"
#include "cuda_check_errors.h"



// 2D float texture
texture<float, cudaTextureType2D, cudaReadModeElementType> texRef;


template<class floatN>
__device__ floatN cubicTex2D(texture<float, 2, cudaReadModeElementType> tex, float x, float y, float2 comp_shift, uint Xdim, uint Ydim)
{
	// transform the coordinate from [0,extent] to [-0.5, extent-0.5]

	float auxX = (Xdim%2==0) ? 0.5 : (trunc((float)(Xdim/2.))/Xdim);
	float auxY = (Ydim%2==0) ? 0.5 : (trunc((float)(Ydim/2.))/Ydim);
	const float2 coord_grid = make_float2(x - auxX + comp_shift.x, y - auxY + comp_shift.y);
	//const float2 coord_grid = make_float2(x - 0.5f + comp_shift.x, y - 0.5f + comp_shift.y);

	const float2 index = floor(coord_grid);
	const float2 fraction = coord_grid - index;

	float2 w0, w1, w2, w3;
	bspline_weights(fraction, w0, w1, w2, w3);

	const float2 g0 = w0 + w1;
	const float2 g1 = w2 + w3;
	const float2 auxh0 = (w1 / g0) + index;

	const float2 h0 = (w1 / g0) + make_float2(-1+auxX, -1+auxY) + index;
	const float2 h1 = (w3 / g1) + make_float2(1+auxX, 1+auxY) + index;
	//const float2 h0 = (w1 / g0) - make_float2(0.5f) + index;  //h0 = w1/g0 - 1, move from [-0.5, extent-0.5] to [0, extent]
	//const float2 h1 = (w3 / g1) + make_float2(1.5f) + index;  //h1 = w3/g1 + 1, move from [-0.5, extent-0.5] to [0, extent]


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
interpolate_kernel2D(float* output, uint width, uint height, double* angle, float2 shift)
{
	uint x = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	uint y = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;
	uint i = __umul24(y, width) + x;

	if(x>=width || y>=height){
		return;
	}

	//for(int i=0; i<1000; i++){

	float x0 = (float)x;
	float y0 = (float)y;

	//AJ to compensate the shift in texture memory
	float2 comp_shift = make_float2(0.5f, 0.5f);

	float x1 = (float)angle[0] * x0 + (float)angle[1] * y0 + shift.x;
	float y1 = (float)angle[3] * x0 + (float)angle[4] * y0 + shift.y;

	output[i] = cubicTex2D<float>(texRef, x1, y1, comp_shift, width, height);

	//}

}


cudaPitchedPtr interpolate2D(uint width, uint height, double* angle)
{
	// Prepare the geometry
	float xOrigin = trunc((float)width/2);
	float yOrigin = trunc((float)height/2);
	float x0 = (float)angle[0] * (xOrigin) + (float)angle[1] * (yOrigin) + (float)angle[2];
	float y0 = (float)angle[3] * (xOrigin) + (float)angle[4] * (yOrigin) + (float)angle[5];
	float xShift = xOrigin - x0;
	float yShift = yOrigin - y0;

	// Allocate the output image
	float* output;
	gpuErrchk(cudaMalloc((void**)&output, width * height * sizeof(float)));
	double* d_angle;
	gpuErrchk(cudaMalloc((void**)&d_angle, 9 * sizeof(double)));
	gpuErrchk(cudaMemcpy(d_angle, angle, 9 * sizeof(double), cudaMemcpyHostToDevice));


	// Visit all pixels of the output image and assign their value
	dim3 blockSize(min(PowTwoDivider(width), 16), min(PowTwoDivider(height), 16));
	dim3 gridSize(width / blockSize.x, height / blockSize.y);
	//AJ to obtain a shift with the same direction that in transform_geometry
	float2 shift = make_float2((float)xShift, (float)yShift);

	interpolate_kernel2D<<<gridSize, blockSize>>>(output, width, height, d_angle, shift);

	cudaDeviceSynchronize();
	gpuErrchk(cudaFree(d_angle));

	return make_cudaPitchedPtr(output, width * sizeof(float), width, height);
}
#endif

