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

#ifndef CUDA_INTERPOLATION3D_H
#define CUDA_INTERPOLATION3D_H

#include "cuda_basic_math.h"
#include "cuda_bspline_rotation.h"
#include "cuda_check_errors.h"

// 3D float texture
texture<float, cudaTextureType3D, cudaReadModeElementType> texRefVol;


template<class floatN>
__device__ floatN cubicTex3D(texture<float, 3, cudaReadModeElementType> tex, float3 coord, float3 comp_shift, uint Xdim, uint Ydim, uint Zdim)
{
	// shift the coordinate from [0,extent] to [-0.5, extent-0.5]
	float auxX = (Xdim%2==0) ? 0.5 : (trunc((float)(Xdim/2.))/Xdim);
	float auxY = (Ydim%2==0) ? 0.5 : (trunc((float)(Ydim/2.))/Ydim);
	float auxZ = (Zdim%2==0) ? 0.5 : (trunc((float)(Zdim/2.))/Zdim);
	const float3  center = make_float3(auxX, auxY, auxZ);
	const float3 coord_grid = coord - center + comp_shift;

	//const float3 coord_grid = coord - 0.5f + comp_shift;
	const float3 index = floor(coord_grid);
	const float3 fraction = coord_grid - index;
	float3 w0, w1, w2, w3;
	bspline_weights(fraction, w0, w1, w2, w3);

	const float3 g0 = w0 + w1;
	const float3 g1 = w2 + w3;

	const float3 h0 = (w1 / g0) + make_float3(-1+auxX, -1+auxY, -1+auxZ) + index;
	const float3 h1 = (w3 / g1) + make_float3(1+auxX, 1+auxY, 1+auxZ) + index;
	//const float3 h0 = (w1 / g0) - 0.5f + index;  //h0 = w1/g0 - 1, move from [-0.5, extent-0.5] to [0, extent]
	//const float3 h1 = (w3 / g1) + 1.5f + index;  //h1 = w3/g1 + 1, move from [-0.5, extent-0.5] to [0, extent]

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
interpolate_kernel3D(float* output, uint width, uint height, uint depth, float3 extent, double* angle, float3 shift)
{
	uint x = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	uint y = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;
	uint z = __umul24(blockIdx.z, blockDim.z) + threadIdx.z;
	uint i = (width * height * z) + (width * y) + x;

	if(x>=width || y>=height || z>=depth){
		return;
	}

	float x0 = (float)x;
	float y0 = (float)y;
	float z0 = (float)z;

	//AJ to compensate the shift in texture memory
	float3 comp_shift = make_float3(0.5f, 0.5f, 0.5f);

	float x1 = (float)angle[0] * x0 + (float)angle[1] * y0 + (float)angle[2] * z0 + shift.x;
	float y1 = (float)angle[4] * x0 + (float)angle[5] * y0 + (float)angle[6] * z0 + shift.y;
	float z1 = (float)angle[8] * x0 + (float)angle[9] * y0 + (float)angle[10] * z0 + shift.z;

	output[i] = cubicTex3D<float>(texRefVol, make_float3(x1, y1, z1), comp_shift, width, height, depth);

}


void interpolate3D(uint width, uint height, uint depth, double* angle, float* output)
{
	// Prepare the geometry
	float xOrigin = trunc((float)width/2);
	float yOrigin = trunc((float)height/2);
	float zOrigin = trunc((float)depth/2);

	float x0 = (float)angle[0] * xOrigin + (float)angle[1] * yOrigin + (float)angle[2] * zOrigin + (float)angle[3];
	float y0 = (float)angle[4] * xOrigin + (float)angle[5] * yOrigin + (float)angle[6] * zOrigin + (float)angle[7];
	float z0 = (float)angle[8] * xOrigin + (float)angle[9] * yOrigin + (float)angle[10] * zOrigin + (float)angle[11];

	float xShift = xOrigin - x0;
	float yShift = yOrigin - y0;
	float zShift = zOrigin - z0;

	double* d_angle;
	gpuErrchk(cudaMalloc((void**)&d_angle, 16 * sizeof(double)));
	gpuErrchk(cudaMemcpy(d_angle, angle, 16 * sizeof(double), cudaMemcpyHostToDevice));

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

	//AJ to obtain a shift with the same direction that in transform_geometry
	float3 shift = make_float3((float)xShift, (float)yShift, (float)zShift);
	float3 extent = make_float3((float)width, (float)height, (float)depth);
	interpolate_kernel3D<<<gridSize, blockSize>>>(output, width, height, depth, extent, d_angle, shift);

	cudaDeviceSynchronize();
	gpuErrchk(cudaFree(d_angle));

}

#endif

