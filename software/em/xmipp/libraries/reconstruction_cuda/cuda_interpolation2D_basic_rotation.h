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

#ifndef CUDA_BASIC_INTERPOLATION2D_H
#define CUDA_BASIC_INTERPOLATION2D_H

#include "cuda_basic_math.h"
#include "cuda_bspline_basic_rotation.h"



// 2D float texture
texture<float, cudaTextureType2D, cudaReadModeElementType> texRefBasic;

/*
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
		float v = index.y + y; //AJ + (float(0.5)/(float)Ydim);
		for (float x=-1; x < 2.5f; x++)
		{
			float bsplineXY = bspline(x-fraction.x) * bsplineY;
			float u = index.x + x; //AJ + (float(0.5)/(float)Xdim);
			result += bsplineXY * tex2D(tex, u, v);
		}
	}
	return result;
}

*/


__global__ void
rotate_kernel_normalized_2D(float *output, size_t Xdim, size_t Ydim, double* angle)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    if(x>=Xdim || y>=Ydim){
    	return;
    }

    // Transform coordinates
    float u = x / (float)Xdim;
    float v = y / (float)Ydim;

    u = (Xdim%2==0) ? u-0.5 : u-(trunc((float)(Xdim/2.))/Xdim);
    v = (Ydim%2==0) ? v-0.5 : v-(trunc((float)(Ydim/2.))/Ydim);
    //u -= 0.5f; //- (float(0.5)/(float)Xdim);
    //v -= 0.5f; //+ (float(0.5)/(float)Ydim);

    float desp_u = ((float)angle[2]/(float)Xdim);
    float desp_v = ((float)angle[5]/(float)Ydim);

    //float tu = u * (float)angle[0] + v * (float)angle[1] + desp_u + 0.5f + (float(0.5)/(float)Xdim);
    //float tv = u * (float)angle[3] + v * (float)angle[4] + desp_v + 0.5f + (float(0.5)/(float)Ydim);
    float tu = u * (float)angle[0] + v * (float)angle[1] + desp_u + (float(0.5)/(float)Xdim);
    float tv = u * (float)angle[3] + v * (float)angle[4] + desp_v + (float(0.5)/(float)Ydim);
    tu = (Xdim%2==0) ? tu+0.5 : tu+(trunc((float)(Xdim/2.))/Xdim);
    tv = (Ydim%2==0) ? tv+0.5 : tv+(trunc((float)(Ydim/2.))/Ydim);

    // Read from texture and write to global memory
   	output[y * Xdim + x] = tex2D(texRefBasic, tu, tv);
}

/*
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

    float tu = u * (float)angle[0] + v * (float)angle[1] + 0.5f; //AJ + (float(0.5)/(float)Xdim);
    float tv = u * (float)angle[3] + v * (float)angle[4] + 0.5f; //AJ + (float(0.5)/(float)Ydim);

    tu = tu*(float)Xdim;
    tv = tv*(float)Ydim;

    // Read from texture and write to global memory
   	output[y * Xdim + x] = cubicTex2DSimple(texRefBasic, tu, tv);
}
*/


#endif
