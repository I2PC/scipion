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



