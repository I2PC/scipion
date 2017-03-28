/***************************************************************************
 * Authors:     AUTHOR_NAME (amaia@cnb.csic.es)
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


// float2 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ float2 make_float2(float s);

inline __host__ __device__ float2 make_float2(int2 a);

// addition
inline __host__ __device__ float2 operator+(float2 a, float2 b);
inline __host__ __device__ void operator+=(float2 &a, float2 b);

// subtract
inline __host__ __device__ float2 operator-(float2 a, float2 b);
inline __host__ __device__ void operator-=(float2 &a, float2 b);


// multiply
inline __host__ __device__ float2 operator*(float2 a, float2 b);
inline __host__ __device__ float2 operator*(float2 a, float s);
inline __host__ __device__ float2 operator*(float s, float2 a);
inline __host__ __device__ void operator*=(float2 &a, float s);

// divide
inline __host__ __device__ float2 operator/(float2 a, float2 b);
inline __host__ __device__ float2 operator/(float2 a, float s);
inline __host__ __device__ float2 operator/(float s, float2 a);  //Danny
inline __host__ __device__ void operator/=(float2 &a, float s);

// dot product
inline __host__ __device__ float dot(float2 a, float2 b);

// length
inline __host__ __device__ float length(float2 v);

// normalize
inline __host__ __device__ float2 normalize(float2 v);

// floor
inline __host__ __device__ float2 floor(const float2 v);

// reflect
inline __host__ __device__ float2 reflect(float2 i, float2 n);

inline __device__ __host__ uint UMIN(uint a, uint b);

inline __device__ __host__ uint PowTwoDivider(uint n);

inline __device__ __host__ float2 operator-(float a, float2 b);

inline __device__ __host__ float3 operator-(float a, float3 b);
