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

#ifndef CUDA_BASICMAT_H
#define CUDA_BASICMAT_H

typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char uchar;
typedef signed char schar;

#define Pole (sqrt(3.0f)-2.0f)  //pole for cubic b-spline

// float2 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ float2 make_float2(float s)
{
    return make_float2(s, s);
}
inline __host__ __device__ float2 make_float2(int2 a)
{
    return make_float2(float(a.x), float(a.y));
}

// addition
inline __host__ __device__ float2 operator+(float2 a, float2 b)
{
    return make_float2(a.x + b.x, a.y + b.y);
}
inline __host__ __device__ void operator+=(float2 &a, float2 b)
{
    a.x += b.x; a.y += b.y;
}

// subtract
inline __host__ __device__ float2 operator-(float2 a, float2 b)
{
    return make_float2(a.x - b.x, a.y - b.y);
}
inline __host__ __device__ void operator-=(float2 &a, float2 b)
{
    a.x -= b.x; a.y -= b.y;
}

// multiply
inline __host__ __device__ float2 operator*(float2 a, float2 b)
{
    return make_float2(a.x * b.x, a.y * b.y);
}
inline __host__ __device__ float2 operator*(float2 a, float s)
{
    return make_float2(a.x * s, a.y * s);
}
inline __host__ __device__ float2 operator*(float s, float2 a)
{
    return make_float2(a.x * s, a.y * s);
}
inline __host__ __device__ void operator*=(float2 &a, float s)
{
    a.x *= s; a.y *= s;
}

// divide
inline __host__ __device__ float2 operator/(float2 a, float2 b)
{
    return make_float2(a.x / b.x, a.y / b.y);
}
inline __host__ __device__ float2 operator/(float2 a, float s)
{
    float inv = 1.0f / s;
    return a * inv;
}
inline __host__ __device__ float2 operator/(float s, float2 a)  //Danny
{
	return make_float2(s / a.x, s / a.y);
}
inline __host__ __device__ void operator/=(float2 &a, float s)
{
    float inv = 1.0f / s;
    a *= inv;
}

// dot product
inline __host__ __device__ float dot(float2 a, float2 b)
{
    return a.x * b.x + a.y * b.y;
}

// length
inline __host__ __device__ float length(float2 v)
{
    return sqrtf(dot(v, v));
}

// normalize
inline __host__ __device__ float2 normalize(float2 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}

// floor
inline __host__ __device__ float2 floor(const float2 v)
{
    return make_float2(floor(v.x), floor(v.y));
}

// reflect
inline __host__ __device__ float2 reflect(float2 i, float2 n)
{
	return i - 2.0f * n * dot(n,i);
}

inline __device__ __host__ uint UMIN(uint a, uint b)
{
	return a < b ? a : b;
}

inline __device__ __host__ uint PowTwoDivider(uint n)
{
	if (n == 0) return 0;
	uint divider = 1;
	while ((n & divider) == 0) divider <<= 1;
	return divider;
}

inline __device__ __host__ float2 operator-(float a, float2 b)
{
	return make_float2(a - b.x, a - b.y);
}

inline __device__ __host__ float3 operator-(float a, float3 b)
{
	return make_float3(a - b.x, a - b.y, a - b.z);
}


// float3 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
inline __host__ __device__ float3 make_float3(float s)
{
    return make_float3(s, s, s);
}
inline __host__ __device__ float3 make_float3(float2 a)
{
    return make_float3(a.x, a.y, 0.0f);
}
inline __host__ __device__ float3 make_float3(float2 a, float s)
{
    return make_float3(a.x, a.y, s);
}
inline __host__ __device__ float3 make_float3(float4 a)
{
    return make_float3(a.x, a.y, a.z);  // discards w
}
inline __host__ __device__ float3 make_float3(int3 a)
{
    return make_float3(float(a.x), float(a.y), float(a.z));
}

// min
static __inline__ __host__ __device__ float3 fminf(float3 a, float3 b)
{
	return make_float3(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z));
}

// max
static __inline__ __host__ __device__ float3 fmaxf(float3 a, float3 b)
{
	return make_float3(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z));
}

// addition
inline __host__ __device__ float3 operator+(float3 a, float3 b)
{
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __host__ __device__ float3 operator+(float3 a, float b)
{
    return make_float3(a.x + b, a.y + b, a.z + b);
}
inline __host__ __device__ void operator+=(float3 &a, float3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}

// subtract
inline __host__ __device__ float3 operator-(float3 a, float3 b)
{
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __host__ __device__ float3 operator-(float3 a, float b)
{
    return make_float3(a.x - b, a.y - b, a.z - b);
}
inline __host__ __device__ void operator-=(float3 &a, float3 b)
{
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
}

// multiply
inline __host__ __device__ float3 operator*(float3 a, float3 b)
{
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __host__ __device__ float3 operator*(float3 a, float s)
{
    return make_float3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ float3 operator*(float s, float3 a)
{
    return make_float3(a.x * s, a.y * s, a.z * s);
}
inline __host__ __device__ void operator*=(float3 &a, float s)
{
    a.x *= s; a.y *= s; a.z *= s;
}

// divide
inline __host__ __device__ float3 operator/(float3 a, float3 b)
{
    return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __host__ __device__ float3 operator/(float3 a, float s)
{
    float inv = 1.0f / s;
    return a * inv;
}
inline __host__ __device__ float3 operator/(float s, float3 a)  //Danny
{
//    float inv = 1.0f / s;
//    return a * inv;
	return make_float3(s / a.x, s / a.y, s / a.z);
}
inline __host__ __device__ void operator/=(float3 &a, float s)
{
    float inv = 1.0f / s;
    a *= inv;
}

// dot product
inline __host__ __device__ float dot(float3 a, float3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// cross product
inline __host__ __device__ float3 cross(float3 a, float3 b)
{
    return make_float3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

// length
inline __host__ __device__ float length(float3 v)
{
    return sqrtf(dot(v, v));
}

// normalize
inline __host__ __device__ float3 normalize(float3 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}

// floor
inline __host__ __device__ float3 floor(const float3 v)
{
    return make_float3(floor(v.x), floor(v.y), floor(v.z));
}

// reflect
inline __host__ __device__ float3 reflect(float3 i, float3 n)
{
	return i - 2.0f * n * dot(n,i);
}

#endif
