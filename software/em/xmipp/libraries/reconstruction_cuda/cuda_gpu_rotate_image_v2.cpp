
//Host includes
#include "cuda_gpu_rotate_image_v2.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
//CUDA includes
#include <cuda_runtime.h>

#define Pole (sqrt(3.0f)-2.0f)  //pole for cubic b-spline

// 2D float texture
texture<float, 2, cudaReadModeElementType> texRef;


//CUDA functions

typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char uchar;
typedef signed char schar;

// float2 functions
////////////////////////////////////////////////////////////////////////////////

// additional constructors
__host__ __device__ float2 make_float2(float s)
{
    return make_float2(s, s);
}
__host__ __device__ float2 make_float2(int2 a)
{
    return make_float2(float(a.x), float(a.y));
}

// addition
__host__ __device__ float2 operator+(float2 a, float2 b)
{
    return make_float2(a.x + b.x, a.y + b.y);
}
 __host__ __device__ void operator+=(float2 &a, float2 b)
{
    a.x += b.x; a.y += b.y;
}

// subtract
 __host__ __device__ float2 operator-(float2 a, float2 b)
{
    return make_float2(a.x - b.x, a.y - b.y);
}
 __host__ __device__ void operator-=(float2 &a, float2 b)
{
    a.x -= b.x; a.y -= b.y;
}

// multiply
 __host__ __device__ float2 operator*(float2 a, float2 b)
{
    return make_float2(a.x * b.x, a.y * b.y);
}
 __host__ __device__ float2 operator*(float2 a, float s)
{
    return make_float2(a.x * s, a.y * s);
}
 __host__ __device__ float2 operator*(float s, float2 a)
{
    return make_float2(a.x * s, a.y * s);
}
 __host__ __device__ void operator*=(float2 &a, float s)
{
    a.x *= s; a.y *= s;
}

// divide
 __host__ __device__ float2 operator/(float2 a, float2 b)
{
    return make_float2(a.x / b.x, a.y / b.y);
}
 __host__ __device__ float2 operator/(float2 a, float s)
{
    float inv = 1.0f / s;
    return a * inv;
}
 __host__ __device__ float2 operator/(float s, float2 a)  //Danny
{
//    float inv = 1.0f / s;
//    return a * inv;
	return make_float2(s / a.x, s / a.y);
}
 __host__ __device__ void operator/=(float2 &a, float s)
{
    float inv = 1.0f / s;
    a *= inv;
}

// dot product
 __host__ __device__ float dot(float2 a, float2 b)
{
    return a.x * b.x + a.y * b.y;
}

// length
 __host__ __device__ float length(float2 v)
{
    return sqrtf(dot(v, v));
}

// normalize
 __host__ __device__ float2 normalize(float2 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    return v * invLen;
}

// floor
 __host__ __device__ float2 floor(const float2 v)
{
    return make_float2(floor(v.x), floor(v.y));
}

// reflect
 __host__ __device__ float2 reflect(float2 i, float2 n)
{
	return i - 2.0f * n * dot(n,i);
}

 __device__ __host__ uint UMIN(uint a, uint b)
{
	return a < b ? a : b;
}

 __device__ __host__ uint PowTwoDivider(uint n)
{
	if (n == 0) return 0;
	uint divider = 1;
	while ((n & divider) == 0) divider <<= 1;
	return divider;
}

 __device__ __host__ float2 operator-(float a, float2 b)
{
	return make_float2(a - b.x, a - b.y);
}

 __device__ __host__ float3 operator-(float a, float3 b)
{
	return make_float3(a - b.x, a - b.y, a - b.z);
}






cudaPitchedPtr CopyVolumeHostToDevice(const float* host, uint width, uint height, uint depth)
{
	cudaPitchedPtr device = {0};
	const cudaExtent extent = make_cudaExtent(width * sizeof(float), height, depth);
	cudaMalloc3D(&device, extent);
	cudaMemcpy3DParms p = {0};
	p.srcPtr = make_cudaPitchedPtr((void*)host, width * sizeof(float), width, height);
	p.dstPtr = device;
	p.extent = extent;
	p.kind = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&p);
	return device;
}



__device__ float InitialCausalCoefficient(
	float* c,			// coefficients
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
		c = (float*)((uchar*)c + step);
	}
	return(Sum);
}


__device__ float InitialAntiCausalCoefficient(
	float* c,			// last coefficient
	uint DataLength,	// number of samples or coefficients
	int step)			// element interleave in bytes
{
	// this initialization corresponds to clamping boundaries
	return((Pole / (Pole - 1.0f)) * *c);
}

__device__ void ConvertToInterpolationCoefficients(
	float* coeffs,		// input samples --> output coefficients
	uint DataLength,	// number of samples or coefficients
	int step)			// element interleave in bytes
{
	// compute the overall gain
	const float Lambda = (1.0f - Pole) * (1.0f - 1.0f / Pole);

	// causal initialization
	float* c = coeffs;
	float previous_c;  //cache the previously calculated c rather than look it up again (faster!)
	*c = previous_c = Lambda * InitialCausalCoefficient(c, DataLength, step);
	// causal recursion
	for (uint n = 1; n < DataLength; n++) {
		c = (float*)((uchar*)c + step);
		*c = previous_c = Lambda * *c + Pole * previous_c;
	}
	// anticausal initialization
	*c = previous_c = InitialAntiCausalCoefficient(c, DataLength, step);
	// anticausal recursion
	for (int n = DataLength - 2; 0 <= n; n--) {
		c = (float*)((uchar*)c - step);
		*c = previous_c = Pole * (previous_c - *c);
	}
}





__global__ void SamplesToCoefficients2DX(
	float* image,		// in-place processing
	uint pitch,			// width in bytes
	uint width,			// width of the image
	uint height)		// height of the image
{
	// process lines in x-direction
	const uint y = blockIdx.x * blockDim.x + threadIdx.x;
	float* line = (float*)((uchar*)image + y * pitch);  //direct access

	ConvertToInterpolationCoefficients(line, width, sizeof(float));
}

__global__ void SamplesToCoefficients2DY(
	float* image,		// in-place processing
	uint pitch,			// width in bytes
	uint width,			// width of the image
	uint height)		// height of the image
{
	// process lines in x-direction
	const uint x = blockIdx.x * blockDim.x + threadIdx.x;
	float* line = image + x;  //direct access

	ConvertToInterpolationCoefficients(line, height, pitch);
}



void CubicBSplinePrefilter2DTimer(float* image, uint pitch, uint width, uint height)
{

	dim3 dimBlockX(min(PowTwoDivider(height), 64));
	dim3 dimGridX(height / dimBlockX.x);
	SamplesToCoefficients2DX<<<dimGridX, dimBlockX>>>(image, pitch, width, height);

	dim3 dimBlockY(min(PowTwoDivider(width), 64));
	dim3 dimGridY(width / dimBlockY.x);
	SamplesToCoefficients2DY<<<dimGridY, dimBlockY>>>(image, pitch, width, height);

}





__device__ void bspline_weights(float2 fraction, float2& w0, float2& w1, float2& w2, float2& w3)
{
	const float2 one_frac = 1.0f - fraction;
	const float2 squared = fraction * fraction;
	const float2 one_sqd = one_frac * one_frac;

	w0 = 1.0f/6.0f * one_sqd * one_frac;
	w1 = 2.0f/3.0f - 0.5f * squared * (2.0f-fraction);
	w2 = 2.0f/3.0f - 0.5f * one_sqd * (2.0f-one_frac);
	w3 = 1.0f/6.0f * squared * fraction;
}

__device__ float cubicTex2D(texture<float, 2, cudaReadModeElementType> tex, float x, float y)
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
	float tex00 = tex2D(tex, h0.x, h0.y);
	float tex10 = tex2D(tex, h1.x, h0.y);
	float tex01 = tex2D(tex, h0.x, h1.y);
	float tex11 = tex2D(tex, h1.x, h1.y);

	// weigh along the y-direction
	tex00 = g0.y * tex00 + g1.y * tex01;
	tex10 = g0.y * tex10 + g1.y * tex11;

	// weigh along the x-direction
	return (g0.x * tex00 + g1.x * tex10);
}


__global__ void
interpolate_kernel(float* output, uint width, float2 extent, float2 a, float2 shift)
{
	uint x = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	uint y = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;
	uint i = __umul24(y, width) + x;

	float x0 = (float)x;
	float y0 = (float)y;
	float x1 = a.x * x0 - a.y * y0 + shift.x;
	float y1 = a.x * y0 + a.y * x0 + shift.y;

	output[i] = cubicTex2D(texRef, x1, y1);

}


cudaPitchedPtr interpolate(uint width, uint height, float angle)
{
	// Prepare the geometry
	float2 a = make_float2((float)cos(angle), (float)sin(angle));
	double xOrigin = floor(width/2);
	double yOrigin = floor(height/2);
	double x0 = a.x * (xOrigin) - a.y * (yOrigin);
	double y0 = a.y * (xOrigin) + a.x * (yOrigin);
	double xShift = xOrigin - x0;
	double yShift = yOrigin - y0;

	// Allocate the output image
	float* output;
	cudaMalloc((void**)&output, width * height * sizeof(float));

	// Visit all pixels of the output image and assign their value
	dim3 blockSize(min(PowTwoDivider(width), 16), min(PowTwoDivider(height), 16));
	dim3 gridSize(width / blockSize.x, height / blockSize.y);
	float2 shift = make_float2((float)xShift, (float)yShift);
	float2 extent = make_float2((float)width, (float)height);
	interpolate_kernel<<<gridSize, blockSize>>>(output, width, extent, a, shift);

	return make_cudaPitchedPtr(output, width * sizeof(float), width, height);
}



void CopyVolumeDeviceToHost(float* host, const cudaPitchedPtr device, uint width, uint height, uint depth)
{
	const cudaExtent extent = make_cudaExtent(width * sizeof(float), height, depth);
	cudaMemcpy3DParms p = {0};
	p.srcPtr = device;
	p.dstPtr = make_cudaPitchedPtr((void*)host, width * sizeof(float), width, height);
	p.extent = extent;
	p.kind = cudaMemcpyDeviceToHost;
	cudaMemcpy3D(&p);
	cudaFree(device.ptr);  //free the GPU volume
}



void cuda_rotate_image_v2(float *image, float *rotated_image, size_t Xdim, size_t Ydim, float ang){

	std::cerr  << "Inside CUDA function " << ang << std::endl;

	//CUDA code
	size_t matSize=Xdim*Ydim*sizeof(float);

	//Filtering process (first step)
	struct cudaPitchedPtr bsplineCoeffs, cudaOutput;
	bsplineCoeffs = CopyVolumeHostToDevice(image, (uint)Xdim, (uint)Ydim, 1);
	CubicBSplinePrefilter2DTimer((float*)bsplineCoeffs.ptr, (uint)bsplineCoeffs.pitch, (uint)Xdim, (uint)Ydim);

	// Init texture
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaArray* cuArray;
    cudaMallocArray(&cuArray, &channelDesc, Xdim, Ydim);
    // Copy to device memory some data located at address h_data in host memory
    cudaMemcpy2DToArray(cuArray, 0, 0, bsplineCoeffs.ptr, bsplineCoeffs.pitch, Xdim * sizeof(float), Ydim, cudaMemcpyDeviceToDevice);

    // Bind the array to the texture reference
    cudaBindTextureToArray(texRef, cuArray, channelDesc);

    // Specify texture object parameters
    texRef.filterMode = cudaFilterModeLinear;
    texRef.normalized = false;

    //Interpolation (second step)
    cudaOutput = interpolate(Xdim, Ydim, ang);

    CopyVolumeDeviceToHost(rotated_image, cudaOutput, Xdim, Ydim, 1);

    cudaFree(cuArray);



}
