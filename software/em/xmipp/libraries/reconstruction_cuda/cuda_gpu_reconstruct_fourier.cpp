/***************************************************************************
 *
 * Authors:     David Strelak (davidstrelak@gmail.com)
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

#include <cuda_runtime_api.h>
#include "reconstruction_cuda/cuda_utils.h" // cannot be in header as it includes cuda headers
#include "cuda_gpu_reconstruct_fourier.h"

#if SHARED_BLOB_TABLE
__shared__ float BLOB_TABLE[BLOB_TABLE_SIZE_SQRT];
#endif

#if SHARED_IMG
__shared__ Point3D<float> SHARED_AABB[2];
extern __shared__ float2 IMG[];
#endif

// FIELDS

// Holding streams used for calculation. Present on CPU
cudaStream_t* streams;

// Wrapper to hold pointers to GPU memory (and have it also accessible from CPU)
std::map<int,FRecBufferDataGPUWrapper*> wrappers;

// Holding blob coefficient table. Present on GPU
float* devBlobTableSqrt = NULL;

__device__ __constant__ int cMaxVolumeIndexX = 0;
__device__ __constant__ int cMaxVolumeIndexYZ = 0;
__device__ __constant__ float cBlobRadius = 0.f;
__device__ __constant__ float cBlobAlpha = 0.f;
__device__ __constant__ float cIw0 = 0.f;
__device__ __constant__ float cIDeltaSqrt = 0.f;


__device__
float bessi0(float x)
{
    float y, ax, ans;
    if ((ax = fabsf(x)) < 3.75f)
    {
        y = x / 3.75f;
        y *= y;
        ans = 1.f + y * (3.5156229f + y * (3.0899424f + y * (1.2067492f
                                          + y * (0.2659732f + y * (0.360768e-1f + y * 0.45813e-2f)))));
    }
    else
    {
        y = 3.75f / ax;
        ans = (expf(ax) / sqrtf(ax)) * (0.39894228f + y * (0.1328592e-1f
                                      + y * (0.225319e-2f + y * (-0.157565e-2f + y * (0.916281e-2f
                                                                + y * (-0.2057706e-1f + y * (0.2635537e-1f + y * (-0.1647633e-1f
                                                                                            + y * 0.392377e-2f))))))));
    }
    return ans;
}

__device__
float bessi1(float x)
{
    float ax, ans;
    float y;
    if ((ax = fabsf(x)) < 3.75f)
    {
        y = x / 3.75f;
        y *= y;
        ans = ax * (0.5f + y * (0.87890594f + y * (0.51498869f + y * (0.15084934f
                               + y * (0.2658733e-1f + y * (0.301532e-2f + y * 0.32411e-3f))))));
    }
    else
    {
        y = 3.75f / ax;
        ans = 0.2282967e-1f + y * (-0.2895312e-1f + y * (0.1787654e-1f
                                  - y * 0.420059e-2f));
        ans = 0.39894228f + y * (-0.3988024e-1f + y * (-0.362018e-2f
                                + y * (0.163801e-2f + y * (-0.1031555e-1f + y * ans))));
        ans *= (expf(ax) / sqrtf(ax));
    }
    return x < 0.0 ? -ans : ans;
}

__device__
float bessi2(float x)
{
    return (x == 0) ? 0 : bessi0(x) - ((2*1) / x) * bessi1(x);
}

__device__
float bessi3(float x)
{
    return (x == 0) ? 0 : bessi1(x) - ((2*2) / x) * bessi2(x);
}

__device__
float bessi4(float x)
{
    return (x == 0) ? 0 : bessi2(x) - ((2*3) / x) * bessi3(x);
}


template<int order>
__device__
float kaiserValue(float r, float a, float alpha)
{
    float rda, rdas, arg, w;

    rda = r / a;
    if (rda <= 1.f)
    {
        rdas = rda * rda;
        arg = alpha * sqrtf(1.f - rdas);
        if (order == 0)
        {
            w = bessi0(arg) / bessi0(alpha);
        }
        else if (order == 1)
        {
            w = sqrtf (1.f - rdas);
            if (alpha != 0.f)
                w *= bessi1(arg) / bessi1(alpha);
        }
        else if (order == 2)
        {
            w = sqrtf (1.f - rdas);
            w = w * w;
            if (alpha != 0.f)
                w *= bessi2(arg) / bessi2(alpha);
        }
        else if (order == 3)
        {
            w = sqrtf (1.f - rdas);
            w = w * w * w;
            if (alpha != 0.f)
                w *= bessi3(arg) / bessi3(alpha);
        }
        else if (order == 4)
        {
            w = sqrtf (1.f - rdas);
            w = w * w * w *w;
            if (alpha != 0.f)
                w *= bessi4(arg) / bessi4(alpha);
        }
        else {
        	printf("order (%d) out of range in kaiser_value(): %s, %d\n", order, __FILE__, __LINE__);
        }
    }
    else
        w = 0.f;

    return w;
}

/**
 * Structure for buffer data on GPU
 * Adds some utility methods for copying data to GPU and device specific code.
 * Variables are stored in the same memory space as holding object, but
 * pointers point to the GPU memory
 */
struct RecFourierBufferDataGPU : public RecFourierBufferData {

	RecFourierBufferDataGPU(RecFourierBufferData* orig) {
		copyMetadata(orig);
		FFTs = CTFs = paddedImages = modulators = NULL;

		// allocate space at GPU
		alloc(orig->FFTs, FFTs, orig); // FFT are always necessary
		alloc(orig->spaces, spaces, orig);
		if ( ! hasFFTs) {
			alloc(orig->paddedImages, paddedImages, orig);
		}
		if (hasCTFs) {
			alloc(orig->CTFs, CTFs, orig);
			alloc(orig->modulators, modulators, orig);
		}
	}

	~RecFourierBufferDataGPU() {
		cudaFree(FFTs);
		cudaFree(CTFs);
		cudaFree(paddedImages);
		cudaFree(modulators);
		cudaFree(spaces);
		gpuErrchk( cudaPeekAtLastError() );

		invalidate();
	}

	/**
	 * Same as in parent struct, but implemented for device
	 */
	__device__
	float* getNthItem(float* array, int itemIndex) {
		if (array == FFTs) return array + (fftSizeX * fftSizeY * itemIndex * 2); // *2 since it's complex
		if (array == CTFs) return array + (fftSizeX * fftSizeY * itemIndex);
		if (array == modulators) return array + (fftSizeX * fftSizeY * itemIndex);
		if (array == paddedImages) return array + (paddedImgSize * paddedImgSize * itemIndex);
		return NULL; // undefined
	}

	/**
	 * Method copies data from buffer. Asynchronous method
	 */
	void copyDataFrom(RecFourierBufferData* orig, int stream) {
		copyMetadata(orig);
		copy(orig->FFTs, FFTs, orig, stream);
		copy(orig->CTFs, CTFs, orig, stream);
		copy(orig->paddedImages, paddedImages, orig, stream);
		copy(orig->modulators, modulators, orig, stream);
		copy(orig->spaces, spaces, orig, stream);
	}

	/**
	 * Same as in parent struct, but implemented for device
	 */
	__device__
	int getNoOfSpaces() {
		return noOfImages * noOfSymmetries;
	}

private:
	/**
	 * Method copies content of srcArray (host) to 'dstArray' (device, must be allocated a priory)
	 * Arrays are supposed to be from the 'orig' buffer, or this object.
	 * Asynchronous method.
	 */
	template<typename T>
	void copy(T* srcArray, T*& dstArray, RecFourierBufferData* orig, int stream) {
		if (NULL != srcArray) {
			size_t bytes = sizeof(T) * orig->getNoOfElements(srcArray);
			cudaHostRegister(srcArray, bytes, 0);
			cudaMemcpyAsync(dstArray, srcArray, bytes, cudaMemcpyHostToDevice, streams[stream]);
			cudaHostUnregister(srcArray);
			gpuErrchk( cudaPeekAtLastError() );
		}
	}

	/**
	 * Method allocates 'dstArray' (device) to the size of 'srcArray' of the 'orig' buffer
	 * Blocking method
	 */
	template<typename T>
	void alloc(T* srcArray, T*& dstArray, RecFourierBufferData* orig) {
		size_t bytes = orig->getMaxByteSize(srcArray);
		cudaMalloc((void **) &dstArray, bytes);
		gpuErrchk( cudaPeekAtLastError() );
	}

	/**
	 * Method copies 'metadata' (fields, not arrays) from the 'orig' buffer
	 * Blocking method
	 */
	void copyMetadata(RecFourierBufferData* orig) {
		hasCTFs = orig->hasCTFs;
		hasFFTs = orig->hasFFTs;
		noOfImages = orig->noOfImages;
		paddedImgSize = orig->paddedImgSize;
		fftSizeX = orig->fftSizeX;
		fftSizeY = orig->fftSizeY;
		maxNoOfImages = orig->maxNoOfImages;
		noOfSymmetries = orig->noOfSymmetries;
	}
};


FRecBufferDataGPUWrapper::FRecBufferDataGPUWrapper(RecFourierBufferData* orig) {
	cpuCopy = new RecFourierBufferDataGPU(orig);
	gpuCopy = NULL;
}

FRecBufferDataGPUWrapper::~FRecBufferDataGPUWrapper() {
	cudaFree(gpuCopy);
	gpuErrchk( cudaPeekAtLastError() );
	delete cpuCopy;
}

void FRecBufferDataGPUWrapper::copyFrom(RecFourierBufferData* orig, int stream) {
	cpuCopy->copyDataFrom(orig, stream);
}

void FRecBufferDataGPUWrapper::copyToDevice(int stream) {
	if (NULL == gpuCopy) {
		cudaMalloc((void **) &gpuCopy, sizeof(RecFourierBufferDataGPU));
		gpuErrchk( cudaPeekAtLastError() );
	}
	cudaMemcpyAsync(gpuCopy, cpuCopy, sizeof(RecFourierBufferDataGPU), cudaMemcpyHostToDevice, streams[stream]);
	gpuErrchk( cudaPeekAtLastError() );
}

float* allocateTempVolumeGPU(float*& ptr, int size, int typeSize) {
	cudaMalloc((void**)&ptr, size * size * size * typeSize);
	cudaMemset(ptr, 0.f, size * size * size * typeSize);
	gpuErrchk( cudaPeekAtLastError() );

	return ptr;
}

void copyTempVolumes(std::complex<float>*** tempVol, float*** tempWeights,
		float* tempVolGPU, float* tempWeightsGPU,
		int size) {
	for (int z = 0; z < size; z++) {
		for (int y = 0; y < size; y++) {
			int index = (z * size * size) + (y * size);
			cudaMemcpy(tempVol[z][y], &tempVolGPU[2 * index], 2 * size * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(tempWeights[z][y] , &tempWeightsGPU[index], size * sizeof(float), cudaMemcpyDeviceToHost);
		}
	}
	gpuErrchk(cudaPeekAtLastError());
}

void releaseTempVolumeGPU(float*& ptr) {
	cudaFree(ptr);
	ptr = NULL;
	gpuErrchk(cudaPeekAtLastError());
}


/** Index to frequency
 *
 * Given an index and a size of the FFT, this function returns the corresponding
 * digital frequency (-1/2 to 1/2)
 */ // FIXME unify with xmipp_fft.h::FFT_IDX2DIGFREQ
__device__
float FFT_IDX2DIGFREQ(int idx, int size) {
	if (size <= 1) return 0;
	return ((idx <= (size / 2)) ? idx : (-size + idx)) / (float)size;
}

/** Returns true if x is in (min, max), i.e. opened, interval */
template <typename T>
__device__
static bool inRange(T x, T min, T max) {
	return (x > min) && (x < max);
}

/** Returns value within the range (included) */
template<typename T, typename U>
__device__
static U clamp(U val, T min, T max) {
	U res = (val > max) ? max : val;
	return (res < min) ? min : res;
}

/**
 * Calculates Z coordinate of the point [x, y] on the plane defined by p0 (origin) and two vectors
 * Returns 'true' if the in point lies within parallelogram
 */
__device__
bool getZ(float x, float y, float& z, const Point3D<float>& a, const Point3D<float>& b, const Point3D<float>& p0) {
	// from parametric eq. of the plane
	float u = ((y-p0.y)*a.x + (p0.x-x)*a.y) / (a.x * b.y - b.x * a.y);
	float t = (-p0.x + x - u*b.x) / (a.x);

	z = p0.z + t*a.z + u*b.z;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

/**
 * Calculates Y coordinate of the point [x, z] on the plane defined by p0 (origin) and two vectors
 * Returns 'true' if the in point lies within parallelogram
 */
__device__
bool getY(float x, float& y, float z, const Point3D<float>& a, const Point3D<float>& b, const Point3D<float>& p0) {
	// from parametric eq. of the plane
	float u = ((z-p0.z)*a.x + (p0.x-x)*a.z) / (a.x * b.z - b.x * a.z);
	float t = (-p0.x + x - u*b.x) / (a.x);

	y = p0.y + t*a.y + u*b.y;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

/**
 * Calculates X coordinate of the point [y, z] on the plane defined by p0 (origin) and two vectors
 * Returns 'true' if the in point lies within parallelogram
 */
__device__
bool getX(float& x, float y, float z, const Point3D<float>& a, const Point3D<float>& b, const Point3D<float>& p0) {
	// from parametric eq. of the plane
	float u = ((z-p0.z)*a.y + (p0.y-y)*a.z) / (a.y * b.z - b.y * a.z);
	float t = (-p0.y + y - u*b.y) / (a.y);

	x = p0.x + t*a.x + u*b.x;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

/** Do 3x3 x 1x3 matrix-vector multiplication */
__device__
void multiply(const float transform[3][3], Point3D<float>& inOut) {
	float tmp0 = transform[0][0] * inOut.x + transform[0][1] * inOut.y + transform[0][2] * inOut.z;
	float tmp1 = transform[1][0] * inOut.x + transform[1][1] * inOut.y + transform[1][2] * inOut.z;
	float tmp2 = transform[2][0] * inOut.x + transform[2][1] * inOut.y + transform[2][2] * inOut.z;
	inOut.x = tmp0;
	inOut.y = tmp1;
	inOut.z = tmp2;
}

/** Compute Axis Aligned Bounding Box of given cuboid */
__device__
void computeAABB(Point3D<float>* AABB, Point3D<float>* cuboid) {
	AABB[0].x = AABB[0].y = AABB[0].z = INFINITY;
	AABB[1].x = AABB[1].y = AABB[1].z = -INFINITY;
	Point3D<float> tmp;
	for (int i = 0; i < 8; i++) {
		tmp = cuboid[i];
		if (AABB[0].x > tmp.x) AABB[0].x = tmp.x;
		if (AABB[0].y > tmp.y) AABB[0].y = tmp.y;
		if (AABB[0].z > tmp.z) AABB[0].z = tmp.z;
		if (AABB[1].x < tmp.x) AABB[1].x = tmp.x;
		if (AABB[1].y < tmp.y) AABB[1].y = tmp.y;
		if (AABB[1].z < tmp.z) AABB[1].z = tmp.z;
	}
	AABB[0].x = ceilf(AABB[0].x);
	AABB[0].y = ceilf(AABB[0].y);
	AABB[0].z = ceilf(AABB[0].z);

	AABB[1].x = floorf(AABB[1].x);
	AABB[1].y = floorf(AABB[1].y);
	AABB[1].z = floorf(AABB[1].z);
}

/**
 * Method will map one voxel from the temporal
 * spaces to the given projection and update temporal spaces
 * using the pixel value of the projection.
 */
template<bool hasCTF>
__device__
void processVoxel(
	float* tempVolumeGPU, float* tempWeightsGPU,
	int x, int y, int z,
	RecFourierBufferDataGPU* const data,
	const RecFourierProjectionTraverseSpace* const space)
{
	Point3D<float> imgPos;
	float wBlob = 1.f;
	float wCTF = 1.f;
	float wModulator = 1.f;
	const float* __restrict__ img = data->getNthItem(data->FFTs, space->projectionIndex);


	float dataWeight = space->weight;
	int xSize = data->fftSizeX;
	int ySize = data->fftSizeY;

	// transform current point to center
	imgPos.x = x - cMaxVolumeIndexX/2;
	imgPos.y = y - cMaxVolumeIndexYZ/2;
	imgPos.z = z - cMaxVolumeIndexYZ/2;
	if (imgPos.x*imgPos.x + imgPos.y*imgPos.y + imgPos.z*imgPos.z > space->maxDistanceSqr) {
		return; // discard iterations that would access pixel with too high frequency
	}
	// rotate around center
	multiply(space->transformInv, imgPos);
	// transform back and round
	// just Y coordinate needs adjusting, since X now matches to picture and Z is irrelevant
	int imgX = clamp((int)(imgPos.x + 0.5f), 0, xSize - 1);
	int imgY = clamp((int)(imgPos.y + 0.5f + cMaxVolumeIndexYZ / 2), 0, ySize - 1);

	int index3D = z * (cMaxVolumeIndexYZ+1) * (cMaxVolumeIndexX+1) + y * (cMaxVolumeIndexX+1) + x;
	int index2D = imgY * xSize + imgX;

	if (hasCTF) {
		const float* __restrict__ CTF = data->getNthItem(data->CTFs, space->projectionIndex);
		const float* __restrict__ modulator = data->getNthItem(data->modulators, space->projectionIndex);
		wCTF = CTF[index2D];
		wModulator = modulator[index2D];
	}

	float weight = wBlob * wModulator * dataWeight;

	 // use atomic as two blocks can write to same voxel
	atomicAdd(&tempVolumeGPU[2*index3D], img[2*index2D] * weight * wCTF);
	atomicAdd(&tempVolumeGPU[2*index3D + 1], img[2*index2D + 1] * weight * wCTF);
	atomicAdd(&tempWeightsGPU[index3D], weight);
}

/**
 * Method will map one voxel from the temporal
 * spaces to the given projection and update temporal spaces
 * using the pixel values of the projection withing the blob distance.
 */
template<bool hasCTF, int blobOrder>
__device__
void processVoxelBlob(
	float* tempVolumeGPU, float *tempWeightsGPU,
	int x, int y, int z,
	RecFourierBufferDataGPU* const data,
	const RecFourierProjectionTraverseSpace* const space,
	const float* blobTableSqrt,
	int imgCacheDim)
{
	Point3D<float> imgPos;
	volatile int xSize = data->fftSizeX;
	volatile int ySize = data->fftSizeY;
	// transform current point to center
	imgPos.x = x - cMaxVolumeIndexX/2;
	imgPos.y = y - cMaxVolumeIndexYZ/2;
	imgPos.z = z - cMaxVolumeIndexYZ/2;
	if ((imgPos.x*imgPos.x + imgPos.y*imgPos.y + imgPos.z*imgPos.z) > space->maxDistanceSqr) {
		return; // discard iterations that would access pixel with too high frequency
	}
	// rotate around center
	multiply(space->transformInv, imgPos);
	// transform back just Y coordinate, since X now matches to picture and Z is irrelevant
	imgPos.y += cMaxVolumeIndexYZ / 2;

	// check that we don't want to collect data from far far away ...
	float radiusSqr = cBlobRadius * cBlobRadius;
	float zSqr = imgPos.z * imgPos.z;
	if (zSqr > radiusSqr) return;

	// create blob bounding box
	int minX = ceilf(imgPos.x - cBlobRadius);
	int maxX = floorf(imgPos.x + cBlobRadius);
	int minY = ceilf(imgPos.y - cBlobRadius);
	int maxY = floorf(imgPos.y + cBlobRadius);
	minX = fmaxf(minX, 0);
	minY = fmaxf(minY, 0);
	maxX = fminf(maxX, xSize-1);
	maxY = fminf(maxY, ySize-1);

	int index3D = z * (cMaxVolumeIndexYZ+1) * (cMaxVolumeIndexX+1) + y * (cMaxVolumeIndexX+1) + x;
	float volReal, volImag, w;
	volReal = volImag = w = 0.f;
#if !SHARED_IMG
	const float* __restrict__ img = data->getNthItem(data->FFTs, space->projectionIndex);
#endif
	float dataWeight = space->weight;

	// ugly spaghetti code, but improves performance by app. 10%
	if (hasCTF) {
		const float* __restrict__ CTF = data->getNthItem(data->CTFs, space->projectionIndex);
		const float* __restrict__ modulator = data->getNthItem(data->modulators, space->projectionIndex);

		// check which pixel in the vicinity should contribute
		for (int i = minY; i <= maxY; i++) {
			float ySqr = (imgPos.y - i) * (imgPos.y - i);
			float yzSqr = ySqr + zSqr;
			if (yzSqr > radiusSqr) continue;
			for (int j = minX; j <= maxX; j++) {
				float xD = imgPos.x - j;
				float distanceSqr = xD*xD + yzSqr;
				if (distanceSqr > radiusSqr) continue;

#if SHARED_IMG
				int index2D = (i - SHARED_AABB[0].y) * imgCacheDim + (j-SHARED_AABB[0].x); // position in img - offset of the AABB
#else
				int index2D = i * xSize + j;
#endif

				float wCTF = CTF[index2D];
				float wModulator = modulator[index2D];
#if PRECOMPUTE_BLOB_VAL
				int aux = (int) ((distanceSqr * cIDeltaSqrt + 0.5f));
	#if SHARED_BLOB_TABLE
				float wBlob = BLOB_TABLE[aux];
	#else
				float wBlob = blobTableSqrt[aux];
	#endif
#else
				float wBlob = kaiserValue<blobOrder>(sqrtf(distanceSqr),cBlobRadius, cBlobAlpha);
#endif
				float weight = wBlob * wModulator * dataWeight;
				w += weight;
#if SHARED_IMG
				volReal += IMG[index2D].x * weight * wCTF;
				volImag += IMG[index2D].y * weight * wCTF;
#else
				volReal += img[2*index2D] * weight * wCTF;
				volImag += img[2*index2D + 1] * weight * wCTF;
#endif
			}
		}
	} else {
		// check which pixel in the vicinity should contribute
		for (int i = minY; i <= maxY; i++) {
			float ySqr = (imgPos.y - i) * (imgPos.y - i);
			float yzSqr = ySqr + zSqr;
			if (yzSqr > radiusSqr) continue;
			for (int j = minX; j <= maxX; j++) {
				float xD = imgPos.x - j;
				float distanceSqr = xD*xD + yzSqr;
				if (distanceSqr > radiusSqr) continue;

#if SHARED_IMG
				int index2D = (i - SHARED_AABB[0].y) * imgCacheDim + (j-SHARED_AABB[0].x); // position in img - offset of the AABB
#else
				int index2D = i * xSize + j;
#endif

#if PRECOMPUTE_BLOB_VAL
				int aux = (int) ((distanceSqr * cIDeltaSqrt + 0.5f));
#if SHARED_BLOB_TABLE
				float wBlob = BLOB_TABLE[aux];
#else
				float wBlob = blobTableSqrt[aux];
#endif
#else
				float wBlob = kaiserValue<blobOrder>(sqrtf(distanceSqr),cBlobRadius, cBlobAlpha);
#endif
				float weight = wBlob * dataWeight;
				w += weight;
#if SHARED_IMG
				volReal += IMG[index2D].x * weight;
				volImag += IMG[index2D].y * weight;
#else
				volReal += img[2*index2D] * weight;
				volImag += img[2*index2D + 1] * weight;
#endif
			}
		}
	}
	// use atomic as two blocks can write to same voxel
	atomicAdd(&tempVolumeGPU[2*index3D], volReal);
	atomicAdd(&tempVolumeGPU[2*index3D + 1], volImag);
	atomicAdd(&tempWeightsGPU[index3D], w);
}

/**
  * Method will process one projection image and add result to temporal
  * spaces.
  */
template<bool useFast, bool hasCTF, int blobOrder>
__device__
void processProjection(
	float* tempVolumeGPU, float *tempWeightsGPU,
	RecFourierBufferDataGPU* const data,
	const RecFourierProjectionTraverseSpace* const tSpace,
	const float* devBlobTableSqrt,
	int imgCacheDim)
{
	// map thread to each (2D) voxel
	volatile int idx = blockIdx.x*blockDim.x + threadIdx.x;
	volatile int idy = blockIdx.y*blockDim.y + threadIdx.y;

	if (tSpace->XY == tSpace->dir) { // iterate XY plane
		if (idy >= tSpace->minY && idy <= tSpace->maxY) {
			if (idx >= tSpace->minX && idx <= tSpace->maxX) {
				if (useFast) {
					float hitZ;
					if (getZ(idx, idy, hitZ, tSpace->u, tSpace->v, tSpace->bottomOrigin)) {
						int z = (int)(hitZ + 0.5f); // rounding
						processVoxel<hasCTF>(tempVolumeGPU, tempWeightsGPU, idx, idy, z, data, tSpace);
					}
				} else {
					float z1, z2;
					bool hit1 = getZ(idx, idy, z1, tSpace->u, tSpace->v, tSpace->bottomOrigin); // lower plane
					bool hit2 = getZ(idx, idy, z2, tSpace->u, tSpace->v, tSpace->topOrigin); // upper plane
					if (hit1 || hit2) {
						z1 = clamp(z1, 0, cMaxVolumeIndexYZ);
						z2 = clamp(z2, 0, cMaxVolumeIndexYZ);
						int lower = floorf(fminf(z1, z2));
						int upper = ceilf(fmaxf(z1, z2));
						for (int z = lower; z <= upper; z++) {
							processVoxelBlob<hasCTF, blobOrder>(tempVolumeGPU, tempWeightsGPU, idx, idy, z, data, tSpace, devBlobTableSqrt, imgCacheDim);
						}
					}
				}
			}
		}
	} else if (tSpace->XZ == tSpace->dir) { // iterate XZ plane
		if (idy >= tSpace->minZ && idy <= tSpace->maxZ) { // map z -> y
			if (idx >= tSpace->minX && idx <= tSpace->maxX) {
				if (useFast) {
					float hitY;
					if (getY(idx, hitY, idy, tSpace->u, tSpace->v, tSpace->bottomOrigin)) {
						int y = (int)(hitY + 0.5f); // rounding
						processVoxel<hasCTF>(tempVolumeGPU, tempWeightsGPU, idx, y, idy, data, tSpace);
					}
				} else {
					float y1, y2;
					bool hit1 = getY(idx, y1, idy, tSpace->u, tSpace->v, tSpace->bottomOrigin); // lower plane
					bool hit2 = getY(idx, y2, idy, tSpace->u, tSpace->v, tSpace->topOrigin); // upper plane
					if (hit1 || hit2) {
						y1 = clamp(y1, 0, cMaxVolumeIndexYZ);
						y2 = clamp(y2, 0, cMaxVolumeIndexYZ);
						int lower = floorf(fminf(y1, y2));
						int upper = ceilf(fmaxf(y1, y2));
						for (int y = lower; y <= upper; y++) {
							processVoxelBlob<hasCTF, blobOrder>(tempVolumeGPU, tempWeightsGPU, idx, y, idy, data, tSpace, devBlobTableSqrt, imgCacheDim);
						}
					}
				}
			}
		}
	} else { // iterate YZ plane
		if (idy >= tSpace->minZ && idy <= tSpace->maxZ) { // map z -> y
			if (idx >= tSpace->minY && idx <= tSpace->maxY) { // map y > x
				if (useFast) {
					float hitX;
					if (getX(hitX, idx, idy, tSpace->u, tSpace->v, tSpace->bottomOrigin)) {
						int x = (int)(hitX + 0.5f); // rounding
						processVoxel<hasCTF>(tempVolumeGPU, tempWeightsGPU, x, idx, idy, data, tSpace);
					}
				} else {
					float x1, x2;
					bool hit1 = getX(x1, idx, idy, tSpace->u, tSpace->v, tSpace->bottomOrigin); // lower plane
					bool hit2 = getX(x2, idx, idy, tSpace->u, tSpace->v, tSpace->topOrigin); // upper plane
					if (hit1 || hit2) {
						x1 = clamp(x1, 0, cMaxVolumeIndexX);
						x2 = clamp(x2, 0, cMaxVolumeIndexX);
						int lower = floorf(fminf(x1, x2));
						int upper = ceilf(fmaxf(x1, x2));
						for (int x = lower; x <= upper; x++) {
							processVoxelBlob<hasCTF, blobOrder>(tempVolumeGPU, tempWeightsGPU, x, idx, idy, data, tSpace, devBlobTableSqrt, imgCacheDim);
						}
					}
				}
			}
		}
	}
}

/**
 * Method will rotate box using transformation matrix around center of the
 * working space
 */
__device__
void rotate(Point3D<float>* box, const float transform[3][3]) {
	for (int i = 0; i < 8; i++) {
		Point3D<float> imgPos;
		// transform current point to center
		imgPos.x = box[i].x - cMaxVolumeIndexX/2;
		imgPos.y = box[i].y - cMaxVolumeIndexYZ/2;
		imgPos.z = box[i].z - cMaxVolumeIndexYZ/2;
		// rotate around center
		multiply(transform, imgPos);
		// transform back just Y coordinate, since X now matches to picture and Z is irrelevant
		imgPos.y += cMaxVolumeIndexYZ / 2;

		box[i] = imgPos;
	}
}

/**
 * Method calculates an Axis Aligned Bounding Box in the image space.
 * AABB is guaranteed to be big enough that all threads in the block,
 * while processing the traverse space, will not read image data outside
 * of the AABB
 */
__device__
void calculateAABB(const RecFourierProjectionTraverseSpace* tSpace, const RecFourierBufferDataGPU* buffer, Point3D<float>* dest) {
	Point3D<float> box[8];
	// calculate AABB for the whole working block
	if (tSpace->XY == tSpace->dir) { // iterate XY plane
		box[0].x = box[3].x = box[4].x = box[7].x = blockIdx.x*blockDim.x - cBlobRadius;
		box[1].x = box[2].x = box[5].x = box[6].x = (blockIdx.x+1)*blockDim.x + cBlobRadius - 1.f;

		box[2].y = box[3].y = box[6].y = box[7].y = (blockIdx.y+1)*blockDim.y + cBlobRadius - 1.f;
		box[0].y = box[1].y = box[4].y = box[5].y = blockIdx.y*blockDim.y- cBlobRadius;

		getZ(box[0].x, box[0].y, box[0].z, tSpace->u, tSpace->v, tSpace->bottomOrigin);
		getZ(box[4].x, box[4].y, box[4].z, tSpace->u, tSpace->v, tSpace->topOrigin);

		getZ(box[3].x, box[3].y, box[3].z, tSpace->u, tSpace->v, tSpace->bottomOrigin);
		getZ(box[7].x, box[7].y, box[7].z, tSpace->u, tSpace->v, tSpace->topOrigin);

		getZ(box[2].x, box[2].y, box[2].z, tSpace->u, tSpace->v, tSpace->bottomOrigin);
		getZ(box[6].x, box[6].y, box[6].z, tSpace->u, tSpace->v, tSpace->topOrigin);

		getZ(box[1].x, box[1].y, box[1].z, tSpace->u, tSpace->v, tSpace->bottomOrigin);
		getZ(box[5].x, box[5].y, box[5].z, tSpace->u, tSpace->v, tSpace->topOrigin);
	} else if (tSpace->XZ == tSpace->dir) { // iterate XZ plane
		box[0].x = box[3].x = box[4].x = box[7].x = blockIdx.x*blockDim.x - cBlobRadius;
		box[1].x = box[2].x = box[5].x = box[6].x = (blockIdx.x+1)*blockDim.x + cBlobRadius - 1.f;

		box[2].z = box[3].z = box[6].z = box[7].z = (blockIdx.y+1)*blockDim.y + cBlobRadius - 1.f;
		box[0].z = box[1].z = box[4].z = box[5].z = blockIdx.y*blockDim.y- cBlobRadius;

		getY(box[0].x, box[0].y, box[0].z, tSpace->u, tSpace->v, tSpace->bottomOrigin);
		getY(box[4].x, box[4].y, box[4].z, tSpace->u, tSpace->v, tSpace->topOrigin);

		getY(box[3].x, box[3].y, box[3].z, tSpace->u, tSpace->v, tSpace->bottomOrigin);
		getY(box[7].x, box[7].y, box[7].z, tSpace->u, tSpace->v, tSpace->topOrigin);

		getY(box[2].x, box[2].y, box[2].z, tSpace->u, tSpace->v, tSpace->bottomOrigin);
		getY(box[6].x, box[6].y, box[6].z, tSpace->u, tSpace->v, tSpace->topOrigin);

		getY(box[1].x, box[1].y, box[1].z, tSpace->u, tSpace->v, tSpace->bottomOrigin);
		getY(box[5].x, box[5].y, box[5].z, tSpace->u, tSpace->v, tSpace->topOrigin);
	} else { // iterate YZ plane
		box[0].y = box[3].y = box[4].y = box[7].y = blockIdx.x*blockDim.x - cBlobRadius;
		box[1].y = box[2].y = box[5].y = box[6].y = (blockIdx.x+1)*blockDim.x + cBlobRadius - 1.f;

		box[2].z = box[3].z = box[6].z = box[7].z = (blockIdx.y+1)*blockDim.y + cBlobRadius - 1.f;
		box[0].z = box[1].z = box[4].z = box[5].z = blockIdx.y*blockDim.y- cBlobRadius;

		getX(box[0].x, box[0].y, box[0].z, tSpace->u, tSpace->v, tSpace->bottomOrigin);
		getX(box[4].x, box[4].y, box[4].z, tSpace->u, tSpace->v, tSpace->topOrigin);

		getX(box[3].x, box[3].y, box[3].z, tSpace->u, tSpace->v, tSpace->bottomOrigin);
		getX(box[7].x, box[7].y, box[7].z, tSpace->u, tSpace->v, tSpace->topOrigin);

		getX(box[2].x, box[2].y, box[2].z, tSpace->u, tSpace->v, tSpace->bottomOrigin);
		getX(box[6].x, box[6].y, box[6].z, tSpace->u, tSpace->v, tSpace->topOrigin);

		getX(box[1].x, box[1].y, box[1].z, tSpace->u, tSpace->v, tSpace->bottomOrigin);
		getX(box[5].x, box[5].y, box[5].z, tSpace->u, tSpace->v, tSpace->topOrigin);
	}
	// transform AABB to the image domain
	rotate(box, tSpace->transformInv);
	// AABB is projected on image. Create new AABB that will encompass all vertices
	computeAABB(dest, box);
}

/**
 * Method returns true if AABB lies within the image boundaries
 */
__device__
bool isWithin(Point3D<float>* AABB, int imgXSize, int imgYSize) {
	return (AABB[0].x < imgXSize)
			&& (AABB[1].x >= 0)
			&& (AABB[0].y < imgYSize)
			&& (AABB[1].y >= 0);
}

/**
 * Method will load data from image at position tXindex, tYindex
 * and return them.
 * In case the data lies outside of the image boundaries, zeros (0,0)
 * are returned
 */
__device__
void getImgData(Point3D<float>* AABB,
		int tXindex, int tYindex,
		RecFourierBufferDataGPU* const buffer, int imgIndex,
		float& vReal, float& vImag) {
	int imgXindex = tXindex + AABB[0].x;
	int imgYindex = tYindex + AABB[0].y;
	if ((imgXindex >=0)
			&& (imgXindex < buffer->fftSizeX)
			&& (imgYindex >=0)
			&& (imgYindex < buffer->fftSizeY))	{
		int index = imgYindex * buffer->fftSizeX + imgXindex; // copy data from image
		vReal = buffer->getNthItem(buffer->FFTs, imgIndex)[2*index];
		vImag = buffer->getNthItem(buffer->FFTs, imgIndex)[2*index + 1];

	} else {
		vReal = vImag = 0.f; // out of image bound, so return zero
	}
}

/**
 * Method will copy imgIndex(th) data from buffer
 * to given destination (shared memory).
 * Only data within AABB will be copied.
 * Destination is expected to be continuous array of sufficient
 * size (imgCacheDim^2)
 */
__device__
void copyImgToCache(float2* dest, Point3D<float>* AABB,
		RecFourierBufferDataGPU* const buffer, int imgIndex,
		 int imgCacheDim) {
	for (int y = threadIdx.y; y < imgCacheDim; y += blockDim.y) {
		for (int x = threadIdx.x; x < imgCacheDim; x += blockDim.x) {
			int memIndex = y * imgCacheDim + x;
			getImgData(AABB, x, y, buffer, imgIndex, dest[memIndex].x, dest[memIndex].y);
		}
	}
}

/**
 * Method will use data stored in the buffer and update temporal
 * storages appropriately.
 */
template<bool useFast, bool hasCTF, int blobOrder>
__global__
void processBufferKernel(
		float* tempVolumeGPU, float *tempWeightsGPU,
		RecFourierBufferDataGPU* buffer,
		float* devBlobTableSqrt,
		int imgCacheDim) {
#if SHARED_BLOB_TABLE
	if ( ! useFast) {
		// copy blob table to shared memory
		volatile int id = threadIdx.y*blockDim.x + threadIdx.x;
		volatile int blockSize = blockDim.x * blockDim.y;
		for (int i = id; i < BLOB_TABLE_SIZE_SQRT; i+= blockSize)
			BLOB_TABLE[i] = devBlobTableSqrt[i];
		__syncthreads();
	}
#endif

	for (int i = 0; i < buffer->getNoOfSpaces(); i++) {
		RecFourierProjectionTraverseSpace* space = &buffer->spaces[i];

#if SHARED_IMG
		if ( ! useFast) {
			// make sure that all threads start at the same time
			// as they can come from previous iteration
			__syncthreads();
			if ((threadIdx.x == 0) && (threadIdx.y == 0)) {
				// first thread calculates which part of the image should be shared
				calculateAABB(space, buffer, SHARED_AABB);
			}
			__syncthreads();
			// check if the block will have to copy data from image
			if (isWithin(SHARED_AABB, buffer->fftSizeX, buffer->fftSizeY)) {
				// all threads copy image data to shared memory
				copyImgToCache(IMG, SHARED_AABB, buffer, space->projectionIndex, imgCacheDim);
				__syncthreads();
			} else {
				continue; // whole block can exit, as it's not reading from image
			}
		}
#endif

		processProjection<useFast, hasCTF, blobOrder>(
			tempVolumeGPU, tempWeightsGPU,
			buffer, space,
			devBlobTableSqrt,
			imgCacheDim);
		__syncthreads(); // sync threads to avoid write after read problems
	}
}

/**
 * Method will process the 'paddedFourier' (not shifted, i.e. low frequencies are in corners)
 * in the following way:
 * high frequencies are skipped (replaced by zero (0))
 * space is shifted, so that low frequencies are in the middle of the Y axis
 * resulting space is cropped.
 * Method returns a 2D array with Fourier coefficients, shifted so that low frequencies are
 * in the center of the Y axis (i.e. semicircle)
 */
__global__
void convertImagesKernel(std::complex<float>* iFouriers, int iSizeX, int iSizeY, int iLength,
		 RecFourierBufferDataGPU* oBuffer, float maxResolutionSqr) {
	// assign pixel to thread
	volatile int idx = blockIdx.x*blockDim.x + threadIdx.x;
	volatile int idy = blockIdx.y*blockDim.y + threadIdx.y;

	int halfY = iSizeY / 2;
	float normFactor = iSizeY*iSizeY;
	int oSizeX = oBuffer->fftSizeX;

	// input is an image in Fourier space (not normalized)
	// with low frequencies in the inner corners
	for (int n = 0; n < iLength; n++) {
		float2 freq;
		if ((idy < iSizeY) // for all input lines
				&& (idx < oSizeX)) { // for all output pixels in the line
			// process line only if it can hold sufficiently high frequency, i.e. process only
			// first and last N lines
			if (idy < oSizeX || idy >= (iSizeY - oSizeX)) {
				// check the frequency
				freq.x = FFT_IDX2DIGFREQ(idx, iSizeY);
				freq.y = FFT_IDX2DIGFREQ(idy, iSizeY);
				if ((freq.x * freq.x + freq.y * freq.y) > maxResolutionSqr) {
					continue;
				}
				// do the shift (lower line will move up, upper down)
				int newY = (idy < halfY) ? (idy + oSizeX) : (idy - iSizeY + oSizeX);
				int oIndex = newY*oSizeX + idx;

				int iIndex = n*iSizeY*iSizeX + idy*iSizeX + idx;
				float* iValue = (float*)&(iFouriers[iIndex]);

				// copy data and perform normalization
				oBuffer->getNthItem(oBuffer->FFTs, n)[2*oIndex] = iValue[0] / normFactor;
				oBuffer->getNthItem(oBuffer->FFTs, n)[2*oIndex + 1] = iValue[1] / normFactor;
			}
		}
	}
}

/**
 * Method takes padded input pictures, performs FFT and convert resulting images
 * as necessary for the algorithm.
 * Asynchronous method.
 */
void convertImages(
		FRecBufferDataGPUWrapper* wrapper,
		float maxResolutionSqr,
		int streamIndex) {

	cudaStream_t stream = streams[streamIndex];

	RecFourierBufferDataGPU* hostBuffer = wrapper->cpuCopy;
	// store to proper structure
	GpuMultidimArrayAtGpu<float> imagesGPU(
			hostBuffer->paddedImgSize, hostBuffer->paddedImgSize, 1, hostBuffer->noOfImages, hostBuffer->paddedImages);
	// perform FFT
	GpuMultidimArrayAtGpu<std::complex<float> > resultingFFT;
	mycufftHandle myhandle;
	imagesGPU.fft(resultingFFT, myhandle);
	myhandle.clear(); // release unnecessary memory
	imagesGPU.d_data = NULL; // unbind the data

	// now we have performed FFTs of the input images
	// buffers have to be updated accordingly
	hostBuffer->hasFFTs = true;
	cudaMemsetAsync(hostBuffer->FFTs, 0.f, hostBuffer->getFFTsByteSize(), stream); // clear it, as kernel writes only to some parts
	wrapper->copyToDevice(streamIndex);
	gpuErrchk( cudaPeekAtLastError() );

	// run kernel, one thread for each pixel of input FFT
	dim3 dimBlock(BLOCK_DIM, BLOCK_DIM);
	dim3 dimGrid(ceil(resultingFFT.Xdim/(float)dimBlock.x), ceil(resultingFFT.Ydim/(float)dimBlock.y));
	convertImagesKernel<<<dimGrid, dimBlock, 0, stream>>>(
			resultingFFT.d_data, resultingFFT.Xdim, resultingFFT.Ydim, resultingFFT.Ndim,
			wrapper->gpuCopy, maxResolutionSqr);
	// now we have converted input images to FFTs in the required format
}

void waitForGPU() {
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
}

void createStreams(int count) {
	streams = new cudaStream_t[count];
	for (int i = 0; i < count; i++) {
		cudaStreamCreate(&streams[i]);
	}
}

void deleteStreams(int count) {
	for (int i = 0; i < count; i++) {
		cudaStreamDestroy(streams[i]);
	}
	delete[] streams;
}

void allocateWrapper(RecFourierBufferData* buffer, int streamIndex) {
	wrappers[streamIndex] = new FRecBufferDataGPUWrapper(buffer);
}

void copyBlobTable(float* blobTableSqrt, int blobTableSize) {
	cudaMalloc((void **) &devBlobTableSqrt, blobTableSize*sizeof(float));
	cudaMemcpy(devBlobTableSqrt, blobTableSqrt, blobTableSize*sizeof(float), cudaMemcpyHostToDevice);
	gpuErrchk( cudaPeekAtLastError() );
}

void releaseBlobTable() {
	cudaFree(devBlobTableSqrt);
	gpuErrchk( cudaPeekAtLastError() );
}

void releaseWrapper(int streamIndex) {
	delete wrappers[streamIndex];
}

void copyConstants(
		int maxVolIndexX, int maxVolIndexYZ,
		float blobRadius, float blobAlpha,
		float iDeltaSqrt, float iw0) {
	cudaMemcpyToSymbol(cMaxVolumeIndexX, &maxVolIndexX,sizeof(maxVolIndexX));
	cudaMemcpyToSymbol(cMaxVolumeIndexYZ, &maxVolIndexYZ,sizeof(maxVolIndexYZ));
	cudaMemcpyToSymbol(cBlobRadius, &blobRadius,sizeof(blobRadius));
	cudaMemcpyToSymbol(cBlobAlpha, &blobAlpha,sizeof(blobAlpha));
	cudaMemcpyToSymbol(cIw0, &iw0,sizeof(iw0));
	cudaMemcpyToSymbol(cIDeltaSqrt, &iDeltaSqrt,sizeof(iDeltaSqrt));
	gpuErrchk( cudaPeekAtLastError() );
}

/**
 * Method will use data stored in the buffer and update temporal
 * storages appropriately.
 * Acuall calculation is done asynchronously, but 'buffer' can be reused
 * once the method returns.
 */
template<int blobOrder>
void processBufferGPU_(float* tempVolumeGPU, float* tempWeightsGPU,
		RecFourierBufferData* buffer,
		float blobRadius, int maxVolIndexYZ, bool useFast,
		float maxResolutionSqr, int streamIndex) {

	cudaStream_t stream = streams[streamIndex];

	// copy all data to gpu
	FRecBufferDataGPUWrapper* wrapper = wrappers[streamIndex];
	wrapper->copyFrom(buffer, streamIndex);
	wrapper->copyToDevice(streamIndex);

	// process input data if necessary
	if ( ! wrapper->cpuCopy->hasFFTs) {
		convertImages(wrapper, maxResolutionSqr, streamIndex);
	}
	// now wait till all necessary data are loaded to GPU (so that host can continue in work)
	cudaStreamSynchronize(stream);

	// enqueue kernel and return control
	int size2D = maxVolIndexYZ + 1;
	int imgCacheDim = ceil(sqrt(2.f) * sqrt(3.f) *(BLOCK_DIM + 2*blobRadius));
	dim3 dimBlock(BLOCK_DIM, BLOCK_DIM);
	dim3 dimGrid(ceil(size2D/(float)dimBlock.x),ceil(size2D/(float)dimBlock.y));

	// by using templates, we can save some registers, especially for 'fast' version
	if (useFast && buffer->hasCTFs) {
		processBufferKernel<true, true, blobOrder><<<dimGrid, dimBlock, 0, stream>>>(
			tempVolumeGPU, tempWeightsGPU,
			wrapper->gpuCopy,
			devBlobTableSqrt,
			imgCacheDim);
		   return;
   }
   if (useFast && !buffer->hasCTFs) {
	   processBufferKernel<true, false, blobOrder><<<dimGrid, dimBlock, 0, stream>>>(
				tempVolumeGPU, tempWeightsGPU,
				wrapper->gpuCopy,
				devBlobTableSqrt,
				imgCacheDim);
	   return;
   }
   // if making copy of the image in shared memory, allocate enough space
   int sharedMemSize = SHARED_IMG ? (imgCacheDim*imgCacheDim*sizeof(float2)) : 0;
   if (!useFast && buffer->hasCTFs) {
	   processBufferKernel<false, true, blobOrder><<<dimGrid, dimBlock, sharedMemSize, stream>>>(
			tempVolumeGPU, tempWeightsGPU,
			wrapper->gpuCopy,
			devBlobTableSqrt,
			imgCacheDim);
	   return;
   }
   if (!useFast && !buffer->hasCTFs) {
	   processBufferKernel<false, false, blobOrder><<<dimGrid, dimBlock, sharedMemSize, stream>>>(
			tempVolumeGPU, tempWeightsGPU,
			wrapper->gpuCopy,
			devBlobTableSqrt,
			imgCacheDim);
	   return;
   }
}

void processBufferGPU(float* tempVolumeGPU, float* tempWeightsGPU,
		RecFourierBufferData* buffer,
		float blobRadius, int maxVolIndexYZ, bool useFast,
		float maxResolutionSqr, int streamIndex, int blobOrder) {
	switch ((int)(blobOrder)) {
	case 0:
		processBufferGPU_<0>(tempVolumeGPU, tempWeightsGPU,
			buffer,
			blobRadius, maxVolIndexYZ, useFast,
			maxResolutionSqr,
			streamIndex);
		break;
	case 1:
		processBufferGPU_<1>(tempVolumeGPU, tempWeightsGPU,
				buffer,
				blobRadius, maxVolIndexYZ, useFast,
				maxResolutionSqr,
				streamIndex);
		break;
	case 2:
		processBufferGPU_<2>(tempVolumeGPU, tempWeightsGPU,
				buffer,
				blobRadius, maxVolIndexYZ, useFast,
				maxResolutionSqr,
				streamIndex);
		break;
	case 3:
		processBufferGPU_<3>(tempVolumeGPU, tempWeightsGPU,
				buffer,
				blobRadius, maxVolIndexYZ, useFast,
				maxResolutionSqr,
				streamIndex);
		break;
	case 4:
		processBufferGPU_<4>(tempVolumeGPU, tempWeightsGPU,
				buffer,
				blobRadius, maxVolIndexYZ, useFast,
				maxResolutionSqr,
				streamIndex);
		break;
	default:
		REPORT_ERROR(ERR_VALUE_INCORRECT, "m out of range [0..4] in kaiser_value()");
	}
}

