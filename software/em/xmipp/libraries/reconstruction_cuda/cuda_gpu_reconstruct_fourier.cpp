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
#include "reconstruction_cuda/cuda_basic_math.h"

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
__device__ __constant__ float cOneOverBlobRadiusSqr = 0.f;
__device__ __constant__ float cBlobAlpha = 0.f;
__device__ __constant__ float cIw0 = 0.f;
__device__ __constant__ float cIDeltaSqrt = 0.f;
__device__ __constant__ float cOneOverBessiOrderAlpha = 0.f;

__device__
float bessi0Fast(float x) { // X must be <= 15
	// stable rational minimax approximations to the modified bessel functions, blair, edwards
	// from table 5
	float x2 = x*x;
	float num = -0.8436825781374849e-19f; // p11
	num = fmaf(num, x2, -0.93466495199548700e-17f); // p10
	num = fmaf(num, x2, -0.15716375332511895e-13f); // p09
	num = fmaf(num, x2, -0.42520971595532318e-11f); // p08
	num = fmaf(num, x2, -0.13704363824102120e-8f);  // p07
	num = fmaf(num, x2, -0.28508770483148419e-6f);  // p06
	num = fmaf(num, x2, -0.44322160233346062e-4f);  // p05
	num = fmaf(num, x2, -0.46703811755736946e-2f);  // p04
	num = fmaf(num, x2, -0.31112484643702141e-0f);  // p03
	num = fmaf(num, x2, -0.11512633616429962e+2f);  // p02
	num = fmaf(num, x2, -0.18720283332732112e+3f);  // p01
	num = fmaf(num, x2, -0.75281108169006924e+3f);  // p00

	float den = 1.f; // q01
	den = fmaf(den, x2, -0.75281109410939403e+3f); // q00

	return num/den;
}

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
        ans = (expf(ax) * rsqrtf(ax)) * (0.39894228f + y * (0.1328592e-1f
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
        ans *= (expf(ax) * rsqrtf(ax));
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
float kaiserValue(float r, float a)
{
    float rda, rdas, arg, w;

    rda = r / a;
    if (rda <= 1.f)
    {
        rdas = rda * rda;
        arg = cBlobAlpha * sqrtf(1.f - rdas);
        if (order == 0)
        {
            w = bessi0(arg) * cOneOverBessiOrderAlpha;
        }
        else if (order == 1)
        {
            w = sqrtf (1.f - rdas);
			w *= bessi1(arg) * cOneOverBessiOrderAlpha;
        }
        else if (order == 2)
        {
            w = sqrtf (1.f - rdas);
            w = w * w;
			w *= bessi2(arg) * cOneOverBessiOrderAlpha;
        }
        else if (order == 3)
        {
            w = sqrtf (1.f - rdas);
            w = w * w * w;
			w *= bessi3(arg) * cOneOverBessiOrderAlpha;
        }
        else if (order == 4)
        {
            w = sqrtf (1.f - rdas);
            w = w * w * w *w;
			w *= bessi4(arg) * cOneOverBessiOrderAlpha;
        }
        else {
        	printf("order (%d) out of range in kaiser_value(): %s, %d\n", order, __FILE__, __LINE__);
        }
    }
    else
        w = 0.f;

    return w;
}

__device__
float kaiserValueFast(float distSqr) {
	float arg = cBlobAlpha * sqrtf(1.f - (distSqr * cOneOverBlobRadiusSqr)); // alpha * sqrt(1-(dist/blobRadius^2))
	return bessi0Fast(arg) * cOneOverBessiOrderAlpha * cIw0;
}

/**
 * Structure for buffer data on GPU
 * Adds some utility methods for copying data to GPU and device specific code.
 * Variables are stored in the same memory space as holding object, but
 * pointers point to the GPU memory
 */
struct RecFourierBufferDataGPU : public RecFourierBufferData {
private: // private to prevent unintended initialization
	RecFourierBufferDataGPU(RecFourierBufferData* orig) {};
	~RecFourierBufferDataGPU() {};
public:

	void create(RecFourierBufferData* orig) {
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

	void destroy() {
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
			cudaMemcpyAsync(dstArray, srcArray, bytes, cudaMemcpyHostToDevice, streams[stream]);
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
	void* ptr;
	cudaMallocHost(&ptr, sizeof(RecFourierBufferDataGPU)); // allocate page-locked
	cpuCopy = (RecFourierBufferDataGPU*)ptr;
	cpuCopy->create(orig);
	gpuCopy = NULL;
}

FRecBufferDataGPUWrapper::~FRecBufferDataGPUWrapper() {
	cudaFree(gpuCopy);
	gpuErrchk( cudaPeekAtLastError() );
	cpuCopy->destroy();
	cudaFreeHost(cpuCopy);
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
	size_t bytes = (size_t)size * size * size * typeSize;
	cudaMalloc((void**)&ptr, bytes);
	cudaMemset(ptr, 0.f, bytes);
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

/**
 * Calculates Z coordinate of the point [x, y] on the plane defined by p0 (origin) and normal
 */
__device__
float getZ(float x, float y, const Point3D<float>& n, const Point3D<float>& p0) {
	// from a(x-x0)+b(y-y0)+c(z-z0)=0
	return (-n.x*(x-p0.x)-n.y*(y-p0.y))/n.z + p0.z;
}

/**
 * Calculates Y coordinate of the point [x, z] on the plane defined by p0 (origin) and normal
 */
__device__
float getY(float x, float z, const Point3D<float>& n, const Point3D<float>& p0){
	// from a(x-x0)+b(y-y0)+c(z-z0)=0
	return (-n.x*(x-p0.x)-n.z*(z-p0.z))/n.y + p0.y;
}


/**
 * Calculates X coordinate of the point [y, z] on the plane defined by p0 (origin) and normal
 */
__device__
float getX(float y, float z, const Point3D<float>& n, const Point3D<float>& p0){
	// from a(x-x0)+b(y-y0)+c(z-z0)=0
	return (-n.y*(y-p0.y)-n.z*(z-p0.z))/n.x + p0.x;
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

		box[0].z = getZ(box[0].x, box[0].y, tSpace->unitNormal, tSpace->bottomOrigin);
		box[4].z = getZ(box[4].x, box[4].y, tSpace->unitNormal, tSpace->topOrigin);

		box[3].z = getZ(box[3].x, box[3].y, tSpace->unitNormal, tSpace->bottomOrigin);
		box[7].z = getZ(box[7].x, box[7].y, tSpace->unitNormal, tSpace->topOrigin);

		box[2].z = getZ(box[2].x, box[2].y, tSpace->unitNormal, tSpace->bottomOrigin);
		box[6].z = getZ(box[6].x, box[6].y, tSpace->unitNormal, tSpace->topOrigin);

		box[1].z = getZ(box[1].x, box[1].y, tSpace->unitNormal, tSpace->bottomOrigin);
		box[5].z = getZ(box[5].x, box[5].y, tSpace->unitNormal, tSpace->topOrigin);
	} else if (tSpace->XZ == tSpace->dir) { // iterate XZ plane
		box[0].x = box[3].x = box[4].x = box[7].x = blockIdx.x*blockDim.x - cBlobRadius;
		box[1].x = box[2].x = box[5].x = box[6].x = (blockIdx.x+1)*blockDim.x + cBlobRadius - 1.f;

		box[2].z = box[3].z = box[6].z = box[7].z = (blockIdx.y+1)*blockDim.y + cBlobRadius - 1.f;
		box[0].z = box[1].z = box[4].z = box[5].z = blockIdx.y*blockDim.y- cBlobRadius;

		box[0].y = getY(box[0].x, box[0].z, tSpace->unitNormal, tSpace->bottomOrigin);
		box[4].y = getY(box[4].x, box[4].z, tSpace->unitNormal, tSpace->topOrigin);

		box[3].y = getY(box[3].x, box[3].z, tSpace->unitNormal, tSpace->bottomOrigin);
		box[7].y = getY(box[7].x, box[7].z, tSpace->unitNormal, tSpace->topOrigin);

		box[2].y = getY(box[2].x, box[2].z, tSpace->unitNormal, tSpace->bottomOrigin);
		box[6].y = getY(box[6].x, box[6].z, tSpace->unitNormal, tSpace->topOrigin);

		box[1].y = getY(box[1].x, box[1].z, tSpace->unitNormal, tSpace->bottomOrigin);
		box[5].y = getY(box[5].x, box[5].z, tSpace->unitNormal, tSpace->topOrigin);
	} else { // iterate YZ plane
		box[0].y = box[3].y = box[4].y = box[7].y = blockIdx.x*blockDim.x - cBlobRadius;
		box[1].y = box[2].y = box[5].y = box[6].y = (blockIdx.x+1)*blockDim.x + cBlobRadius - 1.f;

		box[2].z = box[3].z = box[6].z = box[7].z = (blockIdx.y+1)*blockDim.y + cBlobRadius - 1.f;
		box[0].z = box[1].z = box[4].z = box[5].z = blockIdx.y*blockDim.y- cBlobRadius;

		box[0].x = getX(box[0].y, box[0].z, tSpace->unitNormal, tSpace->bottomOrigin);
		box[4].x = getX(box[4].y, box[4].z, tSpace->unitNormal, tSpace->topOrigin);

		box[3].x = getX(box[3].y, box[3].z, tSpace->unitNormal, tSpace->bottomOrigin);
		box[7].x = getX(box[7].y, box[7].z, tSpace->unitNormal, tSpace->topOrigin);

		box[2].x = getX(box[2].y, box[2].z, tSpace->unitNormal, tSpace->bottomOrigin);
		box[6].x = getX(box[6].y, box[6].z, tSpace->unitNormal, tSpace->topOrigin);

		box[1].x = getX(box[1].y, box[1].z, tSpace->unitNormal, tSpace->bottomOrigin);
		box[5].x = getX(box[5].y, box[5].z, tSpace->unitNormal, tSpace->topOrigin);
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

__device__
void processPixel(
	float* tempVolumeGPU, float* tempWeightsGPU,
	RecFourierBufferDataGPU* const data,
	const RecFourierProjectionTraverseSpace* const space)
{
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y*blockDim.y + threadIdx.y;
//	if(x >= data->fftSizeX  || y >= data->fftSizeY) {
//		return;
//	}

	Point3D<float> imgPos;
	float wBlob = 1.f;
	float wCTF = 1.f;
	float wModulator = 1.f;
	const float* __restrict__ img = data->getNthItem(data->FFTs, space->projectionIndex);

	float dataWeight = space->weight;
	int xSize = data->fftSizeX;
	int ySize = data->fftSizeY;
	int index2D = y * xSize + x;

	if (img[2*index2D] == 0.f && img[2*index2D + 1] == 0.f) {
		return;
	}

	// transform current point to center
//	imgPos.x = x - cMaxVolumeIndexX/2;
//	imgPos.y = y ;
//	imgPos.z = z - cMaxVolumeIndexYZ/2;
	imgPos.x = x;
	imgPos.y = y - cMaxVolumeIndexYZ/2;
	imgPos.z = 0.f;
	// rotate around center
	multiply(space->transform, imgPos);
//	if (imgPos.x*imgPos.x + imgPos.y*imgPos.y + imgPos.z*imgPos.z > space->maxDistanceSqr) {
//		return; // discard iterations that would access pixel with too high frequency
//	}
	// transform back and round
	// just Y coordinate needs adjusting, since X now matches to picture and Z is irrelevant
	int volX = clamp((int)(imgPos.x + cMaxVolumeIndexX/2 + 0.5f), 0, cMaxVolumeIndexX);
	int volY = clamp((int)(imgPos.y + cMaxVolumeIndexYZ/2 + 0.5f), 0, cMaxVolumeIndexYZ);
	int volZ = clamp((int)(imgPos.z + cMaxVolumeIndexYZ/2 + 0.5f), 0, cMaxVolumeIndexYZ);

	int index3D = volZ * (cMaxVolumeIndexYZ+1) * (cMaxVolumeIndexX+1) + volY * (cMaxVolumeIndexX+1) + volX;

//	if (hasCTF) {
//		const float* __restrict__ CTF = data->getNthItem(data->CTFs, space->projectionIndex);
//		const float* __restrict__ modulator = data->getNthItem(data->modulators, space->projectionIndex);
//		wCTF = CTF[index2D];
//		wModulator = modulator[index2D];
//	}

	float weight = wBlob * wModulator * dataWeight;

	 // use atomic as two blocks can write to same voxel
	atomicAdd(&tempVolumeGPU[2*index3D], img[2*index2D] * weight * wCTF);
	atomicAdd(&tempVolumeGPU[2*index3D + 1], img[2*index2D + 1] * weight * wCTF);
	atomicAdd(&tempWeightsGPU[index3D], weight);
}

__device__
void processPixelBlob(
	float* tempVolumeGPU, float* tempWeightsGPU,
	RecFourierBufferDataGPU* const data,
	const RecFourierProjectionTraverseSpace* const space,
	const float* blobTableSqrt)
{
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y*blockDim.y + threadIdx.y;

	Point3D<float> imgPos;
	float wBlob = 1.f;
	float wCTF = 1.f;
	float wModulator = 1.f;
	const float* __restrict__ img = data->getNthItem(data->FFTs, space->projectionIndex);

	float dataWeight = space->weight;
	int xSize = data->fftSizeX;
	int ySize = data->fftSizeY;
	int index2D = y * xSize + x;

	if (img[2*index2D] == 0.f && img[2*index2D + 1] == 0.f) {
		return;
	}

	imgPos.x = x;
	imgPos.y = y - cMaxVolumeIndexYZ/2;
	imgPos.z = 0.f;
	// rotate around center
	multiply(space->transform, imgPos);
	// transform to center of the volume
	float volX = imgPos.x + cMaxVolumeIndexX/2;
	float volY = imgPos.y + cMaxVolumeIndexYZ/2;
	float volZ = imgPos.z + cMaxVolumeIndexYZ/2;

	// create blob bounding box
	int minX = ceilf(volX - cBlobRadius);
	int maxX = floorf(volX + cBlobRadius);
	int minY = ceilf(volY - cBlobRadius);
	int maxY = floorf(volY + cBlobRadius);
	int minZ = ceilf(volZ - cBlobRadius);
	int maxZ = floorf(volZ + cBlobRadius);
	minX = fmaxf(minX, 0);
	minY = fmaxf(minY, 0);
	minZ = fmaxf(minZ, 0);
	maxX = fminf(maxX, cMaxVolumeIndexX);
	maxY = fminf(maxY, cMaxVolumeIndexYZ);
	maxZ = fminf(maxZ, cMaxVolumeIndexYZ);
	float radiusSqr = cBlobRadius * cBlobRadius;
	for (int l = minZ; l <= maxZ; l++) {
		float zSqr = (volZ - l) * (volZ - l);
		for (int i = minY; i <= maxY; i++) {
			float ySqr = (volY - i) * (volY - i);
			float yzSqr = ySqr + zSqr;
			if (yzSqr > radiusSqr) continue;
			for (int j = minX; j <= maxX; j++) {
				float xD = volX - j;
				float distanceSqr = xD*xD + yzSqr;
				if (distanceSqr > radiusSqr) continue;

				int index3D = l * (cMaxVolumeIndexYZ+1) * (cMaxVolumeIndexX+1) + i * (cMaxVolumeIndexX+1) + j;

				int aux = (int) ((distanceSqr * cIDeltaSqrt + 0.5f));

				float wBlob = blobTableSqrt[aux];

				float weight = wBlob * dataWeight;

//				if (index3D < 0 || index3D >= (cMaxVolumeIndexX * cMaxVolumeIndexYZ * cMaxVolumeIndexYZ)) {
//					printf("prusvih, X %d %d Y %d %d Z %d %d, ted %d %d %d\n",
//							minX, maxX, minY, maxY, minZ, maxZ, j, i, l );
//				}

				atomicAdd(&tempVolumeGPU[2*index3D], img[2*index2D] * weight * wCTF);
				atomicAdd(&tempVolumeGPU[2*index3D + 1], img[2*index2D + 1] * weight * wCTF);
				atomicAdd(&tempWeightsGPU[index3D], weight);
			}
		}
	}



//	if (hasCTF) {
//		const float* __restrict__ CTF = data->getNthItem(data->CTFs, space->projectionIndex);
//		const float* __restrict__ modulator = data->getNthItem(data->modulators, space->projectionIndex);
//		wCTF = CTF[index2D];
//		wModulator = modulator[index2D];
//	}

//	float weight = wBlob * wModulator * dataWeight;
//
//	 // use atomic as two blocks can write to same voxel
//	atomicAdd(&tempVolumeGPU[2*index3D], img[2*index2D] * weight * wCTF);
//	atomicAdd(&tempVolumeGPU[2*index3D + 1], img[2*index2D + 1] * weight * wCTF);
//	atomicAdd(&tempWeightsGPU[index3D], weight);
}




/**
 * Method will use data stored in the buffer and update temporal
 * storages appropriately.
 */
//template<bool useFast>
__global__
void processBufferKernelInverse(
		float* tempVolumeGPU, float *tempWeightsGPU,
		RecFourierBufferDataGPU* buffer,
		float* devBlobTableSqrt,
		int imgCacheDim,bool useFast) {

	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y*blockDim.y + threadIdx.y;

	if(x >= buffer->fftSizeX  || y >= buffer->fftSizeY) {
		return;
	}

#if SHARED_BLOB_TABLE
	if ( ! useFast) {
		// copy blob table to shared memory
		int id = threadIdx.y*blockDim.x + threadIdx.x;
		int blockSize = blockDim.x * blockDim.y;
		for (int i = id; i < BLOB_TABLE_SIZE_SQRT; i+= blockSize)
			BLOB_TABLE[i] = devBlobTableSqrt[i];
		__syncthreads();
	}
#endif


	for (int i = 0; i < buffer->getNoOfSpaces(); i++) {
		RecFourierProjectionTraverseSpace* space = &buffer->spaces[i];

		if (useFast) {
		processPixel(
			tempVolumeGPU, tempWeightsGPU,
			buffer, space);
		} else {
			processPixelBlob(tempVolumeGPU, tempWeightsGPU,
					buffer, space,
					devBlobTableSqrt);
		}
		__syncthreads(); // sync threads to avoid write after read problems
	}
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


void pinMemory(RecFourierBufferData* buffer) {
	hostRegister(buffer->CTFs, buffer->getMaxByteSize(buffer->CTFs));
	hostRegister(buffer->FFTs, buffer->getMaxByteSize(buffer->FFTs));
	hostRegister(buffer->paddedImages, buffer->getMaxByteSize(buffer->paddedImages));
	hostRegister(buffer->modulators, buffer->getMaxByteSize(buffer->modulators));
	hostRegister(buffer->spaces, buffer->getMaxByteSize(buffer->spaces));
	hostRegister(buffer, sizeof(*buffer));
}

void unpinMemory(RecFourierBufferData* buffer) {
	hostUnregister(buffer->CTFs);
	hostUnregister(buffer->FFTs);
	hostUnregister(buffer->paddedImages);
	hostUnregister(buffer->modulators);
	hostUnregister(buffer->spaces);
	hostUnregister(buffer);
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
		float iDeltaSqrt, float iw0, float oneOverBessiOrderAlpha) {
	cudaMemcpyToSymbol(cMaxVolumeIndexX, &maxVolIndexX,sizeof(maxVolIndexX));
	cudaMemcpyToSymbol(cMaxVolumeIndexYZ, &maxVolIndexYZ,sizeof(maxVolIndexYZ));
	cudaMemcpyToSymbol(cBlobRadius, &blobRadius,sizeof(blobRadius));
	cudaMemcpyToSymbol(cBlobAlpha, &blobAlpha,sizeof(blobAlpha));
	cudaMemcpyToSymbol(cIw0, &iw0,sizeof(iw0));
	cudaMemcpyToSymbol(cIDeltaSqrt, &iDeltaSqrt,sizeof(iDeltaSqrt));
	cudaMemcpyToSymbol(cOneOverBessiOrderAlpha, &oneOverBessiOrderAlpha,sizeof(oneOverBessiOrderAlpha));
	float oneOverBlobRadiusSqr = 1.f / (blobRadius * blobRadius);
	cudaMemcpyToSymbol(cOneOverBlobRadiusSqr, &oneOverBlobRadiusSqr,sizeof(oneOverBlobRadiusSqr));
	gpuErrchk( cudaPeekAtLastError() );
}


/**
 * Method will use data stored in the buffer and update temporal
 * storages appropriately.
 * Acuall calculation is done asynchronously, but 'buffer' can be reused
 * once the method returns.
 */
void processBufferGPUInverse(float* tempVolumeGPU, float* tempWeightsGPU,
		RecFourierBufferData* buffer,
		float blobRadius, int maxVolIndexYZ, bool useFast,
		float maxResolutionSqr, int streamIndex) {
	cudaStream_t stream = streams[streamIndex];

		// copy all data to gpu
		FRecBufferDataGPUWrapper* wrapper = wrappers[streamIndex];
		wrapper->copyFrom(buffer, streamIndex);
		wrapper->copyToDevice(streamIndex);

		// process input data if necessary
		// FIXME unsupported
		// now wait till all necessary data are loaded to GPU (so that host can continue in work)
		cudaStreamSynchronize(stream);

		// enqueue kernel and return control
		int size2D = maxVolIndexYZ + 1;
		int imgCacheDim = ceil(sqrt(2.f) * sqrt(3.f) *(BLOCK_DIM + 2*blobRadius));
		dim3 dimBlock(BLOCK_DIM, BLOCK_DIM);
		dim3 dimGrid(ceil(size2D/(float)dimBlock.x),ceil(size2D/(float)dimBlock.y), GRID_DIM_Z);

		// by using templates, we can save some registers, especially for 'fast' version
		processBufferKernelInverse<<<dimGrid, dimBlock, 0, stream>>>(
				tempVolumeGPU, tempWeightsGPU,
				wrapper->gpuCopy,
				devBlobTableSqrt,
				imgCacheDim, useFast);
}




/**
 * Method will use data stored in the buffer and update temporal
 * storages appropriately.
 * Actual calculation is done asynchronously, but 'buffer' can be reused
 * once the method returns.
 */
template<int blobOrder, bool useFastKaiser>
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
	// FIXME unsupported
	// now wait till all necessary data are loaded to GPU (so that host can continue in work)
	cudaStreamSynchronize(stream);

	// enqueue kernel and return control
	int size2D = maxVolIndexYZ + 1;
	int imgCacheDim = ceil(sqrt(2.f) * sqrt(3.f) *(BLOCK_DIM + 2*blobRadius));
	dim3 dimBlock(BLOCK_DIM, BLOCK_DIM);
	dim3 dimGrid(ceil(size2D/(float)dimBlock.x),ceil(size2D/(float)dimBlock.y), GRID_DIM_Z);

	// by using templates, we can save some registers, especially for 'fast' version
	if (useFast && buffer->hasCTFs) {
//		processBufferKernel<true, true, blobOrder,useFastKaiser><<<dimGrid, dimBlock, 0, stream>>>(
//			tempVolumeGPU, tempWeightsGPU,
//			wrapper->gpuCopy,
//			devBlobTableSqrt,
//			imgCacheDim);
		   return;
   }
   if (useFast && !buffer->hasCTFs) {
//	   processBufferKernel<true, false, blobOrder,useFastKaiser><<<dimGrid, dimBlock, 0, stream>>>(
//				tempVolumeGPU, tempWeightsGPU,
//				wrapper->gpuCopy,
//				devBlobTableSqrt,
//				imgCacheDim);
	   return;
   }
   // if making copy of the image in shared memory, allocate enough space
   int sharedMemSize = SHARED_IMG ? (imgCacheDim*imgCacheDim*sizeof(float2)) : 0;
   if (!useFast && buffer->hasCTFs) {
//	   processBufferKernel<false, true, blobOrder,useFastKaiser><<<dimGrid, dimBlock, sharedMemSize, stream>>>(
//			tempVolumeGPU, tempWeightsGPU,
//			wrapper->gpuCopy,
//			devBlobTableSqrt,
//			imgCacheDim);
	   return;
   }
   if (!useFast && !buffer->hasCTFs) {
//	   processBufferKernel<false, false, blobOrder,useFastKaiser><<<dimGrid, dimBlock, sharedMemSize, stream>>>(
//			tempVolumeGPU, tempWeightsGPU,
//			wrapper->gpuCopy,
//			devBlobTableSqrt,
//			imgCacheDim);
	   return;
   }
}

void setDevice(int device) {
	cudaSetDevice(device);
	gpuErrchk( cudaPeekAtLastError() );
}

void processBufferGPU(float* tempVolumeGPU, float* tempWeightsGPU,
		RecFourierBufferData* buffer,
		float blobRadius, int maxVolIndexYZ, bool useFast,
		float maxResolutionSqr, int streamIndex, int blobOrder, float blobAlpha) {
	switch (blobOrder) {
	case 0:
		if (blobAlpha <= 15.0) {
			processBufferGPU_<0, true>(tempVolumeGPU, tempWeightsGPU,
				buffer,
				blobRadius, maxVolIndexYZ, useFast,
				maxResolutionSqr,
				streamIndex);
		} else {
			processBufferGPU_<0, false>(tempVolumeGPU, tempWeightsGPU,
				buffer,
				blobRadius, maxVolIndexYZ, useFast,
				maxResolutionSqr,
				streamIndex);
		}
		break;
	case 1:
		processBufferGPU_<1, false>(tempVolumeGPU, tempWeightsGPU,
				buffer,
				blobRadius, maxVolIndexYZ, useFast,
				maxResolutionSqr,
				streamIndex);
		break;
	case 2:
		processBufferGPU_<2, false>(tempVolumeGPU, tempWeightsGPU,
				buffer,
				blobRadius, maxVolIndexYZ, useFast,
				maxResolutionSqr,
				streamIndex);
		break;
	case 3:
		processBufferGPU_<3, false>(tempVolumeGPU, tempWeightsGPU,
				buffer,
				blobRadius, maxVolIndexYZ, useFast,
				maxResolutionSqr,
				streamIndex);
		break;
	case 4:
		processBufferGPU_<4, false>(tempVolumeGPU, tempWeightsGPU,
				buffer,
				blobRadius, maxVolIndexYZ, useFast,
				maxResolutionSqr,
				streamIndex);
		break;
	default:
		REPORT_ERROR(ERR_VALUE_INCORRECT, "m out of range [0..4] in kaiser_value()");
	}
}

