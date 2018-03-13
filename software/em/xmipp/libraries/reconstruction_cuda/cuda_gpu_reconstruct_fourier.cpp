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
__device__ __constant__ float cIDeltaSqrt = 0.f;

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

__device__
void processPixel(
	float* tempVolumeGPU, float* tempWeightsGPU,
	RecFourierBufferDataGPU* const data,
	const RecFourierProjectionTraverseSpace* const space,
	int x, int y)
{
	Point3D<float> imgPos;
	const float* __restrict__ img = data->getNthItem(data->FFTs, space->projectionIndex);

	float dataWeight = space->weight;
	int index2D = y * data->fftSizeX + x;

	if (img[2*index2D] == 0.f && img[2*index2D + 1] == 0.f) {
		return;
	}

	// transform current point to center
	imgPos.x = x;
	imgPos.y = y - cMaxVolumeIndexYZ/2;
	imgPos.z = 0.f;
	// rotate around center
	multiply(space->transform, imgPos);
	// transform back and round
	// just Y coordinate needs adjusting, since X now matches to picture and Z is irrelevant
	int volX = clamp((int)(imgPos.x + cMaxVolumeIndexX/2 + 0.5f), 0, cMaxVolumeIndexX);
	int volY = clamp((int)(imgPos.y + cMaxVolumeIndexYZ/2 + 0.5f), 0, cMaxVolumeIndexYZ);
	int volZ = clamp((int)(imgPos.z + cMaxVolumeIndexYZ/2 + 0.5f), 0, cMaxVolumeIndexYZ);

	int index3D = volZ * (cMaxVolumeIndexYZ+1) * (cMaxVolumeIndexX+1) + volY * (cMaxVolumeIndexX+1) + volX;

	 // use atomic as two blocks can write to same voxel
	atomicAdd(&tempVolumeGPU[2*index3D], img[2*index2D] * dataWeight);
	atomicAdd(&tempVolumeGPU[2*index3D + 1], img[2*index2D + 1] * dataWeight);
	atomicAdd(&tempWeightsGPU[index3D], dataWeight);
}

__device__
void processPixelBlob(
	float* tempVolumeGPU, float* tempWeightsGPU,
	RecFourierBufferDataGPU* const data,
	const RecFourierProjectionTraverseSpace* const space,
	const float* blobTableSqrt,
	int x, int y)
{
	Point3D<float> imgPos;
	const float* __restrict__ img = data->getNthItem(data->FFTs, space->projectionIndex);

	float dataWeight = space->weight;
	int index2D = y * data->fftSizeX + x;

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

#if SHARED_BLOB_TABLE
				float wBlob = BLOB_TABLE[aux];
#else
				float wBlob = blobTableSqrt[aux];
#endif

				float weight = wBlob * dataWeight;

				atomicAdd(&tempVolumeGPU[2*index3D], img[2*index2D] * weight );
				atomicAdd(&tempVolumeGPU[2*index3D + 1], img[2*index2D + 1] * weight );
				atomicAdd(&tempWeightsGPU[index3D], weight);
			}
		}
	}
}

/**
 * Method will use data stored in the buffer and update temporal
 * storages appropriately.
 */
template<bool useFast>
__global__
void processBufferKernelInverse(
		float* tempVolumeGPU, float *tempWeightsGPU,
		RecFourierBufferDataGPU* buffer,
		float* devBlobTableSqrt,
		int imgCacheDim) {
#if TILE > 1
	int id = threadIdx.y * blockDim.x + threadIdx.x;
	int tidX = threadIdx.x % TILE + (id / (blockDim.y * TILE)) * TILE;
	int tidY = (id / TILE) % blockDim.y;
	int idx = blockIdx.x*blockDim.x + tidX;
	int idy = blockIdx.y*blockDim.y + tidY;
#else
	// map thread to each (2D) voxel
	volatile int idx = blockIdx.x*blockDim.x + threadIdx.x;
	volatile int idy = blockIdx.y*blockDim.y + threadIdx.y;
#endif


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

	if(idx >= buffer->fftSizeX  || idy >= buffer->fftSizeY) {
		return;
	}

	for (int i = blockIdx.z; i < buffer->getNoOfSpaces(); i += gridDim.z) {
		RecFourierProjectionTraverseSpace* space = &buffer->spaces[i];

		if (useFast) {
			processPixel(
				tempVolumeGPU, tempWeightsGPU,
				buffer, space, idx, idy);
		} else {
			processPixelBlob(tempVolumeGPU, tempWeightsGPU,
				buffer, space,
				devBlobTableSqrt, idx, idy);
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
		float blobRadius,
		float iDeltaSqrt) {
	cudaMemcpyToSymbol(cMaxVolumeIndexX, &maxVolIndexX,sizeof(maxVolIndexX));
	cudaMemcpyToSymbol(cMaxVolumeIndexYZ, &maxVolIndexYZ,sizeof(maxVolIndexYZ));
	cudaMemcpyToSymbol(cBlobRadius, &blobRadius,sizeof(blobRadius));
	cudaMemcpyToSymbol(cIDeltaSqrt, &iDeltaSqrt,sizeof(iDeltaSqrt));
	gpuErrchk( cudaPeekAtLastError() );
}

void setDevice(int device) {
	cudaSetDevice(device);
	gpuErrchk( cudaPeekAtLastError() );
}


/**
 * Method will use data stored in the buffer and update temporal
 * storages appropriately.
 * Actual calculation is done asynchronously, but 'buffer' can be reused
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
		if (useFast) {
			processBufferKernelInverse<true><<<dimGrid, dimBlock, 0, stream>>>(
					tempVolumeGPU, tempWeightsGPU,
					wrapper->gpuCopy,
					devBlobTableSqrt,
					imgCacheDim);
		} else {
			processBufferKernelInverse<false><<<dimGrid, dimBlock, 0, stream>>>(
					tempVolumeGPU, tempWeightsGPU,
					wrapper->gpuCopy,
					devBlobTableSqrt,
					imgCacheDim);
		}
}
