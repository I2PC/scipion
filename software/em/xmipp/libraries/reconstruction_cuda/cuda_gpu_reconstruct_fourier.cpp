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
//#include "reconstruct_fourier.cu"

__device__ __constant__ int cMaxVolumeIndexX = 0;
__device__ __constant__ int cMaxVolumeIndexYZ = 0;
__device__ __constant__ float cBlobRadius = 0.f;
__device__ __constant__ float cBlobAlpha = 0.f;
__device__ __constant__ float cIw0 = 0.f;
__device__ __constant__ float cIDeltaSqrt = 0.f;

float* devBlobTableSqrt = NULL;

// Holding streams used for calculation. Present on CPU
cudaStream_t* streams;

// Wrapper to hold pointers to GPU memory (and have it also accessible from CPU)
std::map<int,FRecBufferDataGPUWrapper*> wrappers;
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
	__device__ __host__
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
//	cudaMalloc((void **) &devBlobTableSqrt, blobTableSize*sizeof(float));
//	cudaMemcpy(devBlobTableSqrt, blobTableSqrt, blobTableSize*sizeof(float), cudaMemcpyHostToDevice);
//	gpuErrchk( cudaPeekAtLastError() );
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
template<int>
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


//	test(size2D,
//			std::pow(size2D, 3), wrapper->cpuCopy->noOfImages,
//			tempVolumeGPU, tempWeightsGPU,
//			wrapper->cpuCopy->spaces, wrapper->cpuCopy->getNoOfSpaces(),
//			wrapper->cpuCopy->FFTs, wrapper->cpuCopy->CTFs, wrapper->cpuCopy->modulators,
//			wrapper->cpuCopy->fftSizeX, wrapper->cpuCopy->fftSizeY,
//			devBlobTableSqrt,
//			imgCacheDim);


	// by using templates, we can save some registers, especially for 'fast' version
//		processBufferKernel<<<dimGrid, dimBlock, 0, stream>>>(
//			tempVolumeGPU, tempWeightsGPU,
//			wrapper->cpuCopy->spaces, wrapper->cpuCopy->getNoOfSpaces(),
//			wrapper->cpuCopy->FFTs, wrapper->cpuCopy->CTFs, wrapper->cpuCopy->modulators,
//			wrapper->cpuCopy->fftSizeX, wrapper->cpuCopy->fftSizeY,
//			devBlobTableSqrt,
//			imgCacheDim);
}



void processBufferGPU(float* tempVolumeGPU, float* tempWeightsGPU,
		RecFourierBufferData* buffer,
		float blobRadius, int maxVolIndexYZ, bool useFast,
		float maxResolutionSqr, int stream, int blobOrder){
	int streamIndex = 0;
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

