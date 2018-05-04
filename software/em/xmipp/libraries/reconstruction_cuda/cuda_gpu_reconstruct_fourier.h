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

#ifndef XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_
#define XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_

#include <map>

#include <reconstruction/reconstruct_fourier_projection_traverse_space.h>
#include <reconstruction/reconstruct_fourier_buffer_data.h>
#include <reconstruction/reconstruct_fourier_defines.h>
#include <data/xmipp_error.h>

#include <reconstruction_cuda/cuda_xmipp_utils.h>

struct RecFourierBufferDataGPU;

/**
 * This struct simplifies work with pointers on GPU
 * It holds a copy both in CPU and GPU memory space
 */
struct FRecBufferDataGPUWrapper {
	FRecBufferDataGPUWrapper(RecFourierBufferData* orig);
	~FRecBufferDataGPUWrapper();

	/**
	 * Copy the content of the 'orig' buffer to new object
	 * Original buffer can be cleaned after.
	 * Works asynchronously on given stream
	 */
	void copyFrom(RecFourierBufferData* orig, int stream);

	/**
	 * Copy CPU version to GPU, asynchronously on given stream
	 */
	void copyToDevice(int stream);

	// object in CPU memory space (page-locked)
	RecFourierBufferDataGPU* cpuCopy;

	// object in GPU memory space
	RecFourierBufferDataGPU* gpuCopy;
};

/** set device to use for calculation (has to be done by each thread) */
void setDevice(int device);

/**
 * Method will allocate buffer wrapper for given stream
 * 'buffer' is used for size references, i.e. it has to have the same
 * sizes that will be used later, during calculation
 * Blocking operation
 */
void allocateWrapper(RecFourierBufferData* buffer, int streamIndex);

/**
 * Release all resources allocated with the buffer wrapper, for given stream
 * Blocking operation
 */
void releaseWrapper(int streamIndex);

/**
 * Allocate 'count' streams on default GPU
 * Blocking operation
 */
void createStreams(int count);

/**
 * Delete streams allocated on default GPU
 * Blocking operation
 */
void deleteStreams(int count);

/**
 *  Method to allocate 3D array (continuous) of given size^3
 *  Allocated array is cleared (to 0 zero)
 *  Blocking operation
 */
float* allocateTempVolumeGPU(float*& ptr, int size, int typeSize);

/**
 * Release memory at GPU
 * Blocking operation
 */
void releaseTempVolumeGPU(float*& ptr);

/**
 * Method will copy continuous 3D arrays (with side of size) from GPU
 * to non-continuous arrays on CPU
 * Blocking operation
 */
void copyTempVolumes(std::complex<float>*** tempVol, float*** tempWeights,
		float* tempVolGPU, float* tempWeightsGPU,
		int size);

/**
 * Blocking method. Once returns, all operation on default GPU are done
 */
void waitForGPU();

/**
 * Method will allocate space at GPU and copy there content of the table
 * Blocking operation
 */
void copyBlobTable(float* blobTableSqrt, int size);

/**
 * Method will release all resources allocated for the blob table at GPU
 * Blocking operation
 */
void releaseBlobTable();

void pinMemory(RecFourierBufferData* buffer);

void unpinMemory(RecFourierBufferData* buffer);

/**
 * Method will copy constants used for calculation to GPU memory
 * Blocking operation
 */
void copyConstants(
		int maxVolIndexX, int maxVolIndexYZ,
		float blobRadius, float blobAlpha,
		float iDeltaSqrt, float iw0, float oneOverBessiOrderAlpha);

/**
 * Method will copy content of the 'buffer' to GPU and
 * runs the calculation (asynchronously on given stream).
 * Once it returns, 'buffer' object can be reused.
 * See also createStreams(int)
 */
void processBufferGPU(float* tempVolumeGPU, float* tempWeightsGPU,
		RecFourierBufferData* buffer,
		float blobRadius, int maxVolIndexYZ, bool useFast,
		float maxResolutionSqr, int stream, int blobOrder, float blobAlpha);

#endif /* XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_ */
