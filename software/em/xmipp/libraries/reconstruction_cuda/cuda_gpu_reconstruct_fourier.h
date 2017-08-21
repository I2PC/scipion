/*
 * cuda_gpu_reconstruct_fourier.h
 *
 *  Created on: Aug 11, 2017
 *      Author: david
 */

#ifndef XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_
#define XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_

#include <data/projection_data.h>
#include <data/fourier_reconstruction_traverse_space.h>


//static ProjectionData* projData;


//void copyProjectionData(ProjectionData* data);
//static ProjectionDataGPU* copyProjectionData(ProjectionData& data);


/** Method to allocate 3D array (continuous) of given size^3 */
float* allocateGPU(float*& where, int size);
void releaseGPU(float*& where);
void copyTempSpaces(std::complex<float>*** tempVol, float*** tempWeights,
		float* tempVolGPU, float* tempWeightsGPU,
		int size);







void processBufferGPU(
		float* tempVolumeGPU,
				float* tempWeightsGPU,
				ProjectionData* data,int noOfPixel,int bufferSize,
				TraverseSpace* traverseSpaces, MATRIX* transformsInv, int noOfTransforms,
				int maxVolIndexX, int maxVolIndexYZ,
				bool useFast, float blobRadius,
				float iDeltaSqrt,
				float* blobTableSqrt, int blobTableSize);
void getTempSpaces(int size, std::complex<float>***& volume, float***& tempWeights);
void copyBuffer(ProjectionData* data, int size);
void allocateTempSpaces(int size);

void a(std::complex<float>* a);
void b();
void c(std::complex<float>* a);

void moje();

static int num = 20;

#endif /* XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_ */
