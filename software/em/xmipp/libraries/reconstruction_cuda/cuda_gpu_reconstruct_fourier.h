/*
 * cuda_gpu_reconstruct_fourier.h
 *
 *  Created on: Aug 11, 2017
 *      Author: david
 */

#ifndef XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_
#define XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_

#include <reconstruction_adapt_cuda/xmipp_gpu_utils.h>
#include <reconstruction_cuda/cuda_xmipp_utils.h>
#include <map>

#include "data/reconstruct_fourier_projection_traverse_space.h"
#include "data/reconstruct_fourier_buffer_data.h"
#include "data/reconstruct_fourier_defines.h"


//static ProjectionData* projData;



struct FRecBufferDataGPU : public RecFourierBufferData {

	FRecBufferDataGPU(RecFourierBufferData* orig);
	~FRecBufferDataGPU();

	float* getNthItem(float* array, int itemIndex);
	void copyDataFrom(RecFourierBufferData* orig, int stream);
	int getNoOfSpaces();

private:
	template<typename T>
	void copy(T* srcArray, T*& dstArray, RecFourierBufferData* orig, int stream);
	template<typename T>
	void alloc(T* srcArray, T*& dstArray, RecFourierBufferData* orig);
	void copyMetadata(RecFourierBufferData* orig);
};

struct FRecBufferDataGPUWrapper {

	FRecBufferDataGPUWrapper(RecFourierBufferData* orig);
	~FRecBufferDataGPUWrapper();
	void copyFrom(RecFourierBufferData* orig, int stream);
	void copyToDevice(int stream);

	FRecBufferDataGPU* cpuCopy;
	FRecBufferDataGPU* gpuCopy;
};

struct FourierReconstructionData
{
	FourierReconstructionData(int sizeX, int sizeY, int noOfImages, bool erase=true);
	void clean();
	float* getImgOnGPU(int imgIndex);
	float* dataOnGpu = NULL;
	int sizeX = 0;
	int sizeY = 0;
	int noOfImages = 0;

};



struct FourierReconDataWrapper
{
	FourierReconDataWrapper(int sizeX, int sizeY, int noOfImages);
	FourierReconDataWrapper(FourierReconstructionData* cpuCopy);
	~FourierReconDataWrapper();
	FourierReconstructionData* cpuCopy;
	FourierReconstructionData* gpuCopy;
private:
	void copyCpuToGpu();
};
/**
struct ProjectionDataGPU
{
	bool skip;
	float* img;
	float* CTF;
	float* modulator;
	int xSize;
	int ySize;
	int imgIndex;
	float weight;
	float localAInv[3][3];
	float localA[3][3];
public:
	ProjectionDataGPU() {
//		printf("constructor %p\n", this);
			setDefault();
		}
	~ProjectionDataGPU() {
//		printf("destructor %p, img %p\n", this, img);
	}
	ProjectionDataGPU(const ProjectionData& data) {
//		printf("constructor from ProjectionData %p\n", this);
		skip = data.skip;
		if (skip) {
			setDefault();
			return;
		}
//		copy(*data.img, img);
//		if (NULL != data.CTF) {
			img = NULL;
//			copy(*data.CTF, CTF);
//		} else {
			CTF = NULL;
//		}
//		if (NULL != data.modulator) {
//			copy(*data.modulator, modulator);
//		} else {
			modulator = NULL;
//		}
//		xSize = data.img->getXSize();
//		ySize = data.img->getYSize();
		imgIndex = data.imgIndex;
		weight = data.weight;
		data.localAInv.convertTo(localAInv);
		data.localA.convertTo(localA);
	}
	// Remove stored data and set to skip
	void clean();
	void setDefault() {
		skip = true;
		img = NULL;
		CTF = NULL;
		modulator = NULL;
		xSize = 0;
		ySize = 0;
		imgIndex = -1;
		weight = 0;
	}
private:
	template<typename T, typename U>
	void copy(const Array2D<T>& from, U& to);
};
*/


void allocateWrapper(RecFourierBufferData* buffer, int streamIndex);
void releaseWrapper(int streamIndex);


void createStreams(int count);

void deleteStreams(int count);

/** Method to allocate 3D array (continuous) of given size^3 */
float* allocateGPU(float*& where, int size, int typeSize);
void releaseGPU(float*& where);
void copyTempSpaces(std::complex<float>*** tempVol, float*** tempWeights,
		float* tempVolGPU, float* tempWeightsGPU,
		int size);



FourierReconDataWrapper* prepareBuffer(GpuMultidimArrayAtGpu<float>& ffts,
		int sizeX, int sizeY, int paddedImgSize, float maxResolutionSqr, int, float*&);

void waitForGPU();

void copyBlobTable(float* blobTableSqrt, int blobTableSize);

void releaseBlobTable();

void copyConstants(
		int maxVolIndexX, int maxVolIndexYZ,
		bool useFast, float blobRadius,
		float iDeltaSqrt);

void processBufferGPU(float* tempVolumeGPU, float* tempWeightsGPU,
		RecFourierBufferData* buffer,
		float blobRadius, int maxVolIndexYZ,
		float maxResolutionSqr, int stream);
void getTempSpaces(int size, std::complex<float>***& volume, float***& tempWeights);
// void copyBuffer(ProjectionData* data, int size);
void allocateTempSpaces(int size);

void a(std::complex<float>* a);
void b();
void c(std::complex<float>* a);

void moje();

static int num = 20;

#endif /* XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_ */
