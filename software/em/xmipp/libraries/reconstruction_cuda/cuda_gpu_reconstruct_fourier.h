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
#include <reconstruction_adapt_cuda/xmipp_gpu_utils.h>
#include <reconstruction_cuda/cuda_xmipp_utils.h>



//static ProjectionData* projData;

struct FRecBufferData
{
	FRecBufferData(bool hasFFTs, bool hasCTFs,
			int fftSizeX, int fftSizeY,	int paddedImgSize,
			int maxNoOfImages, int noOfSymmetries) :
			hasFFTs(hasFFTs), hasCTFs(hasCTFs),
			fftSizeX(fftSizeX), fftSizeY(fftSizeY), paddedImgSize(paddedImgSize) {
		if (hasFFTs) {
			paddedImages = NULL;
			FFTs = new float[fftSizeX * fftSizeY * maxNoOfImages * 2]; // *2 since it's complex
		} else {
			FFTs = NULL;
			paddedImages = new float[paddedImgSize * paddedImgSize * maxNoOfImages];
		}

		if (hasCTFs) {
			CTFs = new float[fftSizeX * fftSizeY * maxNoOfImages];
			modulators = new float[fftSizeX * fftSizeY * maxNoOfImages];
		} else {
			CTFs = modulators = NULL;
		}

		spaces = new TraverseSpace[maxNoOfImages * noOfSymmetries];
		noOfImages = noOfSpaces = 0;
	};

	~FRecBufferData() {
		delete[] FFTs;
		delete[] CTFs;
		delete[] paddedImages;
		delete[] modulators;
		FFTs = CTFs = paddedImages = modulators = NULL;

		delete[] spaces;
		spaces = NULL;
	}

	int getPaddedImgSize() {
		return paddedImgSize * paddedImgSize;
	}

	int getPaddedImgByteSize() {
		return getPaddedImgSize() * sizeof(float);
	}

	float* getNthItem(float* array, int itemIndex) {
		if (array == FFTs) return array + (fftSizeX * fftSizeY * itemIndex * 2); // *2 since it's complex
		if (array == CTFs) return array + (fftSizeX * fftSizeY * itemIndex);
		if (array == modulators) return array + (fftSizeX * fftSizeY * itemIndex);
		if (array == paddedImages) return array + (paddedImgSize * paddedImgSize * itemIndex);
		return NULL; // undefined
	}

	int getNoOfElements(float* array) {
		if (array == FFTs) return (fftSizeX * fftSizeY * noOfImages * 2); // *2 since it's complex
		if (array == CTFs) return (fftSizeX * fftSizeY * noOfImages);
		if (array == modulators) return (fftSizeX * fftSizeY * noOfImages);
		if (array == paddedImages) return (paddedImgSize * paddedImgSize * noOfImages);
		return -1; // undefined
	}

	float* FFTs;
	float* CTFs;
	float* paddedImages;
	float* modulators;
	TraverseSpace* spaces;
	int noOfSpaces;
	int fftSizeX;
	int fftSizeY;
	int paddedImgSize;
	int noOfImages;
	bool hasFFTs;
	bool hasCTFs;
};

struct FRecBufferDataGPU : FRecBufferData {
	float* getNthItem(float* array, int itemIndex);
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
	/** Remove stored data and set to skip */
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


/** Method to allocate 3D array (continuous) of given size^3 */
float* allocateGPU(float*& where, int size, int typeSize);
void releaseGPU(float*& where);
void copyTempSpaces(std::complex<float>*** tempVol, float*** tempWeights,
		float* tempVolGPU, float* tempWeightsGPU,
		int size);



FourierReconDataWrapper* prepareBuffer(GpuMultidimArrayAtGpu<float>& ffts,
		int sizeX, int sizeY, int paddedImgSize, float maxResolutionSqr, int, float*&);



void processBufferGPU(float* tempVolumeGPU, float* tempWeightsGPU,
		float* paddedImages, float* readyFFTs, int noOfImages,
		TraverseSpace* traverseSpaces, int noOfTransforms,
		int maxVolIndexX, int maxVolIndexYZ,
		bool useFast, float blobRadius,
		float iDeltaSqrt,
		float* blobTableSqrt, int blobTableSize,
		int paddedImgSize, float maxResolutionSqr);
void getTempSpaces(int size, std::complex<float>***& volume, float***& tempWeights);
void copyBuffer(ProjectionData* data, int size);
void allocateTempSpaces(int size);

void a(std::complex<float>* a);
void b();
void c(std::complex<float>* a);

void moje();

static int num = 20;

#endif /* XMIPP_LIBRARIES_RECONSTRUCTION_CUDA_CUDA_GPU_RECONSTRUCT_FOURIER_H_ */
