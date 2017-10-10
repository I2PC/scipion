/*
 * fourier_reconstruction_buffer.h
 *
 *  Created on: Oct 9, 2017
 *      Author: david
 */

#ifndef XMIPP_LIBRARIES_DATA_FOURIER_RECONSTRUCTION_BUFFER_H_
#define XMIPP_LIBRARIES_DATA_FOURIER_RECONSTRUCTION_BUFFER_H_

#include "fourier_reconstruction_traverse_space.h"

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

	int getNoOfElements(TraverseSpace* array) {
		return noOfSpaces;
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

	void setDefault() {
		hasCTFs = hasFFTs = false;
		noOfImages = noOfSpaces = fftSizeX = fftSizeY = paddedImgSize = -1;

		FFTs = CTFs = paddedImages = modulators = NULL;
		spaces = NULL;
	}

protected:
	FRecBufferData() { setDefault(); };
};


#endif /* XMIPP_LIBRARIES_DATA_FOURIER_RECONSTRUCTION_BUFFER_H_ */
