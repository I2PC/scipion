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

#ifndef XMIPP_LIBRARIES_DATA_RECONSTRUCT_FOURIER_BUFFER_DATA_H_
#define XMIPP_LIBRARIES_DATA_RECONSTRUCT_FOURIER_BUFFER_DATA_H_

#include "reconstruct_fourier_projection_traverse_space.h"

/**
 * Struct holding data for calculating Fourier Reconstruction
 *
 * Contract:
 * 'FFTs' holds 'right side of the centered FFT', i.e. only unique values, with high frequencies in the corners
 * if 'hasFFTs', 'FFTs' contains 'noOfImages' of size 'fftSizeX' * 'fftSizeY'
 * if not 'hasFFTs', 'paddedImages' contains 'noOfImages' of size 'paddedImgSize' * 'paddedImgSize'
 * 'spaces' holds 'noOfImages' * 'noOfSymmetries', which are 'describing' the FFTs/paddedImages
 * 'noOfImages' <= 'maxNoOfImages'
 */
struct RecFourierBufferData
{
	/**
	 * Buffer object constructor.
	 * If 'hasFFTs' is true, 'FFTs' array is allocated
	 * If 'hasCTFs' is true, 'CTFs' and 'modulators' arrays are allocated
	 */
	RecFourierBufferData(bool hasFFTs, bool hasCTFs,
			int fftSizeX, int fftSizeY,	int paddedImgSize,
			int maxNoOfImages, int noOfSymmetries) :
			hasFFTs(hasFFTs), hasCTFs(hasCTFs),
			fftSizeX(fftSizeX), fftSizeY(fftSizeY), paddedImgSize(paddedImgSize),
			maxNoOfImages(maxNoOfImages),noOfSymmetries(noOfSymmetries), noOfImages(0) {
		spaces = new RecFourierProjectionTraverseSpace[maxNoOfImages * noOfSymmetries];
		if (hasFFTs) {
			paddedImages = NULL;
			FFTs = new float[fftSizeX * fftSizeY * maxNoOfImages * 2](); // *2 since it's complex
		} else {
			FFTs = NULL;
			paddedImages = new float[paddedImgSize * paddedImgSize * maxNoOfImages]();
		}

		if (hasCTFs) {
			CTFs = new float[fftSizeX * fftSizeY * maxNoOfImages]();
			modulators = new float[fftSizeX * fftSizeY * maxNoOfImages]();
		} else {
			CTFs = modulators = NULL;
		}
	};

	~RecFourierBufferData() {
		delete[] FFTs;
		delete[] CTFs;
		delete[] paddedImages;
		delete[] modulators;
		FFTs = CTFs = paddedImages = modulators = NULL;

		delete[] spaces;
		spaces = NULL;
	}

	/**
	 * returns no of pixels in one image
	 */
	int getPaddedImgSize() {
		return paddedImgSize * paddedImgSize;
	}

	/**
	 * returns size of one image in bytes
	 */
	int getPaddedImgByteSize() {
		return getPaddedImgSize() * sizeof(float);
	}

	/**
	 * returns size of one FFT in bytes
	 */
	int getFFTsByteSize() {
		return getNoOfElements(FFTs) * sizeof(float);
	}

	/**
	 * returns pointer to the beginning of data for Nth image data in the array OR
	 * NULL if 'array' is not from this buffer
	 * no boundaries check
	 */
	float* getNthItem(float* array, int itemIndex) {
		if (array == FFTs) return array + (fftSizeX * fftSizeY * itemIndex * 2); // *2 since it's complex
		if (array == CTFs) return array + (fftSizeX * fftSizeY * itemIndex);
		if (array == modulators) return array + (fftSizeX * fftSizeY * itemIndex);
		if (array == paddedImages) return array + (paddedImgSize * paddedImgSize * itemIndex);
		return NULL; // undefined
	}

	/**
	 * returns no of objects currently stored in the array OR
	 * -1 if 'array' is not from this buffer
	 */
	int getNoOfElements(float* array) {
		if (array == FFTs) return (fftSizeX * fftSizeY * noOfImages * 2); // *2 since it's complex
		if (array == CTFs) return (fftSizeX * fftSizeY * noOfImages);
		if (array == modulators) return (fftSizeX * fftSizeY * noOfImages);
		if (array == paddedImages) return (paddedImgSize * paddedImgSize * noOfImages);
		return -1; // undefined
	}

	/**
	 * returns no of objects currently stored in the array OR
	 * -1 if 'array' is not from this buffer
	 */
	int getNoOfElements(RecFourierProjectionTraverseSpace* array) {
		return noOfImages * noOfSymmetries;
	}

	/**
	 * returns max size of the 'array' in bytes OR
	 * -1 if 'array' is not from this buffer
	 */
	int getMaxByteSize(float* array) {
		if (array == FFTs) return (fftSizeX * fftSizeY * maxNoOfImages * 2 * sizeof(float)); // *2 since it's complex
		if (array == CTFs) return (fftSizeX * fftSizeY * maxNoOfImages * sizeof(float));
		if (array == modulators) return (fftSizeX * fftSizeY * maxNoOfImages * sizeof(float));
		if (array == paddedImages) return (paddedImgSize * paddedImgSize * maxNoOfImages * sizeof(float));
		return -1; // undefined
	}

	/**
	 * returns max size of the 'array' in bytes OR
	 * -1 if 'array' is not from this buffer
	 */
	int getMaxByteSize(RecFourierProjectionTraverseSpace* array) {
		if (array == spaces) return (noOfSymmetries * maxNoOfImages * sizeof(RecFourierProjectionTraverseSpace));
		return -1; // undefined
	}

	// holds 'right side of the centered FFT', i.e. only unique values, with high frequencies in the corners
	float* FFTs;
	// holds CTF correction for each pixel of FFT
	float* CTFs;
	// holds correction for each pixel of FFT
	float* modulators;
	/**
	 * holds pictures that need to be converted to FFT.
	 * Images are padded, with data shifted to corners (i.e there is 'empty' cross in the middle
	 */
	float* paddedImages;
	// holds spaces for each FFT
	RecFourierProjectionTraverseSpace* spaces;

	int fftSizeX;
	int fftSizeY;
	int paddedImgSize; // paddedImage is squared
	int noOfImages; // currently stored in buffer
	int maxNoOfImages; // that can be stored in the buffer
	int noOfSymmetries;
	/**
	 * if 'true', 'FFTs' contains 'noOfImages' of size 'fftSizeX' x 'fftSizeY'
	 * otherwise 'paddedImages' contains 'noOfImages' of size 'paddedImgSize' * 'paddedImgSize'
	 */
	bool hasFFTs;
	bool hasCTFs; // if 'true', CTF correction data are present in 'CTFs' and 'modulators'

	/**
	 * set values to default (invalid) values
	 * does not allocate
	 * does not release allocated memory
	 */
	void invalidate() {
		hasCTFs = hasFFTs = false;
		noOfImages = fftSizeX = fftSizeY = paddedImgSize
				= maxNoOfImages = noOfSymmetries = -1;

		FFTs = CTFs = paddedImages = modulators = NULL;
		spaces = NULL;
	}

protected:
	RecFourierBufferData() { invalidate(); };
};


#endif /* XMIPP_LIBRARIES_DATA_RECONSTRUCT_FOURIER_BUFFER_DATA_H_ */
