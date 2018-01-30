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

#ifndef ARRAY_2D_H_
#define ARRAY_2D_H_

// struct is used for compatibility with OpenCL / C for Cuda
/* Instance of this struct represent a 2D array of arbitrary type.
 * Data are being stored as a dynamic 2D array (i.e. not continuous block of memory).
 * Access to elements will be the fastest if traversing in Y -> X order. */
template<typename T>
struct Array2D {
public:
	/* Empty constructor */
	Array2D() : xSize(0), ySize(0), data(NULL) {};
	/* Constructor, allocates the data immediately */
	Array2D(int xSize, int ySize) :
			xSize(xSize), ySize(ySize) {
		data = new T*[ySize];
		for (int y = 0; y < ySize; y++) {
			data[y] = new T[xSize];
			for (int x = 0; x < xSize; x++) {
				data[y][x] = (T) 0;
			}
		}
	#if DEBUG
		std::cout << "Array2D created (" << getXSize() << "x" << getYSize()
			<< ") at " << data << std::endl;
	#endif
	}

	~Array2D() {
		for (int y = 0; y < ySize; y++) {
			delete[] data[y];
		}
		delete[] data;
	#if DEBUG
		std::cout << "Array2D deleted (" << getXSize() << "x" << getYSize()
			<< ") at " << data << std::endl;
	#endif
		data = 0;
	}

	/* Method to access elements of the array */
	T& operator()(int x, int y) const {
	#if DEBUG
		if (0 > x || x >= xSize)
			std::cout << "Array2D " << data << " X=" << x << " out of range [0.." << getXSize() << ")\n";
		if (0 > y || y >= ySize)
			std::cout << "Array2D " << data << " Y=" << y << " out of range [0.." << getYSize() << ")\n";
		std::cout << std::flush;
	#endif
		return data[y][x];
	}

	int getXSize() const {
		return xSize;
	}

	int getYSize() const {
		return ySize;
	}

	bool inRange(int x, int y) const {
		return inRangeX(x) && inRangeY(y);
	}

	bool inRangeX(int x) const {
		return (x >= 0) && (x < xSize);
	}

	bool inRangeY(int y) const {
		return (y >= 0) && (y < ySize);
	}

	T* getRow(int y) const {
		return data[y];
	}
private:
	int xSize;
	int ySize;
	T** data;
};

#endif /* ARRAY_2D_H_ */

