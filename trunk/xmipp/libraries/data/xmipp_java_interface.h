/***************************************************************************
 *
 * Author: Juanjo Vega (juanjo.vega@gmail.com)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medicine
 * Univ. of California, Los Angeles.
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

#ifndef XMIPP_JAVA_INTERFACE
#define XMIPP_JAVA_INTERFACE

#include <string>
#include <iostream>
#include <unistd.h>
#include <error.h>
#include "image.h"
#include "filename.h"
#include "metadata.h"
#include "image.h"
#include "xmipp_java_interface.h"

using namespace std;

std::string readImage(Image<double> &image, const std::string &filename, int readData, int nimage = 0) {
	std::string msg = "";

	try {
		image.read(filename, readData, nimage);
	} catch (XmippError xe) {
		std::cerr << xe;
		msg = filename + ": " + xe.getDefaultMessage();
	} catch (std::exception& e) {
		std::cerr << e.what();
		msg = filename + ": " + e.what();
	} catch (...) {
		std::cerr << "Unhandled exception";
		msg = filename + ": " + "Unhandled exception";
	}

	return msg;
}

std::string readImagePreview(Image<double> &image, const std::string &filename,
		int w, int h, int slice = 0, int nimage = 0) {
	std::string msg = "";

	try {
		image.readPreview(filename, w, h, slice, nimage);

		image().printShape();
	} catch (XmippError xe) {
		std::cerr << xe;
		msg = filename + ": " + xe.getDefaultMessage();
	} catch (std::exception& e) {
		std::cerr << e.what();
		msg = filename + ": " + e.what();
	} catch (...) {
		std::cerr << "Unhandled exception";
		msg = filename + ": " + "Unhandled exception";
	}

	return msg;
}

std::string readMetaData(MetaData &metadata, const std::string &filename) {
	std::string msg = "";

	try {
		metadata.read(filename);
	} catch (XmippError xe) {
		std::cerr << xe;
		msg = filename + ": " + xe.getDefaultMessage();
	} catch (std::exception& e) {
		std::cerr << e.what();
		msg = filename + ": " + e.what();
	} catch (...) {
		std::cerr << "Unhandled exception";
		msg = filename + ": " + "Unhandled exception";
	}

	return msg;
}

std::string writeImage(Image<double> &image, const std::string &filename) {
	std::string msg = "";

	try {
		image.write(filename);
	} catch (XmippError xe) {
		std::cerr << xe;
		msg = filename + ": " + xe.getDefaultMessage();
	} catch (std::exception& e) {
		std::cerr << e.what();
		msg = filename + ": " + e.what();
	} catch (...) {
		std::cerr << "Unhandled exception";
		msg = filename + ": " + "Unhandled exception";
	}

	return msg;
}

/*
 std::string getStrFromValue(MetaData &metadata, const MDLabel MDlabel, std::string &strOut){
 std::string msg = "";

 try {
 std::string field;

 strOut[0] = metadata.getStrFromValue(MDlabel, field);
 } catch (XmippError xe) {
 std::cerr << xe;
 msg = xe.getDefaultMessage();
 } catch (std::exception& e) {
 std::cerr << e.what();
 msg = e.what();
 } catch (...) {
 std::cerr << "Unhandled exception";
 msg = "Unhandled exception";
 }

 return msg;
 }*/

double *getMatrixArray(MultidimArray<double> &array) {
	std::cout << MULTIDIM_ARRAY(array) << std::endl;

	return MULTIDIM_ARRAY(array);
}

double *getImageArray(Image<double> &image) {
	image.data.printShape();
	return getMatrixArray(image.data);
}

double getImageVoxel_TMP(Image<double> &image, int k, int i, int j) {
	std::string msg = "";

	try {
		return image.data.getVoxel(k, i, j);
	} catch (XmippError xe) {
		std::cerr << xe;
		msg = xe.getDefaultMessage();
	} catch (std::exception& e) {
		std::cerr << e.what();
		msg = e.what();
	} catch (...) {
		std::cerr << "Unhandled exception";
		msg = "Unhandled exception";
	}

	return 0;
}

double getMatrixVoxel_TMP(MultidimArray<double> &matrix, int k, int i, int j) {
	std::string msg = "";

	try {
		return matrix.getVoxel(k, i, j);
	} catch (XmippError xe) {
		std::cerr << xe;
		msg = xe.getDefaultMessage();
	} catch (std::exception& e) {
		std::cerr << e.what();
		msg = e.what();
	} catch (...) {
		std::cerr << "Unhandled exception";
		msg = "Unhandled exception";
	}

	std::cout << "MSG: [" << msg << "]" << std::endl;

	return 0;
}

/********************** DEBUGGING ****************************/
std::string testLibraryLink() {
	std::string msg = "Xmipp_Java_Interface TEST";

	std::cout << msg << std::endl;

	return msg;
}

void printData(Image<double> &image){
	image.data.printShape();
}

int getWidth(Image<double> &image){
	int w, h, d = 0;
	long unsigned int n = 0;

	image.getDimensions(w, h, d, n);

	return w;
}

#endif

