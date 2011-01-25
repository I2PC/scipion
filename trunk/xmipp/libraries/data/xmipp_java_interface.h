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

using namespace std;

std::string readImageHeader(Image<double> &image, const std::string &filename) {
	std::string msg = "";

	try {
		FileName fn(filename);

		image.read2(fn, false);
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

std::string readImage(Image<double> &image, const std::string &filename) {
	std::string msg = "";

	try {
		FileName fn(filename);

		image.read2(fn);
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

std::string readImagePreview(Image<double> &image, const std::string &filename, int w, int h, int slice) {
	std::string msg = "";

	try {
		FileName fn(filename);

		image.readPreview(fn, w, h, slice);
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
		FileName fn(filename);

		metadata.read(fn);
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
		FileName fn(filename);

		image.write(fn);
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

double getImageVoxel(Image<double> &image, int k, int i, int j) {
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

double getMatrixVoxel(MultidimArray<double> &matrix, int k, int i, int j) {
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
// Sometimes java annoys you with "Unsatisfied Link" exceptions. This simple functions
// can be used for testing purposes.
std::string testLibraryLink() {
	std::string msg = "Xmipp_Java_Interface TEST";

	std::cout << msg << std::endl;

	return msg;
}

#endif

