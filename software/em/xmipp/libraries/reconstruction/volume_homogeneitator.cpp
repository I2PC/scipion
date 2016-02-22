/***************************************************************************
 * Authors:     Mohsen Kazemi (mkazemi@cnb.csic.es)
 *
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

#include "volume_homogeneitator.h"

void ProgVolumeHomogeneitator::readParams()
{
    fnVol = getParam("-i");
    fnRef = getParam("-ref");
    fnSetOfImgIn = getParam("-img");
    fnSetOfImgOut = getParam("-o");
}

void ProgVolumeHomogeneitator::defineParams()
{
    addUsageLine("Corrects the set of input images of the input volume with respect to the reference map and using optical flow algorithm");
    addParamsLine(" -i <selfile>         : Selfile with input volume");
    addParamsLine(" -ref <selfile>       : Selfile with input reference volume");
    addParamsLine(" -img <selfile>       : Selfile metadata with input aligned images");
    addParamsLine(" [-o <rootname=\"\">] : Output stk rootname");
    addParamsLine("                      : rootname.xmd contains the list of corrected images with respect to the reference map");
}

// Converts a XMIPP MultidimArray to OpenCV matrix
void ProgVolumeHomogeneitator::xmipp2Opencv(const MultidimArray<double> &xmippArray, cv::Mat &opencvMat)
{
    int h = YSIZE(xmippArray);
    int w = XSIZE(xmippArray);
    opencvMat = cv::Mat::zeros(h, w,CV_32FC1);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
    opencvMat.at<float>(i,j) = DIRECT_A2D_ELEM(xmippArray,i,j);
}

// Converts an OpenCV matrix to XMIPP MultidimArray
void ProgVolumeHomogeneitator::opencv2Xmipp(const cv::Mat &opencvMat, MultidimArray<double> &xmippArray)
{
    int h = opencvMat.rows;
    int w = opencvMat.cols;
    xmippArray.initZeros(h, w);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
    DIRECT_A2D_ELEM(xmippArray,i,j) = opencvMat.at<float>(i,j);
}

void ProgVolumeHomogeneitator::run()
{
	Image<double> inV;
	Image<double> refV;
	MetaData setOfImgIn, setOfImgOut;
	double rot, tilt, psi;
	Projection projIn, projRef;

	inV.read(fnVol);
	inV().setXmippOrigin();

	refV.read(fnRef);
	refV().setXmippOrigin();

	setOfImgIn.read(fnSetOfImgIn);
	setOfImgIn.removeDisabled();

	FOR_ALL_OBJECTS_IN_METADATA (setOfImgIn)
	{
		setOfImgIn.getValue(MDL_ANGLE_ROT, rot, __iter.objId);
		setOfImgIn.getValue(MDL_ANGLE_TILT, tilt, __iter.objId);
		setOfImgIn.getValue(MDL_ANGLE_PSI, psi, __iter.objId);
		projectVolume(inV(), projIn, YSIZE(inV()), XSIZE(inV()), rot, tilt, psi);
		projectVolume(refV(), projRef, YSIZE(refV()), XSIZE(refV()), rot, tilt, psi);
PROJECTIONS ARE READY TO USE IN OF!!!!!!!!!!!

	}
}


