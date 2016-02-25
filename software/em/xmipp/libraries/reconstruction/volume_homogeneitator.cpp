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
    fnVol         = getParam("-i");
    fnRef         = getParam("-ref");
    fnSetOfImgIn  = getParam("-img");
    fnSetOfImgOut = getParam("-o");
    winSize       = getIntParam("--winSize");
    //for test
    fnTestOut1    = getParam("-o1");
    fnTestOut2    = getParam("-o2");
}

void ProgVolumeHomogeneitator::defineParams()
{
    addUsageLine("Corrects the set of input images of the input volume with respect to the reference map and using optical flow algorithm");
    addParamsLine("  -i <selfile>         : Selfile with input volume");
    addParamsLine("  -ref <selfile>       : Selfile with input reference volume");
    addParamsLine("  -img <selfile>       : Selfile metadata with input aligned images");
    addParamsLine(" [-o <rootname=\"\">]  : Output rootname");
    addParamsLine("                       : rootname contains the list of corrected images and their angles with respect to the reference map");
    addParamsLine("                       : NOTE: you do not need to define any extension for the output root name. Program will automatically create one .stk and one .xmd");
    addParamsLine(" [--winSize <int=50>]  : window size for optical flow algorithm");
    //for test
    addParamsLine(" [-o1 <rootname=\"\">]  : Output test .stk input reprojection rootname");
    addParamsLine(" [-o2 <rootname=\"\">]  : Output test .stk reference reprojection rootname");
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

// Converts an OpenCV float matrix to an OpenCV Uint8 matrix
void ProgVolumeHomogeneitator::convert2Uint8(cv::Mat opencvDoubleMat, cv::Mat &opencvUintMat)
{
    double min,max;
    cv::minMaxLoc(opencvDoubleMat, &min, &max);
    opencvDoubleMat.convertTo(opencvUintMat, CV_8U, 255.0/(max - min), -min * 255.0/(max - min));
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
	FileName fn_proj;
	Image<double> inV;
	Image<double> refV;
	MetaData setOfImgIn, setOfImgOut;
	double rot, tilt, psi;
	Projection projIn, projRef, projCorr;
	cv::Mat ProjIn, ProjRef, ProjIn8, ProjRef8;
	cv::Mat flow, ProjCorr;
	cv::Mat planes[]={flow, flow};
	size_t objId;
	int count = 0;
	int projIdx = FIRST_IMAGE;
    FileName stackName = fnSetOfImgOut.removeAllExtensions() + ".stk";
    FileName mdName = fnSetOfImgOut.removeAllExtensions() + ".xmd";

	inV.read(fnVol);
	inV().setXmippOrigin();

	refV.read(fnRef);
	refV().setXmippOrigin();

	setOfImgIn.read(fnSetOfImgIn);
	setOfImgIn.removeDisabled();

	fn_proj = stackName;
	createEmptyFile(fn_proj, XSIZE(inV()), YSIZE(inV()), 1, setOfImgIn.size(), true, WRITE_OVERWRITE);

	//calculating progress time
	init_progress_bar(setOfImgIn.size());
	std::cerr<<"calculating and applying OF algorithm for each new reprojection......\n";

	FOR_ALL_OBJECTS_IN_METADATA (setOfImgIn)
	{
		//get the coordinate of each image for reprojection purpose
		setOfImgIn.getValue(MDL_ANGLE_ROT, rot, __iter.objId);
		setOfImgIn.getValue(MDL_ANGLE_TILT, tilt, __iter.objId);
		setOfImgIn.getValue(MDL_ANGLE_PSI, psi, __iter.objId);

		//reprojection from input and reference volumes to calculate optical flow (OF)
		projectVolume(inV(), projIn, YSIZE(inV()), XSIZE(inV()), rot, tilt, psi);
		projectVolume(refV(), projRef, YSIZE(refV()), XSIZE(refV()), rot, tilt, psi);

		//for test
		projIn.write(fnTestOut1, ALL_IMAGES, true, WRITE_APPEND);
		projRef.write(fnTestOut2, ALL_IMAGES, true, WRITE_APPEND);

		//preparing input data to use for OF algorithm
		xmipp2Opencv(projIn(), ProjIn);
		xmipp2Opencv(projRef(), ProjRef);
		convert2Uint8(ProjIn,ProjIn8);
		convert2Uint8(ProjRef,ProjRef8);

		//OF algorithm
		cv::calcOpticalFlowFarneback(ProjRef8, ProjIn8, flow, 0.5, 6, winSize, 1, 5, 1.1, 0);
		cv::split(flow, planes);

		//for test
		//if (count == 0)
		//{
		//	std::cout<<"count: "<<count<<std::endl;
		//	std::cout<<"planes[0]"<<planes[0]<<"\n"<<std::endl;
		//	std::cout<<"planes[1]"<<planes[1]<<"\n"<<std::endl;
		//}

		for( int row = 0; row < planes[0].rows; row++ )
		    for( int col = 0; col < planes[0].cols; col++ )
		    {
		        planes[0].at<float>(row,col) += (float)col;
		        planes[1].at<float>(row,col) += (float)row;
		    }

		//applying the flow on the input image to use as the corrected one
		cv::remap(ProjIn, ProjCorr, planes[0], planes[1], cv::INTER_CUBIC);

		//preparing output data to use for xmipp
		opencv2Xmipp(ProjCorr, projCorr());

		//preparing output metadata included corrected images and angles
		fn_proj.compose(projIdx, stackName);
		objId = setOfImgOut.addObject();
		setOfImgOut.setValue(MDL_IMAGE, fn_proj, objId);
		setOfImgOut.setValue(MDL_ANGLE_ROT, rot, objId);
		setOfImgOut.setValue(MDL_ANGLE_TILT, tilt, objId);
		setOfImgOut.setValue(MDL_ANGLE_PSI, psi, objId);
		projCorr.write(fn_proj, projIdx, true, WRITE_OVERWRITE);
		projIdx++;

		//calculating progress time
		count++;
		if (count%100 == 0)
			progress_bar(count);
	}
	setOfImgOut.write(mdName);
	progress_bar(setOfImgIn.size());
}


